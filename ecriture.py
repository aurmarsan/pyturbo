try: from paraview import vtk 
except: import vtk
from fonctions_basiques import *
from objets import ObjetPyturbo

try: import vtkIOXMLPython as vtkIOPython
except: 
    try: import vtkIOPython as vtkIOPython
    except: import vtk as vtkIOPython

#__________________________________________________________________________________________
class EcritureV3D(ObjetPyturbo):
    """ecriture de fichiers v3d
    indiquer l'acces_fichier avec l'extension voulue
    pour multiblock, inserer '*' dans acces_fichier qui sera remplacee par le numero de bloc
    
    POUR SAUVEGARDER LE MAILLAGE
        indiquer ecrire_maillage = True
    
    liste des attributs
        * format_binaire    :   0 si formatte, 1 si binaire
        * fmt_fichier       :   "bin" ou "fmt"
        * endian            :   big ou little
        * acces_fichier            :   str pour monobloc
                                [str, str] pour multiblocs
        * donnees_aux_cellules :    0 si donnees au noeuds, 1 si aux cellules
    """
    
    #_________________________________________________________________________________
    def __init__(self, input=None, acces_fichier=None, fmt_fichier='bin', precision='i4r8',
            endian='big', 
            ecrire_maillage=False,
            donnees_aux_cellules = False, 
            numbloc=1):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def get(self, nom_attribut):
        """definit ce qui sera affiche par print
        
        """
        return getattr(self, nom_attribut)
    #__________________________________________________________________________________

    #__________________________________________________________________________________
    def ecriture_monobloc(self, acces_fichier, vtkDataObject, numbloc):
        """ecriture d'un seul bloc
        
        """
        dimensions_bloc = vtkDataObject.GetDimensions()
        # ecriture des coordonnees des points du maillage
        if self.ecrire_maillage:
            data = numpy_support.vtk_to_numpy(vtkDataObject.GetPoints().GetData())
            dict_numpy_arrays = {
                'x': data[:, 0], 
                'y': data[:, 1], 
                'z': data[:, 2]
                }
        # ecriture des donnees aux points
        elif not self.donnees_aux_cellules:
            dict_numpy_arrays = {}
            for numarray in range(vtkDataObject.GetPointData().GetNumberOfArrays()):
                nom_array = vtkDataObject.GetPointData().GetArrayName(numarray)
                data = numpy_support.vtk_to_numpy(vtkDataObject.GetPointData().GetArray(numarray))
                if len(data.shape) == 2:
                    for composante in range(data.shape[1]):
                        dict_numpy_arrays[nom_array + '_{0}'.format(composante)] = \
                            data[:, composante]
                else:
                    dict_numpy_arrays[nom_array] = \
                            data
        # ecriture des donnees aux cellules
        else:
            dict_numpy_arrays = {}
            dimensions_bloc = list([
                dimensions_bloc[k] - 1 if dimensions_bloc[k] > 1 else dimensions_bloc[k]
                for k in range(3)]
                )
            for numarray in range(vtkDataObject.GetCellData().GetNumberOfArrays()):
                nom_array = vtkDataObject.GetCellData().GetArrayName(numarray)
                data = numpy_support.vtk_to_numpy(vtkDataObject.GetCellData().GetArray(numarray))
                if len(data.shape) == 2:
                    for composante in range(data.shape[1]):
                        dict_numpy_arrays[nom_array + '_{0}'.format(composante)] = \
                            data[:, composante]
                else:
                    dict_numpy_arrays[nom_array] = \
                            data
        
        ecrire_v3d(acces_fichier = acces_fichier, dict_numpy_arrays = dict_numpy_arrays, 
            dimensions = dimensions_bloc, 
            numbloc = numbloc, fmt_fichier = self.fmt_fichier, precision = self.precision, 
            endian=self.endian, type_maillage = self.ecrire_maillage)
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def update(self):
        """execute l'ecriture des fichiers
        
        """
        # ecriture monobloc
        if isinstance(self.input, vtk.vtkStructuredGrid):
            self.ecriture_monobloc(acces_fichier = self.acces_fichier, 
                vtkDataObject = self.input, numbloc = self.numbloc)
        
        # ecriture multiblocs
        elif isinstance(self.input, vtk.vtkMultiBlockDataSet):
            if not'*' in self.acces_fichier:
                raise IOError, "input est un fichier multiblock, mais il n'y a pas * dans l'acces_fichier indique"
            for numbloc in get_numeros_blocs_non_vides(self.input):
                self.ecriture_monobloc(
                    acces_fichier = self.acces_fichier.replace('*', str(numbloc)), 
                    vtkDataObject = self.input.GetBlock(numbloc), numbloc = numbloc)
        # non supporte
        else:
            raise IOError, "l'objet a ecrire doit etre soit un vtkStructuredGrid, soit un vtkMultiBlockDataSet"
    #_________________________________________________________________
    
    #_________________________________________________________________
    def ecrire(self):
        """effectue la sauvegarde
        
        """
        self.update()
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class EcritureVTK(ObjetPyturbo):
    """sauvegarder au format vtk
    multiblockdataset ou structuredgrid ou polydata
    
    donner acces_fichier sans extension .vtm, .vtp, .vts
    
    """
    #_________________________________________________________________
    def __init__(self, input=None, acces_fichier=None):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_________________________________________________________________

    #_________________________________________________________________
    def update(self):
        # verifications prealables
        if self.acces_fichier is None:
            raise IOError, "indiquez l'acces_fichier d'abord"
        if self.input is None:
            raise IOError, "input is None"
        # cas multiblockdataset
        if isinstance(self.input, vtk.vtkMultiBlockDataSet):
            w = vtkIOPython.vtkXMLMultiBlockDataWriter()
            w.SetFileName(self.acces_fichier + '.vtm')
        elif isinstance(self.input, vtk.vtkPolyData):
            w = vtkIOPython.vtkXMLPolyDataWriter()
            w.SetFileName(self.acces_fichier + '.vtp')
        elif isinstance(self.input, vtk.vtkStructuredGrid):
            w = vtkIOPython.vtkXMLStructuredGridWriter()
            w.SetFileName(self.acces_fichier + '.vts')
        elif isinstance(self.input, vtk.vtkUnstructuredGrid):
            w = vtkIOPython.vtkXMLUnstructuredGridWriter()
            w.SetFileName(self.acces_fichier + '.vtu')
        vtk_set_input(w, self.input)
        print "Ecriture {0}".format(w.GetFileName())
        # w.Update()
        w.Write()
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def ecrire(self):
        """effectue la sauvegarde
        
        """
        self.update()
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class EcritureNumpy(ObjetPyturbo):
    """ecriture de fichiers binaires numpy
    indiquer l'acces_fichier. Il sera complete automatiquement avec le numero de bloc
    et le nom adequat de la variable. 
    
    POUR SAUVEGARDER LE MAILLAGE en plus des donnees
        indiquer ecrire_maillage = True
    
    liste des attributs
        * ecrire_maillage          :   pour ecrire le maillage
        * acces_fichier            :   str pour monobloc
                                [str, str] pour multiblocs
    """
    
    #_________________________________________________________________________________
    def __init__(self, input=None, acces_fichier=None, ecrire_maillage=False):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def get(self, nom_attribut):
        """definit ce qui sera affiche par print
        
        """
        return getattr(self, nom_attribut)
    #__________________________________________________________________________________

    #__________________________________________________________________________________
    def ecriture_maillage_monobloc(self, vtkDataObject):
        """ecriture de la geometrie d'un monobloc
        
        """
        dict_numpy_arrays = {}
        # ecriture des coordonnees des points du maillage
        # if self.ecrire_maillage:
        xyz = numpy_support.vtk_to_numpy(vtkDataObject.GetPoints().GetData())
        dict_numpy_arrays['xyz'] = xyz
        
        if isinstance(vtkDataObject, vtk.vtkPolyData):
            polys = get_vtk_array_as_numpy_array(vtkDataObject, 'polys')
            dict_numpy_arrays['polys'] = polys
            dict_numpy_arrays['nb_polys'] = vtkDataObject.GetNumberOfPolys()
        else:
            cells = get_vtk_array_as_numpy_array(vtkDataObject, 'cells')
            dict_numpy_arrays['cells'] = cells
            
            cellstypes = get_vtk_array_as_numpy_array(vtkDataObject, 'cellstypes')
            dict_numpy_arrays['cellstypes'] = cellstypes
            
            cellslocations = get_vtk_array_as_numpy_array(vtkDataObject, 'cellslocations')
            dict_numpy_arrays['cellslocations'] = cellslocations
            
        for nom_array in dict_numpy_arrays:
            print "Ecriture {0}_{1}".format(self.acces_fichier, nom_array)
            numpy.save("{0}_{1}".format(self.acces_fichier, nom_array), dict_numpy_arrays[nom_array])
        
        return 0
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def ecriture_donnees_monobloc(self, vtkDataObject):
        """ Ecriture des donnees aux points et aux cellules d'un vtkDataObject
        
        """
        dict_numpy_arrays = {}
        for nom_array in get_noms_arrays_presents(vtkDataObject, 'points'):
            data = get_vtk_array_as_numpy_array(vtkDataObject, nom_array, loc = 'points')
            dict_numpy_arrays[nom_array + '_nodes'] = data
        for nom_array in get_noms_arrays_presents(vtkDataObject, 'cells'):
            data = get_vtk_array_as_numpy_array(vtkDataObject, nom_array, loc = 'cells')
            dict_numpy_arrays[nom_array + '_cells'] = data
        
        for nom_array in dict_numpy_arrays:
            print "Ecriture {0}_{1}".format(self.acces_fichier, nom_array)
            numpy.save("{0}_{1}".format(self.acces_fichier, nom_array), 
                dict_numpy_arrays[nom_array])
        
        return 0
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def update(self):
        """execute l'ecriture des fichiers
        
        """
        # ecriture multiblocs
        if isinstance(self.input, vtk.vtkMultiBlockDataSet):
            for numbloc in get_numeros_blocs_non_vides(self.input):
                EcritureNumpy(input = self.input.GetBlock(numbloc), 
                    acces_fichier = "{0}_bloc{1}".format(self.acces_fichier, numbloc), 
                    ecrire_maillage = self.ecrire_maillage).ecrire()
        # ecriture monobloc
        else:
            # ecriture du maillage
            if self.ecrire_maillage:
                self.ecriture_maillage_monobloc(vtkDataObject = self.input)
            
            # ecriture des donnees
            self.ecriture_donnees_monobloc(self.input)
    #_________________________________________________________________
    
    #_________________________________________________________________
    def ecrire(self):
        """effectue la sauvegarde
        
        """
        self.update()
    #_________________________________________________________________
#__________________________________________________________________________________________
