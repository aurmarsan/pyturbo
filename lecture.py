import numpy
import struct
import glob
import sys
from scipy.io import netcdf


try: 
    from paraview import vtk 
    from paraview.vtk import vtkFiltersExtraction
except: 
    import vtk
    import vtk as vtkFiltersExtraction
try: from paraview import numpy_support 
except: from vtk.util import numpy_support



try :
    import vtkIOXMLPython
    import vtkIOEnSightPython
except:
    try: import vtkIOPython
    except: pass

from fonctions_basiques import *
from fonctions_basiques import rotation
from calculs import *
from objets import ObjetPyturbo

#__________________________________________________________________________________________
class LecteurV3D(ObjetPyturbo):
    """lecture de fichiers v3d
    retourne un objet VTK
    
    le numero du bloc est celui indique dans le nom du fichier (a la place de *)
    
    liste des attributs
        * format_binaire    :   0 si formatte, 1 si binaire
        * fmt_binaire       :   si binaire, i4r8 par defaut
        * endian            :   big ou little
        * acces_fichier            :   utiliser "*" pour lecture multi-blocs
            possibilite d'indiquer ['mai*', 'data*'] 
            pour lire les fichiers solutions dans la foulee 
        * decalage_numbloc_lus  :   utile si les numeros des blocs ecrits dans les 
            fichiers V3D sont decales par rapport a ce qu'ils devraient etre
            (decalage = 1 pour la chaine elsA v2 TM par exemple)
            !!! Voir remarque plus bas sur les numeros des blocs !!! 
        * assembler_vecteurs :  determine si les trois composantes d'un vecteur
            doivent etre assemblees. Dans ce cas 
                - ['rou' 'rov' 'row'] sont assembles en 'momentum'
                - ['nom_0' 'nom_1' 'nom_2'] sont assembles en 'nom'
                    pour tout nom quelconque
        * si traitement_non_numeros est True, 
            alors dans le cas ou le string lu a la place de "*" n'est pas un numero, 
            le block "sans numero" lu est rajoute a la fin
                    pour tout nom quelconque
        * si seulement les numeros est du type [1, 2, 3]
            alors lecture seulement des fichiers pour lesquels "*" correspond a un des numeros
            indiques dans la liste.
    
    POUR LA LECTURE MULTIBLOCS, C'EST LE NUMERO DANS LE NOM DE FICHIER A LA PLACE DE * QUI EST UTILISE ! 
    Voir la fonction update dans le cas multibloc
    
    Tout est par contre cable pour pouvoir utiliser le numero du bloc lu dans le fichier V3D
    
    """
    
    #_________________________________________________________________________________
    def __init__(self, fmt_fichier='bin', precision='i4r8',
            endian='big', acces_fichier=None, 
            compter_saut_de_ligne=False, 
            decalage_numbloc_lus=0, 
            assembler_vecteurs=True, 
            traitement_non_numeros=False, 
            seulement_les_numeros=None):
        attributs = locals().copy()
        del attributs['self']
        self._mettre_a_jour = True
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
        
        #sortie du lecteur
        self.output = None
    #__________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def set(self, nom_attribut, valeur):
        """fonction set specifique
        
        gere la variable locale _changement
        qui sert lorsque l'on appelle la sortie
        a savoir s'il faut recalculer
        """
        setattr(self, nom_attribut, valeur)
        if nom_attribut != '_mettre_a_jour':
            self._mettre_a_jour = True
    #_____________________________________________________________________________________
    
    #__________________________________________________________________________________
    def reset_output(self):
        """definit ce qui sera affiche par print
        
        """
        self.output = None
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def importer_maillage(self, vtkDataObject):
        """importation d'un maillage existant
        
        les donnees deja presentes aux points et aux noeuds sont conservees
        
        """
        self.output = vtk_new_shallowcopy(vtkDataObject)
        
        #nettoyage pour ne pas importer les numeros non demandes
        if self.seulement_les_numeros is not None:
            for numbloc in get_numeros_blocs_non_vides(self.output):
                if numbloc not in self.seulement_les_numeros:
                    self.output.SetBlock(numbloc, None)
            
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def ajouter_au_bloc(self, VTKbloc, vtkArray):
        """ajoute un vtkArray contenant des donnees au bloc VTK
        soit dans les PointsData soit dans les CellData
        en fonction du nombre de tuples du vtkArray"""
        NumberOfTuples = vtkArray.GetNumberOfTuples()
        if NumberOfTuples == VTKbloc.GetNumberOfPoints():
            VTKbloc.GetPointData().AddArray(vtkArray)
        elif NumberOfTuples == VTKbloc.GetNumberOfCells():
            VTKbloc.GetCellData().AddArray(vtkArray)
        else:
            raise ValueError, "le nombre de tuples ne correspond ni \
            au nombre de points, ni au nombre de cellules"
        return 0
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def lire_monobloc(self, acces_fichier, numbloc=None):
        """lecture d'un seul bloc
        
        si numbloc est indique, alors le bloc est stocke comme numbloc 
        et le numero de bloc lu dans le v3d est ignore
        
        """
        
        dictionnaire_lecture = lire_v3d(acces_fichier = acces_fichier, fmt_fichier = self.fmt_fichier, 
            endian= self.endian, precision = self.precision,
            compter_saut_de_ligne = self.compter_saut_de_ligne)
        if numbloc is None:
            numbloc = dictionnaire_lecture['numbloc'] + self.decalage_numbloc_lus
        
        try:
            bloc_temp = vtk_new_shallowcopy(self.output if isinstance(self.output, vtk.vtkStructuredGrid)
                else self.output.GetBlock(numbloc))
        except:
            bloc_temp = vtk.vtkStructuredGrid()
        
        # liste des donnees traitees
        donnees_traitees = []
        
        # GEOMETRIE
        # si des coordonnees (x,y,z) sont comprises dans le fichier lu
        if dictionnaire_lecture['data'].has_key('x') \
                and dictionnaire_lecture['data'].has_key('y') \
                and dictionnaire_lecture['data'].has_key('z'):
            numpyArrayCoords = numpy.vstack((dictionnaire_lecture['data']['x'],
                                        dictionnaire_lecture['data']['y'],
                                        dictionnaire_lecture['data']['z']))\
                                    .transpose(1,0)
            vtkArray = numpy_support.numpy_to_vtk(
                                numpy.ascontiguousarray(numpyArrayCoords)
                                , deep = 1)
            points = vtk.vtkPoints()
            points.SetData(vtkArray)
            bloc_temp.SetDimensions(dictionnaire_lecture['dims'][0], 
                                        dictionnaire_lecture['dims'][1], 
                                        dictionnaire_lecture['dims'][2])
            bloc_temp.SetPoints(points)
            donnees_traitees += ['x', 'y', 'z']
        
        # DONNEES
        if not bloc_temp is None and bloc_temp.GetNumberOfPoints() != 0 :
            # si il y a rou rov row
            if dictionnaire_lecture['data'].has_key('rou') \
                    and dictionnaire_lecture['data'].has_key('rov') \
                    and dictionnaire_lecture['data'].has_key('row') \
                    and self.assembler_vecteurs is True:
                momentum = numpy.vstack((dictionnaire_lecture['data']['rou'],
                                        dictionnaire_lecture['data']['rov'],
                                        dictionnaire_lecture['data']['row'])).transpose(1, 0)
                vtkArray = numpy_support.numpy_to_vtk(
                                    numpy.ascontiguousarray(momentum)
                                    , deep = 1)
                vtkArray.SetName('momentum')
                self.ajouter_au_bloc(bloc_temp, vtkArray)
                donnees_traitees += ['rou', 'rov', 'row']
            # si il y a trois donnees (et seulement trois) du type *_0 *_1 *_2
            liste_temp = [nom[:-2] for nom in dictionnaire_lecture['data'].keys()]
            for key in liste_temp:
                if numpy.all([key + '_{0}'.format(composante) in dictionnaire_lecture['data'].keys()
                        for composante in range(3)]) and \
                        numpy.all([key + '_{0}'.format(composante) 
                        not in donnees_traitees for composante in range(3)]) \
                        and self.assembler_vecteurs is True:
                    vector = numpy.vstack((
                                dictionnaire_lecture['data'][key + '_0'],
                                dictionnaire_lecture['data'][key + '_1'],
                                dictionnaire_lecture['data'][key + '_2'])).transpose(1, 0)
                    vtkArray = numpy_support.numpy_to_vtk(
                                    numpy.ascontiguousarray(vector)
                                    , deep = 1)
                    vtkArray.SetName(key)
                    self.ajouter_au_bloc(bloc_temp, vtkArray)
                    donnees_traitees += [key + '_%s'%(composante) for composante in range(3)]
            for key in dictionnaire_lecture['data']:
                if (key not in donnees_traitees):
                    vtkArray = numpy_support.numpy_to_vtk(
                                        dictionnaire_lecture['data'][key]
                                        , deep = 1)
                    vtkArray.SetName(key)
                    self.ajouter_au_bloc(bloc_temp, vtkArray)
                    donnees_traitees += [key]
        else:
            print "## IGNORE -- LA GEOMETRIE N'EST PAS CHARGEE"
            bloc_temp = None
        
        if isinstance(self.output, vtk.vtkStructuredGrid):
            self.output = bloc_temp
        else:
            self.output.SetBlock(numbloc, bloc_temp)
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def update(self):
        """execute la lecture des fichiers
        
        """
        for acces_fichier in self.acces_fichier if isinstance(self.acces_fichier, list) \
                else [self.acces_fichier]:
            # lecture monobloc
            if not '*' in acces_fichier:
                if self.output is None: self.output = vtk.vtkStructuredGrid()
                self.lire_monobloc(acces_fichier = acces_fichier)
            
            # lecture multiblocs
            else:
                if self.output is None: self.output = vtk.vtkMultiBlockDataSet()
                fichiers_non_numero = []
                for fich in glob.glob(acces_fichier):
                    debut_acces_fichier = acces_fichier.split('*')[0]
                    fin_acces_fichier = acces_fichier.split('*')[1]
                    try:
                        numbloc = int(
                            fich[len(debut_acces_fichier):-len(fin_acces_fichier)]
                            ) if len(fin_acces_fichier) != 0 else int(
                            fich[len(debut_acces_fichier):]
                            )
                    except:
                        if self.traitement_non_numeros is True:
                            fichiers_non_numero.append(fich)
                    else:
                        if self.seulement_les_numeros is None or numbloc in self.seulement_les_numeros:
                            self.lire_monobloc(acces_fichier = fich, numbloc = numbloc)
                for fich in fichiers_non_numero:
                    debut_acces_fichier = acces_fichier.split('*')[0]
                    fin_acces_fichier = acces_fichier.split('*')[1]
                    numbloc = self.output.GetNumberOfBlocks()
                    self.lire_monobloc(acces_fichier = fich, numbloc = numbloc)
        self._mettre_a_jour = False
    #__________________________________________________________________________________
    
    #__________________________________________________________________________________
    def get_output(self):
        """execute la lecture des fichiers
        
        """
        if self._mettre_a_jour:
            self.update()
        return self.output
    #__________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
class LecteurVTK(ObjetPyturbo):
    """lecture au format vtk
    multiblockdataset ou structuredgrid ou polydata
    
    en fonction de l'extension .vtm, .vtp, .vts
    
    possibilite de lire plusieurs .vtp ou .vts et creer un multiblock a partir de ceux-ci
    en indiquant une etoile dans acces_fichier
    ex : acces_fichier = maillage/h*.vts retourne un multiblock 
    
    """
    #_________________________________________________________________
    def __init__(self, acces_fichier=None):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
        
        #Initialisation des attributs
        self._mettre_a_jour = True
        self.output = None
    #_________________________________________________________________
    
    #_____________________________________________________________________________________
    def set(self, nom_attribut, valeur):
        """fonction set specifique
        
        gere la variable locale _changement qui sert lorsque l'on appelle la sortie 
        a savoir s'il faut relire
        """
        setattr(self, nom_attribut, valeur)
        if nom_attribut != '_mettre_a_jour':
            self._mettre_a_jour = True
    #_____________________________________________________________________________________

    #_________________________________________________________________
    def update(self):
        print "lecture {0}".format(self.acces_fichier)
        # cas ou acces_fichier contient une etoile
        if "*" in self.acces_fichier:
			maillage = vtk.vtkMultiBlockDataSet()
			for fichier in glob.glob(self.acces_fichier):
				print fichier
				bloc = LecteurVTK(
					acces_fichier = fichier
					).get_output()
				numero_bloc = int(
					fichier.split(self.acces_fichier.split('*')[0])[1].split(
					self.acces_fichier.split('*')[1])[0]
					)
				maillage.SetBlock(numero_bloc, bloc)
			self.output = maillage
			self._mettre_a_jour = False
		
        #appel des lecteurs vtk standards
        else:
            try: 
                w = vtkIOXMLPython.vtkXMLGenericDataObjectReader() 
            except: 
                try: 
                    w = vtkIOPython.vtkXMLGenericDataObjectReader()
                except:
                    if self.acces_fichier[-4:] == ".vtm":
                        w = vtk.vtkXMLMultiBlockDataReader()
                        w.SetFileName(self.acces_fichier)
                    elif self.acces_fichier[-4:] == ".vtp":
                        w = vtk.vtkXMLPolyDataReader()
                        w.SetFileName(self.acces_fichier)
                    elif self.acces_fichier[-4:] == ".vts":
                        w = vtk.vtkXMLStructuredGridReader()
                        w.SetFileName(self.acces_fichier)
                    elif self.acces_fichier[-4:] == ".vtu":
                        w = vtk.vtkXMLUnstructuredGridReader()
                        w.SetFileName(self.acces_fichier)
                    else:
                        w = vtk.vtkDataSetReader()
                        w.SetFileName(self.acces_fichier)()
            w.SetFileName(self.acces_fichier)
            w.Update()
            self.output = w.GetOutput()
            self._mettre_a_jour = False
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne l'objet vtk lu
        
        """
        if self._mettre_a_jour:
            self.update()
        return self.output
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class LecteurEnSight(ObjetPyturbo):
    """lecture au format EnSight
    Pour l'instant seulement les .case
    
    Si <volume_seulement> est True, alors seuls les blocs tri-dimensionnels du multibloc sont conserves
    Si renommer_pour_pyturbo est True, alors les variables sont renommees pour correspondre a pyturbo
    Si un fichier out CFX est donne, alors on va y chercher le omega
    """
    #_________________________________________________________________
    def __init__(self, acces_fichier=None, volume_seulement=True):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
        
        #Initialisation des attributs
        self.output = None
        self._create_vtkreader()
        self.noms_variables_noeuds = None
        self.noms_variables_cellules = None
    #_________________________________________________________________
    
    #_____________________________________________________________________________________
    def set(self, nom_attribut, valeur):
        """fonction set specifique
        
        gere la variable locale _changement qui sert lorsque l'on appelle la sortie 
        a savoir s'il faut relire
        """
        setattr(self, nom_attribut, valeur)
        
        #self.noms_domaines = []
    #_____________________________________________________________________________________

    #_________________________________________________________________
    def get_noms_variables_aux_noeuds(self):
        """ retourne une liste contenant les noms des variables aux noeuds
        """
        self.noms_variables_noeuds = []
        for k in range(self.vtkreader.GetNumberOfPointArrays()):
            self.noms_variables_noeuds.append(self.vtkreader.GetPointArrayName(k))
            
        return self.noms_variables_noeuds
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_noms_variables_aux_cellules(self):
        """ retourne une liste contenant les noms des variables aux noeuds
        """
        self.noms_variables_cellules = []
        for k in range(self.vtkreader.GetNumberOfCellArrays()):
            self.noms_variables_cellules.append(self.vtkreader.GetCellArrayName(k))
        
        return self.noms_variables_cellules
    #_________________________________________________________________
    
    #_________________________________________________________________
    def set_noms_variables_aux_noeuds_a_lire(self, liste_noms):
        """pour n'indiquer que certaines variables a lire
        
        """
        for nom in self.get_noms_variables_aux_noeuds():
            if nom not in liste_noms:
                self.vtkreader.SetPointArrayStatus(nom, 0)
            else:
                self.vtkreader.SetPointArrayStatus(nom, 1)
    #_________________________________________________________________
    
    #_________________________________________________________________
    def set_noms_variables_aux_cellules_a_lire(self, liste_noms):
        """pour n'indiquer que certaines variables a lire
        
        """
        for nom in self.get_noms_variables_aux_cellules():
            if nom not in liste_noms:
                self.vtkreader.SetCellArrayStatus(nom, 0)
            else:
                self.vtkreader.SetCellArrayStatus(nom, 1)
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_noms_variables_aux_noeuds_a_lire(self):
        """pour n'indiquer que certaines variables a lire
        
        """
        liste = []
        for nom in self.get_noms_variables_aux_noeuds():
            if self.vtkreader.GetPointArrayStatus(nom) == 1:
                liste.append(nom)
        return liste
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_noms_variables_aux_cellules_a_lire(self):
        """pour n'indiquer que certaines variables a lire
        
        """
        liste = []
        for nom in self.get_noms_variables_aux_cellules():
            if self.vtkreader.GetCellArrayStatus(nom) == 1:
                liste.append(nom)
        return liste
    #_________________________________________________________________
    
    #_________________________________________________________________
    def _create_vtkreader(self):
        try: self.vtkreader = vtkIOEnSightPython.vtkGenericEnSightReader() 
        except: 
            try: self.vtkreader = vtkIOPython.vtkGenericEnSightReader()
            except: self.vtkreader = vtk.vtkGenericEnSightReader()
        self.vtkreader.SetCaseFileName(self.acces_fichier)
        self.vtkreader.UpdateInformation()
        
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne l'objet vtk lu
        
        """
        print "lecture {0}".format(self.acces_fichier)
        print "variables aux noeuds {0}".format(self.get_noms_variables_aux_noeuds_a_lire())
        print "variables aux cellules {0}".format(self.get_noms_variables_aux_cellules_a_lire())
        self.vtkreader.Update()
        multibloc = self.vtkreader.GetOutput()
        
        if self.volume_seulement:
            #dans ce cas on ne garde que les blocs dont la premiere cellule est tridimensionnelle
            a_supprimer = []
            for numbloc in get_numeros_blocs_non_vides(self.vtkreader.GetOutput()):
                if multibloc.GetBlock(numbloc).GetCell(0).GetCellDimension() != 3:
                    a_supprimer.append(numbloc)
            a_supprimer.sort()
            for numbloc in a_supprimer[::-1]:
                multibloc.RemoveBlock(numbloc)
        
        self.output = multibloc
        
        return self.output
    #_________________________________________________________________    
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class LecteurNumpy(ObjetPyturbo):
    """lecture d'un numpy binaire contenant des donnees. 
    
    Si un fichier xyz est donne, alors un fichier de connectivite doit aussi etre donne. 
    (extension a faire : si un fichier de connectivite n'est pas donne, alors lecture d'un bloc structure
    et demander les dimensions)
    cells_file doit etre un numpy array dont le premier element est le nombre de cellules, puis les 
    suivants decrivent les cellules. 
    Un bloc non-structure est cree avec ca. 
    Si les noms contiennent une *, alors la lecture passe en mode multibloc, et * est remplace par le numero du bloc. 
    La fonction retourne alors un multibloc non-structure
    
    Si arrays_files est indique, alors indiquer egalement arrays_noms pour indiquer les noms sous lesquels doivent etre ajoutes
    l'array. ET : 
        - soit la geometrie est aussi lue
        - soit un vtkDataObject est donnee en entree du lecteur. Il peut etre structure ou non, peu importe. 
    L'array doit contenir l'ensemble des donnees aux points ou aux cellules, pour tout le multibloc ! Dans 
    l'ordre des points ou cellules, compte-tenu de l'ordre des blocs. 
    
    arrays_files et arrays_noms doivent etre des LISTES !
    
    seulement_les_numeros permet de ne lire que certains blocs dans le cas d'un multibloc
    
    pour un polydata (surfacique), donner polys_file a la place des cellules. 
    
    """
    #_________________________________________________________________
    def __init__(self, 
            arrays_files=None, 
            arrays_noms=None, 
            xyz_file=None, 
            cells_file=None, 
            cellstypes_file=None, 
            cellslocations_file=None, 
            polys_file=None,
            nb_polys_file=None,
            vtkDataObject=None,
            seulement_les_numeros=None):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_________________________________________________________________

    #_________________________________________________________________
    def update(self):
        # mode multibloc avec lecture de la geometrie
        if self.xyz_file is not None and "*" in self.xyz_file:
            self.output = vtk.vtkMultiBlockDataSet()
            
            for xyz_file in glob.glob(self.xyz_file):
                # on determine le numero du bloc
                debut_acces_fichier = self.xyz_file.split('*')[0]
                fin_acces_fichier = self.xyz_file.split('*')[1]
                numero_bloc = xyz_file.replace(debut_acces_fichier, '').replace(fin_acces_fichier, '')
                if self.cells_file is not None:
                    cells_file = self.cells_file.replace("*", numero_bloc)
                    cellstypes_file = self.cellstypes_file.replace("*", numero_bloc)
                    cellslocations_file = self.cellslocations_file.replace("*", numero_bloc)
                    polys_file = None
                    nb_polys_file = None
                else:
                    cells_file = None
                    cellstypes_file = None
                    cellslocations_file = None
                    polys_file = self.polys_file.replace("*", numero_bloc)
                    nb_polys_file = self.nb_polys_file.replace("*", numero_bloc)
                
                if self.seulement_les_numeros is None or int(numero_bloc) in self.seulement_les_numeros:
                    # on lit le bloc et on l'ajoute
                    self.output.SetBlock(
                        int(numero_bloc), 
                        LecteurNumpy(arrays_files = [array_file.replace('*', numero_bloc) for array_file in self.arrays_files]
                                if self.arrays_files is not None else None, 
                            arrays_noms = self.arrays_noms, 
                            xyz_file = xyz_file, 
                            cells_file = cells_file, 
                            cellstypes_file = cellstypes_file, 
                            cellslocations_file = cellslocations_file, 
                            polys_file = polys_file,
                            nb_polys_file = nb_polys_file,
                            vtkDataObject = None
                            ).get_output()
                            )
            return 0
        
        # mode multibloc lecture des donnees seulement
        elif self.arrays_files is not None and "*" in self.arrays_files[0]:
            self.output = vtk.vtkMultiBlockDataSet()
            
            for array_file in glob.glob(self.arrays_files[0]):
                # on determine le numero du bloc
                debut_acces_fichier = self.arrays_files[0].split('*')[0]
                fin_acces_fichier = self.arrays_files[0].split('*')[1]
                numero_bloc = array_file.replace(debut_acces_fichier, '').replace(fin_acces_fichier, '')
                
                # on verifie que c'est bien seulement un numero, sinon on saute
                try:
                    int(numero_bloc)
                except:
                    pass
                else:
                    # on lit le bloc et on l'ajoute, si la geometrie est presente
                    if int(numero_bloc) in get_numeros_blocs_non_vides(self.vtkDataObject) \
                            and (self.seulement_les_numeros is None or int(numero_bloc) in self.seulement_les_numeros):
                        self.output.SetBlock(
                            int(numero_bloc), 
                            LecteurNumpy(arrays_files = [array_file.replace('*', numero_bloc) for array_file in self.arrays_files], 
                                arrays_noms = self.arrays_noms, 
                                xyz_file = None, 
                                cells_file = None, 
                                cellstypes_file = None, 
                                cellslocations_file = None, 
                                polys_file = None,
                                nb_polys_file = None, 
                                vtkDataObject = self.vtkDataObject.GetBlock(int(numero_bloc))
                                ).get_output()
                                )
            return 0
            
        # MONOBLOC
        # verifications prealables
        if self.xyz_file is None and self.vtkDataObject is None:
            raise IOError, "pas de geometrie indiquee en entree !"
        
        if self.xyz_file is not None and self.vtkDataObject is not None:
            raise IOError, "il faut soit indiquer un maillage <vtkDataObject> ou une geometrie a lire <xyz_file>"
        
        if self.arrays_files is not None and self.arrays_noms is None:
            raise IOError, "indiquez le nom de l'array a utiliser pour le stockage dans l'objet VTK"
        
        # Geometrie
        if self.xyz_file is not None:
            if self.polys_file is None and (self.cells_file is None or self.cellstypes_file is None or self.cellslocations_file is None):
                raise IOError, "L'arbre de connectivite, le type des cellules et leur localisation dans l'array des cellules doit aussi etre indique."
                
            print 'Lecture {0}'.format(self.xyz_file)
            self.vtkDataObject = self.create_geometrie_monobloc(self.xyz_file, self.cells_file, 
                self.cellstypes_file, self.cellslocations_file, self.polys_file, self.nb_polys_file)
        
        # chargement des donnees aux points et aux cellules
        self.output = vtk_new_shallowcopy(self.vtkDataObject)
        if self.arrays_files is not None:
            if self.arrays_noms is None or len(self.arrays_files) != len(self.arrays_noms):
                raise IOError, 'Indiquer les noms des arrays'
        
            for indice in range(len(self.arrays_files)):
                array_file = self.arrays_files[indice]
                nom_array = self.arrays_noms[indice]
                print 'Lecture {0}'.format(array_file)
                numpy_array = numpy.load(array_file)
                self.output = ajouter_numpy_array_as_vtk_array(self.output, numpy_array, nom_array)
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def create_geometrie_monobloc(self, xyz_file, cells_file, cellstypes_file, cellslocations_file, polys_file, nb_polys_file):
        """fonction de creation d'une geometrie. Soit 3D, soit 2D PolyData
        """
        # Lecture du fichier geometrie et de l'arbre de connectivite
        xyz = numpy.load(xyz_file)
        if cells_file is not None:
            connect = numpy.load(cells_file)
            cellstypes = numpy.load(cellstypes_file)
            cellslocations = numpy.load(cellslocations_file)
            
            bloc = create_bloc_non_structure_from_numpy_array(
                xyz, connect, cellstypes, cellslocations)
        
        elif polys_file is not None:
            polys = numpy.load(polys_file)
            nb_polys = numpy.load(nb_polys_file)
            bloc = create_polydata_from_numpy_array(
                xyz, polys, nb_polys)
        else:
            raise IOError, "la connectivite n'est pas indiquee"
        
        return bloc
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne le resultat de la lecture
        
        """
        # lecture
        self.update()
        # output
        return self.output
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class LecteurFVUns(ObjetPyturbo):
    """lecture d'un fichier au format FieldView non-structure
    
    Le lecteur le plus lent, restreint, et pourri du module :)
    
    SURFACIQUE
    UNE SEULE VARIABLE, SUPPOSEE ETRE LA PRESSION
    
    """
    #_________________________________________________________________
    def __init__(self, acces_fichier):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_________________________________________________________________

    #_________________________________________________________________
    def update(self):
        print "Lecture {0}".format(self.acces_fichier)
        # Geometrie
        f = file(self.acces_fichier, 'r')
        data = f.read()
        xyz = data.split('Nodes')[1].split('Boundary Faces')[0]
        xyz = xyz[xyz.find('\n') + 1:]
        xyz = xyz[xyz.find('\n') + 1:]
        xyz = numpy.fromstring(xyz, sep = " ", dtype=numpy.float64).reshape(-1, 3)
        
        connect = data.split('Boundary Faces')[1].split('Elements')[0]
        connect = connect[connect.find('\n') + 1:]
        connect = connect[connect.find('\n') + 1:]
        connect = numpy.fromstring(connect, sep = " ", dtype=numpy.int)
        connect_propre = numpy.empty(0, dtype = int)
        k = 1
        liste_nb_element = numpy.empty(0, dtype = int)
        cellslocations = numpy.empty(0)
        while k < connect.size:
            cellslocations = numpy.r_[cellslocations, numpy.sum(liste_nb_element) + liste_nb_element.size]
            nb_element = connect[k]
            liste_nb_element = numpy.r_[liste_nb_element, nb_element]
            connect_propre = numpy.r_[connect_propre, connect[k], connect[k+1 : k + nb_element + 1] - 1]
            k += nb_element + 2
        connect = connect_propre
        
        # connect_propre = []
        # k = 1
        # liste_nb_element = []
        # cellslocations = []
        # while k < connect.size:
            # cellslocations += [numpy.sum(liste_nb_element) + len(liste_nb_element)]
            # nb_element = connect[k]
            # liste_nb_element += [nb_element]
            # connect_propre += [connect[k]] + (connect[k+1 : k + nb_element + 1] - 1).tolist()
            # k += nb_element + 2
        # connect = connect_propre
        # 
        # connect = numpy.asarray(connect)
        # liste_nb_element = numpy.asarray(connect)
        # cellslocations = numpy.asarray(cellslocations)
        
        cellstypes = numpy.zeros(liste_nb_element.size)
        dict_vtk_cell_type = {
            vtk.VTK_TRIANGLE: 3,
            vtk.VTK_QUAD: 4
            }
        for vtk_type in dict_vtk_cell_type:
            cellstypes[numpy.where(liste_nb_element == dict_vtk_cell_type[vtk_type])] = vtk_type
        
        self.vtkDataObject = self.create_geometrie(xyz, connect, 
                cellstypes, cellslocations)
        
        # chargement de l'array de pression
        array = data.split('Variables')[1]
        array = numpy.fromstring(array, sep = " ", dtype=numpy.float64)
        
        self.output = vtk_new_shallowcopy(self.vtkDataObject)
        self.output = ajouter_numpy_array_as_vtk_array(self.output, array, 'Pressure')
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def create_geometrie(self, xyz, connect, cellstypes, cellslocations):
        # Lecture du fichier geometrie et de l'arbre de connectivite
        bloc = create_bloc_non_structure_from_numpy_array(
            xyz, connect, cellstypes, cellslocations)
        return bloc
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne le resultat de la lecture
        
        """
        # lecture
        self.update()
        # output
        return self.output
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class LecteurPlot3D(ObjetPyturbo):
    """lecture d'un fichier au format Plot3D
    
    Le lecteur le plus lent, restreint, et pourri du module :)
    
    """
    #_________________________________________________________________
    def __init__(self, file_geom=None, file_data=None, geom_is_binary=1, res_is_binary=1):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_________________________________________________________________

    #_________________________________________________________________
    def update(self):
        print "Lecture {0}".format(self.acces_fichier)
        # Geometrie
        f = file(self.acces_fichier, 'r')
        data = f.read()
        xyz = data.split('Nodes')[1].split('Boundary Faces')[0]
        xyz = xyz[xyz.find('\n') + 1:]
        xyz = xyz[xyz.find('\n') + 1:]
        xyz = numpy.fromstring(xyz, sep = " ", dtype=numpy.float64).reshape(-1, 3)
        
        connect = data.split('Boundary Faces')[1].split('Elements')[0]
        connect = connect[connect.find('\n') + 1:]
        connect = connect[connect.find('\n') + 1:]
        connect = numpy.fromstring(connect, sep = " ", dtype=numpy.int)
        connect_propre = numpy.empty(0, dtype = int)
        k = 1
        liste_nb_element = numpy.empty(0, dtype = int)
        cellslocations = numpy.empty(0)
        while k < connect.size:
            cellslocations = numpy.r_[cellslocations, numpy.sum(liste_nb_element) + liste_nb_element.size]
            nb_element = connect[k]
            liste_nb_element = numpy.r_[liste_nb_element, nb_element]
            connect_propre = numpy.r_[connect_propre, connect[k], connect[k+1 : k + nb_element + 1] - 1]
            k += nb_element + 2
        connect = connect_propre
        
        # connect_propre = []
        # k = 1
        # liste_nb_element = []
        # cellslocations = []
        # while k < connect.size:
            # cellslocations += [numpy.sum(liste_nb_element) + len(liste_nb_element)]
            # nb_element = connect[k]
            # liste_nb_element += [nb_element]
            # connect_propre += [connect[k]] + (connect[k+1 : k + nb_element + 1] - 1).tolist()
            # k += nb_element + 2
        # connect = connect_propre
        # 
        # connect = numpy.asarray(connect)
        # liste_nb_element = numpy.asarray(connect)
        # cellslocations = numpy.asarray(cellslocations)
        
        cellstypes = numpy.zeros(liste_nb_element.size)
        dict_vtk_cell_type = {
            vtk.VTK_TRIANGLE: 3,
            vtk.VTK_QUAD: 4
            }
        for vtk_type in dict_vtk_cell_type:
            cellstypes[numpy.where(liste_nb_element == dict_vtk_cell_type[vtk_type])] = vtk_type
        
        self.vtkDataObject = self.create_geometrie(xyz, connect, 
                cellstypes, cellslocations)
        
        # chargement de l'array de pression
        array = data.split('Variables')[1]
        array = numpy.fromstring(array, sep = " ", dtype=numpy.float64)
        
        self.output = vtk_new_shallowcopy(self.vtkDataObject)
        self.output = ajouter_numpy_array_as_vtk_array(self.output, array, 'Pressure')
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def create_geometrie(self, xyz, connect, cellstypes, cellslocations):
        # Lecture du fichier geometrie et de l'arbre de connectivite
        bloc = create_bloc_non_structure_from_numpy_array(
            xyz, connect, cellstypes, cellslocations)
        return bloc
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne le resultat de la lecture
        
        """
        # lecture
        self.update()
        # output
        return self.output
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class LecteurFNC(ObjetPyturbo):
    """lecture d'un fichier au format FNC (powerflow)
    pour l'instant adapte seulement aux HEX
    """
    #_________________________________________________________________
    def __init__(self, file_path=None):
        attributs = locals().copy()
        del attributs['self']
        # initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
        # initialisation du dictionnaire des parametres
        self.parameters = None
        
    #_________________________________________________________________

    #_________________________________________________________________
    def lire_parametres(self):
        print "Ouverture {0}".format(self.file_path)
        f = netcdf.netcdf_file(self.file_path, 'r')
        
        self.parameters = {}
        # print "Offsets and scales"
        self.parameters['lx_offsets'] = f.variables['lx_offsets'][:]
        self.parameters['lx_scales'] = f.variables['lx_scales'][:]
        self.parameters['offset_coords'] = f.variables['csys'][0, 0:3, 3]

        # print "Scaling factors"
        self.parameters['coeff_dx'] = f.variables['lx_scales'][0]
        self.parameters['CFL_number'] = f.variables['lx_scales'][-1]
        self.parameters['coeff_dt'] = f.variables['lx_scales'][0] / f.variables['lx_scales'][3]
        self.parameters['dt'] = self.parameters['CFL_number'] * self.parameters['coeff_dt'] 
        self.parameters['coeff_press'] = f.variables['lx_scales'][1] * f.variables['lx_scales'][3]**2
        self.parameters['coeff_vel'] = f.variables['lx_scales'][3]
        self.parameters['coeff_density'] = f.variables['lx_scales'][1]
        self.parameters['offset_pressure'] = f.variables['lx_offsets'][4]

        # print "weight factor"
        self.parameters['t_lat'] = 1.0 / 3.0 # lattive temperature
        self.parameters['r_lat'] = 1.0 # lattice r_gas constant
        self.parameters['rTlat'] = self.parameters['r_lat'] * self.parameters['t_lat'] # weight for density to pressure
        self.parameters['irTlat'] = 1.0 / self.parameters['rTlat'] # weight for pressure to density

        # print "Rotation information"
        self.parameters['factor_rotation_angle'] = f.variables['ref_frame_indices'][:] #-1 (rotation sign) or 0
        self.parameters['theta0'] = f.variables['lrf_initial_angular_rotation'][:]
        self.parameters['tmp'] = f.variables['lrf_constant_angular_vel_mag'][:]
        self.parameters['omega'] = self.parameters['tmp'][0] / self.parameters['coeff_dt']
        self.parameters['lrf_n_polyline_vertices'] = f.variables['lrf_n_polyline_vertices'][:]
        self.parameters['lrf_polyline_vertices'] = f.variables['lrf_polyline_vertices'][:]
        self.parameters['lrf_axis_origin'] = f.variables['lrf_axis_origin'][:] # origine des coordonnees pour le LRF
        self.parameters['rotation_origin'] = self.parameters['lrf_axis_origin'][0, :] # copy and drop non necessary dimension
        self.parameters['lrf_axis_direction'] = f.variables['lrf_axis_direction'][:]  # axe de rotation
        self.parameters['rotation_axis'] = self.parameters['lrf_axis_direction'][0, :] # copy and drop non necessary dimension
        self.parameters['lrf_n_points'] = f.variables['lrf_n_points'][:] # nb de points appartenant au LRF
        self.parameters['ref_frame_indices'] = f.variables['ref_frame_indices'][:] # -1 stationnary 0 LRF
        
        # Time info
        self.parameters['start_time'] = f.variables['start_time'][:]
        self.parameters['end_time'] = f.variables['end_time'][:]
        self.parameters['current_time'] = 0.5 * (f.variables['start_time'][:] + f.variables['end_time'][:]) * self.parameters['dt']
        
        # Variables infos
        self.parameters['var_names'] = numpy.array(f.variables['variable_tiny_names'][:].tostring().split('\x00')[:-1])
                
        # print "Initial phase of the mesh"
        self.parameters['n_rot0'] = f.variables['lrf_initial_n_revolutions'][0] * numpy.pi * 2.0
        self.parameters['alpha0'] = f.variables['lrf_initial_angular_rotation'][0] + self.parameters['n_rot0']
        
        self.parameters['axis_index'] = numpy.where(self.parameters['rotation_axis'])[0][0]
        self.parameters['sign_rotation'] = self.parameters['rotation_axis'][self.parameters['axis_index']]
        self.parameters['init_angle'] = self.parameters['alpha0'] * self.parameters['sign_rotation']
        
        # fermeture du fichier
        f.close()
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_infos_temporelles(self):
        """ affiche les informations temporelles du fichier fnc
        les parametres du fichier FNC doivent avoir ete lus avant en utilisant la fonction self.read_parameters()
        """
        print 'Pas de temps : ', self.parameters['dt']
        print 'Indices temporels min : ', self.parameters['start_time']
        print 'Indices temporels max : ', self.parameters['end_time']
        print 'Temps moyens [s] : ', self.parameters['current_time']
        print 'Nombre de pas de temps : ', self.parameters['current_time'].size
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_infos_variables(self):
        """Retourne les informations sur les variables disponibles
        """
        print 'Variables disponibles : ', self.parameters['var_names']
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def lire_maillage(self):
        """ Lecture de la geometrie. self.maillage est un multibloc. 
        Block 0 = rotor
        Block 1 = stator
        """
        print "Ouverture {0}".format(self.file_path)
        f = netcdf.netcdf_file(self.file_path, 'r')
        
        print 'Chargement de la geometrie'
        # read coordinates
        element_coords = f.variables['coords']
        voxel_scales = f.variables['voxel_scales']
        nelm = element_coords.shape[0]
        ndims = element_coords.shape[1]
        
        # "Spatial offset"
        # element_coords[:, :] =  element_coords[:, :] + self.parameters['offset_coords'][None, :]
        element_coords =  element_coords[:, :] + f.variables['csys'][0, 0:3, 3]
        
        # voxel size
        dx = 2**(1.0 + voxel_scales[:])
        
        # compute voxels vertices coordinates
        vertices_coords = numpy.zeros((nelm, 8, ndims))
        vertices_coords[:, 0, :] = element_coords[:, :]
        vertices_coords[:, 1, :] = element_coords[:, :] + dx[:, None] * numpy.array([1, 0, 0])[None, :]
        vertices_coords[:, 2, :] = element_coords[:, :] + dx[:, None] * numpy.array([1, 1, 0])[None, :]
        vertices_coords[:, 3, :] = element_coords[:, :] + dx[:, None] * numpy.array([0, 1, 0])[None, :]
        vertices_coords[:, 4, :] = element_coords[:, :] + dx[:, None] * numpy.array([0, 0, 1])[None, :]
        vertices_coords[:, 5, :] = element_coords[:, :] + dx[:, None] * numpy.array([1, 0, 1])[None, :]
        vertices_coords[:, 6, :] = element_coords[:, :] + dx[:, None] * numpy.array([1, 1, 1])[None, :]
        vertices_coords[:, 7, :] = element_coords[:, :] + dx[:, None] * numpy.array([0, 1, 1])[None, :]
        
        ###################
        # CREATE_ROTOR
        # ind_elm_rotor = numpy.where(self.parameters['ref_frame_indices'] == 0)[0]
        ind_elm_rotor = numpy.where(f.variables['ref_frame_indices'][:] == 0)[0]
        vertices_coords_rotor = vertices_coords[ind_elm_rotor].reshape((-1, 3))
        
        # remove duplicated vertices and get associated tree connectivity
        coords_view_rotor = numpy.ascontiguousarray(vertices_coords_rotor).view(
            numpy.dtype((numpy.void, vertices_coords_rotor.dtype.itemsize * 3)))
        _, idx_rotor, inv_rotor = numpy.unique(coords_view_rotor, return_index=True, return_inverse=True)
        
        # unique_vertices_coords = vertices_coords[idx] * self.parameters['coeff_dx']
        unique_vertices_coords_rotor = vertices_coords_rotor[idx_rotor] * f.variables['lx_scales'][0]
        connectivity_rotor = numpy.arange(ind_elm_rotor.size * 8).reshape((ind_elm_rotor.size, 8))
        connectivity_rotor = inv_rotor[connectivity_rotor]
        connectivity_rotor = numpy.concatenate(
            ((8 * numpy.ones(connectivity_rotor.shape[0], dtype=numpy.int))[:, None], connectivity_rotor), 
            axis = -1)
        
        # VTK rotor
        vtk_rotor = create_bloc_non_structure_from_numpy_array(
            coords = unique_vertices_coords_rotor, 
            cells = connectivity_rotor.ravel(), 
            cellstypes = vtk.vtkHexahedron().GetCellType() * numpy.ones(connectivity_rotor.shape[0]), 
            cellslocations = numpy.arange(connectivity_rotor.shape[0]) * 9)
        
        ###################
        # STATOR
        # ind_elm_stator = numpy.where(self.parameters['ref_frame_indices'] == 0)[0]
        ind_elm_stator = numpy.where(f.variables['ref_frame_indices'][:] == -1)[0]
        vertices_coords_stator = vertices_coords[ind_elm_stator].reshape((-1, 3))
        
        # remove duplicated vertices and get associated tree connectivity
        coords_view_stator = numpy.ascontiguousarray(vertices_coords_stator).view(
            numpy.dtype((numpy.void, vertices_coords_stator.dtype.itemsize * 3)))
        _, idx_stator, inv_stator = numpy.unique(coords_view_stator, return_index=True, return_inverse=True)
        
        # unique_vertices_coords = vertices_coords[idx] * self.parameters['coeff_dx']
        unique_vertices_coords_stator = vertices_coords_stator[idx_stator] * f.variables['lx_scales'][0]
        connectivity_stator = numpy.arange(ind_elm_stator.size * 8).reshape((ind_elm_stator.size, 8))
        connectivity_stator = inv_stator[connectivity_stator]
        connectivity_stator = numpy.concatenate(
            ((8 * numpy.ones(connectivity_stator.shape[0], dtype=numpy.int))[:, None], connectivity_stator), 
            axis = -1)
        
        # VTK stator
        vtk_stator = create_bloc_non_structure_from_numpy_array(
            coords = unique_vertices_coords_stator, 
            cells = connectivity_stator.ravel(), 
            cellstypes = vtk.vtkHexahedron().GetCellType() * numpy.ones(connectivity_stator.shape[0]), 
            cellslocations = numpy.arange(connectivity_stator.shape[0]) * 9)
        
        ####################
        # ASSEMBLAGE MULTIBLOCK ET STOCKAGE
        maillage = vtk.vtkMultiBlockDataSet()
        maillage.SetBlock(0, vtk_rotor)
        maillage.SetBlock(1, vtk_stator)
        
        self.maillage = maillage
        
        # fermeture du fichier
        f.close()
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_maillage(self):
        """ retoure le maillage sans donnees
        """
        return self.maillage
    #_________________________________________________________________
    
    #_________________________________________________________________
    def importer_maillage(self, vtkMultiBlockDataSet):
        """ permet d'importer un maillage. 
        vtkMultiBlockDataSet doit avoir deux blocs : 1-rotor, 2-stator. 
        """
        print 'IMPORT du maillage'
        self.maillage = vtkMultiBlockDataSet
        return 0
    #_________________________________________________________________
    
    
    #_________________________________________________________________
    def lire_datas(self, ind_temps, noms_vars=None, rotation_rotor=True, axe=2, numbloc=None):
        """ Lecture des datas au temps indique par ind_time
            - ind_temps est donne l'index du temps a lire
            - noms_vars est une liste des variables a lire. Si elle est None, toutes les variables disponibles sont lues
        """
        # initialisation de output
        self.output = vtk_new_shallowcopy(self.maillage)
        if numbloc is not None:
            self.output.SetBlock(0 if numbloc is 1 else 1, None)
        
        print "Ouverture {0}".format(self.file_path)
        f = netcdf.netcdf_file(self.file_path, 'r')
        print "reading variable data"
        
        # indice des voxels du rotor et du stator
        ind_elm_rotor = numpy.where(f.variables['ref_frame_indices'][:] == 0)[0]
        ind_elm_stator = numpy.where(f.variables['ref_frame_indices'][:] == -1)[0]
        
        for nom_var in self.parameters['var_names'] if noms_vars is None else noms_vars:
            if nom_var in self.parameters['var_names']:
                print "Chargement et redimensionnement de {0} a l'index temporel {1}".format(nom_var, ind_temps)
                ind_var = int(numpy.where(self.parameters['var_names'] == nom_var)[0])
                data_var = numpy.array(f.variables['measurements'][ind_temps, ind_var]).ravel()
                
                # redimensionnement
                if nom_var is 'p':
                    data_var = (data_var + self.parameters['offset_pressure']) * self.parameters['coeff_press']
                if nom_var[0] == 'v':
                    data_var = data_var * self.parameters['coeff_vel']
                    
                if (numbloc is None or numbloc == 0) and self.output.GetBlock(0) is not None:
                    self.output.SetBlock(0, 
                        ajouter_numpy_array_as_vtk_array(self.output.GetBlock(0), data_var[ind_elm_rotor], self.parameters['var_names'][ind_var])
                        )
                if (numbloc is None or numbloc == 1) and self.output.GetBlock(1) is not None:
                    self.output.SetBlock(1, 
                        ajouter_numpy_array_as_vtk_array(self.output.GetBlock(1), data_var[ind_elm_stator], self.parameters['var_names'][ind_var])
                        )
        
        # Si nom_vars est none, on derive la pression ou la densite l'un a partir de celui qui est disponible
        if noms_vars is None or 'p' in noms_vars or 'rho' in noms_vars:
            if 'p' in self.parameters['var_names']:
                print 'Calcul de rho a partir de p'
                ind_var_p = int(numpy.where(self.parameters['var_names'] == 'p')[0])
                data_p = numpy.array(f.variables['measurements'][ind_temps, ind_var_p]).ravel() 
                data_rho = data_p * self.parameters['irTlat'] * self.parameters['coeff_density']
                self.output.SetBlock(0, 
                    ajouter_numpy_array_as_vtk_array(self.output.GetBlock(0), data_rho[ind_elm_rotor], 'rho')
                    )
                self.output.SetBlock(1, 
                    ajouter_numpy_array_as_vtk_array(self.output.GetBlock(1), data_rho[ind_elm_stator], 'rho')
                    )
            # amarsan - ATTENTION -- je suis pas certain que ca s'appelle rho...  
            if 'rho' in self.parameters['var_names']:
                print 'Calcul de p a partir de rho'
                ind_var_rho = int(numpy.where(self.parameters['var_names'] == 'rho')[0])
                data_rho = numpy.array(f.variables['measurements'][ind_temps, ind_var_rho]).ravel() 
                data_p = (data_rho * self.parameters['rTlat'] + self.parameters['offset_pressure']) * self.parameters['coeff_press']
                if (numbloc is None or numbloc == 0) and self.output.GetBlock(0) is not None:
                    self.output.SetBlock(0, 
                        ajouter_numpy_array_as_vtk_array(self.output.GetBlock(0), data_p[ind_elm_rotor], 'p')
                        )
                if (numbloc is None or numbloc == 1) and self.output.GetBlock(1) is not None:
                    self.output.SetBlock(1, 
                        ajouter_numpy_array_as_vtk_array(self.output.GetBlock(1), data_p[ind_elm_stator], 'p')
                        )
        # rotation du rotor pour le mettre a la bonne position
        if (numbloc is None or numbloc == 0) and self.output.GetBlock(0) is not None and rotation_rotor is True:
            print "Rotation du rotor autour de l'axe {1}, omega = {0} rad/s".format(self.parameters['omega'], ['x', 'y', 'z'][axe])
            alpha = self.parameters['omega'] * self.parameters['sign_rotation'] * self.parameters['current_time'][ind_temps] \
                + self.parameters['init_angle']
            self.output.SetBlock(0, rotation(self.output.GetBlock(0), numpy.rad2deg(alpha), axe))
        
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne le resultat de la lecture
        
        """
        return self.output
    #_________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class LecteurSNC(ObjetPyturbo):
    """lecture d'un fichier au format FNC (powerflow)
    pour l'instant adapte seulement aux HEX
    """
    #_________________________________________________________________
    def __init__(self, file_path=None):
        attributs = locals().copy()
        del attributs['self']
        # initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
        # initialisation du dictionnaire des parametres
        self.parameters = None
        
    #_________________________________________________________________

    #_________________________________________________________________
    def lire_parametres(self):
        print "Ouverture {0}".format(self.file_path)
        f = netcdf.netcdf_file(self.file_path, 'r')
        
        self.parameters = {}
        # print "Offsets and scales"
        self.parameters['lx_offsets'] = f.variables['lx_offsets'][:]
        self.parameters['lx_scales'] = f.variables['lx_scales'][:]
        self.parameters['offset_coords'] = f.variables['csys'][0, 0:3, 3]

        # print "Scaling factors"
        self.parameters['coeff_dx'] = f.variables['lx_scales'][0]
        self.parameters['CFL_number'] = f.variables['lx_scales'][-1]
        self.parameters['coeff_dt'] = f.variables['lx_scales'][0] / f.variables['lx_scales'][3]
        self.parameters['dt'] = self.parameters['CFL_number'] * self.parameters['coeff_dt'] 
        self.parameters['coeff_press'] = f.variables['lx_scales'][1] * f.variables['lx_scales'][3]**2
        self.parameters['coeff_vel'] = f.variables['lx_scales'][3]
        self.parameters['coeff_density'] = f.variables['lx_scales'][1]
        self.parameters['offset_pressure'] = f.variables['lx_offsets'][4]

        # print "weight factor"
        self.parameters['t_lat'] = 1.0 / 3.0 # lattive temperature
        self.parameters['r_lat'] = 1.0 # lattice r_gas constant
        self.parameters['rTlat'] = self.parameters['r_lat'] * self.parameters['t_lat'] # weight for density to pressure
        self.parameters['irTlat'] = 1.0 / self.parameters['rTlat'] # weight for pressure to density

        # print "Rotation information"
        if 'ref_frame_indices' in f.variables:
            self.parameters['factor_rotation_angle'] = f.variables['ref_frame_indices'][:] #-1 (rotation sign) or 0
            self.parameters['theta0'] = f.variables['lrf_initial_angular_rotation'][:]
            self.parameters['tmp'] = f.variables['lrf_constant_angular_vel_mag'][:]
            self.parameters['omega'] = self.parameters['tmp'][0] / self.parameters['coeff_dt']
            self.parameters['lrf_n_polyline_vertices'] = f.variables['lrf_n_polyline_vertices'][:]
            self.parameters['lrf_polyline_vertices'] = f.variables['lrf_polyline_vertices'][:]
            self.parameters['lrf_axis_origin'] = f.variables['lrf_axis_origin'][:] # origine des coordonnees pour le LRF
            self.parameters['rotation_origin'] = self.parameters['lrf_axis_origin'][0, :] # copy and drop non necessary dimension
            self.parameters['lrf_axis_direction'] = f.variables['lrf_axis_direction'][:]  # axe de rotation
            self.parameters['rotation_axis'] = self.parameters['lrf_axis_direction'][0, :] # copy and drop non necessary dimension
            self.parameters['lrf_n_points'] = f.variables['lrf_n_points'][:] # nb de points appartenant au LRF
            self.parameters['ref_frame_indices'] = f.variables['ref_frame_indices'][:] # -1 stationnary 0 LRF
            
            # print "Initial phase of the mesh"
            self.parameters['n_rot0'] = f.variables['lrf_initial_n_revolutions'][0] * numpy.pi * 2.0
            self.parameters['alpha0'] = f.variables['lrf_initial_angular_rotation'][0] + self.parameters['n_rot0']
            
            self.parameters['axis_index'] = numpy.where(self.parameters['rotation_axis'])[0][0]
            self.parameters['sign_rotation'] = self.parameters['rotation_axis'][self.parameters['axis_index']]
            self.parameters['init_angle'] = self.parameters['alpha0'] * self.parameters['sign_rotation']
            self.parameters['rotor-stator'] = True
        else:
            self.parameters['rotor-stator'] = False
        
        # Time info
        self.parameters['start_time'] = f.variables['start_time'][:]
        self.parameters['end_time'] = f.variables['end_time'][:]
        self.parameters['current_time'] = 0.5 * (f.variables['start_time'][:] + f.variables['end_time'][:]) * self.parameters['dt']
        
        # Variables infos
        self.parameters['var_names'] = numpy.array(f.variables['variable_tiny_names'][:].tostring().split('\x00')[:-1])
        
        
        # fermeture du fichier
        f.close()
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_infos_temporelles(self):
        """ affiche les informations temporelles du fichier fnc
        les parametres du fichier FNC doivent avoir ete lus avant en utilisant la fonction self.read_parameters()
        """
        print 'Pas de temps : ', self.parameters['dt']
        print 'Indices temporels min : ', self.parameters['start_time']
        print 'Indices temporels max : ', self.parameters['end_time']
        print 'Temps moyens [s] : ', self.parameters['current_time']
        print 'Nombre de pas de temps : ', self.parameters['current_time'].size
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_infos_variables(self):
        """Retourne les informations sur les variables disponibles
        """
        print 'Variables disponibles : ', self.parameters['var_names']
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def lire_maillage(self, numbloc=None):
        """ Lecture de la geometrie de la surface
        Si numbloc = 0 ou numbloc = 1 alors seulement les blocs 0 et 1 sont retournes. 
        bloc 0 : rotor
        bloc 1 : stator
        """
        print "Ouverture {0}".format(self.file_path)
        f = netcdf.netcdf_file(self.file_path, 'r')
        
        print 'Chargement de la geometrie'
        # read coordinates
        coords = f.variables['vertex_coords']
        vertex_list = f.variables['vertex_refs']
        first_vertex = f.variables['first_vertex_refs']
        surfel_area = f.variables['surfel_area']
        surfel_norm = f.variables['surfel_normal']
        face_id = f.variables['face']
        
        # "Spatial offset"
        coords =  (coords[:, :] + f.variables['csys'][0, 0:3, 3]) * f.variables['lx_scales'][0]
        
        
        # Create surface
        nb_polys = f.variables['first_vertex_refs'][:].size
        polys = numpy.insert(f.variables['vertex_refs'][:], 
            f.variables['first_vertex_refs'][:], 
            numpy.r_[f.variables['first_vertex_refs'][1:] - f.variables['first_vertex_refs'][:-1], f.variables['vertex_refs'][:].size - f.variables['first_vertex_refs'][-1]]
            )
        surface = create_polydata_from_numpy_array(coords, 
            polys = polys, 
            nb_polys = nb_polys)
        
        # decoupage dans le cas d'un cas rotor-stator
        if self.parameters['rotor-stator'] is True:
            # amarsan 30 aout 2017 - modif pour surface rotor-stator
            # indice des voxels du rotor et du stator
            ind_elm_rotor = numpy.where(f.variables['ref_frame_indices'][:] == 0)[0]
            ind_elm_stator = numpy.where(f.variables['ref_frame_indices'][:] == -1)[0]
            
            # multibloc. bloc 0 sera rotor. bloc 1 sera stator. 
            maillage = vtk.vtkMultiBlockDataSet()
            # ROTOR
            if numbloc is None or numbloc == 0:
                # amarsan - partie a remplacer quand vtkIdList.SetArray sera disponible pour pouvoir l'indiquer d'un coup avec un numpy array
                vtkIdList = vtk.vtkIdList()
                vtkIdList.SetNumberOfIds(ind_elm_rotor.size)
                for k in xrange(ind_elm_rotor.size):
                    vtkIdList.SetId(k, ind_elm_rotor[k])
                extractor = vtkFiltersExtraction.vtkExtractCells()
                extractor.SetInputData(surface)
                extractor.SetCellList(vtkIdList)
                extractor.Update()
                surface_rotor = extractor.GetOutput()
                surface_rotor.GetCellData().RemoveArray('vtkOriginalCellIds')
                
                p_surface_rotor = vtk.vtkPolyData()
                p_surface_rotor.SetPolys(surface_rotor.GetCells())
                p_surface_rotor.SetPoints(surface_rotor.GetPoints())
                maillage.SetBlock(0, p_surface_rotor)
            # STATOR
            if numbloc is None or numbloc == 1:
                # amarsan - partie a remplacer quand vtkIdList.SetArray sera disponible pour pouvoir l'indiquer d'un coup avec un numpy array
                vtkIdList = vtk.vtkIdList()
                vtkIdList.SetNumberOfIds(ind_elm_stator.size)
                for k in xrange(ind_elm_stator.size):
                    vtkIdList.SetId(k, ind_elm_stator[k])
                extractor = vtkFiltersExtraction.vtkExtractCells()
                extractor.SetInputData(surface)
                extractor.SetCellList(vtkIdList)
                extractor.Update()
                surface_stator = extractor.GetOutput()
                surface_stator.GetCellData().RemoveArray('vtkOriginalCellIds')
                
                p_surface_stator = vtk.vtkPolyData()
                p_surface_stator.SetPolys(surface_stator.GetCells())
                p_surface_stator.SetPoints(surface_stator.GetPoints())
                maillage.SetBlock(1, p_surface_stator)
            else:
                maillage.SetBlock(1, None)
            
            # stockage dans l'objet de sortie
            self.maillage = maillage
        
        # cas stator seul, avec uniquement un seul referentiel. L'objet de sortie est un polydata simple. 
        else:
            self.maillage = surface
        
        # fermeture du fichier
        f.close()
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_decoupage_surface(self):
        """ Fonction qui retourne les index des sous-ensemble des surfaces
        et les noms associes
        """
        f = netcdf.netcdf_file(self.file_path, 'r')
        index = numpy.unique(f.variables['face'][:])
        faces_names = numpy.take(f.variables['face_names'][:].tostring().split('\x00'), index)
        print numpy.c_[index, faces_names]
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_maillage(self):
        """ retoure le maillage sans donnees
        """
        return self.maillage
    #_________________________________________________________________
    
    #_________________________________________________________________
    def importer_maillage(self, vtkMultiBlockDataSet):
        """ permet d'importer un maillage. 
        vtkMultiBlockDataSet
        """
        print 'IMPORT du maillage'
        self.maillage = vtkMultiBlockDataSet
        return 0
    #_________________________________________________________________
    
    
    #_________________________________________________________________
    def lire_datas(self, ind_temps, noms_vars=None, numbloc=None, rotation_rotor=True, axe=2, ajouter_faces_id=True, 
            ajouter_normales=True, ajouter_surfaces=False):
        """ Lecture des datas au temps indique par ind_time
                - ind_temps est donne l'index du temps a lire
                - noms_vars est une liste des variables a lire. Si elle est None, toutes les variables disponibles sont lues
            
            Si numbloc = 0 ou numbloc = 1 alors seulement les blocs 0 et 1 sont stockes. 
            bloc 0 : rotor
            bloc 1 : stator
        """
        # initialisation de output
        output = vtk_new_shallowcopy(self.maillage)
        if numbloc is not None:
            output.SetBlock(0 if numbloc is 1 else 1, None)
        
        print "Ouverture {0}".format(self.file_path)
        f = netcdf.netcdf_file(self.file_path, 'r')
        print "reading variable data"
        
        # amarsan 30 aout 2017 - modif pour surface rotor-stator
        # indice des voxels du rotor et du stator
        if self.parameters['rotor-stator'] is True:
            ind_elm_rotor = numpy.where(f.variables['ref_frame_indices'][:] == 0)[0]
            ind_elm_stator = numpy.where(f.variables['ref_frame_indices'][:] == -1)[0]
        
        #####################
        # fonction d'ajout d'une variable en faisant la distinction des cas rotor-stator et stator seul
        def _ajouter_variable(output, data_var, var_name):
            if self.parameters['rotor-stator'] is True:
                if (numbloc is None or numbloc == 0) and output.GetBlock(0) is not None:
                    output.SetBlock(0, 
                        ajouter_numpy_array_as_vtk_array(output.GetBlock(0), data_var[ind_elm_rotor], var_name)
                        )
                if (numbloc is None or numbloc == 1) and output.GetBlock(1) is not None:
                    output.SetBlock(1, 
                        ajouter_numpy_array_as_vtk_array(output.GetBlock(1), data_var[ind_elm_stator], var_name)
                        )
            else:
                output = ajouter_numpy_array_as_vtk_array(output, data_var, var_name)
            return output
        #####################
        
        for nom_var in self.parameters['var_names'] if noms_vars is None else noms_vars:
            if nom_var in self.parameters['var_names']:
                print "Chargement et redimensionnement de {0} a l'index temporel {1}".format(nom_var, ind_temps)
                ind_var = int(numpy.where(self.parameters['var_names'] == nom_var)[0])
                data_var = numpy.array(f.variables['measurements'][ind_temps, ind_var]).ravel()
                
                # redimensionnement
                if nom_var is 'p':
                    data_var = (data_var + self.parameters['offset_pressure']) * self.parameters['coeff_press']
                if nom_var[0] == 'v':
                    data_var = data_var * self.parameters['coeff_vel']
                
                output = _ajouter_variable(output, data_var, self.parameters['var_names'][ind_var])
        
        # On derive la pression ou la densite l'un a partir de celui qui est disponible
        if noms_vars is None or 'p' in noms_vars or 'rho' in noms_vars:
            if 'p' in self.parameters['var_names']:
                print 'Calcul de rho a partir de p'
                ind_var_p = int(numpy.where(self.parameters['var_names'] == 'p')[0])
                data_p = numpy.asarray(f.variables['measurements'][ind_temps, ind_var_p]).ravel() 
                data_var = data_p * self.parameters['irTlat'] * self.parameters['coeff_density']
                output = _ajouter_variable(output, data_var, 'rho')
            
            if 'rho' in self.parameters['var_names']:
                print 'Calcul de p a partir de rho'
                ind_var_rho = int(numpy.where(self.parameters['var_names'] == 'rho')[0])
                data_rho = numpy.asarray(f.variables['measurements'][ind_temps, ind_var_rho]).ravel() 
                data_var = (data_rho * self.parameters['rTlat'] + self.parameters['offset_pressure']) * self.parameters['coeff_press']
                output = _ajouter_variable(output, data_var, 'p')
        
        # on ajoute les face id si demande
        if ajouter_faces_id is True:
            print 'Ajout des faces id'
            data_var = f.variables['face'][:]
            output = _ajouter_variable(output, data_var, 'face_id')
        
        # on ajoute les normales si demande
        if ajouter_normales is True:
            print 'Ajout des normales'
            data_var = f.variables['surfel_normal'][:]
            output = _ajouter_variable(output, data_var, 'surfel_normals')
        
        # on ajoute les surfaces si demande
        if ajouter_surfaces is True:
            print 'Ajout des surfaces'
            data_var = f.variables['surfel_area'][:]
            output = _ajouter_variable(output, data_var, 'surfel_area')
        
        # rotation du rotor pour correspondre a l'instant demande
        if self.parameters['rotor-stator'] is True:
            if (numbloc is None or numbloc == 0) and output.GetBlock(0) is not None and rotation_rotor is True:
                print "Rotation du rotor autour de l'axe z, omega = {0} rad/s".format(self.parameters['omega'])
                alpha = self.parameters['omega'] * self.parameters['current_time'][ind_temps] + self.parameters['init_angle']
                output.SetBlock(0, rotation(output.GetBlock(0), numpy.rad2deg(alpha), axe))
        
        self.output = output
        return 0
    #_________________________________________________________________
    
    #_________________________________________________________________
    def get_output(self):
        """retourne le resultat de la lecture
        
        """
        return self.output
    #_________________________________________________________________
#__________________________________________________________________________________________
