try: 
    from paraview.vtk import vtkFiltersVerdict
    from paraview.vtk import vtkFiltersGeneral
    from paraview.vtk import vtkCommonTransforms
    from paraview.vtk import vtkFiltersGeometry
    from paraview.vtk import vtkFiltersExtraction
except:
    import vtk as vtkFiltersVerdict
    import vtk as vtkFiltersExtraction
    import vtk as vtkFiltersGeneral
    import vtk as vtkCommonTransforms
    import vtk as vtkFiltersGeometry

try :from paraview import numpy_support
except: from vtk.util import numpy_support

import numpy

from UVParametrizationFilter import UVParametrization as UVParametrisation
from objets import ObjetPyturbo
from calculs import CalculettePyturbo
from fonctions_basiques import *

#__________________________________________________________________________________________
class Extraction(ObjetPyturbo):
    """permet d'extraire une surface quelconque
    
    les surface possibles sont
        - i= ;j= ; k=    si l'entree est compose de vtkStructuredGrid (mono ou multiblock)
        - toute grandeur calculable par une CalculettePyturbo
            utiliser coordx, coordy et coordz
        - toute autre grandeur calculable par calculs.CalculettePyturbo
    
    indiquez la formule dans formule_extraction SANS ESPACES
        coordx+coordy=12. par exemple
    
    imin, imax etc. sont utilisables
    
    #ToDo
    completer la fonction pour pouvoir prendre une inegalite
        
    """
    #_____________________________________________________________________________________
    def __init__(self, input=None, formule_extraction=None, calculer_vecteur_normal=True,
            normals_aux_cellules=False, axe=None):
        #initialisation de la classe parente
        attributs = locals().copy()
        del attributs['self']
        ObjetPyturbo.__init__(self, **attributs)
        # initialisation particuliere
        self._mettre_a_jour = True
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    
    def set(self, nom_attribut, valeur):
        """fonction set specifique
        
        gere la variable locale _changement
        qui sert lorsque l'on appelle la sortie
        a savoir s'il faut regenerer la coupe
        """
        setattr(self, nom_attribut, valeur)
        if nom_attribut != '_mettre_a_jour':
            self._mettre_a_jour = True
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def __couper_bloc__(self, vtkDataSet):
        """retourne un filtre vtk adapte a la coupe desiree
        il suffit ensuite de faire GetOutput() pour obtenir le resultat de la coupe
        ne s'applique PAS a un multiblockdataset
        
        """
        # VERIFICATIONS INITIALES
        if isinstance(vtkDataSet, vtk.vtkMultiBlockDataSet):
            raise IOError, '__couper_bloc__ ne prend PAS de MultiBlockDataSet en entree'
        if self.formule_extraction is None:
            raise IOError, "indiquez d'abord la self.formule_extraction pour l'extraction"
        if ' ' in self.formule_extraction:
            raise IOError, "la formule_extraction doit etre indiquee sans espaces"
        if not '=' in self.formule_extraction:
            raise IOError, "pour l'instant, seules les equations sont supportees comme formule d'extraction"
        else:
            cle_coupe = self.formule_extraction.split('=')[0].strip()
            valeur_coupe = self.formule_extraction.split('=')[1].strip()
        
        # EXECUTION
        if cle_coupe in ['i', 'j', 'k']:
            if not isinstance(vtkDataSet, vtk.vtkStructuredGrid):
                raise IOError, "une coupe i, j ou k est demande, mais l'entree n'est pas un bloc structure"
            filtre_vtk = vtkFiltersExtraction.vtkExtractGrid()
            vtk_set_input(filtre_vtk, vtkDataSet)
            extent_vtkDataSet = list(vtkDataSet.GetExtent())
            exec "valeur_coupe = {0}".format(valeur_coupe.replace(
                'imax', str(extent_vtkDataSet[1])).replace(
                'jmax', str(extent_vtkDataSet[3])).replace(
                'kmax', str(extent_vtkDataSet[5])).replace(
                'imin', str(extent_vtkDataSet[0])).replace(
                'jmin', str(extent_vtkDataSet[2])).replace(
                'kmin', str(extent_vtkDataSet[4])))
            voi = extent_vtkDataSet[0] if cle_coupe != 'i' else valeur_coupe, \
                extent_vtkDataSet[1] if cle_coupe != 'i' else valeur_coupe, \
                extent_vtkDataSet[2] if cle_coupe != 'j' else valeur_coupe, \
                extent_vtkDataSet[3] if cle_coupe != 'j' else valeur_coupe, \
                extent_vtkDataSet[4] if cle_coupe != 'k' else valeur_coupe, \
                extent_vtkDataSet[5] if cle_coupe != 'k' else valeur_coupe
            filtre_vtk.SetVOI(voi)
            filtre_vtk.Update()
            data = filtre_vtk.GetOutput()
        else:
            exec "valeur_coupe = float({0})".format(valeur_coupe)
            calculette = CalculettePyturbo(input = vtkDataSet, axe = self.axe) if self.axe is not None else CalculettePyturbo(input = vtkDataSet)
            calculette.set('a_calculer', cle_coupe)
            a_couper = calculette.get_output()
            a_couper = set_scalaires_actifs(
                input = a_couper, loc = 'points', array_name = cle_coupe)
            filtre_vtk = vtk.vtkContourFilter()
            filtre_vtk.SetComputeNormals(0)
            vtk_set_input(filtre_vtk, a_couper)
            filtre_vtk.SetValue(0, valeur_coupe)
            filtre_vtk.Update()
            data = filtre_vtk.GetOutput()
        
        #CALCUL DU VECTEUR NORMAL
        if self.calculer_vecteur_normal == 1:
            data = calculer_vecteur_normal(data, self.normals_aux_cellules)
        return data
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def update(self):
        """genere la coupe
        
        """
        if self.input is None:
            raise IOError, "indiquez d'abord l'objet vtk en entree"
        
        if isinstance(self.input, vtk.vtkMultiBlockDataSet):
            self.output = vtk_new_instance(self.input)
            for numbloc in get_numeros_blocs_non_vides(self.input):
                extraction_bloc = self.__couper_bloc__(self.input.GetBlock(numbloc))
                if extraction_bloc.GetNumberOfPoints() != 0:
                    self.output.SetBlock(numbloc, extraction_bloc)
        else:
            self.output = self.__couper_bloc__(self.input)
        
        self._mettre_a_jour = False
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_output(self):
        """retourne la sortie de la classe
        
        mise a jour effectuee si necessaire
        
        """
        if self._mettre_a_jour:
            self.update()
        return self.output
    #_____________________________________________________________________________________
#_____________________________________________________________________________________
