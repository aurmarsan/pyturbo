try: from paraview import vtk 
except: import vtk
import numpy

#__________________________________________________________________________________________
class ObjetPyturbo:
    """classe generale
    
    definit des methodes communes a tous les objets manipules

    """
    #_____________________________________________________________________________________
    def __init__(self, **arguments_entree):
        for propriete in arguments_entree:
            if propriete != 'self':
                self.set(propriete, arguments_entree[propriete])
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def set(self, attribut, valeur):
        """ methode d'attribution de la valeur d'un attribut
        
                arguments d'entree:
                    * attribut: nom de l'attribut (caracteres)
                    * valeur:   valeur de l'attribut (type == type(attribut))
                argument de sortie:
                    * 1 si attribut existe
                    * 0 sinon
                    
        """
        setattr(self, attribut, valeur)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get(self, attribut):
        """ methode de restitution d'un attribut
        
                argument d'entree:
                    * attribut: nom de l'attribut (caracteres)
                argument de sortie:
                    * valeur de l'attribut
                    
        """
        return getattr(self, attribut)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def __str__(self):
        """definit ce qui sera affiche par print
        La fonction python callable est utilisee :
            | Return whether the object is callable (i.e., some kind of function).
            | Note that classes are callable, as are instances with a __call__() method.
        
        Les attributs doivent etre formates de la facon suivante :
            - tous les attributs commencent par une minuscule, fonctions comprises
            - seules les classes commencent par une majuscule, dont les classes qui incluses
                en tant qu'attributs d'une autre classe plus grosse. 
        
        Pour les arrays numpy, on affiche simplement shape. 
        Compte-tenu du formatage precedent, ceux-ci sont designes par un nom qui commence par une minuscule. 
        """
        chaine = "Objet {0}".format(self.get('__class__'))
        for arg in dir(self):
            if not callable(self.get(arg)) and arg[0].islower() \
                    and not isinstance(self.get(arg), vtk.vtkObject):
                if "numpy" in str(type(self.get(arg))): 
                    chaine += '\n\t- {0} : shape = {1}'.format(arg, self.get(arg).shape)
                else:
                    chaine += '\n\t- {0} : {1}'.format(arg, self.get(arg))
        for arg in dir(self):
            if not callable(self.get(arg)) and arg[0].isupper() \
                    and not isinstance(self.get(arg), vtk.vtkObject):
                chaine += '\n'
                chaine += '\ncompose de {1}'.format(arg, self.get(arg))
        return chaine
    #__________________________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class RefAero(ObjetPyturbo):
    """ classe RefAero:
        classe derivant de(s) classe(s) objets.ref:
        attributs de la classe superieure:
            * nom                   : nom de l'objet (defaut = "")
            * fmt_reel                  : format des nombres reels (defaut = numpy.float32) 
            * fmt_entier                : format des nombres entiers (defaut = numpy.int32)
            
        attributs de la classe:
            * t_ref                 : temperature de reference
            * p_ref                 : pression de reference
            * l_ref                 : longueur de reference
            * gamma_ref             : rapport des capacites calorifiques de reference
            * r_gaz_ref             : constantes de l'air de reference
            * pr_ref                    : nombre de Prandlt de reference
            * prturb_ref                : nombre de Prandlt turbulent  de reference
    """
    #_____________________________________________________________________________________
    def __init__(self, nom=None, 
            t_ref=288.15,
            p_ref=101325.,
            l_ref=1.e-3,
            gamma_ref=1.4,
            r_gaz_ref=287.04,
            pr_ref=1.,
            prturb_ref=1):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_coefficient_redimensionnement_debit(self):
        """retourne un coefficient de redimensionnement d'une valeur de debit
        adimensionnee
        
        attention, il faut eventuellement la multiplier ensuite par le nombre 
        total de canaux de la roue consideree
        
        """
        coeff = self.p_ref * self.l_ref ** 2 * (
                self.gamma_ref / (self.r_gaz_ref * self.t_ref)
            ) ** (1./2)
        return coeff
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_coefficient_redimensionnement_pression(self):
        """retourne un coefficient de redimensionnement d'une valeur de pression
        adimensionnee
        """
        coeff = self.p_ref * self.gamma_ref
        return coeff
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_coefficient_redimensionnement_vitesse(self):
        """retourne un coefficient de redimensionnement d'une valeur de vitesse
        adimensionnee
        """
        coeff = numpy.sqrt(self.gamma_ref * self.r_gaz_ref * self.t_ref)
        return coeff
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_coefficient_redimensionnement_temps(self):
        """retourne un coefficient de redimensionnement d'une valeur de temps
        adimensionnee
        """
        coeff = self.l_ref / self.get_coefficient_redimensionnement_vitesse()
        return coeff
    #_____________________________________________________________________________________
#__________________________________________________________________________________________

#__________________________________________________________________________________________
class Maillage(ObjetPyturbo):
    """classe du maillage
    
    contient le maillage, sans aucune autre donnee que la geometrie
        * numeros_blocs_par_roue :   None par defaut
                            si multibloc, les numeros des blocs par roue
    
    """
    #_____________________________________________________________________________________
    def __init__(self, nom=None, nombre_aubes_par_bloc=None, unite_longueur=1e-3, 
            _vtkDataObject = None):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def __str__(self):
        chaine = ObjetPyturbo.__str__(self)
        chaine += '\n\tvtkDataObject_is_set : {0}'.format(not self._vtkDataObject is None)
        return chaine
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def set_vtkDataObject(self, vtkDataObject):
        if not isinstance(vtkDataObject, vtk.vtkMultiBlockDataSet):
            raise IOError, 'le maillage doit etre un MultiBlockDataSet'
        self.set('_vtkDataObject', vtkDataObject)
    #_____________________________________________________________________________________
    
    
    #_____________________________________________________________________________________
    def get_vtkDataObject(self):
        return self.get('_vtkDataObject')
    #_____________________________________________________________________________________
#__________________________________________________________________________________________
