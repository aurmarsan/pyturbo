try: 
    from paraview import vtk 
except: 
    import vtk
try: 
    from paraview import numpy_support 
except: 
    from vtk.util import numpy_support

try: 
    from paraview.vtk import vtkFiltersGeneral
except:
    pass

import numpy
import copy
from fonctions_basiques import *
from objets import ObjetPyturbo
from objets import RefAero
from UVParametrizationFilter import UVParametrization

#_____________________________________________________________________________________
class CalculetteGenerique(ObjetPyturbo):
    """utilise un vtkArrayCalculator pour effectuer le calcul demande
    
    s'adapte au type de vtkDataObject donne en entree
        - MultiBlockDataSet
        - PolyData
        - StructuredGrid
    
    nom_du_resultat peut etre laisse a None
        auquel cas la formule est utilisee comme nom de l'array resultat
    
    'variables_scalaires' et 'variables_vectorielles' permettent d'indiquer 
        simplement des variables a utiliser telles quelles dans la formule
    
    pour une definition plus precise d'une variable, notamment dans le cas ou le nom
        de la variable dans la formule n'est pas le meme que celui de l'array a utiliser
        utiliser les fonctions ajouter_variable_scalaires et ajouter_variable_vectorielle.
    """

    #_____________________________________________________________________________________
    def __init__(self, input=None, formule=None, nom_du_resultat=None, \
            variables_scalaires=[], variables_vectorielles=[], 
            resultat_en_coordonnees=False):
        # initialisation
        self.input = input
        self.formule = formule
        self.nom_du_resultat = formule if nom_du_resultat is None else nom_du_resultat
        self._mettre_a_jour = True
        self.variables_scalaires = []
        self.variables_vectorielles = []
        self.resultat_en_coordonnees = resultat_en_coordonnees
        
        for variable in variables_scalaires:
            self.ajouter_variable_scalaire(variable, variable)
        for variable in variables_vectorielles:
            self.ajouter_variable_vectorielle(variable, variable)
    #_____________________________________________________________________________________
    
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

    #_____________________________________________________________________________________
    def ajouter_variable_scalaire(self, nom_variable, nom_array, composante=0):
        """ajoute une variable scalaire
        
        nom_variable specifie le nom utilise dans la formule pour faire reference au scalaire
        nom_array est le nom de l'array qui contient le scalaire
        composante specifie la composante de cet array a utiliser"""
        self.variables_scalaires.append((
            nom_variable, nom_array, composante))
        self._mettre_a_jour = True
    #_____________________________________________________________________________________
          
    #_____________________________________________________________________________________  
    def ajouter_variable_vectorielle(self, nom_variable, nom_array, 
            composante_0=0, composante_1=1, composante_2=2):
        """ajoute une variable vectorielle
        
        nom_variable specifie le nom utilise dans la formule pour faire reference au vecteur
        nom_array est le nom de l'array qui contient le vecteur
        composante_ specifient les composante de cet array a utiliser"""
        self.variables_vectorielles.append((
            nom_variable, nom_array, composante_0, composante_1, composante_2))
        self._mettre_a_jour = True
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def __getCalculator__(self):
        c = vtk.vtkArrayCalculator()
        c.AddCoordinateScalarVariable('coordx', 0)
        c.AddCoordinateScalarVariable('coordy', 1)
        c.AddCoordinateScalarVariable('coordz', 2)
        c.AddCoordinateVectorVariable('coords', 0, 1, 2)
        for scalar_description in self.variables_scalaires:
            c.AddScalarVariable(scalar_description[0], 
                    scalar_description[1], scalar_description[2])
        for vector_description in self.variables_vectorielles:
            c.AddVectorVariable(vector_description[0], 
                    vector_description[1], vector_description[2], 
                    vector_description[3], vector_description[4])
        c.SetFunction(self.formule)
        if self.resultat_en_coordonnees:
            c.SetCoordinateResults(1)
        else:
            c.SetResultArrayName(self.nom_du_resultat)
        c.ReplaceInvalidValuesOn()
        c.SetReplacementValue(0.0)
        return c
    #_____________________________________________________________________________________

    #_____________________________________________________________________________________
    def Update(self):
        """execute le calcul"""
        if self.input is None:
            raise IOError, "input n'est pas renseigne"
        # execution du calcul
        if isinstance(self.input, vtk.vtkMultiBlockDataSet):
            self.output = vtk.vtkMultiBlockDataSet()
            for numbloc in get_numeros_blocs_non_vides(self.input):
                c = self.__getCalculator__()
                vtk_set_input(c, self.input.GetBlock(numbloc))
                c.Update()
                self.output.SetBlock(numbloc, 
                    vtk_new_shallowcopy(c.GetOutput()))
        else:
            self.output = vtk_new_instance(self.input)
            c = self.__getCalculator__()
            vtk_set_input(c, self.input)
            c.Update()
            self.output = vtk_new_shallowcopy(c.GetOutput())
        # on indique que la mise a jour a ete effectuee
        self._mettre_a_jour = False
        return 0
    
    #_____________________________________________________________________________________
    def get_output(self):
        if self._mettre_a_jour: 
            self.Update()
        return self.output
    #_____________________________________________________________________________________
#_____________________________________________________________________________________
    

#_____________________________________________________________________________________
class CalculettePyturbo(ObjetPyturbo):
    """calculette qui sait comment calculer les principales grandeurs aerodynamiques
    donner une formule SANS ESPACES
    
    pour le calcul de la majorite des grandeurs, la vitesse de rotation 
    doit etre disponible comme grandeurs stockee aux noeuds
    (a changer a l'avenir, quand les pipeline vtk gereront mieux les FieldData)
    
    Les variables de base en entree sont : ro, roe, momentum, omega. 
    
    les noms des arrays utilises peuvent etre changes par l'utilisateur
    utiliser print !
    
    A AMELIORER POUR POUVOIR UTILISER GetOutputPort en sortie ... 
    --> architecture pipeline. 
    --> difficulte pour les multiblockdataset du setblock en pipeline
    
    ATTENTION ATTENTION 
    Indiquer l'unite du maillage
    utilisee pour faire omega * coordr et passer de vabs a vrel ou inversement
    
    """

    #_____________________________________________________________________________________
    def __init__(self, 
            input=None, a_calculer = None, 
            nom_resultat = None, 
            axe=2, 
            RefAero=RefAero(), \
            unite_maillage = 1e-3,
            momentumRelativeFormulation=True, \
            keepIntermediateVariables=False, 
            # hubFileName = "/home/amarsan/post_doc/data/moyeu_zr",
            # tipFileName = "/home/amarsan/post_doc/data/carter_zr",
            # hubFileName = "/media/FreeAgent GoFlex Drive/DATA_PI4/hub",
            # tipFileName = "/media/FreeAgent GoFlex Drive/DATA_PI4/shroud",
            use_cell_data = False,
            ):
        """fonction d'initialisation
        
        c'est ici qu'est defini le dictionnaire contenant les formules pour le 
        calcul des grandeurs
        
        les noms des arrays a utiliser sont aussi definis
            - vitesse
            - moment cinetique
            - masse volumique
            - vitesse de rotation
            - etc... 
        
        axe doit etre specifie pour pouvoir permettre le calcul de coordr et coordtheta
        0 = x, 1 = y, 2 = z
        
        """
        #initialisation de la classe parente
        attributs = locals().copy()
        del attributs['self']
        ObjetPyturbo.__init__(self, **attributs)
        
        # initialisation particuliere
        self._mettre_a_jour = True
        
        #definition des noms qui vont etre utilises pour le calcul
        # ils peuvent etre changes par l'utilisateur
        self.densityArrayName = 'ro'
        self.totalEnergyPerUnitOfVolumeArrayName = 'roe'
        self.momentumArrayName = 'momentum'
        self.omegaArrayName = 'omega'
        
        self.relativeVelocityArrayName = 'vrel'
        self.absoluteVelocityArrayName = 'vabs'
        self.absoluteCineticEnergyArrayName = 'ecin'
        self.relativeCineticEnergyArrayName = 'ecinrel'
        self.internalEnergyArrayName = 'e_interne'
        self.staticTemperatureArrayName = 'ts'
        self.absoluteTotalTemperatureArrayName = 'tt'
        self.relativeTotalTemperatureArrayName = 'ttrel'
        self.staticPressureArrayName = 'ps'
        self.absoluteTotalPressureArrayName = 'pt'
        self.relativeTotalPressureArrayName = 'ptrel'
        self.absoluteMachNumberArrayName = 'mabs'
        self.relativeMachNumberArrayName = 'mrel'
        self.entropyArrayName = 's'
        self.radialCoordinateArrayName = 'coordr'
        self.angularCoordinateArrayName = 'coordtheta'
        
        self.radialUnitVectorArrayName = 'er'
        self.angularUnitVectorArrayName = 'etheta'
        
        self.rtRelativeAngleArrayName = 'alphaRTrel'
        self.rtAbsoluteAngleArrayName = 'alphaRTabs'
        
        self.xRelativeAngleArrayName = 'alphaXrel'
        self.xAbsoluteAngleArrayName = 'alphaXabs'
        
        self.AbsoluteMeridionalAngleArrayName = 'alpha_m'
        self.RelativeMeridionalAngleArrayName = 'alpha_m_rel'
        
       
        self.dictionnaire_des_formules = {
            'RelativeVelocity': [
                {'omega': self.omegaArrayName, 
                    'ro': self.densityArrayName, 
                    'momentum': self.momentumArrayName},
                'momentum * 1 / ro' if self.momentumRelativeFormulation else 
                    'momentum * 1 / ro + omega * coordz * {0} * jHat - omega * coordy * {0} * kHat'.format(self.unite_maillage) if axe == 0
                    else 'momentum * 1 / ro + omega * coordx * {0} * kHat - omega * coordz * {0} * iHat'.format(self.unite_maillage) if axe == 1
                    else 'momentum * 1 / ro + omega * coordy * {0} * iHat - omega * coordx * {0} * jHat'.format(self.unite_maillage),
                self.relativeVelocityArrayName],
                
            'AbsoluteVelocity': [
                {'omega': self.omegaArrayName, 
                    'ro': self.densityArrayName,
                    'momentum': self.momentumArrayName},
                'momentum * 1 / ro' if not(self.momentumRelativeFormulation) else
                    'momentum * 1 / ro - omega * coordz * {0} * jHat + omega * coordy * {0} * kHat'.format(self.unite_maillage) if axe == 0
                    else 'momentum * 1 / ro - omega * coordx * {0} * kHat + omega * coordz * {0} * iHat'.format(self.unite_maillage) if axe == 1
                    else 'momentum * 1 / ro - omega * coordy * {0} * iHat + omega * coordx * {0} * jHat'.format(self.unite_maillage),
                self.absoluteVelocityArrayName],
                
            'AbsoluteCineticEnergy': [
                {'vabs': self.absoluteVelocityArrayName},
                'vabs . vabs * 1 / 2',
                self.absoluteCineticEnergyArrayName],
                
            'RelativeCineticEnergy': [
                {'vrel': self.relativeVelocityArrayName},
                'vrel . vrel * 1 / 2',
                self.relativeCineticEnergyArrayName],
            
            'InternalEnergy': [
                {'ro': self.densityArrayName,
                    'roEt': self.totalEnergyPerUnitOfVolumeArrayName,
                    'ecin': self.relativeCineticEnergyArrayName if self.momentumRelativeFormulation \
                        else self.absoluteCineticEnergyArrayName},
                'roEt / ro - ecin',
                self.internalEnergyArrayName],
            
            'StaticTemperature': [
                {'e_interne': self.internalEnergyArrayName},
                'e_interne * ({0} - 1) / {1}'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                self.staticTemperatureArrayName],
            
            'AbsoluteTotalTemperature': [
                {'ts': self.staticTemperatureArrayName,
                    'ecin': self.absoluteCineticEnergyArrayName},
                'ts + ecin * ( {0} - 1) / ( {0} * {1} )'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                self.absoluteTotalTemperatureArrayName],
            
            'RelativeTotalTemperature': [
                {'ts': self.staticTemperatureArrayName,
                    'ecinrel': self.relativeCineticEnergyArrayName},
                'ts + ecinrel * ( {0} - 1) / ( {0} * {1} )'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                self.relativeTotalTemperatureArrayName],
            
            'StaticPressure': [
                {'ro': self.densityArrayName,
                    'ts': self.staticTemperatureArrayName},
                'ro * {0} * ts'.format(self.RefAero.r_gaz_ref),
                self.staticPressureArrayName],
            
            'AbsoluteTotalPressure': [
                {'ps': self.staticPressureArrayName,
                    'tt': self.absoluteTotalTemperatureArrayName,
                    'ts': self.staticTemperatureArrayName},
                'ps * (tt / ts) ^ ({0} / ({0} - 1))'.format(self.RefAero.gamma_ref),
                self.absoluteTotalPressureArrayName],
            
            'RelativeTotalPressure': [
                {'ps': self.staticPressureArrayName,
                    'ttrel': self.relativeTotalTemperatureArrayName,
                    'ts': self.staticTemperatureArrayName},
                'ps * (ttrel / ts) ^ ({0} / ({0} - 1))'.format(self.RefAero.gamma_ref),
                self.relativeTotalPressureArrayName],
            
            'AbsoluteMachNumber': [
                {'ts': self.staticTemperatureArrayName,
                    'vabs': self.absoluteVelocityArrayName},
                'mag(vabs) / sqrt({0} * {1} * ts)'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                self.absoluteMachNumberArrayName],
            
            'RelativeMachNumber': [
                {'ts': self.staticTemperatureArrayName, 
                    'vrel': self.relativeVelocityArrayName},
                'mag(vrel) / sqrt({0} * {1} * ts)'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                self.relativeMachNumberArrayName],
            
            'Entropy': [
                {'ts': self.staticTemperatureArrayName,
                    'ps': self.staticPressureArrayName},
                '{0} * {1} / ({0} - 1) * ln(ts / {2}) - {1} * ln(ps / {3})'
                .format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref, self.RefAero.t_ref, self.RefAero.p_ref),
                self.entropyArrayName],
                        'RadialCoordinate': [
                {},
                'sqrt(coordy ^ 2 + coordz ^ 2)' if axe == 0 
                    else 'sqrt(coordx ^ 2 + coordz ^ 2)' if axe == 1 
                    else 'sqrt(coordx ^ 2 + coordy ^ 2)',
                self.radialCoordinateArrayName],
            
            'AngularCoordinate': [
                {'coordr':  self.radialCoordinateArrayName},
                'acos(coordy / coordr) * coordz / abs(coordz) + 2 * acos(-1.0)  * (1 - coordz / abs(coordz))/2.0' if axe == 0
                    else 'acos(coordz / coordr) * coordx / abs(coordx) + 2 * acos(-1.0)  * (1 - coordx / abs(coordx))/2.0' if axe == 1
                    else 'acos(coordx / coordr) * coordy / abs(coordy) + 2 * acos(-1.0)  * (1 - coordy / abs(coordy))/2.0',
                self.angularCoordinateArrayName],
            
            'RadialUnitVector': [
                {'coordtheta':  self.angularCoordinateArrayName},
                'cos(coordtheta) * jHat + sin(coordtheta) * kHat' if axe == 0
                    else 'cos(coordtheta) * kHat + sin(coordtheta) * iHat' if axe == 1
                    else 'cos(coordtheta) * iHat + sin(coordtheta) * jHat',
                self.radialUnitVectorArrayName],
            
            'AngularUnitVector': [
                {'coordtheta':  self.angularCoordinateArrayName},
                '-sin(coordtheta) * jHat + cos(coordtheta) * kHat' if axe == 0
                    else '-sin(coordtheta) * kHat + cos(coordtheta) * iHat' if axe == 1
                    else '-sin(coordtheta) * iHat + cos(coordtheta) * jHat',
                self.angularUnitVectorArrayName],
            
            'YZRelativeAngle': [
                {
                    'coordtheta':  self.angularCoordinateArrayName, 
                    'vrel': self.relativeVelocityArrayName,
                    'er': self.radialUnitVectorArrayName,
                    'etheta': self.angularUnitVectorArrayName},
                'acos( (vrel . er) / mag(vrel) ) * sign(vrel . etheta) * 90.0 / acos(0.0)',
                self.rtRelativeAngleArrayName],
            
            'YZAbsoluteAngle': [
                {
                    'coordtheta':  self.angularCoordinateArrayName, 
                    'vabs': self.absoluteVelocityArrayName,
                    'er': self.radialUnitVectorArrayName,
                    'etheta': self.angularUnitVectorArrayName},
                'acos( (vabs . er) / mag(vabs) ) * sign(vabs . etheta) * 90.0 / acos(0.0)',
                self.rtAbsoluteAngleArrayName],
                        
            
            'XRelativeAngle': [
                {
                    'vrel': self.relativeVelocityArrayName,
                    'er': self.radialUnitVectorArrayName},
                'acos( (vrel - (vrel . er) * er) . iHat / mag((vrel - (vrel . er) * er)) ) * sign((vrel - (vrel . er) * er) . etheta) * 90.0 / acos(0.0)',
                self.xRelativeAngleArrayName],
            
            'XAbsoluteAngle': [
                {
                    'vabs': self.absoluteVelocityArrayName,
                    'er': self.radialUnitVectorArrayName,
                    'etheta': self.angularUnitVectorArrayName},
                'acos(((vabs - (vabs . er) * er) . iHat)/ mag(vabs - (vabs . er) * er)) * sign((vabs - (vabs . er) * er) . etheta) * 90.0 / acos(0.0)',
                self.xAbsoluteAngleArrayName],
            
            'XCoordinate': [
                {},
                'coordx',
                'coordx'],
            
            'YCoordinate': [
                {},
                'coordy',
                'coordy'],
            
            'ZCoordinate': [
                {},
                'coordz',
                'coordz'],
            
            'UVParametrization_RelativeMeridionalAbscissa': [
                {},
                'UVParametrization',
                'xm'],
            
            'UVParametrization_hsH': [
                {},
                'UVParametrization',
                'hsH'],

            'gradPs_adv': [
                {
                    'vabs': self.absoluteVelocityArrayName,
                    'grad(ps)': 'grad(' + self.staticPressureArrayName + ')'
                    },
                'grad(ps).vabs/mag(vabs)',
                'gradPs_adv'],
            
            'angle_meridien_absolu': [
                {
                    'vabs': self.absoluteVelocityArrayName,
                    'er': self.radialUnitVectorArrayName, 
                    'etheta': self.angularUnitVectorArrayName
                    },
                'acos( mag(vabs - (vabs.etheta) * etheta) / mag(vabs) ) * sign(vabs . etheta) * 90.0 / acos(0.0)',
                self.AbsoluteMeridionalAngleArrayName],
            
            'angle_meridien_relatif': [
                {
                    'vrel': self.relativeVelocityArrayName,
                    'er': self.radialUnitVectorArrayName, 
                    'etheta': self.angularUnitVectorArrayName
                    },
                'acos( mag(vrel - (vrel.etheta) * etheta) / mag(vrel) ) * sign(vrel . etheta) * 90.0 / acos(0.0)',
                self.RelativeMeridionalAngleArrayName],
            
            #'Q_criterion': [
                #{},
                #'Q_criterion',
                #'Q_criterion'],
            
            }
    #_____________________________________________________________________________________
    
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

    #_____________________________________________________________________________________
    def __input_has_array__(self, ArrayName):
        """retourne True si input a un array ArrayName aux points
        """
        #si c'est un multiblockdataset on verifie que 
        #l'array est present dans tous les blocs
        if self.get('input') is None:
            raise IOError, "indiquez l'objet VTK sur lequel effectuer le calcul"
        if self.use_cell_data == False:
            return bool(self.input.GetPointData().HasArray(ArrayName))
        elif self.use_cell_data == True:
            return bool(self.input.GetCellData().HasArray(ArrayName))
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_output(self):
        if self._mettre_a_jour:
            self.Update()
        return self.output
    #_____________________________________________________________________________________
        
    #_____________________________________________________________________________________
    def __get_what_to_do__(self):
        if hasattr(self, 'a_calculer') is False:
            raise IOError, 'indiquez les variables a calculer'
        to_do = []
        for quantity in self.a_calculer:
            if ' ' in quantity:
                raise IOError, 'Indiquer la formule suivante SANS ESPACES -- {0}'.format(quantity)
            if self.__input_has_array__(quantity) == False:
                if quantity in numpy.asarray(self.dictionnaire_des_formules.values())[:, -1]:
                    index = numpy.where(numpy.asarray(
                        self.dictionnaire_des_formules.values())[:, -1] == quantity)[0]
                    to_do.append(self.dictionnaire_des_formules.values()[index])
                else:
                    previous_variables = dict.fromkeys(get_variables_in_function(quantity))
                    for key in previous_variables.keys():
                        previous_variables[key] = key
                    dict_quantity = [
                        previous_variables, 
                        quantity.replace(' ', ''),
                        quantity.replace(' ', '')]
                    to_do.append(dict_quantity)
        return to_do
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def SimilarInstance(self):
        """cree une instance similaire
        ne copie par nom_resultat
        """
        newCalculator = CalculettePyturbo()
        for arg in dir(self):
            if not callable(self.get(arg)) and (arg[0].islower() or arg[0].isupper()) \
                    and arg != 'input' and arg != 'output':
                setattr(newCalculator, arg, getattr(self, arg))
        newCalculator.set('input', self.input)
        newCalculator.set('nom_resultat', None)
        return newCalculator
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def Update(self):
        
        self.output = vtk_new_instance(self.input)
        
        # traitement recursif du cas multibloc
        if isinstance(self.output, vtk.vtkMultiBlockDataSet):
            for numbloc in get_numeros_blocs_non_vides(self.input):
                calc_bloc = self.SimilarInstance()
                calc_bloc.input = self.input.GetBlock(numbloc)
                self.output.SetBlock(numbloc, calc_bloc.get_output())
            self._mettre_a_jour = False
            return 0
        
        # cas monobloc
        # dans le cas ou une seule variable est demande, 
        # il faut quand meme que a_calculer soit un tuple
        if isinstance(self.a_calculer, str):
            self.a_calculer = [self.a_calculer]
        
        to_do = self.__get_what_to_do__()
        
        variables_to_have = []
        for i in to_do:
            variables_to_have += i[0].values()
        for i in variables_to_have:
            while variables_to_have.count(i) != 1:
                variables_to_have.remove(i)
        if len(variables_to_have) != 0:
            newCalculator = self.SimilarInstance()
            newCalculator.keepIntermediateVariables = True
            newCalculator.set('a_calculer', list(variables_to_have))
            try: 
                self.output = newCalculator.get_output()
            except:
                print "n'arrive pas a obtenir {0}".format(variables_to_have)
                raise IOError, "impossible de derouler le pipe de calcul"
        else:
            self.output = vtk_new_instance(self.input)
            self.output.ShallowCopy(self.input)
            self.output.SetFieldData(self.input.GetFieldData())
        
        for function in to_do:
            # si il y a quelque chose a faire, mais que la formule associee est vide
            # c'est que le calculateur est perdu
            if function[1] == '':
                raise IOError
            
            if function[2] in get_noms_arrays_presents(self.output):
                #si le array a calculer est deja present au noeuds de self.output, c'est pas la peine de 
                #le recalculer
                pass
            elif function[1] == 'UVParametrization':
                raise Exception, "NE PLUS UTILISER CETTE FONCTION POUR LE CALCUL DE hsH MAIS LA NOUVELLE CLASSE PARAMETRISATION"
                self.output = UVParametrization(self.output, 
                    hubFileName = self.hubFileName, tipFileName = self.tipFileName, 
                    axe = self.axe)
            #elif function[1] == 'Q_criterion':
                #self.output = Q_criterion(self.output, self.relativeVelocityArrayName)
            elif len(function[0]) == 1 and function[1] == 'grad(' + function[0].values()[0] + ')':
                current_bloc = self.output
                try:
                    gradient_calculator = vtkFiltersGeneral.vtkGradientFilter()
                    vtk_set_input(gradient_calculator, current_bloc)
                except:
                    gradient_calculator = vtk.vtkGradientFilter()
                    vtk_set_input(gradient_calculator, current_bloc)
                gradient_calculator.SetInputScalars(0, function[0].values()[0])
                gradient_calculator.SetResultArrayName(function[2])
                gradient_calculator.Update()
                
                current_bloc.ShallowCopy(gradient_calculator.GetOutput())
                # if isinstance(self.output, vtk.vtkMultiBlockDataSet):
                    # for numbloc in get_numeros_blocs_non_vides(self.output):
                        # gradient_calculator = vtk.vtkGradientFilter()
                        # gradient_calculator.SetInputData(self.output.GetBlock(numbloc))
                        # gradient_calculator.SetInputScalars(0, function[0].values()[0])
                        # gradient_calculator.SetResultArrayName(function[2])
                        # gradient_calculator.Update()
                        # self.output.SetBlock(numbloc, gradient_calculator.GetOutput())
                # else:
                    # gradient_calculator = vtk.vtkGradientFilter()
                    # gradient_calculator.SetInputData(self.output)
                    # gradient_calculator.SetInputScalars(0, function[0].values()[0])
                    # gradient_calculator.SetResultArrayName(function[2])
                    # gradient_calculator.Update()
                    # self.output = gradient_calculator.GetOutput()
            else:
                current_bloc = self.output
                
                calc = vtk.vtkArrayCalculator()
                if self.use_cell_data:
                    calc.SetAttributeModeToUseCellData()
                vtk_set_input(calc, current_bloc)
                calc.SetFunction(function[1])
                for var_input in function[0].items():
                    if self.use_cell_data == True:
                        if current_bloc.GetCellData().GetArray(var_input[1]).GetNumberOfComponents() == 3:
                            calc.AddVectorVariable(var_input[0], var_input[1], 0, 1, 2)
                        else:
                            calc.AddScalarVariable(var_input[0], var_input[1], 0)
                    elif self.use_cell_data == False:
                        if current_bloc.GetPointData().GetArray(var_input[1]).GetNumberOfComponents() == 3:
                            calc.AddVectorVariable(var_input[0], var_input[1], 0, 1, 2)
                        else:
                            calc.AddScalarVariable(var_input[0], var_input[1], 0)
                calc.AddCoordinateScalarVariable('coordx', 0)
                calc.AddCoordinateScalarVariable('coordy', 1)
                calc.AddCoordinateScalarVariable('coordz', 2)
                calc.SetResultArrayName(function[2])
                calc.ReplaceInvalidValuesOn()
                calc.SetReplacementValue(0.0)
                calc.Update()
                current_bloc.ShallowCopy(calc.GetOutput())
        
        if self.keepIntermediateVariables == False:
            cleanOutput = vtk_new_shallowcopy(self.output)
            list_to_keep = list(self.a_calculer) + \
                get_noms_arrays_presents(self.input, loc = 'points')

            for quantity in get_noms_arrays_presents(cleanOutput, loc = 'points'):
                if not(quantity in list_to_keep):
                    cleanOutput.GetPointData().RemoveArray(quantity)
            self.output = cleanOutput
        
        #si un nom du resultat est donne, on change le nom de l'array au point
        if self.nom_resultat != None:
            #on convertit d'abord nom_result en une liste si ce n'en est pas une
            if not isinstance(self.nom_resultat, list):
                self.nom_resultat = [self.nom_resultat]
            
            #on verifie de nom_resultat et a_calculer font les memes longueurs
            if len(self.nom_resultat) != len(self.a_calculer):
                raise IOError, "il n'y a pas le meme nombre de a_calculer et nom_resultat"

            #cas multibloc
            for k in range(len(self.a_calculer)):
                avant = self.a_calculer[k]
                apres = self.nom_resultat[k]
                self.output.GetPointData().GetArray(avant).SetName(apres)
            
        self._mettre_a_jour = False
        return 0
    #_____________________________________________________________________________________
#_____________________________________________________________________________________

##_____________________________________________________________________________________
#class CalculettePyturbo(ObjetPyturbo):
    #"""calculette qui sait comment calculer les principales grandeurs aerodynamiques
    #donner une formule SANS ESPACES
    
    #pour le calcul de la majorite des grandeurs, la vitesse de rotation 
    #doit etre disponible comme grandeurs stockee aux noeuds
    #(a changer a l'avenir, quand les pipeline vtk gereront mieux les FieldData)
    
    #Les variables de base en entree sont : ro, roe, momentum, omega. 
    
    #les noms des arrays utilises peuvent etre changes par l'utilisateur
    #utiliser print !
    
    #A AMELIORER POUR POUVOIR UTILISER GetOutputPort en sortie ... 
    #--> architecture pipeline. 
    #--> difficulte pour les multiblockdataset du setblock en pipeline
    
    #ATTENTION ATTENTION 
    #Indiquer l'unite du maillage
    #utilisee pour faire omega * coordr et passer de vabs a vrel ou inversement
    
    #"""

    ##_____________________________________________________________________________________
    #def __init__(self, 
            #input=None, a_calculer = None, 
            #nom_resultat = None, 
            #axe=2, 
            #RefAero=RefAero(), \
            #unite_maillage = 1e-3,
            #momentumRelativeFormulation=True, \
            #keepIntermediateVariables=False, 
            ## hubFileName = "/home/amarsan/post_doc/data/moyeu_zr",
            ## tipFileName = "/home/amarsan/post_doc/data/carter_zr",
            ## hubFileName = "/media/FreeAgent GoFlex Drive/DATA_PI4/hub",
            ## tipFileName = "/media/FreeAgent GoFlex Drive/DATA_PI4/shroud",
            #use_cell_data = False,
            #):
        #"""fonction d'initialisation
        
        #c'est ici qu'est defini le dictionnaire contenant les formules pour le 
        #calcul des grandeurs
        
        #les noms des arrays a utiliser sont aussi definis
            #- vitesse
            #- moment cinetique
            #- masse volumique
            #- vitesse de rotation
            #- etc... 
        
        #axe doit etre specifie pour pouvoir permettre le calcul de coordr et coordtheta
        #0 = x, 1 = y, 2 = z
        
        #"""
        ##initialisation de la classe parente
        #attributs = locals().copy()
        #del attributs['self']
        #ObjetPyturbo.__init__(self, **attributs)
        
        ## initialisation particuliere
        #self._mettre_a_jour = True
        
        ##definition des noms qui vont etre utilises pour le calcul
        ## ils peuvent etre changes par l'utilisateur
        #self.densityArrayName = 'ro'
        #self.totalEnergyPerUnitOfVolumeArrayName = 'roe'
        #self.momentumArrayName = 'momentum'
        #self.omegaArrayName = 'omega'
        
        #self.relativeVelocityArrayName = 'vrel'
        #self.absoluteVelocityArrayName = 'vabs'
        #self.absoluteCineticEnergyArrayName = 'ecin'
        #self.relativeCineticEnergyArrayName = 'ecinrel'
        #self.internalEnergyArrayName = 'e_interne'
        #self.staticTemperatureArrayName = 'ts'
        #self.absoluteTotalTemperatureArrayName = 'tt'
        #self.relativeTotalTemperatureArrayName = 'ttrel'
        #self.staticPressureArrayName = 'ps'
        #self.absoluteTotalPressureArrayName = 'pt'
        #self.relativeTotalPressureArrayName = 'ptrel'
        #self.absoluteMachNumberArrayName = 'mabs'
        #self.relativeMachNumberArrayName = 'mrel'
        #self.entropyArrayName = 's'
        #self.radialCoordinateArrayName = 'coordr'
        #self.angularCoordinateArrayName = 'coordtheta'
        
        #self.radialUnitVectorArrayName = 'er'
        #self.angularUnitVectorArrayName = 'etheta'
        
        #self.rtRelativeAngleArrayName = 'alphaRTrel'
        #self.rtAbsoluteAngleArrayName = 'alphaRTabs'
        
        #self.xRelativeAngleArrayName = 'alphaXrel'
        #self.xAbsoluteAngleArrayName = 'alphaXabs'
        
        #self.AbsoluteMeridionalAngleArrayName = 'alpha_m'
        #self.RelativeMeridionalAngleArrayName = 'alpha_m_rel'
        
       
        #self.dictionnaire_des_formules = {
            #'RelativeVelocity': [
                #{'omega': self.omegaArrayName, 
                    #'ro': self.densityArrayName, 
                    #'momentum': self.momentumArrayName},
                #'momentum * 1 / ro' if self.momentumRelativeFormulation else 
                    #'momentum * 1 / ro + omega * coordz * {0} * jHat - omega * coordy * {0} * kHat'.format(self.unite_maillage) if axe == 0
                    #else 'momentum * 1 / ro + omega * coordx * {0} * kHat - omega * coordz * {0} * iHat'.format(self.unite_maillage) if axe == 1
                    #else 'momentum * 1 / ro + omega * coordy * {0} * iHat - omega * coordx * {0} * jHat'.format(self.unite_maillage),
                #self.relativeVelocityArrayName],
                
            #'AbsoluteVelocity': [
                #{'omega': self.omegaArrayName, 
                    #'ro': self.densityArrayName,
                    #'momentum': self.momentumArrayName},
                #'momentum * 1 / ro' if not(self.momentumRelativeFormulation) else
                    #'momentum * 1 / ro - omega * coordz * {0} * jHat + omega * coordy * {0} * kHat'.format(self.unite_maillage) if axe == 0
                    #else 'momentum * 1 / ro - omega * coordx * {0} * kHat + omega * coordz * {0} * iHat'.format(self.unite_maillage) if axe == 1
                    #else 'momentum * 1 / ro - omega * coordy * {0} * iHat + omega * coordx * {0} * jHat'.format(self.unite_maillage),
                #self.absoluteVelocityArrayName],
                
            #'AbsoluteCineticEnergy': [
                #{'vabs': self.absoluteVelocityArrayName},
                #'vabs . vabs * 1 / 2',
                #self.absoluteCineticEnergyArrayName],
                
            #'RelativeCineticEnergy': [
                #{'vrel': self.relativeVelocityArrayName},
                #'vrel . vrel * 1 / 2',
                #self.relativeCineticEnergyArrayName],
            
            #'InternalEnergy': [
                #{'ro': self.densityArrayName,
                    #'roEt': self.totalEnergyPerUnitOfVolumeArrayName,
                    #'ecin': self.relativeCineticEnergyArrayName if self.momentumRelativeFormulation \
                        #else self.absoluteCineticEnergyArrayName},
                #'roEt / ro - ecin',
                #self.internalEnergyArrayName],
            
            #'StaticTemperature': [
                #{'e_interne': self.internalEnergyArrayName},
                #'e_interne * ({0} - 1) / {1}'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                #self.staticTemperatureArrayName],
            
            #'AbsoluteTotalTemperature': [
                #{'ts': self.staticTemperatureArrayName,
                    #'ecin': self.absoluteCineticEnergyArrayName},
                #'ts + ecin * ( {0} - 1) / ( {0} * {1} )'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                #self.absoluteTotalTemperatureArrayName],
            
            #'RelativeTotalTemperature': [
                #{'ts': self.staticTemperatureArrayName,
                    #'ecinrel': self.relativeCineticEnergyArrayName},
                #'ts + ecinrel * ( {0} - 1) / ( {0} * {1} )'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                #self.relativeTotalTemperatureArrayName],
            
            #'StaticPressure': [
                #{'ro': self.densityArrayName,
                    #'ts': self.staticTemperatureArrayName},
                #'ro * {0} * ts'.format(self.RefAero.r_gaz_ref),
                #self.staticPressureArrayName],
            
            #'AbsoluteTotalPressure': [
                #{'ps': self.staticPressureArrayName,
                    #'tt': self.absoluteTotalTemperatureArrayName,
                    #'ts': self.staticTemperatureArrayName},
                #'ps * (tt / ts) ^ ({0} / ({0} - 1))'.format(self.RefAero.gamma_ref),
                #self.absoluteTotalPressureArrayName],
            
            #'RelativeTotalPressure': [
                #{'ps': self.staticPressureArrayName,
                    #'ttrel': self.relativeTotalTemperatureArrayName,
                    #'ts': self.staticTemperatureArrayName},
                #'ps * (ttrel / ts) ^ ({0} / ({0} - 1))'.format(self.RefAero.gamma_ref),
                #self.relativeTotalPressureArrayName],
            
            #'AbsoluteMachNumber': [
                #{'ts': self.staticTemperatureArrayName,
                    #'vabs': self.absoluteVelocityArrayName},
                #'mag(vabs) / sqrt({0} * {1} * ts)'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                #self.absoluteMachNumberArrayName],
            
            #'RelativeMachNumber': [
                #{'ts': self.staticTemperatureArrayName, 
                    #'vrel': self.relativeVelocityArrayName},
                #'mag(vrel) / sqrt({0} * {1} * ts)'.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref),
                #self.relativeMachNumberArrayName],
            
            #'Entropy': [
                #{'ts': self.staticTemperatureArrayName,
                    #'ps': self.staticPressureArrayName},
                #'{0} * {1} / ({0} - 1) * ln(ts / {2}) - {1} * ln(ps / {3})'
                #.format(self.RefAero.gamma_ref, self.RefAero.r_gaz_ref, self.RefAero.t_ref, self.RefAero.p_ref),
                #self.entropyArrayName],
            
            #'RadialCoordinate': [
                #{},
                #'sqrt(coordy ^ 2 + coordz ^ 2)' if axe == 0 
                    #else 'sqrt(coordx ^ 2 + coordz ^ 2)' if axe == 1 
                    #else 'sqrt(coordx ^ 2 + coordy ^ 2)',
                #self.radialCoordinateArrayName],
            
            #'AngularCoordinate': [
                #{'coordr':  self.radialCoordinateArrayName},
                #'acos(coordy / coordr) * coordz / abs(coordz) + 2 * acos(-1.0)  * (1 - coordz / abs(coordz))/2.0' if axe == 0
                    #else 'acos(coordz / coordr) * coordx / abs(coordx) + 2 * acos(-1.0)  * (1 - coordx / abs(coordx))/2.0' if axe == 1
                    #else 'acos(coordx / coordr) * coordy / abs(coordy) + 2 * acos(-1.0)  * (1 - coordy / abs(coordy))/2.0',
                #self.angularCoordinateArrayName],
            
            #'RadialUnitVector': [
                #{'coordtheta':  self.angularCoordinateArrayName},
                #'cos(coordtheta) * jHat + sin(coordtheta) * kHat' if axe == 0
                    #else 'cos(coordtheta) * kHat + sin(coordtheta) * iHat' if axe == 1
                    #else 'cos(coordtheta) * iHat + sin(coordtheta) * jHat',
                #self.radialUnitVectorArrayName],
            
            #'AngularUnitVector': [
                #{'coordtheta':  self.angularCoordinateArrayName},
                #'-sin(coordtheta) * jHat + cos(coordtheta) * kHat' if axe == 0
                    #else '-sin(coordtheta) * kHat + cos(coordtheta) * iHat' if axe == 1
                    #else '-sin(coordtheta) * iHat + cos(coordtheta) * jHat',
                #self.angularUnitVectorArrayName],
            
            #'YZRelativeAngle': [
                #{
                    #'coordtheta':  self.angularCoordinateArrayName, 
                    #'vrel': self.relativeVelocityArrayName,
                    #'er': self.radialUnitVectorArrayName,
                    #'etheta': self.angularUnitVectorArrayName},
                #'acos( (vrel . er) / mag(vrel) ) * sign(vrel . etheta) * 90.0 / acos(0.0)',
                #self.rtRelativeAngleArrayName],
            
            #'YZAbsoluteAngle': [
                #{
                    #'coordtheta':  self.angularCoordinateArrayName, 
                    #'vabs': self.absoluteVelocityArrayName,
                    #'er': self.radialUnitVectorArrayName,
                    #'etheta': self.angularUnitVectorArrayName},
                #'acos( (vabs . er) / mag(vabs) ) * sign(vabs . etheta) * 90.0 / acos(0.0)',
                #self.rtAbsoluteAngleArrayName],
                        
            
            #'XRelativeAngle': [
                #{
                    #'vrel': self.relativeVelocityArrayName,
                    #'er': self.radialUnitVectorArrayName},
                #'acos( (vrel - (vrel . er) * er) . iHat / mag((vrel - (vrel . er) * er)) ) * sign((vrel - (vrel . er) * er) . jHat) * 90.0 / acos(0.0)',
                #self.xRelativeAngleArrayName],
            
            #'XAbsoluteAngle': [
                #{
                    #'vabs': self.absoluteVelocityArrayName,
                    #'er': self.radialUnitVectorArrayName,
                    #'etheta': self.angularUnitVectorArrayName},
                #'acos(((vabs - (vabs . er) * er) . iHat)/ mag(vabs - (vabs . er) * er)) * sign((vabs - (vabs . er) * er) . etheta) * 90.0 / acos(0.0)',
                #self.xAbsoluteAngleArrayName],
            
            #'XCoordinate': [
                #{},
                #'coordx',
                #'coordx'],
            
            #'YCoordinate': [
                #{},
                #'coordy',
                #'coordy'],
            
            #'ZCoordinate': [
                #{},
                #'coordz',
                #'coordz'],
            
            #'UVParametrization_RelativeMeridionalAbscissa': [
                #{},
                #'UVParametrization',
                #'xm'],
            
            #'UVParametrization_hsH': [
                #{},
                #'UVParametrization',
                #'hsH'],

            #'gradPs_adv': [
                #{
                    #'vabs': self.absoluteVelocityArrayName,
                    #'grad(ps)': 'grad(' + self.staticPressureArrayName + ')'
                    #},
                #'grad(ps).vabs/mag(vabs)',
                #'gradPs_adv'],
            
            #'angle_meridien_absolu': [
                #{
                    #'vabs': self.absoluteVelocityArrayName,
                    #'er': self.radialUnitVectorArrayName, 
                    #'etheta': self.angularUnitVectorArrayName
                    #},
                #'acos( mag(vabs - (vabs.etheta) * etheta) / mag(vabs) ) * sign(vabs . etheta) * 90.0 / acos(0.0)',
                #self.AbsoluteMeridionalAngleArrayName],
            
            #'angle_meridien_relatif': [
                #{
                    #'vrel': self.relativeVelocityArrayName,
                    #'er': self.radialUnitVectorArrayName, 
                    #'etheta': self.angularUnitVectorArrayName
                    #},
                #'acos( mag(vrel - (vrel.etheta) * etheta) / mag(vrel) ) * sign(vrel . etheta) * 90.0 / acos(0.0)',
                #self.RelativeMeridionalAngleArrayName],
            
            ##'Q_criterion': [
                ##{},
                ##'Q_criterion',
                ##'Q_criterion'],
            
            #}
    ##_____________________________________________________________________________________
    
    ##_____________________________________________________________________________________
    #def set(self, nom_attribut, valeur):
        #"""fonction set specifique
        
        #gere la variable locale _changement
        #qui sert lorsque l'on appelle la sortie
        #a savoir s'il faut recalculer
        #"""
        #setattr(self, nom_attribut, valeur)
        #if nom_attribut != '_mettre_a_jour':
            #self._mettre_a_jour = True
    ##_____________________________________________________________________________________

    ##_____________________________________________________________________________________
    #def __input_has_array__(self, ArrayName):
        #"""retourne True si input a un array ArrayName aux points
        #pour tous les blocs si Input ets un MultiBlockDataSet"""
        ##si c'est un multiblockdataset on verifie que 
        ##l'array est present dans tous les blocs
        #if self.get('input') is None:
            #raise IOError, "indiquez l'objet VTK sur lequel effectuer le calcul"
        #if isinstance(self.input, vtk.vtkMultiBlockDataSet):
            #for numbloc in get_numeros_blocs_non_vides(self.input):
                #bloc = self.input.GetBlock(numbloc)
                #if self.use_cell_data == False and bloc.GetPointData().HasArray(ArrayName) == 0:
                    #return False
                #if self.use_cell_data == True and bloc.GetCellData().HasArray(ArrayName) == 0:
                    #return False
            #return True
        #else:
            #if self.use_cell_data == False:
                #return bool(self.input.GetPointData().HasArray(ArrayName))
            #elif self.use_cell_data == True:
                #return bool(self.input.GetCellData().HasArray(ArrayName))
    ##_____________________________________________________________________________________
    
    ##_____________________________________________________________________________________
    #def get_output(self):
        #if self._mettre_a_jour:
            #self.Update()
        #return self.output
    ##_____________________________________________________________________________________
        
    ##_____________________________________________________________________________________
    #def __get_what_to_do__(self):
        #if hasattr(self, 'a_calculer') is False:
            #return IOError, 'indiquez les variables a calculer'
        #to_do = []
        #for quantity in self.a_calculer:
            #if self.__input_has_array__(quantity) == False:
                #if quantity in numpy.asarray(self.dictionnaire_des_formules.values())[:, -1]:
                    #index = numpy.where(numpy.asarray(
                        #self.dictionnaire_des_formules.values())[:, -1] == quantity)[0]
                    #to_do.append(self.dictionnaire_des_formules.values()[index])
                #else:
                    #previous_variables = dict.fromkeys(get_variables_in_function(quantity))
                    #for key in previous_variables.keys():
                        #previous_variables[key] = key
                    #dict_quantity = [
                        #previous_variables, 
                        #quantity.replace(' ', ''),
                        #quantity.replace(' ', '')]
                    #to_do.append(dict_quantity)
        #return to_do
    ##_____________________________________________________________________________________
    
    ##_____________________________________________________________________________________
    #def SimilarInstance(self):
        #"""cree une instance similaire
        #ne copie par nom_resultat
        #"""
        #newCalculator = CalculettePyturbo()
        #for arg in dir(self):
            #if not callable(self.get(arg)) and (arg[0].islower() or arg[0].isupper()) \
                    #and arg != 'input' and arg != 'output':
                #setattr(newCalculator, arg, getattr(self, arg))
        #newCalculator.set('input', self.input)
        #newCalculator.set('nom_resultat', None)
        #return newCalculator
    ##_____________________________________________________________________________________
    
    ##_____________________________________________________________________________________
    #def Update(self):
        ## dans le cas ou une seule variable est demande, 
        ## il faut quand meme que a_calculer soit un tuple
        #if isinstance(self.a_calculer, str):
            #self.a_calculer = [self.a_calculer]
        
        #to_do = self.__get_what_to_do__()
        
        #variables_to_have = []
        #for i in to_do:
            #variables_to_have += i[0].values()
        #for i in variables_to_have:
            #while variables_to_have.count(i) != 1:
                #variables_to_have.remove(i)
        #if len(variables_to_have) != 0:
            #newCalculator = self.SimilarInstance()
            #newCalculator.keepIntermediateVariables = True
            #newCalculator.set('a_calculer', list(variables_to_have))
            #try: 
                #self.output = newCalculator.get_output()
            #except:
                #print "n'arrive pas a obtenir {0}".format(variables_to_have)
                #raise IOError, "impossible de derouler le pipe de calcul"
        #else:
            #self.output = vtk_new_instance(self.input)
            #if isinstance(self.output, vtk.vtkMultiBlockDataSet):
                #for numbloc in get_numeros_blocs_non_vides(self.input):
                    #bloc = vtk_new_shallowcopy(self.input.GetBlock(numbloc))
                    #self.output.SetBlock(numbloc, bloc)
            #else:
                #self.output.ShallowCopy(self.input)
            #self.output.SetFieldData(self.input.GetFieldData())
        
        #for function in to_do:
            ## si il y a quelque chose a faire, mais que la formule associee est vide
            ## c'est que le calculateur est perdu
            #if function[1] == '':
                #raise IOError
            
            #if function[2] in get_noms_arrays_presents(self.output):
                ##si le array a calculer est deja present au noeuds de self.output, c'est pas la peine de 
                ##le recalculer
                #pass
            #elif function[1] == 'UVParametrization':
                #raise Exception, "NE PLUS UTILISER CETTE FONCTION POUR LE CALCUL DE hsH MAIS LA NOUVELLE CLASSE PARAMETRISATION"
                #self.output = UVParametrization(self.output, 
                    #hubFileName = self.hubFileName, tipFileName = self.tipFileName, 
                    #axe = self.axe)
            ##elif function[1] == 'Q_criterion':
                ##self.output = Q_criterion(self.output, self.relativeVelocityArrayName)
            #elif len(function[0]) == 1 and function[1] == 'grad(' + function[0].values()[0] + ')':
                #for numbloc in get_numeros_blocs_non_vides(self.output) \
                        #if isinstance(self.output, vtk.vtkMultiBlockDataSet) else [None]:
                    #current_bloc = self.output.GetBlock(numbloc) if numbloc != None else self.output
                    #try:
                        #gradient_calculator = vtkFiltersGeneral.vtkGradientFilter()
                        #vtk_set_input(gradient_calculator, current_bloc)
                    #except:
                        #gradient_calculator = vtk.vtkGradientFilter()
                        #vtk_set_input(gradient_calculator, current_bloc)
                    #gradient_calculator.SetInputScalars(0, function[0].values()[0])
                    #gradient_calculator.SetResultArrayName(function[2])
                    #gradient_calculator.Update()
                    
                    #current_bloc.ShallowCopy(gradient_calculator.GetOutput())
                ## if isinstance(self.output, vtk.vtkMultiBlockDataSet):
                    ## for numbloc in get_numeros_blocs_non_vides(self.output):
                        ## gradient_calculator = vtk.vtkGradientFilter()
                        ## gradient_calculator.SetInputData(self.output.GetBlock(numbloc))
                        ## gradient_calculator.SetInputScalars(0, function[0].values()[0])
                        ## gradient_calculator.SetResultArrayName(function[2])
                        ## gradient_calculator.Update()
                        ## self.output.SetBlock(numbloc, gradient_calculator.GetOutput())
                ## else:
                    ## gradient_calculator = vtk.vtkGradientFilter()
                    ## gradient_calculator.SetInputData(self.output)
                    ## gradient_calculator.SetInputScalars(0, function[0].values()[0])
                    ## gradient_calculator.SetResultArrayName(function[2])
                    ## gradient_calculator.Update()
                    ## self.output = gradient_calculator.GetOutput()
            #else:
                #for numbloc in get_numeros_blocs_non_vides(self.output) \
                        #if isinstance(self.output, vtk.vtkMultiBlockDataSet) else [None]:
                    #current_bloc = self.output.GetBlock(numbloc) if numbloc != None else self.output
                    
                    #calc = vtk.vtkArrayCalculator()
                    #if self.use_cell_data:
                        #calc.SetAttributeModeToUseCellData()
                    #vtk_set_input(calc, current_bloc)
                    #calc.SetFunction(function[1])
                    #for var_input in function[0].items():
                        #if self.use_cell_data == True:
                            #if current_bloc.GetCellData().GetArray(var_input[1]).GetNumberOfComponents() == 3:
                                #calc.AddVectorVariable(var_input[0], var_input[1], 0, 1, 2)
                            #else:
                                #calc.AddScalarVariable(var_input[0], var_input[1], 0)
                        #elif self.use_cell_data == False:
                            #if current_bloc.GetPointData().GetArray(var_input[1]).GetNumberOfComponents() == 3:
                                #calc.AddVectorVariable(var_input[0], var_input[1], 0, 1, 2)
                            #else:
                                #calc.AddScalarVariable(var_input[0], var_input[1], 0)
                    #calc.AddCoordinateScalarVariable('coordx', 0)
                    #calc.AddCoordinateScalarVariable('coordy', 1)
                    #calc.AddCoordinateScalarVariable('coordz', 2)
                    #calc.SetResultArrayName(function[2])
                    #calc.ReplaceInvalidValuesOn()
                    #calc.SetReplacementValue(0.0)
                    #calc.Update()
                    #current_bloc.ShallowCopy(calc.GetOutput())
        
        #if self.keepIntermediateVariables == False:
            #cleanOutput = vtk_new_shallowcopy(self.output)
            #list_to_keep = list(self.a_calculer) + \
                #get_noms_arrays_presents(self.input, loc = 'points')

            #for quantity in get_noms_arrays_presents(cleanOutput, loc = 'points'):
                #if not(quantity in list_to_keep):
                    #if isinstance(cleanOutput, vtk.vtkMultiBlockDataSet):
                        #for numbloc in get_numeros_blocs_non_vides(cleanOutput):
                            #cleanOutput.GetBlock(numbloc).GetPointData().RemoveArray(quantity)
                    #else:
                        #cleanOutput.GetPointData().RemoveArray(quantity)
            #self.output = cleanOutput
        
        ##si un nom du resultat est donne, on change le nom de l'array au point
        #if self.nom_resultat != None:
            ##on convertit d'abord nom_result en une liste si ce n'en est pas une
            #if not isinstance(self.nom_resultat, list):
                #self.nom_resultat = [self.nom_resultat]
            
            ##on verifie de nom_resultat et a_calculer font les memes longueurs
            #if len(self.nom_resultat) != len(self.a_calculer):
                #raise IOError, "il n'y a pas le meme nombre de a_calculer et nom_resultat"

            ##cas multibloc
            #for k in range(len(self.a_calculer)):
                #avant = self.a_calculer[k]
                #apres = self.nom_resultat[k]
                #if isinstance(self.output, vtk.vtkMultiBlockDataSet):
                    #for numbloc in get_numeros_blocs_non_vides(self.output):
                        #self.output.GetBlock(numbloc).GetPointData().GetArray(avant).SetName(apres)
                #else:
                    #self.output.GetPointData().GetArray(avant).SetName(apres)
            
        #self._mettre_a_jour = False
        #return 0
    ##_____________________________________________________________________________________
##_____________________________________________________________________________________
