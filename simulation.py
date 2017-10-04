try: from paraview import vtk 
except: import vtk
try: from paraview import numpy_support 
except: from vtk.util import numpy_support
import numpy
import glob
import copy
from fonctions_basiques import *
import lecture
from objets import ObjetPyturbo, RefAero, Maillage

#__________________________________________________________________________________________
class Simulation(ObjetPyturbo):
    """
    une simulation standard
    indiquer les grandeurs de reference pour le redimensionnement (type TM)
    
    omega_par_blocs doit etre indique en radians par secondes
    
    axe designe l'axe de rotation: 0 pour x, 1 pour y et 2 pour z. 
    """
    #_____________________________________________________________________________________
    def __init__(self, nom=None, RefAero=RefAero(), 
            Maillage = Maillage(), 
            omega_par_blocs=[29210. * numpy.pi/30.] * 18 + [0.0] * 7, 
            timestep=0.0, itime=None, 
            inititer=None, _vtkDataObject=None, 
            axe = 0):
        attributs = locals().copy()
        del attributs['self']
        #initialisation de la classe parente
        ObjetPyturbo.__init__(self, **attributs)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_vtkDataObject(self):
        return self.get('_vtkDataObject')
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def set_vtkDataObject(self, vtkDataObject):
        self.set('_vtkDataObject', vtkDataObject)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def lire_fichier_log(self, acces_fichier_log, lecture_donnees_maillage=True):
        """lecture des donnees contenues dans un fichier log
        
        utilisation de RefAero pour 
            - redimensionner la vitesse de rotation
            - redimensionner itime
            - redimensionner le timestep
        
        verifier avant de lire les donnees du log que RefAero
        sont correctement renseignes
        
        on suppose que les blocs sont definies par numeros croissants dans le log
        le nombre d'aubes est definis comme axis_ang[0] / axis_ang[1]
        """
        
        def grep(key, string):
            try:
                # on recherche la cle entouree de guillemets
                grep_i = string.split("'" + key + "'")[1:]
                grep = []
                for k in grep_i:
                    grep += k.split(key + '"')
                #on prend ce qu'il y a entre la virgule et la parenthese
                for k in range(len(grep)):
                    grep[k] = grep[k][grep[k].find(',') + 1:].split(')')[0].strip()
                if len(grep) == 1:
                    grep = grep[0]
            except:
                print 'echec de la recherche'
                grep = None
            return grep
        
        log = file(acces_fichier_log, 'r')
        log = log.read()
        
        # recherche de inititer
        if not "('iterations'" in log:
            self.set('inititer', int(grep('inititer', log)))
        else:
            self.set('inititer', int(
                grep('iterations', log).split('[')[1].split(',')[0]
                ))
        
        # recherche de omega_adim
        omega_adim = numpy.asarray(grep('omega', log), dtype = float)
        # redimensionnement omega_adim
        vitesse_ref = numpy.sqrt(
            self.RefAero.get('gamma_ref') * self.RefAero.get('r_gaz_ref') * self.RefAero.get('t_ref')
            )
        omega_dim = omega_adim * vitesse_ref / self.RefAero.get('l_ref')
        self.set('omega_par_blocs', omega_dim)
        
        #nombre d'aubes par bloc du maillage
        if lecture_donnees_maillage:
            nombre_aubes_par_bloc = grep('axis_ang', log)
            for k in range(len(nombre_aubes_par_bloc)):
                exec "nombre_aubes_par_bloc[{0}] = {1}".format(k, nombre_aubes_par_bloc[k]) in locals()
                nombre_aubes_par_bloc[k] = int(nombre_aubes_par_bloc[k][0]) / int(nombre_aubes_par_bloc[k][1])
            self.Maillage.set("nombre_aubes_par_bloc", nombre_aubes_par_bloc)
        
        # timestep
        if "'unsteady'" not in grep('time_algo', log):
            self.set('timestep', 0)
        else:
            # chorochronique : deux vitesses de rotation differentes seulement
            liste_omega = []
            for k in self.get('omega_par_blocs'):
                if k not in liste_omega: liste_omega.append(k)
            if len(liste_omega) != 2:
                raise Exception, 'obligatoirement deux vitesses de rotation differentes'
            # chorochronique : deux nombres d'aubes differents seulement
            liste_nombre_aubes = []
            for k in self.Maillage.get('nombre_aubes_par_bloc'):
                if k not in liste_nombre_aubes: liste_nombre_aubes.append(k)
            if len(liste_nombre_aubes) != 2:
                raise Exception, "obligatoirement deux nombres d'aubes differents"
            # calcul du timestep
            print grep('timestep', log)
            timestep_adim = float(grep('timestep', log))
            timestep = timestep_adim * self.RefAero.get('l_ref') / vitesse_ref
            self.set('timestep', timestep)
            ## calcul de la periode
            #t_tour = 2 * numpy.pi / abs(liste_omega[0] - liste_omega[1])
            ## calcul du nqo
            #nqo = t_tour / (numpy.prod(liste_nombre_aubes) * timestep)
            #if abs(nqo - int(round(nqo))) > 1e-2:
                #raise Exception, "le nqo calcule n'est pas un entier : {0}".format(nqo)
            #self.set('nqo', int(round(nqo)))
        
            # recherche de itime_adim
            itime_adim = float(grep('itime', log))
            # redimensionnement itime_adim
            itime_dim = itime_adim * self.RefAero.get('l_ref') / vitesse_ref
            self.set('itime', itime_dim)
        return 0
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def lire_solution(self, acces_fichier, fmt_fichier='bin_v3d', precision='i4r8', \
            endian='big', decalage_numbloc_lus=0, \
            iteration=None, canaux_a_reconstruire=['tout', 'tout'], \
            rotation_temporelle=True, 
            traitement_non_numeros=False, seulement_les_numeros=None,
            assembler_vecteurs = True):
        """lecture de fichiers contenant des donnees
        
        fmt_fichier peut etre
            - vtk           !!! PAS ENCORE POUR CHOROCHRONIQUE !!!
            - bin_v3d
            - fmt_v3d
        
        par defaut, tous les canaux sont reconstruits
        
        voir la classe LecteurV3D pour l'utilite des cles 
            - traitement_non_numeros
            - seulement_les_numeros
        
        rotation_temporelle peut prendre les valeurs suivantes
            - True          --> rotation normale, avec prise en compte de la vitesse
                                de rotation de chacun des blocs
            - False         --> pas de prise en compte des vitesse de rotation
                                les domaines sont fixes quelle que soit le numero d'iteration
            - un entier     --> l'entier indique le numero de la roue qui doit rester fixe
                                dans le repere absolu. L'autre roue tourne en consequence
                                pour maintenir la coherence a l'interface des deux roues
                    0 pour la premiere roue
                    1 pour la deuxieme roue
        
        """
        self._ajouter_omega_aux_noeuds()
        
        if self.get('timestep') == 0:
            if fmt_fichier[4:] == 'v3d':
                r = lecture.LecteurV3D(acces_fichier = acces_fichier, fmt_fichier = fmt_fichier[:3], 
                    precision = precision, endian = endian, 
                    decalage_numbloc_lus = decalage_numbloc_lus, 
                    traitement_non_numeros = traitement_non_numeros,
                    seulement_les_numeros = seulement_les_numeros,
                    assembler_vecteurs = assembler_vecteurs)
                r.importer_maillage(self.Maillage.get('_vtkDataObject'))
                r.update()
                self.set('_vtkDataObject', r.get_output())
            elif fmt_fichier == 'vtk':
                r = lecture.LecteurVTK(acces_fichier = acces_fichier)
                self.set('_vtkDataObject', r.get_output())
            else:
                raise IOError, 'le format {0} indique pour le fichier est inconnu'.format(fmt_fichier)
        
        else:
            # VERIFICATION DES DONNEES
            if iteration == None:
                raise IOError, "indiquer l'iteration a reconstruire"
            # chorochronique : deux vitesses de rotation differentes seulement
            liste_omega = []
            for k in self.get('omega_par_blocs'):
                if k not in liste_omega: liste_omega.append(k)
            if len(liste_omega) != 2:
                raise IOError, 'obligatoirement deux vitesses de rotation differentes'
            # chorochronique : deux nombres d'aubes differents seulement
            liste_nombre_aubes = []
            for k in self.Maillage.get('nombre_aubes_par_bloc'):
                if k not in liste_nombre_aubes: liste_nombre_aubes.append(k)
            if len(liste_nombre_aubes) != 2:
                raise IOError, "obligatoirement deux nombres d'aubes differents"
            nombre_aubes_par_bloc = self.Maillage.get("nombre_aubes_par_bloc")
            # inititer, itime et timestep doivent etre renseignes
            for propriete in ['inititer', 'itime', 'timestep', 'omega_par_blocs']:
                if self.get(propriete) == None:
                    raise IOError, "{0} doit etre renseigne".format(propriete)
                else:
                    exec "{0} = self.get('{0}')".format(propriete)
            omega_par_blocs = self.get("omega_par_blocs")
            iterations_disponibles = []
            for fich in glob.glob('{0}*'.format(acces_fichier)):
                iteration_trouvee = int(fich.split(acces_fichier.split('*')[1])[1])
                if not(iteration_trouvee in iterations_disponibles):
                    iterations_disponibles.append(iteration_trouvee)
            
            self._vtkDataObject = vtk.vtkMultiBlockDataSet()
            # RECONSTRUCTION
            periode_tour = 2 * numpy.pi / abs(liste_omega[0] - liste_omega[1])
            for numbloc in get_numeros_blocs_non_vides(
                    self.Maillage._vtkDataObject):
                if seulement_les_numeros is None or numbloc in seulement_les_numeros:
                    numbloc = numbloc
                    roue_courante = liste_omega.index(omega_par_blocs[numbloc - 1])
                    roue_opposee = roue_courante - 1
                    numeros_canaux = canaux_a_reconstruire[roue_courante]
                    if numeros_canaux is None:
                        numeros_canaux = []
                    if not isinstance(numeros_canaux, list):
                        numeros_canaux = range(nombre_aubes_par_bloc[numbloc - 1])
                    # calcul de l'angle de la rotation du au temps en degres
                    if rotation_temporelle == True and type(rotation_temporelle) == bool:
                        rotation_temps = ((iteration - inititer) * timestep + itime) \
                            * liste_omega[roue_courante] * 180. / numpy.pi
                    elif type(rotation_temporelle) == int:
                        if rotation_temporelle not in [0, 1]:
                            raise IOError, "rotation_temporelle ne peut valoir que True; False; 0; 1;"
                        roue_fixe = rotation_temporelle
                        rotation_temps = ((iteration - inititer) * timestep + itime) \
                            * (liste_omega[roue_courante] - liste_omega[roue_fixe]) \
                            * 180. / numpy.pi
                    else:
                        rotation_temps = 0.0
                    
                    for numero_canal in numeros_canaux:
                        # calcul de l'angle de la rotation du au numero de canal en degres
                        rotation_numero_canal = numero_canal * 360. / nombre_aubes_par_bloc[numbloc - 1]
                        # numero iteration a chercher
                        dephasage_en_iterations = periode_tour / (
                            nombre_aubes_par_bloc[numbloc - 1] * timestep
                            ) * numpy.sign(
                            liste_omega[roue_courante] - liste_omega[roue_opposee]
                            )
                        iteration_canal = iteration + dephasage_en_iterations * numero_canal
                        # recherche parmi les iterations disponibles
                        # compte-tenu de la periodicite temporelle
                        periode_canal_en_iterations = periode_tour / (
                            liste_nombre_aubes[roue_opposee] * timestep
                            )
                        iteration_a_lire = iteration_canal - (
                            (iteration_canal - numpy.max(iterations_disponibles)) // periode_canal_en_iterations
                            + 1) * periode_canal_en_iterations
                        iteration_a_lire = int(round(iteration_a_lire))
                        if not iteration_a_lire in iterations_disponibles:
                            raise Exception, 'no iteration founded for canal {0}\na candidate would be {1}'\
                                .format(numero_canal, iteration_a_lire)
                        
                        if fmt_fichier[4:] == 'v3d':
                            r = lecture.LecteurV3D(fmt_fichier = fmt_fichier[:3], 
                                precision = precision, endian = endian, 
                                decalage_numbloc_lus = decalage_numbloc_lus, 
                                traitement_non_numeros = traitement_non_numeros,
                                #seulement_les_numeros = seulement_les_numeros,
                                assembler_vecteurs = assembler_vecteurs)
                            r.importer_maillage(self.Maillage._vtkDataObject.GetBlock(numbloc))
                            r.set('acces_fichier', '{0}{1}{2}{3}'.format(
                                acces_fichier.split('*')[0], 
                                numbloc,
                                acces_fichier.split('*')[1],
                                iteration_a_lire
                                ))
                            r.update()
                        
                        numero_bloc_canal = (numbloc + numero_canal * (
                            self.Maillage._vtkDataObject.GetNumberOfBlocks() // 10 + 1
                            )* 10) % (
                            (self.Maillage._vtkDataObject.GetNumberOfBlocks() // 10 + 1) * 10
                            * (nombre_aubes_par_bloc[numbloc - 1] + 1)
                            )
                        self._vtkDataObject.SetBlock(numero_bloc_canal, 
                            rotation(
                                r.get_output(), rotation_temps + rotation_numero_canal, self.axe)
                            )
        return 0
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def _ajouter_omega_aux_noeuds(self):
        """ajoute omega aux noeuds du maillage- necessaire aux calculs de grandeurs
        (a transformer en FieldData quand le pipe vtk le passera correctement)
        
        cette fonction est appelee en meme temps que lire_solution
        de maniere a ce que omega soit disponible aux noeuds
        
        """
        if self.get('omega_par_blocs') is None:
            raise IOError, 'omega_par_blocs doit etre defini avant'
            
        if self.get_maillage_vtkDataObject() is None:
            print 'Maillage._vtkDataObject doit etre defini avant'
        elif isinstance(self.get_maillage_vtkDataObject(), vtk.vtkMultiBlockDataSet):
            for numero_bloc_mai in get_numeros_blocs_non_vides(
                    self.get_maillage_vtkDataObject()):
                bloc = self.get_maillage_vtkDataObject().GetBlock(numero_bloc_mai)
                narray = numpy.ones(bloc.GetNumberOfPoints(), dtype = float) \
                    * self.omega_par_blocs[numero_bloc_mai - 1]
                varray = numpy_support.numpy_to_vtk(narray, deep = 1)
                varray.SetName('omega')
                bloc.GetPointData().AddArray(varray)
        else:
            bloc = self.get_maillage_vtkDataObject()
            narray = numpy.ones(bloc.GetNumberOfPoints(), dtype = float) * self.omega_par_blocs
            varray = numpy_support.numpy_to_vtk(narray, deep = 1)
            varray.SetName('omega')
            bloc.GetPointData().AddArray(varray)
        
        #partie qui ajoute omega aux noeuds du vtkDataObject
        #utile quand on lit directement un vtm pour charger la solution
        #mais ne fonctionne pas lorsque la lecture est instationnaire ! et qu'on lit pour la deuxieme fois
        #if self.get_vtkDataObject() is None:
            #pass
        #elif isinstance(self.get_vtkDataObject(), vtk.vtkMultiBlockDataSet):
            #for numero_bloc_mai in get_numeros_blocs_non_vides(
                    #self.get_vtkDataObject()):
                #bloc = self.get_vtkDataObject().GetBlock(numero_bloc_mai)
                #narray = numpy.ones(bloc.GetNumberOfPoints(), dtype = float) \
                    #* self.omega_par_blocs[numero_bloc_mai - 1]
                #varray = numpy_support.numpy_to_vtk(narray, deep = 1)
                #varray.SetName('omega')
                #bloc.GetPointData().AddArray(varray)
        #else:
            #bloc = self.get_vtkDataObject()
            #narray = numpy.ones(bloc.GetNumberOfPoints(), dtype = float) * self.omega_par_blocs
            #varray = numpy_support.numpy_to_vtk(narray, deep = 1)
            #varray.SetName('omega')
            #bloc.GetPointData().AddArray(varray)
        
        return 0
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def redimensionner(self, \
            DensityArrayName='ro', MomentumArrayName='momentum', \
            TotalEnergyPerUnitOfVolumeName='roe'):
        """redimensionnement des grandeurs en fonction des grandeurs 
        de reference de la Simulation.RefAero
        Ne redimensionne pas les coordonnees des points du maillage
        (qui restent donc normalement en mm
        
        """
        def MultiplyVTKArray(vtkArray, factor):
                numarray = numpy_support.vtk_to_numpy(vtkArray)
                numarray *= factor
                del numarray
                #vtkArray.DataChanged()
                #vtkArray.Modified()
        
        _vtkDataObject = self.get('_vtkDataObject')
        for numbloc in get_numeros_blocs_non_vides(_vtkDataObject):
            for data in [_vtkDataObject.GetBlock(numbloc).GetPointData(), 
                    _vtkDataObject.GetBlock(numbloc).GetCellData()]:
                #redimensionnement de ro
                if data.HasArray(DensityArrayName):
                    MultiplyVTKArray(
                        data.GetArray(DensityArrayName)
                        , self.RefAero.p_ref / (self.RefAero.r_gaz_ref * self.RefAero.t_ref)
                        )
                if data.HasArray(MomentumArrayName):
                    #redimensionnement de momentum
                    MultiplyVTKArray(
                        data.GetArray(MomentumArrayName)
                        , self.RefAero.p_ref / (self.RefAero.r_gaz_ref * self.RefAero.t_ref) \
                            * numpy.sqrt(self.RefAero.gamma_ref * self.RefAero.r_gaz_ref * self.RefAero.t_ref)
                        )
                if data.HasArray(TotalEnergyPerUnitOfVolumeName):
                    #redimensionnement de l'energie totale par unite de volume
                    MultiplyVTKArray(
                        data.GetArray(TotalEnergyPerUnitOfVolumeName)
                        , self.RefAero.gamma_ref * self.RefAero.p_ref
                        )
        return 0
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def set_maillage_vtkDataObject(self, vtkDataObject):
        """methode pour renseigner le vtkDataObject du maillage lie a la simulation
        
        c'est en fait un simple appel de self.Maillage.set_vtkDataObject
        mais cela facilite grandement l'ecriture de script utilisateur
        
        """
        self.Maillage.set_vtkDataObject(vtkDataObject)
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_maillage_vtkDataObject(self):
        """methode pour renseigner le vtkDataObject du maillage lie a la simulation
        
        c'est en fait un simple appel de self.Maillage.set_vtkDataObject
        mais cela facilite grandement l'ecriture de script utilisateur
        
        """
        return self.Maillage.get_vtkDataObject()
    #_____________________________________________________________________________________
#__________________________________________________________________________________________
