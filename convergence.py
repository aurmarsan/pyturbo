from fonctions_basiques import *
from objets import ObjetPyturbo

class FrontiereChoro(ObjetPyturbo):
    """
    classe qui permet l'analyse des fichiers chorochroniques
    
    Les coefficients choro sont stockes en tant que self.coeffs_a et  self.coeffs_b
    shape = (nb_harmoniques, nb_var, dims_0, dims_1)
    dans l'ordre i, j, k intuitif (et pas l'ordre k, j, i de stockage VTK)
    
    
    Pour initialisation, donner :
        - maillage
        - liste_acces_fichiers, 
        - liste_numblocs
        - liste_frontieres, 
        - liste_directions=[],
        - dict_fenetres_blocs = {}
    
    
    si une liste de dossiers est donnees, alors les fichiers choro indiques par liste_acces_fichiers
    sont successivement lus dans tous les dossiers indiques
    en retour, coeffs_a et coeffs_b sont des arrays dont la permiere dimensions est egale aux nombre
    de dossiers indiques.
    shape = (n_calculs, n_harmoniques, n_grandeurs, dimensions_frontieres)
    """
    
    #_____________________________________________________________________________________
    def __init__(self, 
            maillage, 
            liste_acces_fichiers, 
            liste_numblocs, 
            liste_frontieres, 
            liste_directions=[],
            dict_fenetres_blocs={},
            liste_dossiers = None, 
            ):
        #Initialisation de la classe parente
        attributs = locals().copy()
        del attributs['self']
        ObjetPyturbo.__init__(self, **attributs)
        
        #Coefficients de Fourier
        coeffs_a = None
        coeffs_b = None
        
        #lecture des fichiers chorochroniques si le dictionnaire de description des frontieres a ete renseigne
        self.update()
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________    
    def update(self):
        """fonction de chargement des donnees depuis les fichiers chorochroniques"""
        if self.liste_dossiers is None:
            coeffs_a, coeffs_b = self.__lire_fichiers_choro__(
                maillage = self.maillage,
                liste_acces_fichiers = self.liste_acces_fichiers,
                liste_numblocs = self.liste_numblocs,
                liste_frontieres = self.liste_frontieres,
                liste_directions = self.liste_directions,
                dict_fenetres_blocs = self.dict_fenetres_blocs
                )
            self.coeffs_a = coeffs_a
            self.coeffs_b = coeffs_b
        else:
            coeffs_a_tot = None
            coeffs_b_tot = None
            for doss in self.liste_dossiers:
                liste_acces_fichiers = [doss + "/" + acces_fichier for 
                    acces_fichier in self.liste_acces_fichiers]
                coeffs_a, coeffs_b = self.__lire_fichiers_choro__(
                    maillage = self.maillage,
                    liste_acces_fichiers = liste_acces_fichiers,
                    liste_numblocs = self.liste_numblocs,
                    liste_frontieres = self.liste_frontieres,
                    liste_directions = self.liste_directions,
                    dict_fenetres_blocs = self.dict_fenetres_blocs
                    )
                if coeffs_a_tot is None:
                    coeffs_a_tot = coeffs_a.reshape((1,) + coeffs_a.shape)
                    coeffs_b_tot = coeffs_b.reshape((1,) + coeffs_b.shape)
                else:
                    coeffs_a_tot = numpy.concatenate(
                        (coeffs_a_tot,
                        coeffs_a.reshape((1,) + coeffs_a.shape)
                        ),
                        axis = 0)
                    coeffs_b_tot = numpy.concatenate(
                        (coeffs_b_tot,
                        coeffs_b.reshape((1,) + coeffs_b.shape)
                        ),
                        axis = 0)
            self.coeffs_a = coeffs_a_tot
            self.coeffs_b = coeffs_b_tot
        return 0
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def __lire_fichiers_choro__(self, 
            maillage, 
            liste_acces_fichiers, 
            liste_numblocs, 
            liste_frontieres, 
            liste_directions=[],
            dict_fenetres_blocs = {}):
        """Cette fonction lit les fichiers binaires v3d issus des calculs chorochroniques
        d'elsA. 
            - maillage doit etre le maillage
            - liste_acces_fichiers donne les chemins vers les fichiers chrochroniques a lire
            - liste_numblocs donne la liste des blocs concernes par les fichiers chorochroniques, 
                dans le meme ordre
            - liste_frontieres indique quelles sont les frontieres des blocs concernees
                les cles peuvent etre "imin", "imax", "jmin", "jmax", "kmin", "kmax"
            - liste_directions est une liste de taille N-1 ou N est la taille des listes precedentes
                elle indique comment coller les frontieres des blocs pour reconstruire la frontiere complete
                en fait, c'est directement l'argument axe de numpy.append
                    donc si direction = 0, les coeff_choro sont ajoutes en ligne, un espece de i de la frontiere
                         si direction = 1, les coeff_choro sont ajoutes en colonne, un espece de j de la frontiere
            - dict_fenetres_blocs
                sert dans le cas d'une frontiere qui ne repose pas sur l'ensemble d'une fenetre jmin (par exemple)
                les indices doivent etre donnes en mode elsA, i.e. commencent a 1
                en fait, il est possible de recopier directement la definition de la fenetre depuis la carte de lancement elsA
        
        En sortie, cette fonction retourne 
            - la matrice des a_n en chacun des points de la frontiere 
                dimension N x M x N_harmoniques, ou N et M sont les nombres de points de la frontiere dans les deux directions 
            - la matrice des b_n en chacun des points de la frontiere 
                dimension N x M x N_harmoniques, ou N et M sont les nombres de points de la frontiere dans les deux directions 
        
        Pour l'appel de lire_v3d, les parametres par defaut sont conserves
        """
        
        #Verification des donnes d'entree
        if len(liste_acces_fichiers) != len(liste_numblocs):
            raise IOError, "liste_acces_fichiers et liste_numblocs doivent faire la meme longueur"
        if len(liste_acces_fichiers) != len(liste_frontieres):
            raise IOError, "liste_acces_fichiers et liste_frontieres doivent faire la meme longueur"
        if len(liste_acces_fichiers) != (len(liste_directions) + 1) and liste_directions != []:
            raise IOError, "len(liste_directions) doit etre egal a (liste_acces_fichiers) - 1"
        
        
        #Lecture des fichiers chorochroniques, et mises sous les bonnes dimensions
        #i.e. nb_variables, nk, nj, ni.
        coeffs_choro = []
        for k in range(len(liste_acces_fichiers)):
            fich = liste_acces_fichiers[k]
            
            #cas ou le fichier couvre plusieurs blocs (stage_choro)
            if isinstance(liste_numblocs[k], list):
                #lecture de toutes les dimensions
                liste_dims = None
                for numbloc in liste_numblocs[k][:-1]:
                    if numbloc not in dict_fenetres_blocs:
                        dims_bloc = list(maillage.GetBlock(numbloc).GetDimensions())
                    else:
                        dims_fenetre = dict_fenetres_blocs[numbloc]
                        dims_bloc = [
                            dims_fenetre[1] - dims_fenetre[0] + 1, 
                            dims_fenetre[3] - dims_fenetre[2] + 1, 
                            dims_fenetre[5] - dims_fenetre[4] + 1
                            ]
                    if liste_dims is None:
                        liste_dims = dims_bloc
                    else:
                        liste_dims = numpy.c_[liste_dims, dims_bloc]
                liste_dims = liste_dims.transpose()
                
                #somme des dimensions
                if liste_numblocs[k][-1][0] is 'i':
                    dims = [numpy.sum(liste_dims[:, 0] - 1) + 1, liste_dims[0, 1], liste_dims[0, 2]] 
                elif liste_numblocs[k][-1][0] is 'j':
                    dims = [liste_dims[0, 0], numpy.sum(liste_dims[:, 1] - 1) + 1, liste_dims[0, 2]] 
                elif liste_numblocs[k][-1][0] is 'k':
                    dims = [liste_dims[0, 0], liste_dims[0, 1], numpy.sum(liste_dims[:, 2] - 1) + 1] 
            #cas simple ou le fichier n'est que sur un bloc
            elif liste_numblocs[k] not in dict_fenetres_blocs:
                dims = maillage.GetBlock(liste_numblocs[k]).GetDimensions()
            else:
                dims_fenetre = dict_fenetres_blocs[liste_numblocs[k]]
                dims = [
                    dims_fenetre[1] - dims_fenetre[0] + 1, 
                    dims_fenetre[3] - dims_fenetre[2] + 1, 
                    dims_fenetre[5] - dims_fenetre[4] + 1
                    ]
            frontiere = liste_frontieres[k]
            #lecture du fichier V3D
            coeffs_choro.append(lire_v3d(acces_fichier = fich)['data'])
            
            #reformatage - sachant que les coeff choro sont aux centres des faces
            if frontiere[0] is "i":
                dims_frontiere = (dims[2] - 1, dims[1] - 1)
            elif frontiere[0] is "j":
                dims_frontiere = (dims[2] - 1, dims[0] - 1)
            elif frontiere[0] is "k":
                dims_frontiere = (dims[1] - 1, dims[0] - 1)
            else:
                raise IOError, "frontiere {0} pas comprise".format(frontiere)
            for harmonique in coeffs_choro[-1].keys():
                coeffs_choro[-1][harmonique] = coeffs_choro[-1][harmonique].reshape(
                    (7, ) + dims_frontiere).transpose(0, 2, 1)
                    
        #assemblage des coefficients choro complet pour recreer la frontiere totale
        coeffs_choro_global = coeffs_choro[0]
        for numbloc in range(1, len(coeffs_choro)):
            choro_bloc = coeffs_choro[numbloc]
            axe = liste_directions[numbloc -1]
            
            for key in coeffs_choro_global.keys():
                coeffs_choro_global[key] = numpy.append(
                    coeffs_choro_global[key],
                    choro_bloc[key],
                    axis = 1 + axe
                    )
        
        rangs = []
        for k in coeffs_choro_global.keys():
            if k[0] is "a":
                rangs.append(int(k[1:]))
        rangs.sort()
        
        coeffs_a = numpy.array(
            [coeffs_choro_global["a{0}".format(rang)] for rang in rangs]
            )
        coeffs_b = numpy.array(
            [coeffs_choro_global["b{0}".format(rang)] for rang in rangs]
            )
        return coeffs_a, coeffs_b
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_amplitudes(self):
        """
        fonction de calcul des amplitudes de chacune des harmoniques, pour chacune des variables
        amplitudes = numpy.sqrt(self.coeffs_a ** 2 + self.coeffs_b ** 2)
        """
        amplitudes = numpy.sqrt(self.coeffs_a ** 2 + self.coeffs_b ** 2)
        return amplitudes
    #_____________________________________________________________________________________
    
    #_____________________________________________________________________________________
    def get_phases(self):
        """
        fonction de calcul des phases de chacune des harmoniques, pour chacune des variables
        amplitudes = numpy.arctan(self.coeffs_b / self.coeffs_a)
        """
        amplitudes = numpy.arctan(self.coeffs_b / self.coeffs_a)
        return amplitudes
    #_____________________________________________________________________________________
    
    
