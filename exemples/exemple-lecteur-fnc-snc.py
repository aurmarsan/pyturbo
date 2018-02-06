import numpy, scipy
import pyturbo as p


######################################
# pour lire un fichier surface SNC
path_snc = "indiquer.snc"
r = p.LecteurSNC(file_path = path_snc)
r.lire_parametres()
r.get_infos_temporelles() # ca affiche des trucs
r.get_infos_variables() # ca affiche d'autres trucs

# on peut lire le maillage (ca prend du temps parce qu'il faut reconstruire l'arbre de connectivite entre les voxels)
r.lire_maillage()
# dans ce cas, c'est une bonne idee de l'enregistrer au format VTK en utilisant le writer VTK
# l'extension .vtm .vtp .vtk est ajoutee automatiquement
p.EcritureVTK(r.get_maillage(), 'lendroit-ou-enregistrer').ecrire()

# ainsi, si on a deja lu le maillage, on peut simplement l'importer et ca va plus vite
# penser a inclure l'extension .vtm .vtp .vtk dans le nom de fichier
r.importer_maillage(p.LecteurVTK('lendroit-ou-enregistrer.vtm').get_output())

# enfin on peut lire les donnees qu'on veut
# le reader separe le domaine en deux blocs: 0 le rotor, 1 le stator. On peut aussi laisse numbloc=None pour lire rotor et stator. 
# ici on lit 'p' a la frame temporelle 10. On tourne le rotor pour le placer en position. On ajouter les faces_id qui peuvent
# etre bien pratique pour afficher seulement des morceaux de la surface (hub, blade, etc...)
# on ajoute les normales mais pas les surfaces des cellules. 
r.lire_datas(11, ['p'], rotation_rotor=True, numbloc=0, ajouter_faces_id=True, ajouter_normales=True, ajouter_surfaces=False)
# on recupere la sortie du reader
data = r.get_output()

######################################
# pour lire un fichier FNC, c'est pas mal pareil. 
path_fnc = "indiquer.fnc"
r = p.LecteurFNC(file_path = path_fnc)
r.lire_parametres()
r.get_infos_temporelles()
r.get_infos_variables()

# ici on va dire qu'on a deja lu et sauvegarde le maillage, alors on ne fait que l'importer.
r.importer_maillage(p.LecteurVTK('tu-connais-le-chemin').get_output())

# puis on lit les donnees de vitesse, on place le rotor en position. 
# on lit le stator et le rotor (numbloc = None)
r.lire_datas(11, ['vx', 'vy', 'vz'], 
    numbloc = None, 
    rotation_rotor = True)

volume = r.get_output()



######################################
# ici un exemple un peu particulier. 
# on veut lire une surface, et la tourner du meme angle que le domaine du rotor (volume) est tourne a un certain instant. 
# C'est utile pour la visualisation quand les temps physique d'extraction de la surface et du volume ne sont pas les memes
r = p.LecteurSNC(file_path = path_snc)
r.lire_parametres()
r.get_infos_temporelles()
r.get_infos_variables()
r.importer_maillage(p.LecteurVTK('{0}/VTK-maillages/{1}.vtm'
    .format(path_snc[:path_snc.rfind('/')], path_snc[path_snc.rfind('/') + 1 : ])
    ).get_output())
# on lit seulement la surface rotor, sans tourner le rotor. On recupere ainsi le maillage a la position t=0. 
r.lire_datas(last_frame_index, [], rotation_rotor=False, numbloc=0, ajouter_faces_id=False, ajouter_normales=False, ajouter_surfaces=False)
surface = r.get_output()

# on va chercher dans le fichier FNC les informations qui nous permettent de calculer de combien il faut tourner la surface
# on ne lit pas les donnees ! Seulement les informations. 
rvol = p.LecteurFNC(file_path = path_fnc)
rvol.lire_parametres()
rvol.get_infos_temporelles()
rvol.get_infos_variables()

# Puis on tourne la surface de l'angle qui va bien. Et voila. 
surface = p.rotation(surface, numpy.rad2deg(
    rvol.parameters['omega'] * rvol.parameters['sign_rotation'] * rvol.parameters['current_time'][ind] + rvol.parameters['init_angle']
    ))

