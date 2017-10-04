import numpy, scipy
import pyturbo as p
from matplotlib import pyplot as plt

#####################################################
# LECTURE FNC
#####################################################
path_fnc = '/storage/amarsan/MIGMIG/Meas_Cyl_Midpsan.fnc'

r = p.LecteurFNC(file_path = path_fnc)
r.lire_parametres()
r.get_infos_temporelles()
r.get_infos_variables()

# lecture du maillage
# r.lire_maillage()
# p.EcritureVTK(input = r.maillage, acces_fichier = path_fnc[:-4]).ecrire()
r.importer_maillage(p.LecteurVTK(path_fnc[:-3] + 'vtm').get_output())


#####################################################
# LECTURE FNC DATAS
#####################################################
# lecture des donnees
ind_temps = 112
r.lire_datas(ind_temps, noms_vars = None, rotation_rotor = True, axe = 1, numbloc=None)

# recuperation de l'objet vtk en sortie du reader
volume = r.get_output()



#####################################################################
# CALCUL du rayon aux centres des cellules et ajout dans l'objet VTK
#####################################################################
p.get_noms_arrays_presents(volume, loc='c')

coords_centers = p.get_vtk_array_as_numpy_array(volume, 'centers')
pressure = p.get_vtk_array_as_numpy_array(volume, 'p')

# calcul du rayon
radius = numpy.linalg.norm(coords_centers[:, [0,2]], axis=-1)
# ajoutdu rayon dans l'objet VTK
volume = p.ajouter_numpy_array_as_vtk_array(volume, radius, 'rayon')

# sauvegarder si desire
# p.EcritureVTK(input = volume, acces_fichier = 'sauvegarde').ecrire()

#####################################################################
# CALCUL du rayon aux points et ajout dans l'objet VTK
#####################################################################
coords_points = p.get_vtk_array_as_numpy_array(volume, 'coords')
radius = numpy.linalg.norm(coords_points[:, [0,2]], axis=-1)
volume = p.ajouter_numpy_array_as_vtk_array(volume, radius, 'rayon')


#####################################################################
# PROBE avec une ligne a rayon constant
# La ligne est cree manuellement
#####################################################################
coordr = 0.0673
coordy = -0.132389520735594
from vtk.util import numpy_support
# generation du polydata pour le probe
theta = numpy.deg2rad(numpy.arange(-180, 180, 0.5))
coords = numpy.c_[
    coordr * numpy.sin(theta),
    numpy.ones(theta.size) * coordy, 
    coordr * numpy.cos(theta)
    ]
segments = numpy.concatenate(
    (2 * numpy.ones(theta.size)[:, None,], numpy.arange(theta.size)[:, None,], numpy.roll(numpy.arange(theta.size), -1)[:, None,]), 
    axis = 1).ravel()
ligne_probe = p.create_polydata_from_numpy_array(coords, polys = segments, nb_polys = theta.size, isline=True)


# Probe avec le polydata
ligne_meas = p.VTKProbe(ligne_probe, volume, tolerance = 0.1)

# Trace
p.get_noms_arrays_presents(ligne_meas, loc='points')
pressure = p.get_vtk_array_as_numpy_array(ligne_meas, 'p')

plt.figure()
plt.plot(coordr * theta, pressure, '-x')
plt.xlabel('r-theta [m]')
plt.ylabel('pressure [Pa]')
plt.title('Mig Mig')
plt.grid(True)
# plt.xticks()
plt.show()

# import vtk
# from vtk.util import numpy_support
# filtre = vtk.vtkProbeFilter()
# filtre.SetInputData(ligne_probe)
# filtre.SetSourceData(volume.GetBlock(1))
# filtre.ComputeToleranceOff()
# filtre.SetTolerance(0.1)
# filtre.Update()
# result = filtre.GetOutput()

# plt.figure()
# plt.plot(p.get_vtk_array_as_numpy_array(result, 'vtkValidPointMask'))
# plt.show()


# volume_points = p.interpoler_cellules_aux_points(volume)
# filtre = vtk.vtkProbeFilter()
# filtre.SetInputData(ligne_probe)
# filtre.SetSourceData(volume_points.GetBlock(1))
# filtre.ComputeToleranceOff()
# filtre.SetTolerance(0.1)
# filtre.Update()
# result_2 = filtre.GetOutput()

#####################################################################
# Coupes a y constant puis a rayon constant
# Pour avoir la meme ligne que precedememnt, mais avec des coupes
#####################################################################
plan = p.Extraction(input = volume, formule_extraction = 'coordy={0}'.format(coordy)).get_output()
ligne = p.Extraction(input = plan, formule_extraction = 'coordr={0}'.format(coordr), axe=1).get_output()



# trace 
pressure = p.get_vtk_array_as_numpy_array(ligne, 'p')
coords = p.get_vtk_array_as_numpy_array(ligne, 'centers')
theta = numpy.arctan2(coords[:, 0], coords[:, 2])

argsort = numpy.argsort(theta)
theta = theta[argsort]
pressure = pressure[argsort]

plt.figure(1)
plt.plot(theta, pressure, '-x', label = 'ligne coupes')

# trace 
pressure = p.get_vtk_array_as_numpy_array(ligne_meas, 'p')
coords = p.get_vtk_array_as_numpy_array(ligne_meas, 'coords')
theta = numpy.arctan2(coords[:, 0], coords[:, 2])

argsort = numpy.argsort(theta)
theta = theta[argsort]
pressure = pressure[argsort]

plt.figure(1)
plt.plot(theta, pressure, '-x', label = 'ligne a la main')

plt.xlabel('theta [rad]')
plt.ylabel('pressure [Pa]')
plt.title('Mig Mig')
plt.grid(True)
# plt.xticks()
plt.legend(loc = 'upper left')
plt.show()

