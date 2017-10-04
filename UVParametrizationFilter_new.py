try: from paraview import vtk 
except: import vtk
import numpy
try: from paraview import numpy_support 
except: from vtk.util import numpy_support

#######################################################
############ TENTATIVE DE SPHERE INSCRITES ############
#######################################################

#__________________________________________________________________________________________
def trouver_zero_fonction_croissante(valeur_1, valeur_2, fonction_a_annuler, precision=1e-6):
    """recherche le zero de la fonction dans l'intervalle [valeur_1, valeur_2]
    fonction est supposee etre CROISSANTE dans l'intervalle [valeur_1, valeur_2]
    
    """
    milieu = (valeur_1 + valeur_2) / 2.0
    if numpy.abs(valeur_1 - valeur_2) < precision:
        return milieu
    
    f_milieu = fonction_a_annuler(milieu)
    if f_milieu < 0:
        return trouver_zero_fonction_croissante(milieu, valeur_2, fonction_a_annuler, precision)
    else:
        return trouver_zero_fonction_croissante(valeur_1, milieu, fonction_a_annuler, precision)
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def trouver_centres_spheres_inscrites(surface_1, surface_2, 
        precision=1e-6, pas_initial=10):
    """fonction qui recherche les centres des spheres inscrites entre les surface_1 et surface_2
    
    le vecteur 'Normals' doit etre present aux noeuds de surface_1
    
    la recherche se fait par trouver_zero_fonction_croissante, pour chacun des points de surface_1
    on cherche un point sur la ligne orthogonale a surface_1, equidistant de surface_1 et surface_2
        --> centre de la sphere inscrite entre surface_1 et surface_2
            passant par le point de surface_1 considere
    
    surface_1 et surface_2 doivent etre des vtkPolyData
    """
    
    # initialisation des outils vtk de calcul des distances
    calcul_distance_1 = vtk.vtkKdTreePointLocator()
    calcul_distance_1.SetDataSet(surface_1)
    calcul_distance_1.BuildLocator()
    calcul_distance_2 = vtk.vtkKdTreePointLocator()
    calcul_distance_2.SetDataSet(surface_2)
    calcul_distance_2.BuildLocator()
    vtkMath = vtk.vtkMath()
    
    # lecture des donnees de la surface
    liste_normals_surface_1 = numpy_support.vtk_to_numpy(surface_1.GetPointData().GetArray('Normals'))
    liste_points_surface_1 = numpy_support.vtk_to_numpy(surface_1.GetPoints().GetData())
    # ACTION
    liste_centres = numpy.empty((0,3))
    
    print "progression ",
    printed = False
    for numero_point_surface_1 in range(0, surface_1.GetNumberOfPoints()):
        progress = ((numero_point_surface_1 + 1) * 100) // surface_1.GetNumberOfPoints()
        if progress % 10 == 0 and printed is False:
            print "{0}%".format(progress), 
            printed = True
        elif progress % 10 != 0:
            printed = False
        """le vecteur normal et les coordonnees de point definissent une ligne
        sur laquelle chercher le centre de la sphere"""
        normal = liste_normals_surface_1[numero_point_surface_1]
        point = liste_points_surface_1[numero_point_surface_1]
        
        def fonction_a_annuler(x):
            centre = point + x * normal
            p1 = surface_1.GetPoint(calcul_distance_1.FindClosestPoint(centre))
            p2 = surface_2.GetPoint(calcul_distance_2.FindClosestPoint(centre))
            dist_1 = vtkMath.Distance2BetweenPoints(p1, centre)
            dist_2 = vtkMath.Distance2BetweenPoints(p2, centre)
            nv_ecart = dist_1 ** (1. / 2) - dist_2 ** (1. / 2)
            return nv_ecart
        
        """ recherche du centre de la sphere inscrite par trouver_zero_fonction_croissante"""
        # il est possible que le pas initial indique soit trop petit
        pas = pas_initial
        while fonction_a_annuler(pas) < 0:
            pas *= 2.
        # recherche par dichotomie
        x_solution = trouver_zero_fonction_croissante(0.0, pas, fonction_a_annuler, precision)
        centre = point + x_solution * normal
        liste_centres = numpy.r_[liste_centres, centre.reshape((1, 3))]
    print "termine"
    
    # creation d'un object vtkPolyData
    vtkPoints_mean = vtk.vtkPoints()
    vtkPoints_mean.SetData(numpy_support.numpy_to_vtk(liste_centres, deep = 1))
    
    surface_moyenne = vtk.vtkPolyData()
    surface_moyenne.SetPoints(vtkPoints_mean)
    surface_moyenne.SetPolys(surface_1.GetPolys())
    surface_moyenne.Update()

    return surface_moyenne
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def convertir_surface_en_polydata(surface, nettoyage_tolerance=5e-2):
    """fonction de conversion d'une surface en polydata
    avec fusion des blocs si surface est un multiblock
    """
    # si surface est un multiblock, il faut les regrouper en un seul polydata
    if isinstance(surface, vtk.vtkMultiBlockDataSet):
        appendFilter = vtk.vtkAppendPolyData()
        appendFilter.MergePointsOn()
        iter = surface.NewIterator()
        iter.UnRegister(None)
        iter.InitTraversal()
        while not iter.IsDoneWithTraversal():
            curInput = iter.GetCurrentDataObject()
            if not isinstance(curInput, vtk.vtkPolyData):
                conversion = vtk.vtkGeometryFilter()
                conversion.SetInputData(curInput)
                conversion.Update()
                appendFilter.AddInput(conversion.GetOutput())
            else:
                appendFilter.AddInput(curInput)
            iter.GoToNextItem()
        appendFilter.Update()
        surface = appendFilter.GetOutput()
    
    # si surface n'est pas un multiblock, il faut simplemement le convertir en vtkPolyData
    elif not isinstance(surface, vtk.vtkPolyData):
        conversion = vtk.vtkGeometryFilter()
        conversion.SetInputData(surface)
        conversion.Update()
        surface = conversion.GetOutput()
    
    #LA PARTIE SUIVANTE N'EST PAS NECESSAIRE LORSQUE L'ON UTILISE MERGEPOINTSON DE L'APPENDFILTER
    # nettoyage des points confondus (aux frontieres entre les blocs)
    #nettoyage = vtk.vtkCleanPolyData()
    #nettoyage.SetInputData(surface)
    #nettoyage.Update()
    #surface = nettoyage.GetOutput()
    
    return surface
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def extension_maillage_structure(input, direction, longueur = 10):
    """fonction d'extrapolation d'un maillage structure
    
    input doit etre un vtkStructuredGrid
    direction est un entier:
        - 0     extension dans la direction 'i'
        - 1     extension dans la direction 'j'
        - 2     extension dans la direction 'k'
    """
    
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def calculer_surface_moyenne(surface_1, surface_2, precision=1e-6, pas_initial=10):
    """fonction qui calcule une surface moyenne par la methode des spheres inscrites
    extension de la methode des cercles inscrits
    """
    # conversion en polydata et fusion des blocs
    surface_1 = convertir_surface_en_polydata(surface_1)
    surface_2 = convertir_surface_en_polydata(surface_2)
    
    #calcul du vecteur Normals
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(surface_1)
    normals.Update()
    surface_1_modif = normals.GetOutput()
    
    # derafinement de la surface_1
    #derafiner = vtk.vtkCleanPolyData()
    #derafiner.SetInputData(surface_1_modif)
    #derafiner.SetToleranceIsAbsolute(1)
    #derafiner.SetAbsoluteTolerance(1)
    #derafiner.Update()
    #surface_1_modif = derafiner.GetOutput()
    
    # calcul de la surface moyenne
    surface_moyenne = trouver_centres_spheres_inscrites(surface_1_modif, surface_2, 
        precision = precision, pas_initial = pas_initial)
    
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(surface_moyenne)
    normals.ConsistencyOn()
    normals.Update()
    surface_moyenne = normals.GetOutput()
    
    coords = numpy_support.vtk_to_numpy(surface_moyenne.GetPoints().GetData())
    normals = numpy_support.vtk_to_numpy(surface_moyenne.GetPointData().GetArray('Normals'))
    
    inf = coords - normals * 1e3
    sup = coords + normals * 1e3
    
    source = vtk.vtkLineSource()
    source.SetPoint1(inf[0])
    source.SetPoint2(sup[0])
    # return source.GetOutput()
    
    return surface_moyenne
#__________________________________________________________________________________________
