try: from paraview import vtk 
except: import vtk

try :
    from paraview import numpy_support
except:
    from vtk.util import numpy_support

try: 
    from paraview.vtk import vtkFiltersModeling
except:
    import vtk as vtkFiltersModeling

from fonctions_basiques import *
from extractions import *
import numpy

#__________________________________________________________________________________________
def ecrire_interface_file(acces_fichier, dimensions, indice_paraview):
    """ fonction pour ecrire les interface file a destination d'elsA
    en vue de faire du prelevement par conditions limites
    Les indices sont comme dans Paraview - i.e. l'indice du premier point est 0. 
    ecriture de fichiers au format fmt_v3d
    """
    dimensions = numpy.array(dimensions)
    narray = range(1, numpy.prod(dimensions - 1) + 1)
    narray_fente = []
    
    indice_paraview = numpy.asarray(indice_paraview)
    if len(indice_paraview.shape) == 1:
        for i in range(indice_paraview[0], indice_paraview[1]):
            for j in range(indice_paraview[2], indice_paraview[3]):
                narray_fente.append(float((i + 1) + j * (dimensions[0] - 1)))
    elif len(indice_paraview.shape) == 2:
        for fente in range(indice_paraview.shape[0]):
            indices_fente = indice_paraview[fente]
            for i in range(indices_fente[0], indices_fente[1]):
                for j in range(indices_fente[2], indices_fente[3]):
                    narray_fente.append(float((i + 1) + j * (dimensions[0] - 1)))
    else:
        raise IOError, 'Pas compriiis : indice_paraview'
    narray_fente = numpy.asarray(narray_fente)
    
    for num in narray_fente:
        narray.remove(num)
    narray = numpy.asarray(narray)
    
    ecrire_v3d(acces_fichier = 
        acces_fichier + '_fente.bnd',
        dict_numpy_arrays = {'interface_index': numpy.asarray(narray_fente)},
        numbloc = 1, 
        fmt_fichier = 'fmt', 
        dimensions = (narray_fente.size, 1, 1)
        )
    ecrire_v3d(acces_fichier = 
        acces_fichier + '_paroi.bnd',
        dict_numpy_arrays = {'interface_index': numpy.asarray(narray)},
        numbloc = 1, 
        fmt_fichier = 'fmt', 
        dimensions = (narray.size, 1, 1)
        )
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def lire_interface_file(input, acces_fichier, type_fichier = "v3d", \
        fmt_fichier= "fmt", endian= "big" , \
        precision = 'i4r8', nom_array_interface = "interface_index"):
    """ fonction de lecture des interface file 
    
    la surface doit etre indique en entree
    1 est attribue aux cellules indiquee dans le fichier interface file, 0 aux autres
    
    le fichier interface_file doit etre donne a format v3d
    """
    output = vtk_new_shallowcopy(input)
    
    if type_fichier == 'v3d':
        data = lire_v3d(acces_fichier = acces_fichier, fmt_fichier = fmt_fichier,
            endian = endian, precision = precision)
    elif type_fichier == 'tp':
        data = lire_fichier_tecplot(acces_fichier = acces_fichier)
    else:
        raise IOError, 'format de fichier non implemente'
    
    cell_interf = numpy.zeros(output.GetNumberOfCells())
    cell_interf[numpy.asarray(data['data'][nom_array_interface], dtype = int) - 1] = 1
    
    cell_interf = numpy_support.numpy_to_vtk(cell_interf, deep = 1)
    cell_interf.SetName(nom_array_interface)
    
    output.GetCellData().AddArray(cell_interf)
    return output
#__________________________________________________________________________________________

##____________________________________________________________________________
#def trouver_col(volume, nb_aubes, bloc_aubage=None, coupe="coordx=65", formule_extraction_surface_aubage="j=jmin", 
        #surface_aubage=None, 
        #axe=2):
    #"""POUR DIFFUSEUR RADIAL SEULEMENT - extension a programmer
    #fonction qui retourne un plan correspondant au col
    
    
    #nb_aubes est le nombre d'aubes correspondant 
    #a la grille consideree. Ce nombre est utilise pour determiner la position
    #de l'aube voisine.
    
    #coupe est la coupe a effectuer pour obtenir un profil 2d a partir du 3d
    #et determiner la ligne du col, support du plan au col
    
    #Si surface_aubage est donnee, sous la forme d'un objet VTK, alors elle est directement utilisee
    #Typiquement utilise pour les maillages non structure
    
    #axe indique l'axe de rotation : 0 pour x -- 1 pour y -- 2 pour z
    #"""
    #if surface_aubage is None:
        #surface = Extraction(input = bloc_aubage, 
            #formule_extraction = formule_extraction_surface_aubage,
            #calculer_vecteur_normal = 1).get_output()
    #else:
        #surface = surface_aubage
    
    ##on extrait une ligne sur le profil. La recherche du col est faite en 2D
    #ligne = Extraction(input = surface, formule_extraction = coupe).get_output()
    #coords = get_vtk_array_as_numpy_array(input = ligne, nom_array = 'coords')
    
    ##Calcul du rayon pour chacun des points sur le profil
    #if axe == 0:
        #coordr = numpy.sqrt(coords[:, 1] ** 2 + coords[:, 2] ** 2)
    #elif axe == 1:
        #coordr = numpy.sqrt(coords[:, 0] ** 2 + coords[:, 2] ** 2)
    #elif axe == 2:
        #coordr = numpy.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    #else:
        #raise IOError, "pyturbo pas comprendre axe"
    
    ##on cherche le bord d'attaque avec le minimum du rayon
    ##amarsan : a verifier pour turbine radiale
    #coords_ba = coords[numpy.argmin(coordr), :]
    
    ##calcul des coordonnees du bord d'attaque de l'aubage voisin
    #angle_rotation = 2 * numpy.pi / nb_aubes
    #if axe == 0:
        #coords_ba = [
            #coords_ba[0], 
            #coords_ba[1] * numpy.cos(angle_rotation) - coords_ba[2] * numpy.sin(angle_rotation), 
            #coords_ba[2] * numpy.cos(angle_rotation) + coords_ba[1] * numpy.sin(angle_rotation)
            #]
    #elif axe == 1:
        #coords_ba = [
            #coords_ba[0] * numpy.cos(angle_rotation) + coords_ba[2] * numpy.sin(angle_rotation),
            #coords_ba[1],
            #coords_ba[2] * numpy.cos(angle_rotation) - coords_ba[0] * numpy.sin(angle_rotation), 
            #]
    #else: #a ce moment la, l'axe est forcement 0, 1 ou 2. Plus d'erreur possible. 
        #coords_ba = [
            #coords_ba[0] * numpy.cos(angle_rotation) - coords_ba[1] * numpy.sin(angle_rotation), 
            #coords_ba[1] * numpy.cos(angle_rotation) + coords_ba[0] * numpy.sin(angle_rotation), 
            #coords_ba[2]
            #]
    
    ##calcul de la distance entre chacun des points de la ligne et le bord d'attaque de l'aubage voisin
    #dist_ba = numpy.sqrt(numpy.sum((coords - coords_ba) ** 2, axis = 1))
    
    ##on prend le minimum de cette distance
    #coords_col = coords[dist_ba.argmin(), :]
    
    ##calcul du vecteur normal au plan
    #vect_col = coords_ba - coords_col
    #normal = numpy.cross(
        #[axe == 0, axe == 1, axe == 2], vect_col) / numpy.linalg.norm(vect_col)
    
    ##generation du plan de coupe
    #plan = vtk.vtkPlane()
    #plan.SetOrigin(coords_ba)
    #plan.SetNormal(normal)
    
    ##coupe du volume duplique (pour avoir un canal complet)
    #coupe = vtk.vtkCutter()
    #coupe.SetCutFunction(plan)
    #volume_duplique = dupliquer_canal(volume, angle_rotation * 180. / numpy.pi, axe = axe)
    
    #plan = appliquer_sur_multibloc(coupe, volume_duplique)
    
    #plan = merge_multibloc(plan)
    
    ##calcul pour exclure les exterieurs du plan, et ne garder que le col
    ##on fait le produit scalaire avec le segment jusqu'au bord d'attaque et on prend la zone entre 0 et 1
    #coords = get_vtk_array_as_numpy_array(plan, 'coords')
    #data = numpy.dot(coords[:, 1:] - coords_col[1:], 
        #vect_col[1:]) / numpy.linalg.norm(vect_col) ** 2
    #data = numpy_support.numpy_to_vtk(data, deep = 1)
    #data.SetName("dist")
    #plan.GetPointData().AddArray(data)
    
    #plan = set_scalaires_actifs(plan, "dist")
    
    #select = vtk.vtkThreshold()
    #vtk_set_input(select, plan)
    #select.ThresholdBetween(1e-4, 1 - 1e-4)
    #select.Update()
    #col = select.GetOutput()
    #return col
##____________________________________________________________________________

#____________________________________________________________________________
def trouver_col(volume, surface_aubage, nb_aubes, 
        # coupe="coordz=0.005", 
        coupe="coordx=64.7500961", 
        axe=2):
    """Fonction qui retourne un plan correspondant au col
    Pour un profil 2D seulement. 
    Extension a programmer. 
    
    - volume est le volume du canal (un seul canal), sous la forme d'un objet VTK
    - surface_aubage est la surface de l'aubage dans ce canal, sous la forme d'un objet VTK
    - nb_aubes est le nombre d'aubes de la roue consideree. Il est utilise pour dupliquer le canal
    et determiner la position de l'aube voisine, par rotation autour de l'axe. 
    
    - coupe est la coupe a effectuer pour obtenir un profil 2d a partir du 3d
    et determiner la ligne du col, support du plan au col
    
    axe indique l'axe de rotation : 0 pour x -- 1 pour y -- 2 pour z
    """
    
    #on extrait une ligne sur le profil. La recherche du col est faite en 2D
    profil = Extraction(input = surface_aubage, formule_extraction = coupe).get_output()
    ligne = get_vtk_array_as_numpy_array(input = profil, nom_array = 'coords')
    
    #calcul des coordonnees de la ligne de l'aubage voisin, par rotation autour de l'axe
    angle_rotation = 360. / nb_aubes
    profil_voisin = rotation(profil, angle_rotation, axe=axe)
    ligne_voisine = get_vtk_array_as_numpy_array(input = profil_voisin, nom_array = 'coords')
    
    #calcul des distances point a point entre les deux lignes. 
    #En notant n le nombre de point de la ligne, on obtient alors une matrice (n, n)
    """
    ligne = [0 0 0 0 ...
             1 1 1 1 ...
             ...
            ]
    ligne_voisine = [0 1 2 ...
                     0 1 2 ...
                     0 1 2 ...
                     ...
                     ]
    """
    ligne = numpy.repeat(ligne, ligne.shape[0], axis=-1).reshape(ligne.shape + (ligne.shape[0],)).transpose(0, 2, 1)
    ligne_voisine = numpy.repeat(ligne_voisine, ligne_voisine.shape[0], axis=-1).reshape(
        ligne_voisine.shape + (ligne_voisine.shape[0],)).transpose(2, 0, 1)
    distances = numpy.sqrt(numpy.sum((ligne - ligne_voisine) ** 2, axis = -1))
    
    where_min = numpy.where(distances == numpy.min(distances))
    pt_ligne = ligne[where_min[0], 0].ravel()
    pt_ligne_voisine = ligne_voisine[0, where_min[1]].ravel()
    
    #calcul du vecteur normal au plan
    vect_col = pt_ligne_voisine - pt_ligne
    normal = numpy.cross(
        [axe == 0, axe == 1, axe == 2], vect_col) / numpy.linalg.norm(vect_col)
    
    #generation du plan de coupe
    plan = vtk.vtkPlane()
    plan.SetOrigin(pt_ligne)
    plan.SetNormal(normal)
    
    #coupe du volume duplique (pour avoir un canal complet)
    coupe = vtk.vtkCutter()
    coupe.SetCutFunction(plan)
    volume_duplique = dupliquer_canal(volume, 360. / nb_aubes, axe = axe)
    
    plan = appliquer_sur_multibloc(coupe, volume_duplique)
    
    plan = merge_multibloc(plan)
    
    #calcul pour exclure les exterieurs du plan, et ne garder que le col
    #on fait le produit scalaire avec le segment jusqu'au bord d'attaque et on prend la zone entre 0 et 1
    coords = get_vtk_array_as_numpy_array(plan, 'coords')
    data = numpy.dot(coords - pt_ligne, 
        vect_col) / numpy.linalg.norm(vect_col) ** 2
    data = numpy_support.numpy_to_vtk(data, deep = 1)
    data.SetName("dist")
    plan.GetPointData().AddArray(data)
    
    plan = set_scalaires_actifs(plan, "dist")
    
    select = vtk.vtkThreshold()
    vtk_set_input(select, plan)
    select.ThresholdBetween(0.0, 1.0)
    select.Update()
    col = select.GetOutput()
    col = convertir_en_polydata(col)
    return col
#____________________________________________________________________________


#____________________________________________________________________________
def dupliquer_canal(input, angle_periodicite, nb_canaux=2, axe=2):
    """
    duplique pour mettre plusieurs canaux cote-a-cote
    s'applique a un multibloc ou a un monobloc
    
    angle_periodicite en degres
    ### PEUT ETRE UNE LISTE pour un multibloc ###
    
    Axe designe l'axe de rotation. 
    0 pour x, 1 pour y, 2 pour z. 
    *****
    pour un multibloc
    possibilite d'indiquer des listes d'angle de periodicite et de nb_canaux
    dans ce cas, pour le bloc numero N, la valeur de l'angle utilisee est le N-ieme angle dans la liste
    idem pour nb_canaux
    """
    #cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        multibloc_duplique = vtk_new_shallowcopy(input)
        nb_blocs = multibloc_duplique.GetNumberOfBlocks()
        
        if isinstance(angle_periodicite, list):
            if not(isinstance(nb_canaux, list)):
                raise IOError, "lorsque une liste d'angles est donnee en entree, le nb_canaux doit aussi etre une liste"
            for numbloc in get_numeros_blocs_non_vides(multibloc_duplique):
                angle_rotation = angle_periodicite[numbloc]
                nb_reconstruire = nb_canaux[numbloc]
                for num_canal in range(2, nb_reconstruire + 1):
                    multibloc_duplique.SetBlock(numbloc + (num_canal - 1) * (nb_blocs // 10 + 1) * 10,
                        rotation(multibloc_duplique.GetBlock(numbloc), (num_canal - 1) * angle_rotation, axe)
                        )
        else:
            for numbloc in get_numeros_blocs_non_vides(multibloc_duplique):
                for num_canal in range(2, nb_canaux + 1):
                    multibloc_duplique.SetBlock(numbloc + (num_canal - 1) * (nb_blocs // 10 + 1) * 10,
                        rotation(multibloc_duplique.GetBlock(numbloc), (num_canal - 1) * angle_periodicite, axe)
                        )
    
    #monobloc
    else:
        if isinstance(angle_periodicite, list):
            raise IOError, "impossible de donner plusieurs angles de periodicite pour un monobloc... !"
        
        multibloc_duplique = vtk.vtkMultiBlockDataSet()
        multibloc_duplique.SetBlock(1, input)
        for num_canal in range(2, nb_canaux + 1):
            multibloc_duplique.SetBlock(num_canal, 
                rotation(input, (num_canal - 1) * angle_periodicite, axe))
    return multibloc_duplique
#____________________________________________________________________________

#____________________________________________________________________________
def dupliquer_canaux(input, angle_periodicite, num_canaux=[0, 1], axe=2):
    """
    duplique pour mettre plusieurs canaux cote-a-cote
    s'applique a un multibloc ou a un monobloc
    
    angle_periodicite en degres
    ### PEUT ETRE UNE LISTE pour un multibloc ###
    
    Axe designe l'axe de rotation. 
    0 pour x, 1 pour y, 2 pour z. 
    *****
    pour un multibloc
    possibilite d'indiquer des listes d'angle de periodicite et de nb_canaux
    dans ce cas, pour le bloc numero N, la valeur de l'angle utilisee est le N-ieme angle dans la liste
    et num_canaux doit etre une liste de listes contenant les numeros de canaux a reconstruire, pour chaque bloc
    
    ******* ATTENTION : les numeros des canaux doivent etre POSITIFS (ou nuls) ******* 
    
    """
    #cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        multibloc_duplique = vtk_new_instance(input)
        nb_blocs = input.GetNumberOfBlocks()
        
        if isinstance(angle_periodicite, list):
            if not(isinstance(num_canaux[0], list)):
                raise IOError, "lorsque une liste d'angles est donnee en entree, le num_canaux doit etre une liste de listes contenant les numeros de canaux a reconstruire, pour chaque bloc"
            for numbloc in get_numeros_blocs_non_vides(input):
                angle_rotation = angle_periodicite[numbloc]
                num_a_reconstruire = num_canaux[numbloc]
                for num_canal in num_a_reconstruire:
                    multibloc_duplique.SetBlock(numbloc + num_canal * (nb_blocs // 10 + 1) * 10,
                        rotation(input.GetBlock(numbloc), num_canal * angle_rotation, axe)
                        )
        else:
            for numbloc in get_numeros_blocs_non_vides(input):
                for num_canal in num_canaux:
                    multibloc_duplique.SetBlock(numbloc + num_canal * (nb_blocs // 10 + 1) * 10,
                        rotation(input.GetBlock(numbloc), num_canal * angle_periodicite, axe)
                        )
    
    #monobloc
    else:
        if isinstance(angle_periodicite, list):
            raise IOError, "impossible de donner plusieurs angles de periodicite pour un monobloc... !"
        
        multibloc_duplique = vtk.vtkMultiBlockDataSet()
        multibloc_duplique.SetBlock(1, input)
        for num_canal in num_canaux:
            multibloc_duplique.SetBlock(num_canal, 
                rotation(input, num_canal * angle_periodicite, axe))
    return multibloc_duplique
#____________________________________________________________________________

#____________________________________________________________________________
def hacher_diffuseur(volume, surface, coupe_2, coupe_1="coordx=65", retourner_normal_plan=False):
    """teste sur DIFFUSEUR RADIAL seulement
    fonction qui retourne un plan perpendiculaire a surface, generalement celle de l'aubage
    
    coupe_1 et coupe_2 sont les Extractions effectuees pour trouver le point sur surface
        qui servira de centre pour le plan
        typiquement coupe_1 = "coordx=65", coupe_2 = "coordr=158"
    
    la normale au plan est alors determinee par un produit scalaire entre 
        le vecteur iHat du repere et le vecteur normal a la surface en ce point
    
    surface est converti en polydata si besoin
    """
    surface = convertir_en_polydata(surface, calculer_vecteur_normal = True)
    
    ligne = Extraction(input = surface, formule_extraction = coupe_1).get_output()
    point = Extraction(input = ligne, formule_extraction = coupe_2,
        ).get_output()
    if point.GetNumberOfPoints() == 0:
        raise IOError, "impossible de trouver le point central du plan"
    normal = get_vtk_array_as_numpy_array(input = point, nom_array = 'Normals')[0, :]
    coords = get_vtk_array_as_numpy_array(input = point, nom_array = 'coords')[0, :]
    
    normal_plan = numpy.cross([1, 0, 0], normal)
    
    plan = vtk.vtkPlane()
    plan.SetOrigin(coords)
    plan.SetNormal(normal_plan)
    
    coupe = vtk.vtkCutter()
    coupe.SetCutFunction(plan)
    
    if not isinstance(volume, vtk.vtkMultiBlockDataSet):
        vtk_set_input(coupe, volume)
        coupe.Update()
        output = coupe.GetOutput()
    else:
        output = appliquer_sur_multibloc(coupe, volume)
    
    if retourner_normal_plan:
        return output, normal_plan
    else:
        return output
#____________________________________________________________________________

#__________________________________________________________________________
def moyenne_azimutale_hsH(vtkobject, liste_hsh, quantite, moyenne_debit=False, \
        nom_array_hsh = "hsH"):
    """fonction qui effectue une moyenne azimutale a differentes hauteurs de veine
    En entree:
        - vtkobject
            l'objet VTK auquel appliquer la fonction indifferemment mono-bloc ou multi-blocs
        - liste_hsh
            contient la lste des hauteurs auxquelles realiser la moyenne azimutale
        - quantite
            la quantite a moyenner
        - moyenne_debit
            False ou True en fonction de si la moyenne azimutale doit etre ponderee 
            par l'extension azimutale (r*theta) ou le debit
        - nom_array_hsh
            par defaut "hsH", mais peut etre change pour s'adapter a d'autres cas
    """
    #extraction du plan
    plan, surface_plan = calculer_surfaces_cellules(vtkobject, retourner_surface_totale = True)
    
    #calcul de la quantite a moyenner
    if moyenne_debit:
        calc = CalculettePyturbo(input = plan, 
                a_calculer = (quantite + "*abs(momentum.Normals)", "abs(momentum.Normals)", 
                ))
    else:
        calc = CalculettePyturbo(input = plan, 
                a_calculer = (quantite, "coordtheta*coordr")
                )
    
    #calcul de la moyenne pour chacune des hauteurs demandees
    valeurs_moyennes = []
    for hsH in liste_hsh:
        #extraction de la ligne a hauteur constante
        ligne = Extraction(input = calc.get_output(),
            formule_extraction = 'hsH={0}'.format(hsH)
                ).get_output()
        #MOYENNE SUR LES LIGNES
        if moyenne_debit:
            valeurs_moyennes.append(
                integrer_sur_la_surface(input = ligne, 
                    array_a_integrer = quantite + "*abs(momentum.Normals)")
                    / integrer_sur_la_surface(input = ligne, 
                    array_a_integrer = "abs(momentum.Normals)")
                )
        else:
            valeurs_moyennes.append(
                integrer_sur_la_surface(input = ligne, array_a_integrer = quantite, 
                    array_poids = "coordtheta*coordr")
                / integrer_sur_la_surface(input = ligne, array_a_integrer = 1, 
                    array_poids = "coordtheta*coordr")
                )
    return valeurs_moyennes
#__________________________________________________________________________

#__________________________________________________________________________
def moyenne_surfacique(vtkobject, quantite, moyenne_debit=False, array_momentum='momentum'):
    """fonction qui calcule la moyenne surfacique d'une grandeur
    En entree:
        - vtkobject
            l'objet VTK auquel appliquer la fonction indifferemment mono-bloc ou multi-blocs
        - quantite
            la quantite a moyenner
        - moyenne_debit
            False ou True en fonction de si la moyenne azimutale doit etre ponderee 
            par l'extension azimutale (r*theta) ou le debit
            Dans ce cas, il faut que momentum = ro*Velocity soit disponible aux noeuds !
    """
    #extraction du plan
    plan, surface_plan = calculer_surfaces_cellules(vtkobject, retourner_surface_totale = True)
    
    #calcul de la quantite a moyenner
    if moyenne_debit:
        plan = CalculettePyturbo(input = plan, 
                a_calculer = (quantite + "*abs({0}.Normals)".format(array_momentum), "abs({0}.Normals)".format(array_momentum))
                ).get_output()
    else:
        plan = CalculettePyturbo(input = plan, a_calculer = quantite).get_output()
    
    #calcul de la moyenne sur la surface
    if moyenne_debit:
        valeurs_moyenne = \
            integrer_sur_la_surface(input = plan, array_a_integrer = quantite + "*abs({0}.Normals)".format(array_momentum)) \
            / integrer_sur_la_surface(input = plan, array_a_integrer = "abs({0}.Normals)".format(array_momentum))
            
    else:
        valeurs_moyenne = \
            integrer_sur_la_surface(input = plan, array_a_integrer = quantite) \
            / surface_plan
            
    return valeurs_moyenne
#__________________________________________________________________________

#__________________________________________________________________________
def integration_surfacique(vtkobject, quantite):
    """fonction qui calcule l'integration surfacique d'une grandeur. Peut etre typiquement
    utilisee pour calculer le debit a travers une fonction 
        integration_surfacique(plan, "momentum.Normals")
    En entree:
        - vtkobject
            l'objet VTK auquel appliquer la fonction indifferemment mono-bloc ou multi-blocs
        - quantite
            la quantite a moyenner
    """
    #extraction du plan
    plan = calculer_surfaces_cellules(vtkobject)
    
    #calcul de la quantite a integrer
    plan = CalculettePyturbo(input = plan, a_calculer = quantite).get_output()
    
    #integration
    resultat_integration = integrer_sur_la_surface(input = plan, 
        array_a_integrer = quantite, array_poids = 'CellSurface')
        
    return resultat_integration
#__________________________________________________________________________


#____________________________________________________________________________
def moyenne_sur_epaisseur(paroi, data, liste_ep):
    """Fonction qui realise la moyenne des donnees presentes aux points
    sur plusieurs nappes decalees par rapport a la premiere le long du vecteur 'Normals'
    
    utilise fonctions_basiques.VTKProbe pour prober mono/multi-blocs
    si paroi est un vtkStructuredGrid, alors il est transforme en vtkPolyData en utilisant
    le filtre vtkGeometryFilter, avant d'entrer dans le vtkProbleFilter
    """
    paroi.GetPointData().SetActiveVectors('Normals')
    for nom_array in get_noms_arrays_presents(paroi, 'points'):
        if nom_array != 'Normals':
            paroi.GetPointData().RemoveArray(nom_array)
    for nom_array in get_noms_arrays_presents(paroi, 'cellules'):
        paroi.GetCellData().RemoveArray(nom_array)
    
    liste_noms_arrays = get_noms_arrays_presents(data, loc = 'points')
    dict_data = dict.fromkeys(liste_noms_arrays)
    
    f = vtk.vtkGeometryFilter()
    vtk_set_input(f, paroi)
    f.Update()
    paroi = f.GetOutput()
    
    for ep in liste_ep:
        print "epaisseur ", ep
        warp = vtk.vtkWarpVector()
        vtk_set_input(warp, paroi)
        warp.SetScaleFactor(ep)
        warp.Update()
        warp = warp.GetOutput()
       
        warp = VTKProbe(input = warp, source = data)
        for array in liste_noms_arrays:
            if dict_data[array] is None:
                dict_data[array] = get_vtk_array_as_numpy_array(warp, array, True).reshape(
                    warp.GetPointData().GetArray(array).GetNumberOfTuples(), 
                    warp.GetPointData().GetArray(array).GetNumberOfComponents())
            else:
                dict_data[array] += get_vtk_array_as_numpy_array(warp, array, True).reshape(
                    warp.GetPointData().GetArray(array).GetNumberOfTuples(), 
                    warp.GetPointData().GetArray(array).GetNumberOfComponents())
    for array in dict_data:
        dict_data[array] /= len(liste_ep)
        varray = numpy_support.numpy_to_vtk(dict_data[array], deep = 1)
        varray.SetName(array)
        paroi.GetPointData().AddArray(varray)
    
    return paroi
#____________________________________________________________________________

#____________________________________________________________________________
def calculer_champ_meridien(vtkDataObject, 
        maillage_regulier=None,
        retourner_maillage_regulier=False, 
        hubFileName = "/media/FreeAgent GoFlex Drive/DATA_PI4/hub",
        tipFileName = "/media/FreeAgent GoFlex Drive/DATA_PI4/shroud", 
        stepSize=2.0, relativeExtension=0.1, 
        numberOfPoints_xm = 75, numberOfPoints_h = 10, 
        dtheta = 2 * numpy.pi / (21 * 10)
        ):
    """Fonction qui retourne un plan, contenant le champ meridien
    Le champ est le resultat d'une moyenne azimutale ponderee simple
    """
    #Creation du maillage regulier
    if maillage_regulier is None:
        from UVParametrizationFilter import CreateSpline
        spline, numberOfPoints = CreateSpline(hubFileName, tipFileName, stepSize,
            relativeExtension=relativeExtension)
        
        coords = get_vtk_array_as_numpy_array(vtkDataObject, "coords")
        coordx = coords[:, 0]
        coordy = coords[:, 1]
        coordz = coords[:, 2]
        coordtheta = numpy.arctan2(coordz, coordy)
        
        min_theta = coordtheta.min()
        max_theta = coordtheta.max()
        
        ntheta = int((max_theta - min_theta) // dtheta)
        print "Dimensions du grid ", numberOfPoints_xm, numberOfPoints_h, ntheta
        
        points = vtk.vtkPoints()
        uvMin = -relativeExtension
        uvMax = 1.0 + relativeExtension
        for _v in numpy.linspace(uvMin, uvMax, numberOfPoints_xm):
            print _v
            for _u in numpy.linspace(uvMin, uvMax, numberOfPoints_h):
                for theta in numpy.linspace(min_theta, max_theta, ntheta):
                    pnt = spline(_u, _v)
                    x = pnt[0]
                    y = pnt[1] * numpy.cos(theta)
                    z = pnt[1] * numpy.sin(theta)
                    points.InsertNextPoint(x, y, z)
    
        grid = vtk.vtkStructuredGrid()
        grid.SetDimensions(ntheta, numberOfPoints_h, numberOfPoints_xm)
        grid.SetPoints(points)
        
        #on retourne le maillage_regulier si demande
        if retourner_maillage_regulier == 1:
            return grid
    else:
        grid = maillage_regulier
    
    
    
    #Probe et moyenne azimutale, pour finalement retourner le plan
    grid = VTKProbe(input = grid, source = vtkDataObject)
    
    #pour le plan a retourner, on choisit arbitrairement le plan i=0 du maillage regulier
    plan = Extraction(input = grid, formule_extraction='i=0').get_output()
    
    #moyenne de chacune des grandeurs le long des lignes a x, r, constants
    for nom_grandeur in get_noms_arrays_presents(vtkDataObject):
        grandeur = get_vtk_array_as_numpy_array(grid, nom_grandeur)
        if grandeur.size == grid.GetNumberOfPoints():
            grandeur = grandeur.reshape(grid.GetDimensions()[::-1])
        else:
            grandeur = grandeur.reshape(grid.GetDimensions()[::-1] + (3,))
        grandeur = numpy.ma.masked_less(grandeur, 0.1) #peut-etre a modifier amarsan
        grandeur = numpy.mean(grandeur, axis=2)
        grandeur = grandeur.ravel()
        plan = ajouter_numpy_array_as_vtk_array(plan, grandeur, nom_grandeur)
    
    return plan
#____________________________________________________________________________

#____________________________________________________________________________
def soustraction_vtk(vtk1, vtk2):
    """soustraction de deux objets vtk. vtk1 - vtk2
    Les geometries doivent etre les memes. 
    """
    #somme de tous les champs
    output = vtk_new_shallowcopy(vtk1)
    for var in get_noms_arrays_presents(output):
        voutput = get_vtk_array_as_numpy_array(output, var) - get_vtk_array_as_numpy_array(vtk2, var)
        output = ajouter_numpy_array_as_vtk_array(output, voutput, var)
    return output
#____________________________________________________________________________
#____________________________________________________________________________
def addition_vtk(vtk1, vtk2):
    """addition de deux objets vtk. vtk1 + vtk2
    Les geometries doivent etre les memes. 
    
    Utiliser la fonction somme_vtk pour ajouter plus que deux objets vtk. 
    """
    #somme de tous les champs
    output = vtk_new_shallowcopy(vtk1)
    for var in get_noms_arrays_presents(output):
        voutput = get_vtk_array_as_numpy_array(output, var) + get_vtk_array_as_numpy_array(vtk2, var)
        output = ajouter_numpy_array_as_vtk_array(output, voutput, var)
    return output
#____________________________________________________________________________

#____________________________________________________________________________
def somme_vtk(liste_vtk):
    """somme de plusieurs objets vtk. 
    Les geometries doivent etre les memes. 
    """
    #somme de tous les champs
    output = vtk_new_shallowcopy(liste_vtk[0])
    for data in liste_vtk[1:]:
        for var in get_noms_arrays_presents(output):
            voutput = get_vtk_array_as_numpy_array(output, var) + get_vtk_array_as_numpy_array(data, var)
            output = ajouter_numpy_array_as_vtk_array(output, voutput, var)
    return output
#____________________________________________________________________________
#____________________________________________________________________________
def moyenne_vtk(liste_vtk):
    """moyenne de plusieurs objets vtk. 
    Les geometries doivent etre les memes. 
    """
    #somme de tous les champs
    output = vtk_new_shallowcopy(liste_vtk[0])
    nb = 1.0
    for data in liste_vtk[1:]:
        nb += 1.
        for var in get_noms_arrays_presents(output):
            voutput = (get_vtk_array_as_numpy_array(output, var) * (nb - 1.) + 
                get_vtk_array_as_numpy_array(data, var)) / nb
            output = ajouter_numpy_array_as_vtk_array(output, voutput, var)
    return output
#____________________________________________________________________________

#____________________________________________________________________________
def extruder(input, longueur, nb_step=None):
    """Extrusion. Repose sur vtkLinearExtrusionFilter. Voir l'aide de VTK pour plus de details. 
    extrusion dans la direction du vecteur Normals seulement. 
    D'autres options sont disponibles dans le filtre VTK mais ne sont pas interfacees ici. Voir l'aide VTK. 
    
    Normals doit etre present. 
    
    input doit etre un vtkPolyData
    ou un multibloc compose de vtkPolyData
    """
    # cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, extruder(input.GetBlock(numbloc), longueur, nb_step))
        return output
    
    # cas monobloc
    if nb_step is None:
        f = vtkFiltersModeling.vtkLinearExtrusionFilter()
        vtk_set_input(f, input)
        f.SetExtrusionTypeToNormalExtrusion()
        f.SetScaleFactor(longueur)
        f.Update()
        output = f.GetOutput()
    else:
        f = vtk.vtkAppendPolyData()
        dl = numpy.linspace(0, longueur, nb_step)
        numbloc = -1
        for k in range(nb_step - 1):
            numbloc += 1
            data = decaler_paroi(input, dl[k])
            data = extruder(data, dl[k + 1] - dl[k], None)
            # m.SetBlock(numbloc, data)
            try:
                f.AddInputData(data)
            except:
                f.AddInput(data)
        # output = merge_multibloc(m)
        f.Update()
        output = f.GetOutput()
    return output
#____________________________________________________________________________
##____________________________________________________________________________
#def extrusion(input, longueur, discretisation):
    #"""Extrusion. Repose sur vtkLinearExtrusionFilter. Voir l'aide de VTK pour plus de details. 
    #extrusion dans la direction du vecteur Normals seulement. 
    #D'autres options sont disponibles dans le filtre VTK mais ne sont pas interfacees ici. Voir l'aide VTK. 
    
    #Normals doit etre present. 
    
    #input doit etre un vtkPolyData
    #ou un multibloc compose de vtkPolyData
    #"""
    ## cas multibloc
    #if isinstance(input, vtk.vtkMultiBlockDataSet):
        #output = vtk_new_instance(input)
        #for numbloc in get_numeros_blocs_non_vides(input):
            #output.SetBlock(numbloc, extrusion(input.GetBlock(numbloc), longueur, discretisation))
        #return output
    
    ## cas monobloc
    #coords_init = get_vtk_array_as_numpy_array(input, 'coords', 1)[None]
    #normals = get_vtk_array_as_numpy_array(input, 'Normals')[None]
    
    #coords = coords_init
    #for dx in numpy.linspace(0, longueur, discretisation)[1:]:
        #coords = numpy.r_[coords, coords_init + dx * normals]
    
    #coords = coords[None]
    #coords = coords.transpose(2, 1, 0, 3)
    #output = create_bloc_structure_from_numpy_array(coords)
    
    #return output
##____________________________________________________________________________
