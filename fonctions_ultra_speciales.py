#ce fichier regroupe des fonctions absolument pas pensees pour etre generales
#mais utilisees au cours du post-doc pour exploiter les calculs CFX
try: from paraview import vtk 
except: import vtk

import numpy

from fonctions_basiques import *
from calculs import *
from ecriture import *

#__________________________________________________________________________________________
def lire_fichier_out_cfx(acces_fichier):
    """Fonction extremement speciale. 
    Lecture d'un fichier .out CFX et va ensuite rajouter omega aux blocs
    identitifes par leurs noms de domaine
    
    """
    #Lecture des omega par noms de domaine dans le fichier .out
    f = file(acces_fichier)
    info = f.read()
    info = info[:info.rfind('END')]
    
    #lecture des parametres
    parametres = info[info.find('EXPRESSIONS:'):]
    parametres = parametres[:parametres.find('END')]
    parametres = parametres.split('\n')[1:]
    for param in parametres:
        if '[' in param and ']' in param:
            param = param[:param.find('[')]
            param = param.strip()
            exec param
    
    #creation du dictionnaire des informations, pour chacun des domaines et chacunes des frontieres
    domaines = info.split('DOMAIN: ')[1:]
    for k in domaines:
        k = k[:k.rfind('END')]
    
    dict_domaines = {}
    for dom in domaines:
        nom = dom[:dom.find('\n')]
    
    return omega
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def reconstruire_serie_fourier(input, variable, time, periode, dephasage=0, suppression_coeff_fourier=True, liste_rangs=None):
    """fonction qui reconstruit une serie de Fourier
    Les coefficient doivent etre stockes au points. 
    Le nom de la variable doit contenir <variable>
    Le coefficient est determine par la fin du nom de la variable
        A0, A1, ..., A15, ... , B0, ..., B17, ... 
        
    <input> est un objet VTK, qui peut etre un multibloc
    
    la phase est calculee avec 
        phase = k * 2 * numpy.pi * (time + dephasage) / periode
    ou k est le range de l'harmonique qui est en train d'etre calculee
    """
    ####################################
    # cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, reconstruire_serie_fourier(input.GetBlock(numbloc), variable, time, periode, dephasage, suppression_coeff_fourier, liste_rangs))
    
    ####################################
    # cas monobloc - on etablit la liste des coefficients de Fourier presents
    else:
        # on se base sur le coefficient A0
        base_variable_name = None
        for var in get_noms_arrays_presents(input, 'points'):
            if var[-2:] == 'A0' and variable in var:
                base_variable_name = var[:-2]
                break
        if base_variable_name is None:
            raise Exception, "n'a pas reussi a determiner le nom generique de la variable..."
        
        # on cherche tous les coefficients disponibles
        liste_coeff = []
        for var in get_noms_arrays_presents(input, 'points'):
            if base_variable_name in var:
                liste_coeff.append(var)
        
        # on ecrit la formule
        formule = None
        for coeff in liste_coeff:
            rang = coeff.split(base_variable_name)[1]
            if rang[0] == 'A':
                op = 'cos'
            elif rang[0] == 'B':
                op = 'sin'
            else:
                raise Exception, 'Pas reussi a determiner a quoi correspond le coeff {0}'.format(coeff)
            
            try:
                rang = int(rang[1:])
            except:
                raise Exception, 'Pas reussi a determiner le rang du coefficient {0}'.format(coeff)
                
            phase = rang * 2 * numpy.pi * (time + dephasage) / periode
            
            if liste_rangs is None or rang in liste_rangs:
                if formule is None:
                    formule = '{0}*{1}({2})'.format(coeff, op, phase)
                else:
                    formule += '+{0}*{1}({2})'.format(coeff, op, phase)
        
        print 'reconstruction de la serie de Fourier \n {0} ...'.format(formule[:150])
        
        # si les coefficients sont de types scalaires, alors c'est un scalaire que l'on calcule
        if input.GetPointData().GetArray(liste_coeff[0]).GetNumberOfComponents() == 1:
            output = CalculetteGenerique(input = input, 
                formule = formule, 
                nom_du_resultat = variable,
                variables_scalaires = liste_coeff
                ).get_output()
        # sinon c'est un vecteur
        elif input.GetPointData().GetArray(liste_coeff[0]).GetNumberOfComponents() == 3:
            output = CalculetteGenerique(input = input, 
                formule = formule, 
                nom_du_resultat = variable,
                variables_vectorielles = liste_coeff
                ).get_output()
        else:
            raise Exception, 'pas compris si la variable etait un scalaire ou un vecteur' 
        
        if suppression_coeff_fourier is True:
            print 'suppression des arrays {0} ...'.format(liste_coeff[:5])
            for coeff_fourier in liste_coeff:
                output.GetPointData().RemoveArray(coeff_fourier)
        
    return output
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def applatir_surface(input, axe=2, etheta_min=None, etheta_max=None):
    """Fonction qui ramene une surface dans le plan (iHat, jHat) en supprimant la composante tangentielle. 
    Traite les vecteurs correctement de facon a conserver la topologie parietale.
    
    Indiquer l'axe de rotation: 
    - 0 pour iHat
    - 1 pour jHat
    - 2 pour kHat
    
    si etheta_min etheta_max sont indiques (ou seulement l'un des deux) alors 
    threshold sur normals.etheta pour selectionner l'une des faces de l'aubage
    """
    # cas multibloc - recursif
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, applatir_surface(input.GetBlock(numbloc), axe, etheta_min, etheta_max))
    
        return output
    
    
    # cas monobloc
    # calcul des vecteurs tangents a la surface qui forment un repere
    surface = convertir_en_polydata(input)
    surface = calculer_vecteur_normal(surface)
    normals = get_vtk_array_as_numpy_array(surface, 'Normals')
    surface = CalculettePyturbo(surface, ['er', 'coordr'], axe = axe).get_output()
    er = get_vtk_array_as_numpy_array(surface, 'er')
    v1 = numpy.cross(normals, er)
    
    v1 = v1 / numpy.sqrt(numpy.sum(v1**2, axis = 1))[:, None]
    v2 = numpy.cross(normals, v1)
    
    surface = ajouter_numpy_array_as_vtk_array(surface, v1, 'v1')
    surface = ajouter_numpy_array_as_vtk_array(surface, v2, 'v2')
    
    if etheta_min is not None or etheta_max is not None:
        surface = CalculettePyturbo(surface, 'etheta.Normals').get_output()
        surface = VTKThreshold(surface, 'etheta.Normals', etheta_min, etheta_max)
        surface = convertir_en_polydata(surface)
    
    # changement de repere
    coords = get_vtk_array_as_numpy_array(surface, 'coords')
    coordr = get_vtk_array_as_numpy_array(surface, 'coordr')
    coords = numpy.concatenate(
        (coords[:, axe], coordr, numpy.zeros(coordr.size)), 
        axis = None).transpose()
    
    vtkArray = numpy_support.numpy_to_vtk(numpy.ascontiguousarray(coords), deep = 1)
    points = vtk.vtkPoints()
    points.SetData(vtkArray)
    
    # initialisation de output
    output = vtk_new_shallowcopy(surface)
    output.SetPoints(points)
    # on doit recalculer le vecteur normal pour le rendu graphique
    output = calculer_vecteur_normal(output)
    
    # on traite les vecteurs aux points
    for array in get_noms_arrays_presents(output):
        if output.GetPointData().GetArray(array).GetNumberOfComponents() == 3:
            output = CalculettePyturbo(output, 
                a_calculer = "({0}.v1)*((v1.{1})*iHat+(v1.er)*jHat)+({0}.v2)*((v2.{1})*iHat+(v2.er)*jHat)"
                    .format(array, ['iHat', 'jHat', 'kHat'][axe]),
                nom_resultat = array).get_output()
    output = calculer_vecteur_normal(output)
    
    return output
#__________________________________________________________________________________________

