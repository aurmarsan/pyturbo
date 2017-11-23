try: from paraview import vtk 
except: import vtk
try: from paraview import numpy_support 
except: from vtk.util import numpy_support

#partie pour s'adapter aux version 5. et 6. de VTK
try: 
    from paraview.vtk import vtkFiltersExtraction
    from paraview.vtk import vtkFiltersVerdict
    from paraview.vtk import vtkFiltersGeneral
    from paraview.vtk import vtkCommonTransforms
    from paraview.vtk import vtkFiltersGeometry
except:
    import vtk as vtkFiltersExtraction
    import vtk as vtkFiltersVerdict
    import vtk as vtkFiltersGeneral
    import vtk as vtkCommonTransforms
    import vtk as vtkFiltersGeometry

from scipy import interpolate
import struct
import numpy
import numpy.ma as ma
import copy
import os, glob, sys


#________________________________________________________________________________
def vtk_set_input(filtre, input):
    """fonction de merde simplement pour la compatibilite
    a supprimer plus tard...
    """
    try:
        filtre.SetInputData(input)
    except:
        filtre.SetInput(input)
    return 0
#________________________________________________________________________________

#________________________________________________________________________________
def lire_v3d(acces_fichier, fmt_fichier= "bin", endian= "big" , \
        precision = 'i4r8', compter_saut_de_ligne=False):
    """fonction de lecture d'un fichier v3d
    
    Retourne un dictionnaire contenant
        - 'dims' : dimensions du bloc v3d lu
        - 'data' : dictionnaire contenant les donnees lues
            les numpy arrays ne sont pas mis en forme
    
    si utilisation de numpy.reshape, attention a l'ordre d'empilement 
            ce doit etre .reshape(dim_k, dim_j, dim_i)
    
    compter_saut_de_ligne ne devrait normalement pas etre utilise... Mais peut
        toujours servir si jamais ca plante sur un cas
    """
    
    #verification de l'existence du fichier v3d
    #et ouverture
    try:
        fichier_v3d = open(acces_fichier, 'rU') if fmt_fichier=='fmt'\
                    else open(acces_fichier, 'rb')
    except IOError:
        raise IOError, "le fichier n'existe pas"
        return 1
    else:
        print 'Lecture %s'%(acces_fichier)
    
    #dictionnaire de lecture binaire
    dictionnaire_fmt_binaire = {"i4":[4,"l"], "i8":[8,"q"], 
                                "r4":[4,"f"], "r8":[8,"d"], 
                                "big":">", "little":"<"}

    data = {}
    if fmt_fichier == "bin":
        #options de lecture pour struct.unpack
        ordre_bin = '>' if endian == 'big' else '<'
        fmt_entier = "l" if precision[:2] == 'i4' else 'q'
        taille_entier = int(precision[1])
        fmt_reel = "f" if precision[-2:] == 'i4' else 'd'
        taille_reel = int(precision[3])
        
        a_lire = struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0]
        nb_vars = struct.unpack(ordre_bin + fmt_entier, fichier_v3d.read(a_lire))[0]
        if struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0] != a_lire: 
            raise Exception, 'erreur au cours de la lecture'
        
        for incr_var in xrange(nb_vars):
            a_lire = struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0]
            nom_var = struct.unpack('20s', fichier_v3d.read(a_lire))[0]
            var=nom_var.replace('va','').replace(' ','')
            if struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0] != a_lire: 
                raise Exception, 'erreur au cours de la lecture'
            
            a_lire = struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0]
            if nom_var[:2] == 'va': 
                info_dims = struct.unpack(ordre_bin + '5' + fmt_entier, 
                    fichier_v3d.read(a_lire))
            else: 
                info_dims = struct.unpack(ordre_bin + '4' + fmt_entier, 
                    fichier_v3d.read(a_lire))
            if struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0] != a_lire: 
                raise Exception, 'erreur au cours de la lecture'
            dims = numpy.asarray(info_dims[-3:]).astype("Int32")
            numbloc = info_dims[1 if nom_var[:2] == 'va' else 0]
            
            a_lire = struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0]
            vars = struct.unpack(ordre_bin + str(int(a_lire / taille_reel)) + fmt_reel, 
                fichier_v3d.read(a_lire))
            vars = numpy.asarray(vars)
            if struct.unpack(ordre_bin + 'l', fichier_v3d.read(4))[0] != a_lire: 
                raise Exception, 'erreur au cours de la lecture'
            
            data[var] = vars
    
    elif fmt_fichier == "fmt":
        nb_vars = int(fichier_v3d.readline().strip())
        for var_iter in xrange(nb_vars):
            tampon = fichier_v3d.readline()
            nom_var = tampon[:20]
            length_format = int(tampon[20:40].split('e')[1].split('.')[0])
            read = fichier_v3d.readline()[:-1]
            dims=[]
            for i in range(len(read) / 6):
                dims.append(int(read[6 * i: 6 * (i + 1)]))
            is_variable = (len(dims) == 5)
            numbloc = dims[1 if is_variable else 0]
            dims = numpy.asarray(dims[2:] if is_variable else dims[1:])
            vars = fichier_v3d.read(length_format * dims.prod() + dims.prod() / 6 + ((dims.prod()%6) != 0))
            # vars = fichier_v3d.read(length_format * dims.prod())
            # vars = vars.replace('\n','')
            vars = vars.replace('\n',' ' if compter_saut_de_ligne else '')
            vars = numpy.fromstring(vars, numpy.dtype('|S' + str(length_format)), count = dims.prod())
            vars = numpy.asarray(vars, dtype=numpy.dtype('float')).reshape(dims.prod(),1)[:,0]
            var=nom_var.replace('va','').replace(' ','')
            data[var] = vars
    fichier_v3d.close()
    OutputDictionnary = {'numbloc': numbloc, 'dims' : dims, 'data' : data}
    return OutputDictionnary
#________________________________________________________________________________

#________________________________________________________________________________
def get_numeros_blocs_non_vides(vtkMultiBlockDataSet):
    """retourne le numero des blocs non-vides d'un MultiBlockDataset(
    """
    list_of_blocks = []
    for numbloc in range(vtkMultiBlockDataSet.GetNumberOfBlocks()):
        if vtkMultiBlockDataSet.GetBlock(numbloc) != None:
            if vtkMultiBlockDataSet.GetBlock(numbloc).GetNumberOfPoints() != 0:
                list_of_blocks.append(numbloc)
    return list_of_blocks
#________________________________________________________________________________

#________________________________________________________________________________
def vtk_new_instance(vtkDataObject):
    """genere un nouvelle instance de la classe
    en faisant attention a ce qu'elle ne soit pas instanciee deux fois
    
    """
    output = vtkDataObject.NewInstance()
    while(output.GetReferenceCount() > 1):
        output.UnRegister(None)
    return output
#________________________________________________________________________________

#________________________________________________________________________________
def vtk_new_shallowcopy(vtkDataObject, mode_multibloc=True):
    """genere un nouvelle instance de la classe
    en faisant attention a ce qu'elle ne soit pas instanciee deux fois
    
    et copie superficielle de l'objet vtk fourni en entree
    
    si mode_multibloc, alors les blocs composant le multibloc sont eux 
    aussi shallow copies
    """
    # creation d'un nouvel objet de la meme classe
    output = vtkDataObject.NewInstance()
    while(output.GetReferenceCount() > 1):
        output.UnRegister(None)
    # ShallowCopy
    if isinstance(vtkDataObject, vtk.vtkMultiBlockDataSet) and mode_multibloc is True:
        for numbloc in get_numeros_blocs_non_vides(vtkDataObject):
            output.SetBlock(numbloc, 
                vtk_new_shallowcopy(vtkDataObject.GetBlock(numbloc))
                )
    else:
        output.ShallowCopy(vtkDataObject)
    return output
#________________________________________________________________________________

#______________________________________________________________________________________________________________________________
def ecrire_v3d(acces_fichier, dict_numpy_arrays, dimensions, numbloc=0, fmt_fichier = "bin",\
    precision="i4r8", endian="big", type_maillage=False):
    """Fonction d'ecriture des fichiers au format Voir3D (v3d).
    
    dict_numpy_arrays est un dictionnaire qui contient les variables a ecrire
    dimensions indique les dimensions qui doivent etre ecrites dans le fichier
        SOUS LA FORME D'UN TUPLE
    
    type_maillage conditionne si le fichier v3d ecrit est du type maillage ou donnees
        (en-tete differents : dans le cas d'une variable, "va " est rajoute 
        devant le nom de la variable au moment de l'ecriture, 
        et le numero de la variable est ecrit avant les dimensions)
    
    """
    ## definition des dictionnaires:
    format_binaire = {"i4" : [4,"l"], "i8" : [8,"q"], "r4" : [4,"f"], "r8" : [8,"d"], "big":">", 
        "little" : "<"}
    
    nb_vars = len(dict_numpy_arrays)
    if nb_vars == 0:
        print 'None to save'
        return 1
    
    ## ouverture du fichier:
    try: 
        fichier_v3d = open(acces_fichier,"w" if fmt_fichier == 'fmt' else "wb")
    except:
        raise IOError, "Le fichier ne peut pas etre ecrit"
    else:
        print 'Ecriture %s'%(acces_fichier)
    
    ## ecriture du fichier:
    if fmt_fichier == "fmt":
        fichier_v3d.write((5*" " + str(nb_vars))[-5:] + "\n")
        incr_var = 0
        for nom_var in dict_numpy_arrays if not type_maillage else ['x', 'y', 'z']:
            numpyArray = dict_numpy_arrays[nom_var]
            
            while len(dimensions) < 3:
                dimensions += (1,)
            numpyArray = numpyArray.ravel()
            
            if not type_maillage:
                var = "va " + nom_var
            else:
                var = nom_var
            
            if type_maillage:
                fichier_v3d.write((var + 20*" ")[:20] + (20*" " + "6e15.7")[-20:] + "\n" )
                fichier_v3d.write((6*" " + str(numbloc))[-6:] +\
                    (6*" " +  str(dimensions[0]))[-6:] +\
                    (6*" " +  str(dimensions[1]))[-6:] +\
                    (6*" " +  str(dimensions[2]))[-6:] + "\n")
            else:
                fichier_v3d.write((var + 20*" ")[:20] +  (20*" " + "6e14.7")[-20:] + "\n")
                fichier_v3d.write((6*" " + str(incr_var + 1))[-6:] +\
                    (6*" " + str(numbloc))[-6:] +\
                    (6*" " +  str(dimensions[0]))[-6:] +\
                    (6*" " +  str(dimensions[1]))[-6:] +\
                    (6*" " +  str(dimensions[2]))[-6:] + "\n")
                incr_var += 1

            incr = 0
            if type_maillage:
                for value in numpyArray:
                    incr += 1
                    fichier_v3d.write(
                        '{0: 15.7E}{1}'.format(value, 
                        "\n" if incr%6 == 0 or incr == numpyArray.size
                        else ""))
            else:
                for value in numpyArray:
                    incr += 1
                    fichier_v3d.write(
                        '{0: 14.7E}{1}'.format(value, 
                        "\n" if incr%6 == 0 or incr == numpyArray.size
                        else ""))
        
    elif fmt_fichier == "bin":
        # options de lecture pour struct.unpack
        ordre_bin = '>' if endian == 'big' else '<'
        fmt_entier = "l" if precision[:2] == 'i4' else 'q'
        taille_entier = int(precision[1])
        fmt_reel = "f" if precision[-2:] == 'i4' else 'd'
        taille_reel = int(precision[3])
        
        fichier_v3d.write(struct.pack(ordre_bin + fmt_entier, taille_entier))
        fichier_v3d.write(struct.pack(ordre_bin + fmt_entier, nb_vars))
        fichier_v3d.write(struct.pack(ordre_bin + fmt_entier, taille_entier))
        
        incr_var = 0
        for nom_var in dict_numpy_arrays if not type_maillage else ['x', 'y', 'z']:
            numpyArray = dict_numpy_arrays[nom_var]
            
            while len(dimensions) < 3:
                dimensions += (1,)
            numpyArray = numpyArray.ravel()
            
            if not type_maillage:
                var = ("va  " + nom_var + 20 * " ")[:20]
            else:
                var = (nom_var + 20 * " ")[:20]
            fichier_v3d.write(struct.pack(ordre_bin + "l", 20))
            fichier_v3d.write(struct.pack("20s", var))
            fichier_v3d.write(struct.pack(ordre_bin + "l" ,20))
            
            if type_maillage: 
                fichier_v3d.write(struct.pack(ordre_bin + "l", 4 * taille_entier))
                fichier_v3d.write(struct.pack(ordre_bin + "4" + fmt_entier, 
                    numbloc, dimensions[0], dimensions[1], dimensions[2])
                    )
                fichier_v3d.write(struct.pack(ordre_bin + "l", 4 * taille_entier))
            else:
                fichier_v3d.write(struct.pack(ordre_bin + "l", 5 * taille_entier))
                fichier_v3d.write(struct.pack(ordre_bin + "5" + fmt_entier,
                    incr_var + 1 , numbloc, dimensions[0], dimensions[1], dimensions[2])
                    )
                fichier_v3d.write(struct.pack(ordre_bin + "l", 5 * taille_entier))
                incr_var += 1
            
            fichier_v3d.write(struct.pack(ordre_bin + "l", numpy.prod(dimensions) * taille_reel))
            for value_point in numpyArray:
                fichier_v3d.write(struct.pack(ordre_bin + fmt_reel, value_point))
            fichier_v3d.write(struct.pack(ordre_bin + "l", numpy.prod(dimensions) * taille_reel))
    fichier_v3d.close()
    return 0
#______________________________________________________________________________________________________________________________

#_____________________________________________________________________________________
def modifier_triedre(input, nv_triedre):
    """manipulation du triedre (i, j, k)
    
    s'applique a un vtkStructuredGrid ou a un MultiBlockDataSet
    """
    
    if isinstance(input, vtk.vtkMultiBlockDataSet):
		maillage_new = vtk.vtkMultiBlockDataSet()
		for numbloc in get_numeros_blocs_non_vides(input):
			bloc_new = modifier_triedre(input = input.GetBlock(numbloc), nv_triedre = nv_triedre)
			maillage_new.SetBlock(numbloc, bloc_new)
		return maillage_new
    
    output = vtk.vtkStructuredGrid()
    # changement des dimensions
    dims = input.GetDimensions()
    output.SetDimensions([dims[abs(nv_triedre[0]) - 1], dims[abs(nv_triedre[1]) - 1], 
                dims[abs(nv_triedre[2]) - 1]])
    # lecture des coordonnees des points
    coords = numpy_support.vtk_to_numpy(input.GetPoints().GetData())
    coords = coords.reshape([dims[2], dims[1] , dims[0], 3]).transpose(2, 1, 0, 3)
    # modification du triedre
    new_coords = copy.deepcopy(numpy.ascontiguousarray(coords))
    new_coords = new_coords[::numpy.sign(nv_triedre[0]), ::numpy.sign(nv_triedre[1]), 
            ::numpy.sign(nv_triedre[2])]
    new_coords = new_coords.transpose(abs(nv_triedre[0]) - 1, abs(nv_triedre[1]) - 1, 
        abs(nv_triedre[2]) - 1, 3)
    # remise dans l'ordre k,j,i pour stockage en vtk puis ravel
    new_coords = new_coords.transpose(2, 1, 0, 3).reshape(new_coords.size / 3, 3)
    # modification des points du inputs
    nv_points = vtk.vtkPoints()
    nv_points.SetData(numpy_support.numpy_to_vtk(
        numpy.ascontiguousarray(new_coords), deep = 1))
    output.SetPoints(nv_points)
    output.Update()
    
    # il faut modifier aussi l'ordre des donnees aux points
    for numarray in range(input.GetPointData().GetNumberOfArrays()):
        nb_composantes = input.GetPointData().GetArray(numarray).GetNumberOfComponents()
        data = numpy_support.vtk_to_numpy(input.GetPointData().GetArray(numarray)).reshape(
            dims + (nb_composantes, )).transpose(2, 1, 0, nb_composantes)
        nv_data = numpy.array(data)
        nv_data = nv_data[::numpy.sign(nv_triedre[0]), ::numpy.sign(nv_triedre[1]), 
            ::numpy.sign(nv_triedre[2])]
        nv_data = nv_data.transpose(abs(nv_triedre[0]) - 1, abs(nv_triedre[1]) - 1, 
            abs(nv_triedre[2]) - 1, nb_composantes)
        # remise dans l'ordre k,j,i pour stockage en vtk puis ravel
        nv_data = nv_data.transpose(2, 1, 0, nb_composantes).reshape(nv_data.size / nb_composantes, 
            nb_composantes)
        vtkArray = numpy_support.numpy_to_vtk(numpy.ascontiguousarray(nv_data), deep = 1)
        vtkArray.SetName(input.GetPointData().GetArrayName(numarray))
        output.GetPointData().AddArray(vtkArray)
    
    # il faut modifier aussi l'ordre des donnees aux cellules
    for numarray in range(input.GetCellData().GetNumberOfArrays()):
        nb_composantes = input.GetCellData().GetArray(numarray).GetNumberOfComponents()
        data = numpy_support.vtk_to_numpy(input.GetCellData().GetArray(numarray)).reshape(
            dims[2] - 1, dims[1] - 1, dims[0] - 1, nb_composantes).transpose(2, 1, 0, nb_composantes)
        nv_data = numpy.array(data)
        nv_data = nv_data[::numpy.sign(nv_triedre[0]), ::numpy.sign(nv_triedre[1]), 
            ::numpy.sign(nv_triedre[2])]
        nv_data = nv_data.transpose(abs(nv_triedre[0]) - 1, abs(nv_triedre[1]) - 1, 
            abs(nv_triedre[2]) - 1, nb_composantes)
        # remise dans l'ordre k,j,i pour stockage en vtk puis ravel
        nv_data = nv_data.transpose(2, 1, 0, nb_composantes).reshape(nv_data.size / nb_composantes, 
            nb_composantes)
        vtkArray = numpy_support.numpy_to_vtk(numpy.ascontiguousarray(nv_data), deep = 1)
        vtkArray.SetName(input.GetCellData().GetArrayName(numarray))
        output.GetCellData().AddArray(vtkArray)
    print 'Modification du triedre faite : ', nv_triedre
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def get_variables_in_function(function):
    """retourne une liste contenant les variables necessaires 
    au calcul de la fonction"""
    op_vect = ['grad', 'div', 'rot', 'lapl']
    op_calc = ['(', ')', '+', '-', '*', '/', '.', '^', 'abs', 'acos', 'asin', 'atan', 
        'ceil', 'cos', 'cosh', 'exp', 'floor', 'log', 'max', 'min', 'mag', 'norm', 'sign', 
        'sin', 'sinh', 'sqrt', 'tan', 'tanh', 'iHat','jHat','kHat', 'arccos', 'arcsin', 'arctan']
    liste_vars = []
    read_g=0
    read_d=1
    while read_g < len(function):
        bloc = function[read_g:read_d]
        if bloc in op_calc or bloc.isdigit() or bloc == ' ' or bloc == '':
            read_g = read_d
        elif bloc in op_vect:
            i=0
            nb_parentheses = 0
            Loop = True
            while Loop:
                if function[read_d + i] == '(':
                    nb_parentheses = nb_parentheses + 1
                if function[read_d + i] == ')':
                    nb_parentheses = nb_parentheses - 1 
                if nb_parentheses == 0:
                    Loop = False
                    ##liste_vars.append(function[read_d + 1:read_d+i])
                    liste_vars.append((bloc + '('+function[read_d + 1:read_d+i]+')').replace(' ', ''))
                    #liste_vars.append((function[read_d + 1:read_d+i]).replace(' ', ''))
                    read_g = read_d + i + 1
                    read_d = read_g
                i = i+1
        elif bloc[-1] in op_calc:
            read_g = read_d
            liste_vars.append(bloc[:-1].replace(' ', ''))
        elif read_d == len(function):
            read_g = read_d
            liste_vars.append(bloc.replace(' ', ''))
        read_d = read_d + 1
    if len(liste_vars) == 1 and liste_vars[0] == function:
        liste_vars[0] = liste_vars[0][liste_vars[0].find('(') + 1 : liste_vars[0].rfind(')')]
    liste_a_parcourir = copy.copy(liste_vars)
    for name in liste_a_parcourir:
        while liste_vars.count(name) > 1:
            liste_vars.remove(name)
    return liste_vars
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def get_noms_arrays_presents(input, loc='points'):
    """retourne la liste des array presents
    si multibloc, alors ne regarde que le premier bloc non vide
    
    si loc est different de points, alors regarde aux cellules
    """
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        names = get_noms_arrays_presents(input.GetBlock(get_numeros_blocs_non_vides(input)[0]), loc = loc)
    else:
        bloc = input
        names = []
        for numarray in range(bloc.GetPointData().GetNumberOfArrays() if loc == 'points'
                else bloc.GetCellData().GetNumberOfArrays()):
            names.append(bloc.GetPointData().GetArray(numarray).GetName() if loc == 'points'
                else bloc.GetCellData().GetArray(numarray).GetName())
    return names
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def renommer_arrays(input, dict_renommer, loc='points'):
    """renomme les arrays. 
    dict_renommer est un dictionnaire {ancien_nom: nouveau_nom, ... }
    
    """
    output = vtk_new_shallowcopy(input)
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, renommer_arrays(input.GetBlock(numbloc), dict_renommer, loc))
    else:
        data = input.GetPointData() if loc == 'points' else input.GetCellData()
        for old_name in dict_renommer:
            data.GetArray(old_name).SetName(dict_renommer[old_name])
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def set_scalaires_actifs(input, array_name, loc = 'points'):
    """definit un champ de scalaires comme actif aux points ou cellules
    mono ou multi blocs
    
    """
    output = vtk_new_shallowcopy(input)
    if not array_name in get_noms_arrays_presents(output, loc='points') \
            + get_noms_arrays_presents(output, loc='cellules'):
        raise IOError, "output n'a pas d'array {0}".format(array_name)
    if isinstance(output, vtk.vtkMultiBlockDataSet):
        for numbloc in get_numeros_blocs_non_vides(output):
            output.SetBlock(numbloc, set_scalaires_actifs(
                input.GetBlock(numbloc), loc = loc, array_name = array_name))
    else:
        output.GetPointData().SetActiveScalars(array_name) if loc == 'points' \
            else output.GetCellData().SetActiveScalars(array_name)
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def VTKProbe(input, source, tolerance=None):
    """probe generique - mono/multi bloc par mono/multi bloc
    
    input contient la GEOMETRIE sur laquelle interpoler les donnees NOUVEAU
    source contient le maillage qui contient les donnees ANCIEN
    
    tolerance permet de faire des trucs. Essayer des valeurs, 0.1 par exemple. 
    Voir la doc du vtkProbeFilter. 
    """
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk.vtkMultiBlockDataSet()
        for numbloc in get_numeros_blocs_non_vides(input):
            print "probe par le bloc {0} de input".format(numbloc)
            output.SetBlock(numbloc, VTKProbe(input = input.GetBlock(numbloc), source = source, tolerance=tolerance))
    elif isinstance(source, vtk.vtkMultiBlockDataSet):
        output = vtk_new_shallowcopy(input)
        dict_data = {}
        
        for numbloc in get_numeros_blocs_non_vides(source):
            print "probe du bloc {0} de la source".format(numbloc)
            bloc = VTKProbe(input, source.GetBlock(numbloc), tolerance=tolerance)
            
            for nom_array in get_noms_arrays_presents(bloc, loc = 'points'):
                array = get_vtk_array_as_numpy_array(bloc, nom_array)
                if dict_data.has_key(nom_array):
                    dict_data[nom_array] += array
                else:
                    dict_data[nom_array] = array
        
        # if len(dict_data.keys()) > 0:
            # argwhere = numpy.argwhere(dict_data['vtkValidPointMask'] == 0)
            # dict_data['vtkValidPointMask'][argwhere] = 1
            # for nom_array in dict_data:
                # if nom_array != 'vtkValidPointMask':
                    # if len(dict_data[nom_array].shape) == 2:
                        # for k in range(dict_data[nom_array].shape[1]):
                            # dict_data[nom_array][:, k] /= dict_data['vtkValidPointMask']
                    # else:
                        # dict_data[nom_array] /= dict_data['vtkValidPointMask']
                    # 
                # vtk_array = numpy_support.numpy_to_vtk(dict_data[nom_array], deep = 1)
                # vtk_array.SetName(nom_array)
                # output.GetPointData().AddArray(vtk_array)
        for nom_array in dict_data:
            vtk_array = numpy_support.numpy_to_vtk(dict_data[nom_array], deep = 1)
            vtk_array.SetName(nom_array)
            output.GetPointData().AddArray(vtk_array)

    else:
        #conversion en polydata -- VTKProbe ne fonctionne pas bien avec les StructuredGrid
        #et il ne respecte pas non plus la structure des objets...
        geom = vtk.vtkPolyData()
        geom.SetPoints(input.GetPoints())
        
        # on cree output et on supprime les arrays qui sont presents dans input
        # pour eviter la confusion
        output = vtk_new_shallowcopy(input)
        for nom_array in get_noms_arrays_presents(input, 'points'):
            output.GetPointData().RemoveArray(nom_array)
        for nom_array in get_noms_arrays_presents(input, 'cellules'):
            output.GetCellData().RemoveArray(nom_array)
        
        #verification de l'intersection des boites de input et source
        bounds_geom = geom.GetBounds()
        bounds_source = source.GetBounds()
        if bounds_geom[0] > bounds_source[1] \
                or bounds_geom[1] < bounds_source[0] \
                or bounds_geom[2] > bounds_source[3] \
                or bounds_geom[3] < bounds_source[2] \
                or bounds_geom[4] > bounds_source[5] \
                or bounds_geom[5] < bounds_source[4]:
            return output
        
        filtre = vtk.vtkProbeFilter()
        vtk_set_input(filtre, geom)
        try:
            filtre.SetSource(source)
        except:
            filtre.SetSourceData(source)
        if tolerance is not None:
            filtre.ComputeToleranceOff()
            filtre.SetTolerance(tolerance)
            print 'TOLERANCE pour vtkProbeFilter reglee manuellement a {0}'.format(tolerance)
        filtre.Update()
        
        #on ajoute ce qu'on a probe a output
        for nom_array in get_noms_arrays_presents(filtre.GetOutput()):
            #print nom_array
            output.GetPointData().AddArray(
                filtre.GetOutput().GetPointData().GetArray(nom_array))
    if numpy.max(get_vtk_array_as_numpy_array(output, 'vtkValidPointMask')) > 1:
        print 'ATTENTION ATTENTION ATTENTION : verifier le resultat. La sonde intersecte plusieurs domaines.'
    return output
#_____________________________________________________________________________________

##_____________________________________________________________________________________
#def interpoler_avec_scipy(source, input, \
        #methode_interpolation='nearest', fill_value=0):
    #"""fonction d'interpolation
        #- d'un multiblock sur un autre, les blocs etant confondus
        #- d'un bloc sur un autre
    #utilise la routine scipy.interpolate.griddata
    
    #- input contient le nouveau maillage
    #- source contient l'ancien maillage, avec les donnees
    #"""
    #output = vtk_new_shallowcopy(input)
    
    #if isinstance(output, vtk.vtkMultiBlockDataSet):
        #for numbloc in get_numeros_blocs_non_vides(output):
            #print "interpolation bloc ", numbloc
            #output.SetBlock(numbloc, 
                #interpoler_avec_scipy(source = source.GetBlock(numbloc), input = output.GetBlock(numbloc))
                #)
    #else:
        ## interpolation des donnees aux noeuds
        #coords_source = numpy_support.vtk_to_numpy(source.GetPoints().GetData())
        #coords_input = numpy_support.vtk_to_numpy(output.GetPoints().GetData())
        #for numarray in range(source.GetPointData().GetNumberOfArrays()):
            #nom_array = source.GetPointData().GetArrayName(numarray)
            #print 'array ', nom_array
            #z_source = numpy_support.vtk_to_numpy(source.GetPointData().GetArray(numarray))
            
            #z_input = interpolate.griddata(
                #coords_source, z_source, 
                #coords_input, 
                #method = 'nearest', 
                #fill_value = 0
                #)
            
            #z_input = numpy_support.numpy_to_vtk(z_input, deep = 1)
            #z_input.SetName(nom_array)
            #output.GetPointData().AddArray(z_input)
        ## interpolation des donnees aux cellules
        #f = vtk.vtkCellCenters()
        #f.SetInputData(source)
        #f.Update()
        #coords_cellules_source = numpy_support.vtk_to_numpy(f.GetOutput().GetPoints().GetData())
        #f = vtk.vtkCellCenters()
        #f.SetInputData(output)
        #f.Update()
        #coords_cellules_input = numpy_support.vtk_to_numpy(f.GetOutput().GetPoints().GetData())
        #for numarray in range(source.GetCellData().GetNumberOfArrays()):
            #nom_array = source.GetCellData().GetArrayName(numarray)
            #print 'array ', nom_array
            #z_source = numpy_support.vtk_to_numpy(source.GetCellData().GetArray(numarray))
            
            #z_input = interpolate.griddata(
                #coords_cellules_source, z_source, 
                #coords_cellules_input, 
                #method = methode_interpolation, 
                #fill_value = fill_value
                #)
            #z_input = numpy_support.numpy_to_vtk(z_input, deep = 1)
            #z_input.SetName(nom_array)
            #output.GetCellData().AddArray(z_input)
    #return output
##_____________________________________________________________________________________

#_____________________________________________________________________________________
def interpoler_avec_scipy(source, input, \
        methode_interpolation='nearest', fill_value=0):
    """fonction d'interpolation -- nouvelle version qui pourrait (?) gerer mieux les multibloc
        - d'un multiblock sur un autre
        - d'un bloc sur un autre
    utilise la routine scipy.interpolate.griddata
    
    - source = ANCIEN contient l'ancien maillage, avec les donnees
    - input = NOUVEAU contient le nouveau maillage
    """
    output = vtk_new_shallowcopy(input)
    
    if isinstance(output, vtk.vtkMultiBlockDataSet):
        for numbloc in get_numeros_blocs_non_vides(output):
            print "interpolation bloc ", numbloc
            output.SetBlock(numbloc, 
                interpoler_avec_scipy(source = source, input = output.GetBlock(numbloc))
                )
    else:
        #a partir de la, output est un bloc et source peut etre un bloc ou un multibloc
        
        # interpolation des donnees aux noeuds
        coords_source = get_vtk_array_as_numpy_array(source, 'coords')
        coords_input = get_vtk_array_as_numpy_array(output, 'coords')
        for nom_array in get_noms_arrays_presents(source, loc='points'):
            print 'array ', nom_array
            z_source = get_vtk_array_as_numpy_array(source, nom_array)
            
            z_input = interpolate.griddata(
                coords_source, z_source, 
                coords_input, 
                method = methode_interpolation, 
                fill_value = 0
                )
            
            z_input = numpy_support.numpy_to_vtk(z_input, deep = 1)
            z_input.SetName(nom_array)
            output.GetPointData().AddArray(z_input)
        # interpolation des donnees aux cellules
        f = vtkFiltersGeneral.vtkCellCenters()
        if isinstance(source, vtk.vtkMultiBlockDataSet):
            source_cells = appliquer_sur_multibloc(f, source)
        else:
            vtk_set_input(f, source)
            f.Update()
            source_cells = f.GetOutput()
        coords_cellules_source = get_vtk_array_as_numpy_array(source_cells, 'coords')
        
        f = vtkFiltersGeneral.vtkCellCenters()
        vtk_set_input(f, output)
        f.Update()
        coords_cellules_input = get_vtk_array_as_numpy_array(f.GetOutput(), 'coords')
        for nom_array in get_noms_arrays_presents(source, loc = 'cellules'):
            print 'array ', nom_array
            z_source = get_vtk_array_as_numpy_array(source, nom_array)
            
            z_input = interpolate.griddata(
                coords_cellules_source, z_source, 
                coords_cellules_input, 
                method = methode_interpolation, 
                fill_value = fill_value
                )
            z_input = numpy_support.numpy_to_vtk(z_input, deep = 1)
            z_input.SetName(nom_array)
            output.GetCellData().AddArray(z_input)
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def get_cell_centers(input):
    """Utilise le filtre vtkCellCenters pour retourner les coordonnees des centres des cellules sous forme de numpy array
    vtkCellCenters is a filter that takes as input any dataset and generates on output points at the center of the cells in the dataset. 
    """
    f = vtkFiltersGeneral.vtkCellCenters()
    vtk_set_input(f, input)
    f.Update()
    mesh_centers = f.GetOutput()
    return get_vtk_array_as_numpy_array(mesh_centers, 'coords')
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def lire_fichier_tecplot(acces_fichier, sep=" "):
    """fonction de lecture d'un fichier au format tecplot ASCII
    
    """
    f = file(acces_fichier, 'r')
    f.readline()
    read = f.readline()
    read = read.split('"')
    #liste des variables contenues dans les donnees
    variables = [read[2*i+1] for i in range(int(len(read)/2))]
    #lecture des dimensions
    ligne = f.readline()
    dims = []
    for cle in ['I=', 'J=', 'K=']:
        if cle in ligne:
            dims.append(int(ligne.split(cle)[1].split(',')[0]))
        else:
            dims.append(1)
    data_out = {}
    data_out['dims'] = dims
    #lecture des donnees
    read = f.read()
    read = read.replace('\n', sep).split(sep)
    while '' in read:
        read.remove('')
    #read contient toutes les valeurs dans l'ordre, sous forme de string
    #Separation des variables
    NITER = len(read)/len(variables)
    
    data_out['data'] = {}
    for var in range(0, len(variables)):
        data_out['data'][variables[var]] = [float(read[i]) for i in range(NITER*var, NITER*(var+1))]
    return data_out
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def get_vtk_array_as_numpy_array(input, nom_array, copy=False, loc=None):
    """fonction qui retourne le array du vtkDataObject donne en input
    sous forme d'un numpy array en utilisant numpy_support
    
    la donnee peut etre indiferemment aux cellules ou aux points
    si input est un vtkMultiBlockDataSet, alors les donnees sont concatenees 
    dans l'ordre des blocs croissants
    
    nom_array peut etre coords, pour avoir les coordonnees des points
    nom_array peut etre centers, pour avoir les coordonnees des centres des cellules
    
    nom_array peut etre polys pour avoir l'arbre de connectivite (pour un vtkPolyData uniquement)
    
    Pour un maillage non-structure :
        - nom_array peut etre cells
        - nom_array peut etre cellstypes
        - nom_array peut etre cellslocations
    
    
    Si loc est None, alors la donnee est recherchee en priorite aux points, puis aux cellules 
        si elle n'a pas ete trouvee aux points.
    Si loc est donnee, alors on ne cherche que aux points ou au cellules
        loc = 'points' pour chercher la donnee aux points uniquement
        loc = <autre chose> pour chercher la donnee aux cellules uniquement
    """
    # cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = None
        for numbloc in get_numeros_blocs_non_vides(input):
            if output is None:
                output = get_vtk_array_as_numpy_array(input.GetBlock(numbloc), nom_array)  
            elif nom_array == 'polys': 
                raise IOError, 'pas disponible pour le cas multibloc'
            elif nom_array == 'cells': 
                raise IOError, 'pas disponible pour le cas multibloc'
            elif nom_array == 'cellstypes': 
                raise IOError, 'pas disponible pour le cas multibloc'
            elif nom_array == 'cellslocations': 
                raise IOError, 'pas disponible pour le cas multibloc'
            else: 
                output = numpy.concatenate((
                    output, get_vtk_array_as_numpy_array(input.GetBlock(numbloc), nom_array)), 
                    axis = 0)
    # cas monobloc
    else:
        if nom_array == 'coords':
            output = numpy_support.vtk_to_numpy(input.GetPoints().GetData())
        elif nom_array == 'centers':
            output = get_cell_centers(input)
        elif nom_array == 'polys':
            output = numpy_support.vtk_to_numpy(input.GetPolys().GetData())
            # if output.size % 4 == 0:
                # output = output.reshape(output.size / 4, 4)[:, 1:]
            # elif output.size % 5 == 0:
                # output = output.reshape(output.size / 5, 5)[:, 1:]
            # else:
                # raise IOError, "type de cellule non reconnu"
        elif nom_array == 'cells':
            output = numpy_support.vtk_to_numpy(input.GetCells().GetData())
        elif nom_array == 'cellstypes':
            output = numpy_support.vtk_to_numpy(input.GetCellTypesArray())
        elif nom_array == 'cellslocations':
            output = numpy_support.vtk_to_numpy(input.GetCellLocationsArray())
        else:
            if loc == None:
                array_aux_points = input.GetPointData().HasArray(nom_array)
                array_aux_cellules = input.GetCellData().HasArray(nom_array)
                if array_aux_points and array_aux_cellules:
                    raise Exception, "conflit : {0} est present aux noeuds ET aux cellules".format(nom_array)
                elif not array_aux_points and not array_aux_cellules:
                    raise Exception, "erreur : {0} n'est present NI aux noeuds NI aux cellules".format(nom_array)
                if array_aux_points:
                    print 'La donnee est aux POINTS'
                    loc_data = input.GetPointData() 
                else:
                    print 'La donnee est aux CELLULES'
                    loc_data = input.GetCellData()
            else:
                loc_data = input.GetPointData() if loc == 'points' else input.GetCellData()
            output = numpy_support.vtk_to_numpy(loc_data.GetArray(nom_array))
    if copy == True:
        output = numpy.array(output)
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def calculer_surfaces_cellules(input, retourner_surface_totale=False):
    """fonction qui ajoute aux cellules la valeur de l'aire de chacune 
    des cellules
    
    mono/multi bloc(s)
    s'appuie sur la classe vtk.vtkCellQuality
    
    """
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk.vtkMultiBlockDataSet()
        aire_surface = 0
        for numbloc in get_numeros_blocs_non_vides(input):
            bloc, aire_surface_bloc = calculer_surfaces_cellules(input.GetBlock(numbloc), 
                retourner_surface_totale = True)
            output.SetBlock(numbloc, bloc)
            aire_surface += aire_surface_bloc
    else:
        if isinstance(input, vtk.vtkStructuredGrid):
            dimensions = list(input.GetDimensions())
            if dimensions.count(1) != 1:
                raise IOError, "le vtkStructuredGrid en entree n'est pas une surface"
        
        f = vtkFiltersVerdict.vtkMeshQuality()
        vtk_set_input(f, input)
        
        f.SetTriangleQualityMeasureToArea()
        f.SetQuadQualityMeasureToArea()
        f.SetTetQualityMeasureToVolume()
        f.SetHexQualityMeasureToVolume()
        f.Update()
        output = f.GetOutput()
        output.GetCellData().GetArray('Quality').SetName('CellSurface')
        
        aire_surface = numpy.sum(
            get_vtk_array_as_numpy_array(output, 'CellSurface'))
    if retourner_surface_totale == True:
        return output, aire_surface
    else:
        return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def integrer_sur_la_surface(input, array_a_integrer, array_poids = 'CellSurface', 
        remplacer_nan=False):
    """fonction qui integre in array sur une surface. Aux cellules ou aux noeuds. 
    mono/multi bloc(s)
    possiblite d'indiquer 1 pour array_a_integrer 
    
    Si remplacer_nan est vrai, alors les nan sont remplaces par zero
    
    array_a_integrer et array_poids doivent deja etre presents aux cellules ou aux noeuds
    
    array aux points et poids aux cellules --> on interpole array aux cellules
    array aux cellules et poids aux points --> on interpole poids aux cellules
    """
    for array in [array_a_integrer, array_poids]:
        if not(array in get_noms_arrays_presents(input, loc = 'cellules')
                or array in get_noms_arrays_presents(input, loc = 'points')) \
                and array != 1:
            raise IOError, '{0} present ni aux noeuds ni aux cellules'.format(array)
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        resultat_integration = 0
        for numbloc in get_numeros_blocs_non_vides(input):
            resultat_integration_bloc = integrer_sur_la_surface(input.GetBlock(numbloc), 
                array_a_integrer, array_poids)
            resultat_integration += resultat_integration_bloc
    else:
        #si array_a_integrer est egal a 1
        if array_a_integrer == 1:
            array_poids = get_vtk_array_as_numpy_array(input, array_poids)
            array_a_integrer = numpy.ones(shape=array_poids.shape)
        #les deux sont situes au meme endroit
        elif array_a_integrer in get_noms_arrays_presents(input, 'points') and \
                array_poids in get_noms_arrays_presents(input, 'points') or \
                array_a_integrer in get_noms_arrays_presents(input, 'cellules') and \
                array_poids in get_noms_arrays_presents(input, 'cellules'):
            array_a_integrer = get_vtk_array_as_numpy_array(input, array_a_integrer)
            array_poids = get_vtk_array_as_numpy_array(input, array_poids)
        #array aux points et poids aux cellules --> on interpole array aux cellules
        elif array_a_integrer in get_noms_arrays_presents(input, 'points') and \
                array_poids in get_noms_arrays_presents(input, 'cellules'):
            array_poids = get_vtk_array_as_numpy_array(input, array_poids)
            f = vtk.vtkPointDataToCellData()
            #f.SetInputData(input)
            vtk_set_input(f, input)
            f.Update()
            array_a_integrer = get_vtk_array_as_numpy_array(f.GetOutput(), array_a_integrer)
        #array aux cellules et poids aux points --> on interpole poids aux cellules
        elif array_a_integrer in get_noms_arrays_presents(input, 'cellules') and \
                array_poids in get_noms_arrays_presents(input, 'points'):
            array_poids = get_vtk_array_as_numpy_array(input, array_poids)
            f = vtk.vtkCellDataToPointData()
            #f.SetInputData(input)
            vtk_set_input(f, input)
            f.Update()
            array_a_integrer = get_vtk_array_as_numpy_array(f.GetOutput(), array_a_integrer)
        
        if remplacer_nan:
            print 'Remplacement dans array_a_integrer de {0} NaN'.format(
                numpy.where(numpy.isnan(array_a_integrer))[0].size)
            array_a_integrer = numpy.nan_to_num(array_a_integrer)
            print 'Remplacement dans array_poids de {0} NaN'.format(
                numpy.where(numpy.isnan(array_poids))[0].size)
            array_poids = numpy.nan_to_num(array_poids)
        
        resultat_integration = numpy.sum(array_a_integrer * array_poids)
    return resultat_integration
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def interpoler_cellules_aux_points(input):
    """interpolation CellDataToPointData
    """
    f = vtk.vtkCellDataToPointData()
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = appliquer_sur_multibloc(f, input)
    else:
        vtk_set_input(f, input)
        f.Update()
        output = f.GetOutput()
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def interpoler_points_aux_cellules(input):
    """interpolation CellDataToPointData
    """
    f = vtk.vtkPointDataToCellData()
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = appliquer_sur_multibloc(f, input)
    else:
        vtk_set_input(f, input)
        f.Update()
        output = f.GetOutput()
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def ajouter_suffixe_aux_arrays(input, suffixe):
    """ajoute un suffixe a tous les arrays de input
    
    sert quand on veut resampler deux dataset ensemble par exemple
    """
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk.vtkMultiBlockDataSet()
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, 
                ajouter_suffixe_aux_arrays(input.GetBlock(numbloc), suffixe)
                )
    else:
        output = vtk_new_shallowcopy(input)
        for nom_array in get_noms_arrays_presents(output, 'points'):
            array = numpy_support.vtk_to_numpy(output.GetPointData().GetArray(nom_array))
            varray = numpy_support.numpy_to_vtk(array, deep = 1)
            varray.SetName(nom_array + suffixe)
            output.GetPointData().AddArray(varray)
            output.GetPointData().RemoveArray(nom_array)
        for nom_array in get_noms_arrays_presents(output, 'cellules'):
            array = get_vtk_array_as_numpy_array(output.GetCellData().GetArray(nom_array))
            varray = numpy_support.numpy_to_vtk(array, deep = 1)
            varray.SetName(nom_array + suffixe)
            output.GetCellData().AddArray(varray)
            output.GetCellData().RemoveArray(nom_array)
    return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def Rotation(input, alpha, axe=2):
    """ancien appel, avec majsucule -- pour compatibilite seulement. NE PAS UTILISER
    """
    return rotation(input, alpha, axe)
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def rotation(input, alpha, axe=2):
    """ Rotation d'un input autour de l'axe - alpha en DEGRES
    L'axe est indique par son numero
    0 pour x, 1 pour y, 2 pour z 
    
    Utile principalement pour dupliquer des canaux
    fait appel a vtkTransform et a vtkArrayCalculator pour faire 
    tourner les vecteurs (arrays a trois composantes) 
    """
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, 
                rotation(input.GetBlock(numbloc), alpha, axe=axe)
                )
        return output
    else:
        #pour eviter que vtkTransform ne tourne un vecteur au points
        #VTK ne permet apparait que l'ajout d'un seul vecteur aux points
        
        
        from time import time
        t1 = time()
        
        input.GetPointData().SetActiveVectors(None)
        
        transform = vtkCommonTransforms.vtkTransform()
        if axe == 0:
            print 'rotation x'
            transform.RotateX(alpha)
        elif axe == 1:
            print 'rotation y'
            transform.RotateY(alpha)
        elif axe == 2:
            print 'rotation z'
            transform.RotateZ(alpha)
        
        transformFilter = vtkFiltersGeneral.vtkTransformFilter()
        
        
        #transformFilter.SetInputData(input)
        vtk_set_input(transformFilter, input)
        
        transformFilter.SetTransform(transform)
        transformFilter.Update()
        output = transformFilter.GetOutput()
        
        NomsVecteurs = []
        for i in range(input.GetPointData().GetNumberOfArrays()):
            if input.GetPointData().GetArray(i).GetNumberOfComponents() == 3:
                NomsVecteurs.append(input.GetPointData().GetArrayName(i))
        
        t2 = time()
        # liste_calculateurs = []
        for vect in NomsVecteurs:
            # calculator = vtk.vtkArrayCalculator()
            # liste_calculateurs.append(calculator)
            # 
            # calculator.SetResultArrayType(vtk.VTK_FLOAT) # 2014.02.14: on force la sortie a etre un float. Sinon ca devient un VTK_DOUBLE et le streamtracer ne fonctionne pas. 
            # #calculator.SetInputData(output)
            # if len(liste_calculateurs) == 1:
                # vtk_set_input(calculator, output)
            # else:
                # calculator.SetInputConnection(liste_calculateurs[-2].GetOutputPort())
            # calculator.AddScalarVariable(vect + "_X", vect, 0)
            # calculator.AddScalarVariable(vect + "_Y", vect, 1)
            # calculator.AddScalarVariable(vect + "_Z", vect, 2)
            # calculator.SetResultArrayName(vect)
            # if axe == 0:
                # calculator.SetFunction('({0}_X)*iHat + (({0}_Y)*cos({1})-({0}_Z)*sin({1}))*jHat + (({0}_Z)*cos({1}) + ({0}_Y)*sin({1}))*kHat'\
                    # .format(vect, alpha * numpy.pi / 180.))
            # elif axe == 1:
                # calculator.SetFunction('({0}_Y)*jHat + (({0}_Z)*cos({1})-({0}_X)*sin({1}))*kHat + (({0}_X)*cos({1}) + ({0}_Z)*sin({1}))*iHat'\
                    # .format(vect, alpha * numpy.pi / 180.))
            # elif axe == 2:
                # calculator.SetFunction('({0}_Z)*kHat + (({0}_X)*cos({1})-({0}_Y)*sin({1}))*iHat + (({0}_Y)*cos({1}) + ({0}_X)*sin({1}))*jHat'\
                    # .format(vect, alpha * numpy.pi / 180.))
            # else:
                # raise IOError, "gni -- moi pas comprendre l'axe de rotation"
        # if len(liste_calculateurs) > 0:
            # calculator = liste_calculateurs[-1]
            # calculator.Update()
            # output = vtk_new_shallowcopy(calculator.GetOutput())
            print 'rotation du vecteur ', vect
            narray = get_vtk_array_as_numpy_array(output, vect)
            narray_new = numpy.ones(narray.shape)
            alpha_rad = numpy.deg2rad(alpha)
            if axe == 0:
                narray_new[:, 0] = (narray[:, 0])
                narray_new[:, 1] = ((narray[:, 1]) * numpy.cos(alpha_rad)-(narray[:, 2]) * numpy.sin(alpha_rad))
                narray_new[:, 2] = ((narray[:, 2]) * numpy.cos(alpha_rad) + (narray[:, 1]) * numpy.sin(alpha_rad))
            elif axe == 1:
                narray_new[:, 1] = (narray[:, 1])
                narray_new[:, 2] = ((narray[:, 2]) * numpy.cos(alpha_rad)-(narray[:, 0]) * numpy.sin(alpha_rad))
                narray_new[:, 0] = ((narray[:, 0]) * numpy.cos(alpha_rad) + (narray[:, 2]) * numpy.sin(alpha_rad))
            elif axe == 2:
                narray_new[:, 2] = (narray[:, 2])
                narray_new[:, 0] = ((narray[:, 0]) * numpy.cos(alpha_rad)-(narray[:, 1]) * numpy.sin(alpha_rad))
                narray_new[:, 1] = ((narray[:, 1]) * numpy.cos(alpha_rad) + (narray[:, 0]) * numpy.sin(alpha_rad))
            else:
                raise IOError, "gni -- moi pas comprendre l'axe de rotation"
            output = ajouter_numpy_array_as_vtk_array(output, narray_new, vect)
        t3 = time()
        # print t2 - t1
        # print t3 - t2
        
        return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def redimensionner(input, coefficient):
    """ redimensionne un maillage, sans toucher aux donnees aux points
    typiquement pour passer de mm en metres ou inversement
    """
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, 
                redimensionner(input.GetBlock(numbloc), coefficient)
                )
        return output
    else:
        #pour eviter que vtkTransform ne tourne un vecteur au points
        #VTK ne permet apparait que l'ajout d'un seul vecteur aux points
        input.GetPointData().SetActiveVectors(None)
        
        transform = vtkCommonTransforms.vtkTransform()
        transform.Scale(coefficient, coefficient, coefficient)
        
        transformFilter = vtkFiltersGeneral.vtkTransformFilter()
        
        
        #transformFilter.SetInputData(input)
        vtk_set_input(transformFilter, input)
        
        transformFilter.SetTransform(transform)
        transformFilter.Update()
        output = transformFilter.GetOutput()
        
        return output
#_____________________________________________________________________________________

#_____________________________________________________________________________________
def ecrire_fichier_colonnes(acces_fichier=None, dictionnaire_donnees=None, ecrasement=False):
    """ fonction d'ecriture d'un fichier formatte contenant des donnees en colonnes
    la premiere ligne contient les noms des variables, les lignes 
    suivantes contiennent les donnees
    
    les donnees contenues dans le dictionnaire peuvent ne pas contenir le 
    meme nombre de donnees
    
    """
    # serie de tests prealables a l'ecriture du fichier
    if acces_fichier is None:
        raise IOError, "acces_fichier n'est pas renseigne"
    if dictionnaire_donnees is None:
        raise IOError, "dictionnaire_donnees n'est pas renseigne"
    if os.path.exists(acces_fichier) and ecrasement is False:
        raise IOError, "un fichier existe deja a l'emplacement indique"
    # ecriture du fichier
    print "Ecriture ", acces_fichier
    fichier = file(acces_fichier, "w")
    # en-tetes
    nombre_de_cles = len(dictionnaire_donnees.keys())
    liste_cles = dictionnaire_donnees.keys()
    for numero_cle in range(nombre_de_cles):
        fichier.write(liste_cles[numero_cle])
        fichier.write('\n' if numero_cle == nombre_de_cles - 1 else '\t')
    # ecriture des donnees
    index = 0
    while True:
        ligne = ''
        compteur = 0
        for cle in liste_cles:
            if not isinstance(dictionnaire_donnees[cle], list):
                valeur = [dictionnaire_donnees[cle]]
            else:
                valeur = dictionnaire_donnees[cle]
            if index < len(valeur):
                compteur += 1
                ligne += str(valeur[index])
            ligne += "\n" if cle == liste_cles[-1] else '\t'
        if compteur != 0:
            fichier.write(ligne)
            index += 1
        else:
            break
    fichier.close()
#_________________________________________________________________

#_____________________________________________________________________________________
def lire_fichier_colonnes(acces_fichier=None, dtype=str, sep="\t"):
    """ fonction de lecture d'un fichier formatte contenant des donnees en colonnes
    la premiere ligne contient les noms des variables
    les lignes suivantes contiennent les donnees
    colonnes separees par une tabulation. 
    
    """
    # serie de tests prealables a l'ecriture du fichier
    if acces_fichier is None:
        raise IOError, "acces_fichier n'est pas renseigne"
    if os.path.exists(acces_fichier) is False:
        raise IOError, "pas de fichier a l'emplacement indique"
    # ecriture du fichier
    fichier = file(acces_fichier, "r")
    bloc = fichier.read()
    bloc = bloc.split('\n')
    cles = bloc[0]
    cles = bloc[0].split(sep)
    
    data = dict.fromkeys(cles)
    for key in data:
        data[key] = []
    
    for ligne in bloc[1:]:
        ligne = ligne.split(sep)
        for num in range(len(ligne)):
            valeur = ligne[num]
            cle = cles[num]
            if valeur != '':
                data[cle].append(valeur)
    for cle in data:
        data[cle] = numpy.array(data[cle], dtype = dtype)
    return data
#_____________________________________________________________________________________

#__________________________________________________________________________________________
def calculer_vecteur_normal(input, normals_aux_cellules=False):
    """fonction qui calcul le vecteur normal
    
    normals_aux_cellules permet d'indiquer si le vecteur normal doit etre calcule
    aux centres des faces ou au noeuds
    
    Dans le cas d'un StructuredGrid, la fonction de calcul de la normale est fait maison
    Dans le cas d'un PolyData, on utilise le filtre VTKPolyDataNormals
    """
    
    #on commence par shallow copier pour ne pas modifier l'objet input
    data = vtk_new_shallowcopy(input)
    
    #le cas structure : il est ecrit a la main parce qu'il n'y a pas de filtre 
    #directement disponible dans vtk pour calculer le vecteur normal en structure
    #une autre solution aurait pu consister a transformer le bloc structure en polydata
    #puis a utiliser vtkPolyDataNormals, mais il ne semble pas alors qu'il soit possible
    #de faire les calculs en double precision
    if isinstance(data, vtk.vtkStructuredGrid):
        coords = get_vtk_array_as_numpy_array(data, "coords")
        #lecture shape_bloc et inversion car stockage k, j, i
        shape_bloc = numpy.asarray(data.GetDimensions())[::-1]
        if not 1 in shape_bloc:
            raise IOError, "le StructuredGrid n'est pas une surface"
        shape_frontiere = numpy.take(
            shape_bloc,
            numpy.argwhere(shape_bloc !=1)
            ).ravel()
        coords = coords.reshape(tuple(shape_frontiere) + (3,))
        
        #Calcul du vecteur normal aux cellules
        if normals_aux_cellules:
            #calcul des diagonales
            vec1 = (coords[1:, 1:] - coords[:-1, :-1])
            vec2 = (coords[1:, :-1] - coords[:-1 ,1:])
            #produit vectoriel
            normals = numpy.cross(vec1, vec2, axis = -1)
            normals = normals.reshape(data.GetNumberOfCells(), 3)
            normals /= numpy.apply_along_axis(numpy.linalg.norm, -1, normals)[:, None]
            Normals = numpy_support.numpy_to_vtk(normals, deep=1)
            Normals.SetName("Normals")
            data.GetCellData().AddArray(Normals)
        #Calcul du vecteur normal aux points
        elif normals_aux_cellules == False:
            
            #___________________________________________________________________________________
            def traitement_frontiere(coords_front):
                """fonction qui calcule le vecteur normal sur une frontiere bien specifique
                en entree, la liste des coordonnees des points, comme si c'etait une frontiere jmin
                coords_front shape doit donc etre (2, nb_points, 3)
                
                Elle sera appelee ensuite de facon a simplifier le code
                """
                normals_frontiere = numpy.zeros((coords_front.shape[1], 3))
                
                normals_frontiere[:-1] = numpy.cross(
                    coords_front[0, 1:] - coords_front[0, :-1], 
                    coords_front[1, :-1] - coords_front[0, :-1],
                    axis = -1)
                normals_frontiere[-1] = numpy.cross(
                    coords_front[1, -1] - coords_front[0, -1],
                    coords_front[0, -2] - coords_front[0, -1],
                    axis = -1)
                normals_frontiere /= numpy.linalg.norm(normals_frontiere)
                return normals_frontiere
            #___________________________________________________________________________________
            
            #calcul des diagonales pour les points au coeur du maillage
            vec1 = (coords[2:, 2:] - coords[:-2, :-2])
            vec1 = vec1 / numpy.linalg.norm(vec1)
            vec2 = (coords[2:, :-2] - coords[:-2 ,2:])
            vec2 = vec2 / numpy.linalg.norm(vec2)
            
            #produit vectoriel
            normals = numpy.zeros(tuple(shape_frontiere) + (3,))
            normals[1:-1, 1:-1] = numpy.cross(vec1, vec2, axis = -1)
            
            #frontiere jmin
            normals[0, :] = traitement_frontiere(coords[:2, :])
            #frontiere jmax
            normals[-1, :] = traitement_frontiere(coords[-2:, :][::-1, ::-1])[::-1]
            #frontiere imin - sans recalcul de normals au coin
            normals[1:-1, 0] = traitement_frontiere(coords[1:-1, :2].transpose(1, 0, 2)[:,::-1])[::-1]
            #frontiere imax
            normals[1:-1, -1] = traitement_frontiere(coords[1:-1, -2:].transpose(1, 0, 2)[::-1])

            #normalisation
            normals = normals.reshape(data.GetNumberOfPoints(), 3)
            normals /= numpy.apply_along_axis(numpy.linalg.norm, -1, normals)[:, None]
            Normals = numpy_support.numpy_to_vtk(normals, deep = 1)
            Normals.SetName("Normals")
            data.GetPointData().AddArray(Normals)
        else:
            raise IOError, normals_aux_cellules + "... C'est quoi ce bidule ? "
    
    #cas d'un polydata : on reutilise simplement le filtre vtkPolyDataNormals
    elif isinstance(input, vtk.vtkPolyData):
        normals = vtk.vtkPolyDataNormals()
        normals.SetComputeCellNormals(normals_aux_cellules)
        normals.SetComputePointNormals(not normals_aux_cellules)
        #normals.SetInputData(data)
        vtk_set_input(normals, data)
        normals.Update()
        if normals_aux_cellules == False:
            data.GetPointData().AddArray(
                normals.GetOutput().GetPointData().GetArray("Normals"))
        else:
            data.GetCellData().AddArray(
                normals.GetOutput().GetCellData().GetArray("Normals"))
    
    #cas d'un multibloc : appel recursif
    elif isinstance(input, vtk.vtkMultiBlockDataSet):
        for numbloc in get_numeros_blocs_non_vides(data):
            data.SetBlock(numbloc, 
                calculer_vecteur_normal(data.GetBlock(numbloc), normals_aux_cellules)
                )
    else:
        raise IOError
    return data
#__________________________________________________________________________________________

#____________________________________________________________________________
def decaler_paroi(paroi, decalage):
    """Fonction qui decale une surface selon le vecteur normal
    
    la surface doit avoir un vecteur normal aux points
    la surface doit etre un polydata
    
    normals doit etre present
    """
    # cas multibloc
    if isinstance(paroi, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(paroi)
        for numbloc in get_numeros_blocs_non_vides(paroi):
            output.SetBlock(numbloc, decaler_paroi(paroi.GetBlock(numbloc), decalage))
        return output
    
    # cas monobloc polydata
    # paroi.GetPointData().SetActiveVectors('Normals')
    # for nom_array in get_noms_arrays_presents(paroi, 'points'):
        # if nom_array != 'Normals':
            # paroi.GetPointData().RemoveArray(nom_array)
    # for nom_array in get_noms_arrays_presents(paroi, 'cellules'):
        # paroi.GetCellData().RemoveArray(nom_array)
    
    #liste_noms_arrays = get_noms_arrays_presents(data, loc = 'points')
    #dict_data = dict.fromkeys(liste_noms_arrays)
    
    #f = vtkFiltersGeometry.vtkGeometryFilter()
    ##f.SetInputData(paroi)
    #vtk_set_input(f, paroi)
    #f.Update()
    #paroi = f.GetOutput()
    
    # warp = vtkFiltersGeneral.vtkWarpVector()
    #warp.SetInputData(paroi)
    # vtk_set_input(warp, paroi)
    # warp.SetScaleFactor(decalage)
    # warp.Update()
    # warp = warp.GetOutput()
    
    coords = get_vtk_array_as_numpy_array(paroi, 'coords')
    normals = get_vtk_array_as_numpy_array(paroi, 'Normals')
    
    new_coords = coords + decalage * normals
    vtk_array = numpy_support.numpy_to_vtk(new_coords, deep=1)
    
    point = vtk.vtkPoints()
    point.SetData(vtk_array)
    output = vtk_new_shallowcopy(paroi)
    output.SetPoints(point)
       
    return output
#____________________________________________________________________________

#____________________________________________________________________________
def appliquer_sur_multibloc(filtre_vtk, multibloc):
    """applique un filtre VTK de maniere iterative sur un multibloc
    
    le filtre_vtk doit etre deja completement configure, sauf Input
    """
    output = vtk.vtkMultiBlockDataSet()
    for numbloc in get_numeros_blocs_non_vides(multibloc):
        vtk_set_input(filtre_vtk, multibloc.GetBlock(numbloc))
        filtre_vtk.Update()
        output.SetBlock(numbloc, vtk_new_shallowcopy(filtre_vtk.GetOutput()))
    return output
#____________________________________________________________________________

#__________________________________________________________________________________________
def convertir_en_polydata(input, calculer_vecteur_normal=True):
    """fonction de conversion en polydata
    avec fusion des blocs si input est un multiblock
    
    utilise les filtres vtkAppendPolyData et vtkGeometryFilter
    """
    # si input est un multiblock, il faut les regrouper en un seul polydata
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        appendFilter = vtk.vtkAppendPolyData()
        for numbloc in get_numeros_blocs_non_vides(input):
            # un test pour s'adapter aux differentes versions de VTK -- changement de noms des methodes
            try:
                appendFilter.AddInput(convertir_en_polydata(input.GetBlock(numbloc)))
            except:
                appendFilter.AddInputData(convertir_en_polydata(input.GetBlock(numbloc)))
        appendFilter.Update()
        input = appendFilter.GetOutput()
    
    # si input n'est pas un multiblock, il faut simplemement le convertir en vtkPolyData
    elif not isinstance(input, vtk.vtkPolyData):
        conversion = vtkFiltersGeometry.vtkGeometryFilter()
        #conversion.SetInputData(input)
        vtk_set_input(conversion, input)
        conversion.Update()
        input = conversion.GetOutput()
    
    if calculer_vecteur_normal == True:
        normals = vtk.vtkPolyDataNormals()
        vtk_set_input(normals, input)
        normals.Update()
        input = normals.GetOutput()
    
    return input
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def merge_multibloc(multibloc, merge_points=True):
    """fonction de conversion d'une multibloc en polydata
    avec fusion des blocs si multibloc est un multiblock
    """
    appendFilter = vtk.vtkAppendFilter()
    if merge_points:
        try:
            appendFilter.MergePointsOn()
        except:
            print 'ATTENTION : La version de VTK ne permet pas de merger les points (vtkAppendFilter.MergePointsOn)'
    else:
        appendFilter.MergePointsOff()
    for numbloc in get_numeros_blocs_non_vides(multibloc):
        try:
            appendFilter.AddInput(convertir_en_polydata(multibloc.GetBlock(numbloc)))
        except:
            appendFilter.AddInputData(convertir_en_polydata(multibloc.GetBlock(numbloc)))
    appendFilter.Update()
    multibloc = appendFilter.GetOutput()
    
    return multibloc
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def ajouter_constante(vtkDataObject, nom_constante, valeur):
    """fonction qui ajoute une constante en chacun des noeuds de l'objet vtk
    Retourne une copie de vtkDataObject (vtk_new_shallowcopy), avec la constante ajoutee
    
    vtkDataObject peut etre indiferemment mono-bloc ou multi-blocs, et consitute de 
    n'importe quel type d'objet vtk qui possede la methode GetPointData()
    
    Typiquement, cette fonction peut etre utilisee pour ajouter une constante omega aux noeuds
    qui est necessaire au calcul des grandeurs
        ajouter_constante(vtkDataObject, 'omega', 0.0)
    
    UNE EXTENSION POURRA ETRE ENVISAGEE POUR PROPOSER D'AJOUTER LA CONSTANTE AUX CENTRES
    """
    if isinstance(vtkDataObject, vtk.vtkMultiBlockDataSet):
        output = vtk.vtkMultiBlockDataSet()
        for numbloc in get_numeros_blocs_non_vides(vtkDataObject):
            output.SetBlock(numbloc,
                ajouter_constante(vtkDataObject.GetBlock(numbloc), nom_constante, valeur)
                )
    else:
        output = vtk_new_shallowcopy(vtkDataObject)
        varray = numpy_support.numpy_to_vtk(
            numpy.ones(vtkDataObject.GetNumberOfPoints()) * valeur, deep = 1)
        varray.SetName(nom_constante)
        output.GetPointData().AddArray(varray)
    return output
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def supprimer_array(vtkDataObject, nom_array, loc='points'):
    """fonction qui supprime un array
    """
    if isinstance(vtkDataObject, vtk.vtkMultiBlockDataSet):
        output = vtk.vtkMultiBlockDataSet()
        for numbloc in get_numeros_blocs_non_vides(vtkDataObject):
            output.SetBlock(numbloc,
                supprimer_array(vtkDataObject.GetBlock(numbloc), nom_array, loc)
                )
    else:
        output = vtk_new_shallowcopy(vtkDataObject)
        if loc == 'points':
            output.GetPointData().RemoveArray(nom_array)
        else:
            output.GetCellData().RemoveArray(nom_array)
    return output
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def ajouter_numpy_array_as_vtk_array(input, numpy_array, nom):
    """fonction qui ajoute un numpy array 
    input peut etre un multibloc
    il faut alors que l'array concatene soit dans le meme ordre que ce qu'il etait 
    en sortie de get_vtk_array_as_numpy_array
    
    typiquement, shape(nk, nj, ni),ravel() pour un multibloc
    
    dans le cas d'une variables vectorielle, le shape doit alors etre
        (nb_points, nb_composantes)
        (nb_cellules, nb_composantes)
    le numpy_array peut aussi etre donne en ligne (operation numpy.ravel), les composantes dans l'ordre
        [vx, vy, vz, vx, vy, ...]
    """
    #conversion en numpy.float64 pour eviter des probleme ensuite avec numpy_to_vtk
    numpy_array = numpy.asarray(numpy_array, dtype=numpy.float64)
    #on determine en fonction de la longueur du numpy array si la donnee doit aller aux points ou aux cellules
    
    #multiblock n'a pas de fonction GetNumberOfCells, alors on somme a la main, en supposant
    #que c'est un multibloc a une seule profondeur
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        nombre_cellules = sum(
            [input.GetBlock(numbloc).GetNumberOfCells() 
                for numbloc in get_numeros_blocs_non_vides(input)]
            )
    else:
        nombre_cellules = input.GetNumberOfCells()
    
    if numpy_array.size == input.GetNumberOfPoints():
        loc = 'points'
        mode_vecteur = 0
    elif numpy_array.size == nombre_cellules:
        loc = 'cells'
        mode_vecteur = 0
    elif numpy_array.size == input.GetNumberOfPoints() * 3:
        loc = 'points'
        mode_vecteur = 1
    elif numpy_array.size == nombre_cellules * 3:
        loc = 'cells'
        mode_vecteur = 1
    else:
        raise Exception, "taille du numpy array ne correspond ni au nombre de points ni au nombre de cellules"
    
    if mode_vecteur == 1:
        numpy_array = numpy_array.reshape(-1, 3)
    
    # print 'nb cells', nombre_cellules
    # print 'narray size', numpy_array.size
    # print 'loc', loc
    # print 'mode_vecteur', mode_vecteur
    
    #print "Ajout de l'array numpy"
    #print 'localisation ', loc
    #print 'mode vecteur ', mode_vecteur
    
    #on cree output, meme instance que input
    output = vtk_new_shallowcopy(input)
    
    #cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        num_blocs = get_numeros_blocs_non_vides(input)
        if loc == "points":
            liste_nb_par_bloc = [input.GetBlock(numbloc).GetNumberOfPoints() 
                for numbloc in num_blocs]
        elif loc == "cells":
            liste_nb_par_bloc = [input.GetBlock(numbloc).GetNumberOfCells() 
                for numbloc in num_blocs]
        
        liste_arg = [
            numpy.arange(
                int(numpy.sum(liste_nb_par_bloc[:rang_bloc])), 
                int(numpy.sum(liste_nb_par_bloc[:rang_bloc]) + liste_nb_par_bloc[rang_bloc])
                )
                for rang_bloc in range(len(num_blocs))
                ]
        
        for rang_bloc in range(len(num_blocs)):
            narray = numpy_array[liste_arg[rang_bloc]]
            numbloc = num_blocs[rang_bloc]
            output.SetBlock(numbloc, 
                ajouter_numpy_array_as_vtk_array(
                    output.GetBlock(numbloc), 
                    narray, 
                    nom)
                )
    
    #monobloc - polydata, structured ou autre
    else:
        numpy_array = numpy.ascontiguousarray(numpy_array)
        varray = numpy_support.numpy_to_vtk(numpy_array, deep = 1)
        varray.SetName(nom)
        if loc == 'points':
            # if len(numpy_array.shape) == 2 and numpy_array.shape[1] == 3:
                # output.GetPointData().SetVectors(varray)
                # output.GetPointData().SetActiveVectors(nom)
            # else:
                # output.GetPointData().AddArray(varray)
            output.GetPointData().AddArray(varray)
            output.GetPointData().Update()
        elif loc == "cells":
            # if len(numpy_array.shape) == 2 and numpy_array.shape[1] == 3:
                # output.GetCellData().SetVectors(varray)
                # output.GetCellData().SetActiveVectors(nom)
            # else:
                # output.GetCellData().AddArray(varray)
            output.GetCellData().AddArray(varray)
            output.GetCellData().Update()
    return output
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def VTKThreshold(input, nom_array, valeur_min=None, valeur_max=None, loc='points', UseContinuousCellRange=False):
    """Fonction Threshold avec nom_array
    multi ou monobloc
    indiquer seulement valeur_min pour threshold by min
    indiquer seulement valeur_max pour threshold by max
    indiquer les deux pour threshold between
    
    loc = 'points', sinon on utilise les donnees aux cellules
    
    amarsan complement 03 novembre 2017
    modification de la fonction pour pouvoir clipper avec une coordonnee
    'coordsX', 'coordsY', 'coordsZ'
    """
    # dans le cas ou on veut Threshold avce une coordonnee spatiale, on ajoute d'abord cette coordonnee spatiale aux points
    if nom_array[:6] == 'coords':
        print 'ajout de la coordonnee {0} aux points'.format(nom_array)
        coords = get_vtk_array_as_numpy_array(input, 'coords')
        input = ajouter_numpy_array_as_vtk_array(input, 
            coords[:, {'coordsX': 0, 'coordsY': 1, 'coordsZ': 2}[nom_array]], 
            nom_array)
    
    #cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, VTKThreshold(input.GetBlock(numbloc), 
                nom_array, valeur_min, valeur_max, loc)
                )
    #cas monobloc
    else:
        if nom_array in get_noms_arrays_presents(input, loc = 'points'):
            input.GetPointData().SetActiveScalars(nom_array)
            input.GetCellData().SetActiveScalars(None)
        elif nom_array in get_noms_arrays_presents(input, loc = 'cells'):
            input.GetCellData().SetActiveScalars(nom_array)
            input.GetPointData().SetActiveScalars(None)
        else:
            raise IOError, "{0} n'est pas present dans input".format(nom_array)
        select = vtk.vtkThreshold()
        select.AllScalarsOff()
        select.UseContinuousCellRangeOn() if UseContinuousCellRange else select.UseContinuousCellRangeOff()
        vtk_set_input(select, input)
        if valeur_min != None and valeur_max != None:
            select.ThresholdBetween(valeur_min, valeur_max)
        elif valeur_min != None:
            select.ThresholdByUpper(valeur_min)
        elif valeur_max != None:
            select.ThresholdByLower(valeur_max)
        else:
            raise IOError, "indiquez valeur_min ou valeur_max"
        select.Update()
        
        output = select.GetOutput()
        
    return output
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def VTKBoxClip(input, xmin, xmax, ymin, ymax, zmin, zmax):
    """ clip avec une boite. 
    """
    print 'Fonction extremement lente, je ne sais pas pourquoi. (amarsan 3 novembre 2017)'
    if 0:
        clipper =  vtk.vtkBoxClipDataSet()
        vtk_set_input(clipper, input)
        clipper.SetBoxClip(xmin, xmax, ymin, ymax, zmin, zmax)
        clipper.Update()
        return clipper.GetClippedOutput()
    return None
#__________________________________________________________________________________________


#__________________________________________________________________________________________
def VTKSubset(input, i_gardes = None, j_gardes = None, k_gardes = None):
    """Fonction Subset avec les indices i, j, k
    multi ou monobloc, ne marche que sur un vtkStructuredGrid
    les valeurs des indices sont celles qu'on mettrait dans Paraview
    (donc attention en particulier au decalage de 1 par rapport a elsA !!)
    
    indice_gardes = [indice_min, indice_max], 3 paires d'indices a donner
    par defaut, si on ne donne pas d'info sur un indice, pas de decoupage dans cette direction
    """
    #cas multibloc
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, VTKThreshold(input.GetBlock(numbloc), 
                nom_array, valeur_min, valeur_max)
                )
    #cas monobloc
    else:
        dim_gardees = list(input.GetExtent())
        if i_gardes != None:
            dim_gardees[0], dim_gardees[1] = i_gardes[0], i_gardes[1]
        if j_gardes != None:
            dim_gardees[2], dim_gardees[3] = j_gardes[0], j_gardes[1]
        if k_gardes != None:
            dim_gardees[4], dim_gardees[5] = k_gardes[0], k_gardes[1]
        extract = vtkFiltersExtraction.vtkExtractGrid()
        vtk_set_input(extract, input)
        extract.SetVOI(dim_gardees)
        extract.Update()
        output = extract.GetOutput() 
    return output
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def modifier_axe(objetVTK, nv_repere):
    """fonction qui modifie les axes (x,y,z) de l'objet VTK donne en entree
    sans modifier la facon dont il est indexe (i,j,k)
    Retourne un objet VTK avec le nouveau repere souhaite
    
    objetVTK peut etre mono-bloc ou multi-blocs
    
    Exemple d'utilisation : on veut intervertir les axes x et z
    dans tous les cas, l'ancien repere est represente par (x,y,z) = (1,2,3)
    On exprime le nouveau repere en fonction de l'ancien :
    x => z , y => y , z => x    :   cela se traduit par nv_repere = (3,2,1)
    
    On peut aussi changer le sens d'un axe. Exemple :
    Pour x => y , y => -z , z => -x : on donnera nv_repere = (2,-3,-1)
    
    FONCTION A COMPLETER ULTERIEUREMENT :
    ca marche actuellement pour un maillage, il faudra par la suite
    l'etendre a d'autres types de donnees comme des champs aeros, par exemple
    """
    
    if isinstance(objetVTK, vtk.vtkMultiBlockDataSet):
        n = vtk.vtkMultiBlockDataSet()
        
        for numbloc in get_numeros_blocs_non_vides(objetVTK):
            n.SetBlock(numbloc,
                modifier_axe(objetVTK.GetBlock(numbloc), nv_repere)
                )
    
        return n
    
    else:
        dims = objetVTK.GetDimensions()
        coords = get_vtk_array_as_numpy_array(objetVTK, "coords")
        coords = coords.reshape(dims[::-1] + (3,))
        coords_new = numpy.concatenate(
        [
        coords[:, :, :, abs(nv_repere[0]) - 1, None] * numpy.sign(nv_repere[0]),
        coords[:, :, :, abs(nv_repere[1]) - 1, None] * numpy.sign(nv_repere[1]),
        coords[:, :, :, abs(nv_repere[2]) - 1, None] * numpy.sign(nv_repere[2])
        ], axis = 3
        )
        coords_new = numpy.array(coords_new)
        coords_new = coords_new.reshape(
        coords_new.size / 3, 3)

        vtk_array = numpy_support.numpy_to_vtk(coords_new, deep=1)

        point = vtk.vtkPoints()
        point.SetData(vtk_array)

        objetVTK = vtk_new_instance(objetVTK)
        objetVTK.SetPoints(point)
        objetVTK.SetDimensions(dims)
    
        return objetVTK
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def get_noms_blocs(multibloc):
    """Fonction qui retourne les noms des blocs d'un multibloc
    sous forme de liste. 
    La position dans la liste correspond au numero du bloc
    """
    liste = []
    for numbloc in range(multibloc.GetNumberOfBlocks()):
        nom_bloc = multibloc.GetMetaData(numbloc).Get(vtk.vtkCompositeDataSet.NAME())
        if nom_bloc is not None:
            liste.append(nom_bloc.replace('\n', '').replace(' ', ''))
        else:
            liste.append(nom_bloc)
    return liste
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def get_bloc_par_nom(multibloc, nom_bloc):
    """Fonction qui retourne les noms des blocs d'un multibloc
    sous forme de liste. 
    La position dans la liste correspond au numero du bloc
    """
    liste_noms = get_noms_blocs(multibloc)
    if not nom_bloc in liste_noms:
        raise Exception, 'Aucun bloc ne porte ce nom. Les noms sont ' + str(liste_noms)
    return multibloc.GetBlock(liste_noms.index(nom_bloc))
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def get_numero_bloc_par_nom(multibloc, nom_bloc):
    """Fonction qui retourne les noms des blocs d'un multibloc
    sous forme de liste. 
    La position dans la liste correspond au numero du bloc
    """
    liste_noms = get_noms_blocs(multibloc)
    if not nom_bloc in liste_noms:
        raise Exception, 'Aucun bloc ne porte ce nom. Les noms sont ' + str(liste_noms)
    return liste_noms.index(nom_bloc)
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def create_bloc_structure_from_numpy_array(coords):
    """Fonction qui retourne un bloc structure a partir d'un numpy array qui contient les coordonnees des points. 
    
    Le numpy array doit avoir les dimensions suivantes (ni, nj, nk, 3)
    """
    #creation du bloc
    bloc = vtk.vtkStructuredGrid()
    #conversion du numpy array en vtk array
    vtkArray = numpy_support.numpy_to_vtk(numpy.ascontiguousarray(
        coords.transpose(2, 1, 0, 3).reshape(coords.size / 3, 3)), deep = 1)
    #creation des points
    points = vtk.vtkPoints()
    points.SetData(vtkArray)
    #remplissage du bloc
    bloc.SetDimensions(coords.shape[0], coords.shape[1], coords.shape[2])
    bloc.SetPoints(points)
    return bloc
#_____________________________________________________________________

#__________________________________________________________________________________________
def create_bloc_non_structure_from_numpy_array(coords, cells, cellstypes, cellslocations):
    """Fonction qui retourne un bloc non-structure a partir de numpy array
    Les entrees a donner sont :
        - coords :          un numpy_array de dimension (nb_points, 3) qui contient les coordonnees des points.
        - cells :           l'arbre de connectivite, qui definit les cellules. 
                            (nb_points_dans_cellule_1, indice_point_1, indice_point_2, indice_point_3, ..., nb_points_dans_cellule2, indice_point_, indice_point_, indice_point_, ...)
        - cellstypes :      qui definit le type des cellules, une par une
        - cellslocations :  qui definit la position de la definition de la cellules dans l'array <cells>
    """
    #creation du bloc
    bloc = vtk.vtkUnstructuredGrid()
    
    #conversion du numpy array en vtk array
    vtkArray = numpy_support.numpy_to_vtk(numpy.ascontiguousarray(coords), deep = 1)
    
    #creation des points
    points = vtk.vtkPoints()
    points.SetData(vtkArray)
    
    #remplissage du bloc
    bloc.SetPoints(points)
    
    # definition des cellules
    vtkCells = vtk.vtkCellArray()
    vtkCells.SetCells(cellstypes.size, numpy_support.numpy_to_vtk(cells, deep = 1, array_type = vtk.vtkIdTypeArray().GetDataType()))
    
    bloc.SetCells(
        numpy_support.numpy_to_vtk(cellstypes, deep = 1, array_type = vtk.vtkUnsignedCharArray().GetDataType()), 
        numpy_support.numpy_to_vtk(cellslocations, deep = 1, array_type = vtk.vtkIdTypeArray().GetDataType()), 
        vtkCells
        )
    # bloc.Update()
    # bloc.UpdateData()
    
    return bloc
#_____________________________________________________________________

#__________________________________________________________________________________________
def create_polydata_from_numpy_array(coords, polys=None, nb_polys=None, isline=False):
    """Fonction qui retourne une surface polydata a partir de numpy array
    Les entrees a donner sont :
        - coords :          un numpy_array de dimension (nb_points, 3) qui contient les coordonnees des points.
        - polys :           les polydata = connectivite 
        - nb_polys :        le nombre de cellules
        - isline :          dans le cas ou c'est une ligne qu'on veut creer, et non pas une surface. 
    """
    #creation du bloc
    bloc = vtk.vtkPolyData()
    
    #conversion du numpy array en vtk array
    vtkArray = numpy_support.numpy_to_vtk(numpy.ascontiguousarray(coords), deep = 1)
    
    #creation des points
    points = vtk.vtkPoints()
    points.SetData(vtkArray)
    
    #remplissage du bloc
    bloc.SetPoints(points)
    
    if polys is not None:
        if nb_polys is None:
            raise IOError, 'Indiquer le nombre total de facettes'
        # definition des polys
        vtkCells = vtk.vtkCellArray()
        vtkCells.SetCells(int(nb_polys), numpy_support.numpy_to_vtk(
            polys.ravel(), deep = 1, array_type = vtk.vtkIdTypeArray().GetDataType()
            ))
        
        bloc.SetPolys(vtkCells) if isline is False else bloc.SetLines(vtkCells)
        
    return bloc
#_____________________________________________________________________

#_____________________________________________________________________
def extraire_surface(input, region_to_extract = None, angle_split=None):
    """fonction qui extrait les surfaces
    
    region_to_extract doit etre une liste d'entiers.
    
    si angle_split = None
    """
    
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, extraire_surface(input.GetBlock(numbloc), region_to_extract, angle_split))
        return output
    
    # extraction des surfaces exterieures
    f = vtkFiltersGeometry.vtkGeometryFilter()
    vtk_set_input(f, input)
    f.Update()
    surface = f.GetOutput()
    # print surface

    # calcul des normales
    f = vtk.vtkPolyDataNormals()
    if angle_split is not None:
        f.SplittingOn()
        f.SetFeatureAngle(angle_split)
    
    vtk_set_input(f, surface)
    f.Update()
    surface = f.GetOutput()
    
    # calcul de l'indice de connectivite
    f = vtk.vtkConnectivityFilter()
    vtk_set_input(f, surface)
    if region_to_extract is not None:
        f.SetExtractionModeToSpecifiedRegions()
        for id in region_to_extract:
            f.AddSpecifiedRegion(id)
    else:
        f.SetExtractionModeToAllRegions()
    f.ColorRegionsOn()
    f.Update()
    surface = f.GetOutput()

    # on passe ce filtre, parce que sinon les cellules sont supprimees, mais il reste les points dans le dataset
    if region_to_extract is not None:
        f = vtkFiltersGeometry.vtkDataSetSurfaceFilter()
        vtk_set_input(f, surface)
        f.Update()
        surface = f.GetOutput()
    
    return surface
#_____________________________________________________________________

#_____________________________________________________________________
def ajouter_types_cellules(input):
    """fonction qui ajoute le type de cellules aux cellules
    
    """
    # cas multibloc - appel recursif
    if isinstance(input, vtk.vtkMultiBlockDataSet):
        output = vtk_new_instance(input)
        for numbloc in get_numeros_blocs_non_vides(input):
            output.SetBlock(numbloc, ajouter_types_cellules(input.GetBlock(numbloc)))
    # cas monobloc
    else:
        cellstypes = get_vtk_array_as_numpy_array(input, 'cellstypes')
        output = ajouter_numpy_array_as_vtk_array(input, cellstypes, 'CellTypes')
    return output
#_____________________________________________________________________

#_____________________________________________________________________
def extraire_blocs(input, liste_numeros):
    """fonction d'extraction des blocs d'un multibloc
    """
    m = vtk.vtkMultiBlockDataSet()
    for numbloc in liste_numeros:
        m.SetBlock(numbloc, input.GetBlock(numbloc))
    return m
#_____________________________________________________________________

