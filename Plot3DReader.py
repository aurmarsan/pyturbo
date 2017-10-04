import os,sys,re
import cStringIO
import numpy as np
from numpy import rec

# Attention, pour les fichiers binaire c'est un lecteur pour des fichiers ecrits avec du fortran.
# Si ecrit avec du python resultat pas garanti car l'ecriture se fait avec des routines C++

### IMPORTANT! pour verifier l'endianness sys.byteorder() => 'little' petit endian
###                                                    => 'big' big endian

def readHeader(fileName,binary):
    """
    Read the header of PLOT3D files
    input  : filename                         (fileName, string)
             type of file, binary or formatted (binary, boolean)
    output : mesh size                       (dims, tuple of integers)
             variables name                 (vars, tuple of strings)
    """
    
    if binary:
        header = readBinHeader(fileName)
    else:
        header = readFmtHeader(fileName)
    return header


## def readBinHeader(fileName):
##   """
##   lecture du header du fichier PLOT3D binaire
##   entree : nom du fichier (fileName, str)
##   sortie : dimensions de la grille (dims, tuple)
##            nombre de variables    (nbVars, int)
##            type de champ maillage ou aero (isMesh, Booleen)
##   """
##   fichier = open(fileName,"rb")
##   s = struct.Struct(">4s3I")
##   firstline = fichier.read(4+3*4)
##   firstline = s.unpack(firstline)
##   deb = firstline[0]
##   dims = [int(dim) for dim in firstline[1:4]]
##   s = struct.Struct(">4s")
##   firstline = fichier.read(4)
##   firstline = s.unpack(firstline)
##   if deb == firstline[-1]:
##       nbVars = 3
##       fichier.close()
##       isMesh = True
##   else : 
##       fichier.seek(-4,1)
##       s = struct.Struct(">I4s")
##       firstline = fichier.read(4+4)
##       firstline = s.unpack(firstline)
##       fichier.close()
##       assert(deb==firstline[-1])
##       nbVars = firstline[0]
##       isMesh = False
##   print "Dimensions of the mesh:",dims
##   print "Number of variables:",nbVars
##   return nbVars,tuple(dims),isMesh



## def unpack(str,strStruct,endian):
##   """
##   Function to unpack binary data
##   inlet  : binary string to unpack
##   outlet : 
##   """
##   data = rec.fromstring(fichier[start_i:end_i],names='beg, data, end',formats='4S, 3I, 4S',byteorder=endian)


def readBinHeader(fileName,endian='big'):
    """
    reads the header of a binary PLOT3D file using the fromstring function of numpy
    input  : filename                  (fileName, str)
             endianness              (endian, 'big'/'small', default value : 'big')
    output : mesh size                (dims, tuple of integers)
             number of variables        (nbVars, int)
             type of field mesh or aero (isMesh, boolean)
    """
    
    try:
        fichier = open(fileName,"rb")
        end_i = 2*np.dtype('S4').itemsize + 3*np.dtype('int32').itemsize # taille de la chaine de 4 caracteres, 3 entiers, 4 caracteres
        data = fichier.read(end_i)
        fichier.close()
        data = rec.fromstring(data,names='beg, dims, end',formats='4S, 3I, 4S',byteorder=endian) # chaine de 4 caracteres, 3 entiers, 4 caracteres, encode selon le type 'big' ou 'small' endian
        assert(data[0]['beg']==data[0]['end'])   ### creer un raise specifique?
        nbVars = 3
        isMesh = True
        return nbVars,tuple(data[0]['dims']),isMesh
    
    except AssertionError:
        fichier = open(fileName,"rb")
        end_i = 2*np.dtype('S4').itemsize + 4*np.dtype('int32').itemsize
        data = fichier.read(end_i)  # chaine de 4 caracteres, 4 entiers, 4 caracteres
        fichier.close()
        data = rec.fromstring(data,names='beg, dims, nbVars, end',formats='4S, 3I, I, 4S',byteorder=endian) # chaine de 4 caracteres, 4 entiers, 4 caracteres, encode selon le type 'big' ou 'small' endian
        try:
            assert(data[0]['beg']==data[0]['end'])
        except AssertionError:
            print "Error : the file is neither a binary PLOT 3D mesh file nor a binary PLOT 3D aero file. Check the file."
        isMesh = False
        return data[0]['nbVars'],tuple(data[0]['dims']),isMesh
            

def readFmtHeader(fileName):
    """
    reads the header of the formatted PLOT3D file
    input  : filename           (fileName, string)
    output : mesh size         (dims, tuple of integers) 
             number of variables (nbVars, integer)
    """
    fichier = open(fileName,"r")
    firstline = fichier.readline().split()
    dims = [int(dim) for dim in firstline[:3]]
    try:
        nbVars = int(firstline[3])
    except IndexError:
        nbVars = 3
    print "Dimensions of the mesh:",dims
    print "Number of variables:",nbVars
    return nbVars,tuple(dims)


def readData(fileName,header,binary):
    """
    Lecture des donnees du fichier PLOT3D
    entree : nom du fichier (fileName, str)
             binaire ou non (binary, bool)
    sortie : data (data, dictionnaire de numpy array)
    """
    if binary:
        data = readBinData(fileName,header)
    else:
        data = readFmtData(fileName,header)
    return data


## def readBinData(fileName,header):
##   """
##   Lecture effective des donnees du fichier binaire PLOT3D
##   entree : nom du fichier      (fileName, str)
##            dimensions          (dims, tuple)
##            nombre de variables (nbVars, int)
##   sortie : data (data, dictionnaire de numpy array))
##   """
##   nbVars,dims,isMesh = header[0],header[1],header[2]
##   fichier = open(fileName,"rb")
##   if (isMesh) : fichier.seek(20,0)
##   else: fichier.seek(24,0)
##   s = struct.Struct(">4x%dd4x"%(nbVars*dims[0]*dims[1]*dims[2]))            
##   dataArray = fichier.read(nbVars*dims[0]*dims[1]*dims[2]*8+2*4)
##   dataArray = s.unpack(dataArray)
##   fichier.close()
##   # mise des dimensions dans le bon ordre et mise des donnees dans un dictionnaire
##   dataArray = np.array(dataArray)
##   dataArray.shape = nbVars,dims[0]*dims[1]*dims[2]

##   data = {}
##   for l in xrange(nbVars):
##       var = "var%d"%l
##       data[var] = dataArray[l,:]
##       assert(dims[0]*dims[1]*dims[2] == len(data[var]))
##       data[var].shape = dims[2],dims[1],dims[0]
##       data[var] = np.swapaxes(data[var],0,2)

##   return data

def readBinData(fileName,header,endian='big'):
    """
    Effective reading of the binary PLOT3D data file using the fromfile function of numpy
    input  : filename     (fileName, str)
             header's data (tuple containing the number of variables, the mesh size and the type of file)
             endianness (endian, 'big'/'small', default value : 'big')
    output : data         (data, dictionnary of numpy arrays of the size of the mesh)
    """
    nbVars,dims,isMesh = header[0],header[1],header[2]
    if (isMesh):
        dataArray = rec.fromfile(fileName,names='beg, data, end', formats='4S,%dd,4S'%(nbVars*dims[0]*dims[1]*dims[2]),offset=20,byteorder=endian)
        assert(dataArray[0]['beg']== dataArray[0]['end'])
    else:
        dataArray = rec.fromfile(fileName,names='beg, data, end', formats='4S,%dd,4S'%(nbVars*dims[0]*dims[1]*dims[2]),offset=24,byteorder=endian)
        assert(dataArray[0]['beg']== dataArray[0]['end'])
    dataArray = dataArray[0]['data']
    dataArray.shape = nbVars,dims[0]*dims[1]*dims[2]
    data = {}
    for l in xrange(nbVars):
        var = "var%d"%l
        data[var] = dataArray[l,:]
        assert(dims[0]*dims[1]*dims[2] == len(data[var]))
        data[var].shape = dims[2],dims[1],dims[0]
        data[var] = np.swapaxes(data[var],0,2)

    return data 
    

def readFmtData(fileName,header):
    """
    Effective reading of the ASCII PLOT3D data file using the loadtxt function of numpy
    input  : filename (fileName, string)
             dimensions of the mesh (dims, tuple of integers)
             number of variables (nbVars, integer)
    output : data (data, dictionnary of numpy array)
    """
    nbVars,dims = header[0],header[1]
    fichier = open(fileName,"r")
    data=fichier.readlines()[1:]
    fichier.close()
    dataArray = np.loadtxt(cStringIO.StringIO(''.join(data).replace('\n','')),dtype='float64')
    del data
    
    dataArray.shape = nbVars,dims[0]*dims[1]*dims[2]
    # mise des dimensions dans le bon ordre et mise des donnees dans un dictionnaire
    data = {}
    for l in xrange(nbVars):
        var = "var%d"%l
        data[var] = dataArray[l,:]
        assert(dims[0]*dims[1]*dims[2] == len(data[var]))
        data[var].shape = dims[2],dims[1],dims[0]
        data[var] = np.swapaxes(data[var],0,2)
    return data


def plt3dReader(fileName,binary):
    """
    reader for the PLOT3D file format
    """
    assert(isinstance(fileName,str))
    assert(os.path.isfile(fileName))

    return readData(fileName,readHeader(fileName,binary),binary)


if __name__ == "__main__":
    fileName = sys.argv[1]
    dict = plt3dReader(fileName,binary=True)
