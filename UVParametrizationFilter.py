from math import atan2
from math import cos
from math import pi
from math import sin

from numpy import asarray
from numpy import atleast_1d
from numpy import atleast_2d
from numpy import cross
from numpy import diff
from numpy import dot
from numpy import linspace
from numpy import loadtxt
from numpy import r_
from numpy import squeeze
from numpy import sqrt
from numpy import transpose
from numpy import zeros
from numpy import zeros_like
from numpy.linalg import norm
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve

try: from paraview import vtk 
except: import vtk

try: 
    from paraview.numpy_support import vtk_to_numpy
    from paraview.numpy_support import numpy_to_vtk
    from paraview.vtk import vtkFloatArray
    from paraview.vtk import vtkPoints
    from paraview.vtk import vtkStructuredGrid
except: 
    from vtk.util.numpy_support import vtk_to_numpy
    from vtk.util.numpy_support import numpy_to_vtk
    from vtk import vtkFloatArray
    from vtk import vtkPoints
    from vtk import vtkStructuredGrid

PRECISION = 1.0e-12


""" ***** IMPORTANT *****
   - ajuster le parametre de smoothing dans parametricunivariatespline
   - ajuster egalement le pas pour le calcul de hsh
   - ajuster axe pour l'axe de rotation. hsH est calcule dans la direction orthogonale
"""

class ParametricUnivariateSpline:
    """une classe pour interpoler une spline, en fonction d'un parametre u. 
    Celui-ci peut etre soit indique. Si il est None, alors c'est l'abscisse meridienne qui est calculee. 
    x doit etre un array (nb_coords, n) ou n est le nombre de points
    
    Quand appelee, comme une fonction, retourne un array (nb_coords, n)
    ou n est alors le nombre de points demandes (le nombre de u)
    
    s est a ajuster 
    """
    def __init__(self, x, u=None, k=3, s=1e-8):
        if u is None:
            du = sqrt((diff(x, axis=1) ** 2).sum(axis=0))
            u = r_[0.0, du.cumsum() / du.sum()]
        self._degree = k
        self._param = u
        self._splines = [UnivariateSpline(u, p, k=k, s=s) for p in x]
        self._smoothing = s
    
    def __call__(self, u, nu=0):
        u = atleast_1d(u)
        pmin = self._param.min()
        pmax = self._param.max()
        
        if nu > 0:
            u = u.clip(pmin, pmax)
            return transpose([spl(u, nu=nu) for spl in self._splines])
        
        xmin = asarray([spl(pmin) for spl in self._splines])
        xmax = asarray([spl(pmax) for spl in self._splines])
        dmin = self.__call__(pmin, nu=1)
        dmax = self.__call__(pmax, nu=1)
        
        #amarsan
        xmin = xmin.ravel()
        xmax = xmax.ravel()
        dmin = dmin.ravel()
        dmax = dmax.ravel()
        
        res = zeros([2] + list(u.shape))
        for idx, v in enumerate(u):
            if v < pmin:
                #res[:, idx] = xmin + (v - pmin) * dmin
                #amarsan
                res[:, idx] = (xmin + (v - pmin) * dmin).squeeze()
            elif v > pmax:
                res[:, idx] = (xmax + (v - pmax) * dmax)
            else:
                res[:, idx] = [float(spl(v)) for spl in self._splines]
        return res.squeeze()
    
    def length(self):
        der1 = self.__call__(self._param, nu=1)
        dl = sqrt((der1 ** 2).sum(axis=1))
        spl = UnivariateSpline(self._param, dl, k=self._degree, s=self._smoothing)
        return spl.integral(self._param.min(), self._param.max())
    
    def curvature(self, u):
        der1 = atleast_2d(self.__call__(u, nu=1))
        der2 = atleast_2d(self.__call__(u, nu=2))
        der1n = asarray(map(norm, der1)).clip(min=PRECISION)
        crv = asarray([abs(cross(d1, d2)) for d1, d2 in zip(der1, der2)])
        return squeeze(crv / der1n) / self.length()
    
    def curvature_radius(self, u):
        return 1.0 / self.curvature(u).clip(min=PRECISION)


class ParametricBivariateSpline:
    def __init__(self, x, u=None, v=None, ku=3, kv=3, s=1e-8):
        if u is None:
            du = sqrt((diff(x, axis=1) ** 2).sum(axis=0)).mean(axis=1)
            u = r_[0.0, du.cumsum() / du.sum()]
        if v is None:
            dv = sqrt((diff(x, axis=2) ** 2).sum(axis=0)).mean(axis=0)
            v = r_[0.0, dv.cumsum() / dv.sum()]
        self._degree = (ku, kv)
        self._splines = [RectBivariateSpline(u, v, p, kx=ku, ky=kv, s=s) for p in x]
        
    def __call__(self, u, v):
        return squeeze([spl(u, v) for spl in self._splines])


def CreateSpline(hubFileName, tipFileName, stepSize=1.0, relativeExtension=0.1,
        recherche_depuis_carter_spline=1):
    """retourne une fonction qui a u,v associe les points moyeu et carter
    prolongee d'une certaine valeur <relativeExtension>. 
    
    Si recherche_depuis_carter_spline == 1 alors recherche de la ligne centrale en utilisant 
    la ligne carter, au lieu de la ligne moyeu. 
    """
    hub = loadtxt(hubFileName)
    tip = loadtxt(tipFileName)
    
    spl_hub = ParametricUnivariateSpline(hub.T)
    spl_tip = ParametricUnivariateSpline(tip.T)
    
    numberOfPoints = 1 + int(max([spl_hub.length(), spl_tip.length()]) // stepSize)
    print "nombre de points le long de l'abscisse meridienne : ", numberOfPoints
    
    v = 0.0
    R = 0.5 * norm(spl_hub(v) - spl_tip(v))
    centers = []
    liste_theta = []
    for uh in linspace(0.0, 1.0, numberOfPoints):
        if recherche_depuis_carter_spline:
            P = spl_tip(uh)
            t = spl_tip(uh, nu=1)
        else:
            P = spl_hub(uh)
            t = spl_hub(uh, nu=1)
        
        #theta = atan2(-t[0], t[1])
        #amarsan
        theta = atan2(-t[:, 0], t[:, 1]) #atan(y, x) = atan y / x
        
        liste_theta.append(theta)
        
        def origin(R):
            O = (P[0] - R * cos(theta),
                 P[1] - R * sin(theta))
            return O
            
        def err(p):
            v, R = p
            O = origin(R)
            Q = spl_tip(v)
            T = spl_tip(v, nu=1)
            #return norm(Q - O) - R, dot(Q - O, T)
            #amarsan
            return norm(Q - O) - R, dot(Q - O, T.squeeze())
        if norm(err([v, R])) > PRECISION:
            v, R = fsolve(err, [v, R])
        centers.append(origin(R))
    
    #from matplotlib import pyplot as plt
    #import sys, numpy
    #data = transpose(centers)
    #plt.plot(data[0], data[1])
    #data_hub = spl_hub(linspace(0.0, 1.0, numberOfPoints))
    #print data_hub.shape
    #plt.plot(data_hub[0, :], data_hub[1, :])
    #data_tip = spl_tip(linspace(0.0, 1.0, numberOfPoints))
    #plt.plot(data_tip[0, :], data_tip[1, :])
    #
    #plt.figure()
    #plt.plot(numpy.cos(liste_theta))
    #plt.plot(numpy.sin(liste_theta))
    #
    #plt.show()
    #sys.exit()
    
    spl_mid = ParametricUnivariateSpline(transpose(centers))
    
    uvMin = -relativeExtension
    uvMax = 1.0 + relativeExtension
    spl_u = linspace(uvMin, uvMax, numberOfPoints)
    spl_v = (uvMin, uvMax)
    
    u = 0.0
    v = 0.0
    points = []
    for um in spl_u:
        O = spl_mid(um)
        t = spl_mid(um, nu=1)
        
        def err(p):
            u, v = p
            P = spl_hub(u)
            Q = spl_tip(v)
            #return dot(P - O, t), dot(Q - O, t))
            #amarsan
            return dot(P - O, t.squeeze()), dot(Q - O, t.squeeze())
        if norm(err([u, v])) > PRECISION:
            u, v = fsolve(err, [u, v])
        pinf = spl_hub(u)
        psup = spl_tip(v)
        #~ dp = relativeExtension * (psup - pinf) / norm(psup - pinf)
        dp = relativeExtension * (psup - pinf)
        points.append([pinf - dp, psup + dp])
    
    return ParametricBivariateSpline(transpose(points), u=spl_v, v=spl_u,
                                     ku=1), numberOfPoints


def parametrize(spline, numberOfPoints=101, relativeExtension=0.1):
    """retourne le multibloc structure qui quadrille le domaine
    pour pouvoir ensuite faire le probe dans le plan xr
    """
    arrayU = vtkFloatArray()
    arrayV = vtkFloatArray()
    
    arrayU.SetName('hsH')
    arrayV.SetName('xm')
    
    points = vtkPoints()
    uvMin = -relativeExtension
    uvMax = 1.0 + relativeExtension
    
    for _v in linspace(uvMin, uvMax, numberOfPoints):
        for _u in (uvMin, uvMax):
            arrayU.InsertNextTuple1(_u)
            arrayV.InsertNextTuple1(_v)
            
            pnt = spline(_u, _v)
            points.InsertNextPoint(pnt[0], pnt[1], 0.0)
    
    grid = vtkStructuredGrid()
    grid.SetExtent(1, 2, 1, numberOfPoints, 1, 1)
    grid.SetPoints(points)
    
    pointData = grid.GetPointData()
    pointData.AddArray(arrayU)
    pointData.AddArray(arrayV)
    return grid


def parametrize_dataset(dataset, grid, axe):
    """L'interpolation est fait dans un plan xrt
    axe indique l'axe de rotation du maillage initial
    axe = 0 pour x
    axe = 1 pour y
    axe = 0 pour z
    """
    xyz = vtk_to_numpy(dataset.GetPoints().GetData())
    
    from numpy import sqrt
    from numpy import zeros_like
    xrt = zeros_like(xyz)
    if axe == 0:
        xrt[:, 0] = xyz[:, 0]
        xrt[:, 1] = sqrt(xyz[:, 1] ** 2 + xyz[:, 2] ** 2)
    elif axe == 1:
        xrt[:, 0] = xyz[:, 1]
        xrt[:, 1] = sqrt(xyz[:, 0] ** 2 + xyz[:, 2] ** 2)
    elif axe == 2:
        xrt[:, 0] = xyz[:, 2]
        xrt[:, 1] = sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
    else:
        raise IOError, "valeur non valide pour axe"
        
    
    points = vtkPoints()
    points.SetData(numpy_to_vtk(xrt, deep=1))
    
    input = vtk.vtkPolyData()
    input.SetPoints(points)
    
    probe = vtk.vtkProbeFilter()
    try:
        probe.SetInputData(input)
        probe.SetSourceData(grid)
    except:
        probe.SetInput(input)
        probe.SetSource(grid)
    probe.Update()
    
    outputPointData = probe.GetOutput().GetPointData()
    pointData = dataset.GetPointData()
    pointData.AddArray(outputPointData.GetArray('hsH'))
    pointData.AddArray(outputPointData.GetArray('xm'))


#if __name__ == "__main__":
def UVParametrization(input, hubFileName, tipFileName, axe, stepSize=1.0, relativeExtension=0.1):
    #~ hubFileName = "E:\Documents and Settings\u000160\Bureau\sketch\hub.dat"
    #~ tipFileName = "E:\Documents and Settings\u000160\Bureau\sketch\shroud.dat"
    #~ stepSize = 2.5
    #~ relativeExtension = 0.1
    
    #hubFileName = self.hubFileName
    #tipFileName = self.tipFileName
    #stepSize = self.stepSize
    #relativeExtension = self.relativeExtension
    
    spline, numberOfPoints = CreateSpline(hubFileName, tipFileName, stepSize,
                                    relativeExtension=relativeExtension)
    grid = parametrize(spline, numberOfPoints=numberOfPoints,
                       relativeExtension=relativeExtension)

    if input.IsA('vtkMultiBlockDataSet'):
        output = vtk.vtkMultiBlockDataSet()
        output.CopyStructure(input)
        iter = input.NewIterator()
        iter.UnRegister(None)
        iter.InitTraversal()
        while not iter.IsDoneWithTraversal():
            curInput = iter.GetCurrentDataObject()
            curOutput = curInput.NewInstance()
            curOutput.UnRegister(None)
            output.SetDataSet(iter, curOutput)
            curOutput.ShallowCopy(curInput)
            print 'parametrize'
            parametrize_dataset(curOutput, grid, axe)
            iter.GoToNextItem()
    else:
        output = input.NewInstance()
        output.UnRegister(None)
        output.ShallowCopy(input)
        parametrize_dataset(output, grid, axe)
    return output
