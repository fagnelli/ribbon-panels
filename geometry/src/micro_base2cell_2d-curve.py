# -*- coding: utf-8 -*-

#==============================================================================#
# Author(s)  : Filippo AGNELLI (LMS / X / CNRS)                                #
#              e-mail: filippo.agnelli@polytechnique.edu                       #
#==============================================================================#            
# Description: Converts an NURBS-Python b-spline object into a mesh. Using     #
#              geomdl to describe b-spline objects and pyvista for the mesh.   #
#              meshio converts the mesh into any desired format.               #
#              Work published in F. Agnelli, M. Tricarico, A. Constantinescu.  #
#              Shape-shifting panel from 3D printed undulated ribbon lattice   #
#              Extreme Mechanics Letters, Elsevier BV, 2020, 42, 101089        #
#==============================================================================#
# Version    : v.2023-02-19 .......................................... pass    #
#==============================================================================#
# Risks      : file and directory may be changed over time                     #
#==============================================================================#

# Options
out=True                                       # set to True to export mesh data
graph=True                             # set to True for graphical visualization

# Loading external modules
import copy
import numpy as np
from geomdl import exchange                           # import & export b-spline
from geomdl import multi                                     # geomdl containers
from geomdl import operations

import matplotlib.pyplot as plt

#==============================================================================#
# Input arguments

listnu=['-0.0','-0.2','-0.4','-0.6','-0.8']

# Input
idirname='../b-spline/'

def ifilename(nu):                                    # generate input file name
    nfile='micro_nu='+nu+'_base_2d-curve.json'
    return nfile
	
# Output
odirname='../b-spline/'

def ofilename(nu):                                   # generate output file name
    nfile='micro_nu='+nu+'_cell_3d-curve.json'
    return nfile

#==============================================================================#
def symmetry(obj, p1, p2, **kwargs):

    """ Translates curves, surface or volumes by the input vector.

    Keyword Arguments:
        * ``inplace``: if False, operation applied to a copy of the object. *Default: False*

    :param obj: input geometry
    :type obj: abstract.SplineGeometry or multi.AbstractContainer
    :param vec: translation vector
    :type vec: list, tuple
    :return: translated geometry object
    """
    # Input validity checks
    if not p1 or not p2 or not isinstance(p1, (tuple, list)) or not isinstance(p2, (tuple, list)):
        print("The input must be a list or a tuple")
        exit()

    # Input validity checks
    if len(p1) != obj.dimension or len(p2) != obj.dimension:
        print("The input vector must have " + str(obj.dimension) + " components")
        exit()
        
    # Keyword arguments
    inplace=kwargs.get('inplace', False)

    if not inplace:
        geom=copy.deepcopy(obj)
    else:
        geom=obj
   
#   symmetry control points
    am=np.linalg.inv(np.array([[p2[1]-p1[1],p1[0]-p2[0]],[p2[0]-p1[0],p2[1]-p1[1]]]))

    for g in geom:
        new_ctrlpts=[]
        for pt in g.ctrlpts:
            b=np.array([(p1[1]-p2[1])*(pt[0]-2*p1[0])+(p1[0]-p2[0])*(2*p1[1]-pt[1]),(p2[0]-p1[0])*pt[0]+(p2[1]-p1[1])*pt[1]])
            temp=np.matmul(am, b)
            new_ctrlpts.append(temp.tolist())
        new_ctrlpts=[cp2+[0] for cp2 in new_ctrlpts]
        g.ctrlpts=new_ctrlpts

    return geom

#==============================================================================#
# Build the unit cell from base wall

dictbase=dict()             # dictionnary containing the b-spline for each shape
dictcell=dict()

for nu in listnu:
    pspline=multi.CurveContainer()
    cross0=multi.CurveContainer()
    cell0=multi.CurveContainer()

#   import elementary pattern	
    bsplcurve=exchange.import_json(idirname+ifilename(nu))[0]
    bsplcurve.ctrlpts=[cp2+[0] for cp2 in bsplcurve.ctrlpts]     # from 2D to 3D
    pspline.append(bsplcurve)

#   construct cross by rotating the base around the center
    l0=bsplcurve
    l1=operations.rotate(l0, 90)
    l2=operations.rotate(l0,180)
    l3=operations.rotate(l0,270)
    cross0.add([l0,l1,l2,l3])
    del l1, l2, l3

#   construct unit cell by symmetrizing the cross
    cross1=symmetry(cross0, (0.5,0.,0.), (0.5,1.,0.))
    cross2=symmetry(cross1, (0.,0.5,0.), (1.,0.5,0.))
    cross3=symmetry(cross0, (0.,0.5,0.), (1.,0.5,0.))
    cell0.add(cross0)
    cell0.add(cross1)
    cell0.add(cross2)
    cell0.add(cross3)

    del cross0, cross1, cross2, cross3
	
    dictbase[nu]=pspline
    dictcell[nu]=cell0

    if out:
        exchange.export_json(dictcell[nu], odirname+ofilename(nu))

#==============================================================================#

# Figure 1.
print("Figure 1\n")

fig=plt.figure(frameon=False)
plt.axis('off')
ax=plt.gca()
ax.set_aspect('equal')

lc=['C0','C1','C2','C3','C4']

for i, nu in enumerate(listnu):
	for curve in dictcell[nu]:
		plt.plot(np.array(curve.evalpts)[:,0], np.array(curve.evalpts)[:,1], lw=3, color=lc[i])
#for curve in dictpspline:    
#    plt.plot(np.array(curve.ctrlpts)[:,0],np.array(curve.ctrlpts)[:,1], 
#             c='grey', lw=1, ls='dashdot', marker='o', mfc='k')
