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
# Version    : v.2023-02-20 .......................................... pass    #
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


idirname='../b-spline/'                                                  # input
ifilename='ribbon_base_3d-surf.json'


odirname='../b-spline/'                                                 # output
ofilename='ribbon_wall_3d-surf.json'

#==============================================================================#
def symmetry3d(obj, p1, p2, **kwargs):

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

        
    # Keyword arguments
    inplace = kwargs.get('inplace', False)

    if not inplace:
        geom = copy.deepcopy(obj)
    else:
        geom = obj
   
# Symmetry control points
    am = np.linalg.inv(np.array([[p2[1]-p1[1],p1[0]-p2[0]],[p2[0]-p1[0],p2[1]-p1[1]]]))

    for g in geom:
        new_ctrlpts = []
        for pt in g.ctrlpts:
            b = np.array([(p1[1]-p2[1])*(pt[0]-2*p1[0])+(p1[0]-p2[0])*(2*p1[1]-pt[1]),(p2[0]-p1[0])*pt[0]+(p2[1]-p1[1])*pt[1]])
            temp = np.matmul(am, b)
            temp = np.append(temp, [pt[2]], axis=0)
            new_ctrlpts.append(temp.tolist())
        g.ctrlpts = new_ctrlpts

    return geom

#==============================================================================#
# Build the unit cell from base wall

pspline=multi.SurfaceContainer()
line0=multi.SurfaceContainer()
cell0=multi.SurfaceContainer()

# import elementary pattern	
bsplcurve=exchange.import_json(idirname+ifilename)[0]
pspline.append(bsplcurve)

# construct cross by rotating the base around the center
l0=bsplcurve
l1=operations.rotate(l0,180)
line0.add([l0,l1])
del l1

#   construct unit cell by symmetrizing the cross
line1=symmetry3d(line0, (0.5,0., 0), (0.5,1., 0))
cell0.add(line0)
cell0.add(line1)

del line0, line1
	
if out:
    exchange.export_json(cell0, odirname+ofilename)
    exchange.export_obj(cell0,"../figures/ribbon_wall_3d-surf.obj")

#==============================================================================#
