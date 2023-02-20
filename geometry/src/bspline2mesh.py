# -*- coding: utf-8 -*-

#==============================================================================#
# Author(s)  : Filippo AGNELLI (LMS / X / CNRS)                                #
#              e-mail: filippo.agnelli@polytechnique.edu                       #
#==============================================================================#            
# Description: Converts an NURBS-Python b-spline object into a mesh. Using     #
#              geomdl to describe b-spline objects and pyvista for the mesh.   #
#              meshio converts the mesh into any desired format                #
#==============================================================================#
# Version    : v.2023-02-19 .......................................... pass    #
#==============================================================================#
# Risks      : file and directory may be changed over time                     #
#==============================================================================#

import numpy as np
import pyvista as pv

def dichotomysolver(bspline, z):
    
    vmax=1; vmin=0;                              # second curvilinear coordinate

    zmax=bspline.evaluate_single([vmax,0])[2]
    zmin=bspline.evaluate_single([vmin,0])[2]
    zsol=0

    # if we are at the boundary the solution is known.
    if   zmax==z: return float(vmax) 
    elif zmin==z: return float(vmin)
        
    # Check that our z in within our interval.
    if (zmax<z) or (zmin>z):
        print('error range')
        return None

    while abs(zsol-z)>1e-10:
        vsol=0.5*(vmax+vmin)
        zsol=bspline.evaluate_single([vsol,0])[2]
        if   zsol>z: vmax=vsol
        elif zsol<z: vmin=vsol
        
    return vsol

#==============================================================================#

def bspline2mesh(bspline, dens):
    m=pv.PolyData()
    
    if str(bspline)=='container':
        for shape in bspline:


    # Case I - B-spline curve
            if str(shape)=='curve':
                p=shape.evaluate_list(np.linspace(0, 1, dens)) # evaluate points
                m+=pv.lines_from_points(p)          # generate lines from points	

    # Case II - B-spline surface
            elif str(shape)=='surface':
                p=[]
                h=shape.bbox[1][2]-shape.bbox[0][2]
                for pz in range(int(h/0.02)+1):
                    u=dichotomysolver(shape, round(0.02*pz,6))
                    coor=np.transpose(np.append([np.linspace(0, 1, dens)], np.full((1, dens),u),axis=0))
                    coor[:, [1, 0]] = coor[:, [0, 1]]
                    p.extend(shape.evaluate_list(coor))
                p=np.around(np.array(p), decimals=4)
                m+=pv.PolyData(p).delaunay_2d(alpha=0.035)   # Delaunay triangulation
            else:
                print('oups')
                exit(1)
                
    elif str(bspline)=='curve':
        p=bspline.evaluate_list(np.linspace(0, 1, dens))       # evaluate points 
        m=pv.lines_from_points(p)                   # generate lines from points	

    return m

#==============================================================================#
