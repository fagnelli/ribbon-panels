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
import numpy as np
from geomdl import exchange                           # import & export b-spline
from geomdl import multi                                     # geomdl containers
import pyvista as pv

from bspline2mesh import bspline2mesh

#==============================================================================#
# Input arguments

# Input
idirname='../b-spline/'
ifilename='ribbon_wall_3d-surf.json'
	
# Output
odirname='../mesh/3d-shell/'
def ofilename(nu,h): 
    nfile='ribbon_wall_nu='+str('{:.2f}'.format(nu[0]))+'-'+str('{:.2f}'.format(nu[1]))+'_h='+str('{:.2f}'.format(h))+'_3d-shell'
    return nfile

nbno=15                                                   # number of mesh nodes

# Height
#=======    

#h=1
#lh=[h]
lh = [0.8, 0.68, 0.56, 0.48, 0.4, 0.34, 0.28, 0.24]
#lh = [round(0.02*i + 0.1,3) for i in range(46)]
#lh = [0.3]



# Domain
#=======
domain = (0.,1.)
ld = [domain]
#ld = [(round(0.05*i,3),round(1-0.05*i,3)) for i in range(10)]

#==============================================================================#
# Main code

for h in lh:
    
    print('Height',h)
    for domain in ld:

# import elementary pattern
        cell0=multi.SurfaceContainer()	
        cell0.add(exchange.import_json(idirname+ifilename))

        for surf in cell0:
            surf.ctrlpts=(np.array(surf.ctrlpts)*np.array([1,1,h])).tolist()

#       convert b-spline to mesh
        cellmesh=bspline2mesh(cell0, nbno) 
	
# export mesh to any Meshio format
        if out:
            pv.save_meshio(odirname+'avs-ucd/'+ofilename(domain,h)+'.avs',
                           cellmesh, file_format="avsucd")
            pv.save_meshio(odirname+'abaqus/'+ofilename(domain,h)+'.inp',
                           cellmesh, file_format="abaqus")
            pv.save_meshio(odirname+'stl/'+ofilename(domain,h)+'.stl',
                           cellmesh, file_format="stl", binary=True)

#==============================================================================#
