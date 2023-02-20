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
from geomdl import exchange                           # import & export b-spline
from geomdl import multi                                     # geomdl containers
import pyvista as pv

from bspline2mesh import bspline2mesh

#==============================================================================#
# Input arguments

listnu=['-0.0','-0.2','-0.4','-0.6','-0.8']

# Input
idirname='../b-spline/'

def ifilename(nu):                                    # generate input file name
    nfile='micro_nu='+nu+'_cell_3d-curve.json'
    return nfile
	
# Output
odirname='../mesh/2d-beam/'

def ofilename(nu):                                   # generate output file name
    nfile='micro_nu='+nu+'_cell_2d-beam'
    return nfile

nbno=15                                                   # number of mesh nodes

#==============================================================================#
# Main code

dictcell=dict()
dictcellmesh=dict()        # dictionnary containing the mesh for each cell shape

for nu in listnu:

    # import elementary pattern
    cell0=multi.CurveContainer()	
    cell0.add(exchange.import_json(idirname+ifilename(nu)))
    dictcell[nu]=cell0

    # convert b-spline to mesh
    dictcellmesh[nu]=bspline2mesh(dictcell[nu], nbno) 
	
	# export mesh to any Meshio format
    if out:
        pv.save_meshio(odirname+'avs-ucd/'+ofilename(nu)+'.avs',
                       dictcellmesh[nu], file_format="avsucd")
        pv.save_meshio(odirname+'abaqus/'+ofilename(nu)+'.inp',
                       dictcellmesh[nu], file_format="abaqus")

#==============================================================================#
