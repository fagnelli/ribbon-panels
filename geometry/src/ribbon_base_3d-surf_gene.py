# -*- coding: utf-8 -*-

#==============================================================================#
# Author(s)  : Filippo AGNELLI (LMS / X / CNRS)                                #
#              e-mail: filippo.agnelli@polytechnique.edu                       #
#==============================================================================#            
# Description: Compute the control points of a b-spline that best fit the      #
#              shape with desired effective Poisson's ratio. Work published in #
#              F. Agnelli, M. Tricarico, A. Constantinescu.                    #
#              Shape-shifting panel from 3D printed undulated ribbon lattice   #
#              Extreme Mechanics Letters, Elsevier BV, 2020, 42, 101089        #
#==============================================================================#
# Version    : v.2023-02-20 .................................... not tested    #
#==============================================================================#
# Risks      : - file and directory may be changed over time                   #
#              - color of the figures                                          #
#==============================================================================#

# Options
out=True                                       # set to True to export mesh data
graph=True                             # set to True for graphical visualization

# Loading external modules
from geomdl import BSpline
from geomdl import exchange                           # import & export b-spline
from geomdl import utilities
from geomdl.visualization import VisMPL

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
ofilename='ribbon_base_3d-surf.json'

#==============================================================================#

#lh=[0, 0.25, 0.5, 0.75, 1];             # 3rd coordinate (simple interpolation)
lh=[0, 0.21274969, 0.509597215, 0.733649002, 1]; # 3rd coordinate (obtained experimentally)

cpsurf=[]

for i, nu in enumerate(listnu):
    bsplcurve=exchange.import_json(idirname+ifilename(nu))[0]
    bsplcurve.ctrlpts=[cp2+[lh[i]] for cp2 in bsplcurve.ctrlpts] # from 2D to 3D
    cpsurf.extend([cpt for cpt in bsplcurve.ctrlpts])

surf=BSpline.Surface()
surf.degree_u=bsplcurve.degree
surf.degree_v=bsplcurve.degree
surf.set_ctrlpts(cpsurf, len(bsplcurve.ctrlpts), len(listnu))
surf.knotvector_u=utilities.generate_knot_vector(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v=utilities.generate_knot_vector(surf.degree_v, surf.ctrlpts_size_v)
surf.delta=(bsplcurve.delta,bsplcurve.delta)

if out:
    exchange.export_json(surf, odirname+ofilename)
    exchange.export_obj(surf,"../figures/ribbon_base_3d-surf.obj")

#==============================================================================#
