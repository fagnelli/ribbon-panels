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
# Version    : v.2023-02-09 .......................................... pass    #
#==============================================================================#
# Risks      : - file and directory may be changed over time                   #
#              - color of the figures                                          #
#==============================================================================#

# Loading external modules

import numpy as np
import matplotlib.pyplot as plt

import cv2                                                    # image processing
from skimage.morphology import skeletonize                # package for skeleton
from geomdl import BSpline
from geomdl import fitting
from geomdl import linalg
from geomdl import utilities

#==============================================================================
# Options

out=True                                       # set to True to export mesh data
graph=True                             # set to True for graphical visualization

#==============================================================================#

plt.close('all')                                        # close existing windows

if graph:
	plt.ion() 
else:
	plt.ioff()    

#==============================================================================#

# Personal Functions
def robustSkeleton(img):
    """
    Function that computes a morphological skeleton, assuming that the picture
    is periodic, which enhance the accuracy at the boundary.
    """

    # PART I - Get a skeleton
    h, w=img.shape
	# enlarge the image by replicating the borders 
    img=cv2.copyMakeBorder(img, 20, 20, 20, 20, cv2.BORDER_REPLICATE)
    skel=skeletonize(~img/np.max(img))                                # skeleton
    skel=255*skel[20:(h+20), 20:(w+20)]                          # crop skeleton

    # Part II - Get the skeleton of a quarter by superimposing all walls
    sk1=skel[0:int(h/2), int(w/2):w]
    sk2=cv2.rotate(sk1, cv2.ROTATE_90_CLOCKWISE)
    sk3=cv2.rotate(sk1, cv2.ROTATE_180)
    sk4=cv2.rotate(sk1, cv2.ROTATE_90_COUNTERCLOCKWISE)

    skel2 = np.maximum(np.maximum(sk1, sk2),np.maximum(sk3, sk4))
   
    return skel, skel2
 
def visImg(dictshape,listshape):
    """
    Function for visualisation.
    """
    for nu in listshape:
        plt.figure()
        ax=plt.gca()
        ax.set_title(r'$\nu^*$ = '+nu)
        ax.set_aspect('equal')
        plt.imshow(dictshape[nu], cmap=plt.cm.gray)

def bSplineFit(cloud, degree, num_cpts,pi,po):
    
    num_dpts=len(cloud)           # corresponds to variable "r" in the algorithm

    # Get keyword arguments
    degree=degree
    use_centripetal=False
    #num_cpts=kwargs.get('ctrlpts_size', num_dpts - 1)
    num_cpts=num_cpts
    dim=len(cloud[0])                                                # dimension
    
    uk=fitting.compute_params_curve(list(cloud), use_centripetal)       # Get uk
    kv=fitting.compute_knot_vector2(degree, num_dpts, num_cpts, uk)             # Compute knot vector
    
    matrix_n=[]                                               # Compute matrix N
    for i in range(1, num_dpts-1):
        m_temp = []
        for j in range(1, num_cpts-1):
            m_temp.append(fitting.helpers.basis_function_one(degree, kv, j, uk[i]))
        matrix_n.append(m_temp)
    
    matrix_nt=linalg.matrix_transpose(matrix_n)                     # compute NT
    matrix_ntn=linalg.matrix_multiply(matrix_nt, matrix_n)  # compute NTN matrix
    matrix_l, matrix_u=linalg.lu_decomposition(matrix_ntn)    # LU-factorization
    
    # Initialize control points array
    ctrlpts = [[0.0 for _ in range(dim)] for _ in range(num_cpts)]
    
    # Fix start and end points
    ctrlpts[0]  = pi                              #ctrlpts[0] = list(points[0])
    ctrlpts[-1] = po                            #ctrlpts[-1] = list(points[-1])
    
    # Compute Rk - Eqn 9.63
    pt0 = cloud[0]                                                                 # Qzero
    ptm = cloud[-1]                                                                # Qm
    rk = []
    
    for i in range(1, num_dpts - 1):
        ptk = cloud[i]
        n0p = fitting.helpers.basis_function_one(degree, kv, 0, uk[i])
        nnp = fitting.helpers.basis_function_one(degree, kv, num_cpts - 1, uk[i])
        elem2 = [c * n0p for c in pt0]
        elem3 = [c * nnp for c in ptm]
        rk.append([a - b - c for a, b, c in zip(ptk, elem2, elem3)])
    
    vector_r = [[0.0 for _ in range(dim)] for _ in range(num_cpts - 2)]            # Compute R - Eqn. 9.67
    for i in range(1, num_cpts - 1):
        ru_tmp = []
        for idx, pt in enumerate(rk):
            ru_tmp.append([p * fitting.helpers.basis_function_one(degree, kv, i, uk[idx + 1]) for p in pt])
        for d in range(dim):
            for idx in range(len(ru_tmp)):
                vector_r[i - 1][d] += ru_tmp[idx][d]
    
    # Compute control points
    for i in range(dim):
        b = [pt[i] for pt in vector_r]
        y = linalg.forward_substitution(matrix_l, b)
        x = linalg.backward_substitution(matrix_u, y)
        for j in range(1, num_cpts - 1):
            ctrlpts[j][i] = x[j - 1]

    return ctrlpts
    
#==============================================================================#
# Import images

#ndir="D:/Documents/Polytechnique/archives/Clausen/"
ndir="/media/fagnelli/Data/Documents/Polytechnique/archives/Clausen/"
img=cv2.imread(ndir+"Fig3_Original.jpg",0)      # read image & convert grayscale

img00=img[43:(43+132), 601:(601+132)]                     # crop to 132 x 132 px
ret,img00=cv2.threshold(img00,200,255,cv2.THRESH_BINARY)              # binarize 
img02=img[44:(44+130), 454:(454+130)]                     # crop to 130 x 130 px
ret,img02=cv2.threshold(img02,200,255,cv2.THRESH_BINARY)              # binarize

img04=img[44:(44+130), 307:(307+130)]                     # crop to 130 x 130 px
ret,img04=cv2.threshold(img04,200,255,cv2.THRESH_BINARY)              # binarize

img06=img[44:(44+130), 160:(160+130)]                     # crop to 130 x 130 px
ret,img06=cv2.threshold(img06,200,255,cv2.THRESH_BINARY)              # binarize

img08=img[44:(44+130), 12:(12+130)]                       # crop to 130 x 130 px
ret,img08=cv2.threshold(img08,200,255,cv2.THRESH_BINARY)              # binarize

dictimg={'-0.0':img00,'-0.2':img02,'-0.4':img04,'-0.6':img06,'-0.8':img08}
listnu=list(dictimg.keys()) 
del ret, img00, img02, img04, img06, img08

#==============================================================================#

# Skeleton
dictskel = {}; dictsket = {}

for nu in dictimg:
    skel, sket=robustSkeleton(dictimg[nu])
    dictskel[nu]=skel
    dictsket[nu]=sket

del skel, sket

#==============================================================================

# Visualisation
listvis=[listnu[2]]

visImg(dictimg,listvis)
visImg(dictskel,listvis)
visImg(dictsket,listvis)

#==============================================================================

# Cloud of points
dictcloud={}

for nu in dictimg:
    indices = np.where(dictsket[nu]!= [0])
    cloud = np.stack((indices[0], indices[1]), axis=-1)/dictskel[nu].shape[0]

    if (nu =='-0.8') or (nu =='-0.6'):
        cloud = cloud[np.logical_and(cloud[:,0] < 0.25, cloud[:,1] < 0.25)]
    else:
        cloud = cloud[np.logical_and(cloud[:,0] < cloud[:,1], cloud[:,1] < 0.5 - cloud[:,0])]
    
    dictcloud[nu] = cloud

#=============================================================================#
#       B-spline fit                                                          #
#=============================================================================#
    
# Compute the b-spline curve
degree=3
num_cpts=5

ctrlpts=dict()

for nu in dictimg:

    cp=bSplineFit(dictcloud[nu], degree, num_cpts, [0.,0.1],[0.25,0.25])
#    ctrlpts[nu] = cp
#    print('nu=', nu, '; ctrlpts=', cp)

#==============================================================================

# Manual construction of the control point - Trial and error
ctrlpts00=[[0., 0.1], [0.04, 0.1], [0.03, 0.32],   [0.18, 0.26],  [0.25, 0.25]] # nu=-0.0
ctrlpts02=[[0., 0.1], [0.06, 0.1], [0.076, 0.262], [0.15, 0.25],  [0.25, 0.25]] # nu=-0.2
ctrlpts04=[[0., 0.1], [0.06, 0.1], [0.117, 0.16],  [0.14, 0.25],  [0.25, 0.25]] # nu=-0.4
ctrlpts06=[[0., 0.1], [0.06, 0.1], [0.117, 0.16],  [0.213, 0.22], [0.25, 0.25]] # nu=-0.6
ctrlpts08=[[0., 0.1], [0.08, 0.1], [0.16, 0.12],   [0.24, 0.16],  [0.25, 0.25]] # nu=-0.8

ctrlpts={'-0.0':ctrlpts00,
		 '-0.2':ctrlpts02,
		 '-0.4':ctrlpts04,
		 '-0.6':ctrlpts06,
		 '-0.8':ctrlpts08}
del ctrlpts00, ctrlpts02, ctrlpts04, ctrlpts06, ctrlpts08

if out:
    for nu in ctrlpts:
        ctrlpts[nu].reverse()
        np.savetxt('./b-spline/nu='+str(nu)+'_ctrlpts_2d-curve.csv', np.array(ctrlpts[nu]), fmt='%.3f', delimiter=',')

#==============================================================================#
        
# Generate B-spline curve
dictbspline = dict()

for nu in listvis:

    dictbspline[nu]=BSpline.Curve()
    dictbspline[nu].degree=degree
    dictbspline[nu].ctrlpts=ctrlpts[nu]
    dictbspline[nu].knotvector=utilities.generate_knot_vector(degree, dictbspline[nu].ctrlpts_size)
    dictbspline[nu].delta=0.02

    plt.figure()
    plt.xlim(0, 0.25)
    plt.ylim(0.05, 0.3)
    plt.gca().set_aspect('equal')
    plt.scatter(dictcloud[nu][:,0],dictcloud[nu][:,1], c='C1')
    plt.plot(np.array(dictbspline[nu].evalpts)[:,0],np.array(dictbspline[nu].evalpts)[:,1],'k-',label='curve')
    plt.plot(np.array(dictbspline[nu].ctrlpts)[:,0],np.array(dictbspline[nu].ctrlpts)[:,1],'bo', ls='dashdot',label='control points')
    plt.legend(loc=3, fontsize='small', fancybox=True)
    
#==============================================================================
