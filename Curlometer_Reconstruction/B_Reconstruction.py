# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 12:20:02 2021

@author: Teddy Broeren
"""
#from math import atan2
#import h5py
import numpy as np
import numpy.matlib
import scipy.linalg as la
from itertools import combinations
from datetime import datetime
import matplotlib.pyplot as plt

startTime = datetime.now()
def find_polys(loc,n):
    # find_poly creates a list of all n length combinations of numbers of 0 through N
    N = np.size(loc,0)
    temp = list(range(N))
    return list(combinations(temp,n))

def get_positions():
    # read .txt files containing locations of spacecraft in HelioSwarm mission design
    length = np.zeros(9)
    for name in range(9):
        data = np.loadtxt("positions/n%s_clean.txt" % name,dtype = str)
        length[name] = np.size(data[:,4])
    
    L = np.min(length).astype(int)
    times = data[:L,0:4]
    positions = np.zeros([9,L,3])
    for name in range(9):
        data = np.loadtxt("positions/n%s_clean.txt" % name,dtype = str)    
        positions[name,:,:] = data[:L,4:7].astype(float)
    return [positions, times]

def calc_RLEP(r): 
    # function calc_RLEP takes in array, each row of which represents a point in 3d space
    # find the number of satellites as the number of rows of r
    N = np.size(r[:,0])
    # find the mesocentre rb
    rb = np.mean(r, axis=0)
    # calculate the Volumetric Tensor R
    R = np.zeros([3,3])
    for i in range(N):
        R += np.outer(r[i,:]-rb, r[i,:]-rb)/N
    # find the eigenvalues of R as value in lambdas
    temp = la.eig(R)
    lambdas = temp[0]
    # find semiaxes of quasi-ellipsoid a,b,c
    # check if eigenvalues are real
    if any(np.imag(lambdas) != np.zeros(3)):
        raise ValueError('Eigenvalue has imaginary component')
    lambdas_real = np.real(lambdas)
    [c,b,a] = np.sqrt( np.sort(lambdas_real) )
    # calculate L,E,P
    L = 2*a
    E = 1 - b/a
    P = 1 - c/b
    return [R,L,E,P]

def CurlometerO1(r,B,x0):
    # First order curlometer method
    # r is 4X3 array of s/c absolute positions
    # B is 4x3 array of s/c measured B values
    # x0 is 3x1 vector (location to calculate J)
    
    # recenter problem to desired location
    r = r - np.vstack((x0,x0,x0,x0))
    
    # set constants vector of measured B values
    b = np.zeros([12,1])
    b[0:4,0], b[4:8,0], b[8:12,0] = B[:,0], B[:,1], B[:,2]
    
    # set linking matrix of positions
    H = np.zeros([4,4])
    H[:,0] = np.ones([4])
    H[:,1:4] = r
    A = np.zeros([12,12])
    A[0:4,0:4],A[4:8,4:8],A[8:12,8:12] = H,H,H
    
    # solve linear system Ax=b
    x = np.linalg.solve(A, b)
    
    # extract divergence, mag field at x0, and J
    mu_0 = 1
    J = (1/mu_0)*np.array([x[10]-x[7],x[3]-x[9],x[5]-x[2]]).T
    divB = x[1] + x[6] + x[11]
    B_x0 = np.array([x[0],x[4],x[8]]).T
        
    return divB, B_x0, J

'''
The two required conditions for using a tetrahedron to reconstruct B at a given point k are:
    1. The tetrahedron must have a shape parameter less than chi_thres
            \chi := \sqrt{E^2 + P^2} < chi_thres
            
    2. The reconstructed point k must be near the tetrahedron's barycenter 
           |k - r_0| < L_coeff*L
    where L is the characteristic size of the selected tetrahedron
'''
# %% Set initial conditions
hour = 205              # hour of HS configuration to select [94, 144, or 205]
L_coeff = 1             # radius of reconstruction around each tetrahedron's barycenter
chi_thres = 1.0         # shape threshold for using a tetrahedron in reconstruction (chi := sqrt{E^2 + P^2})

# %% load data
B_analytical = np.load('example_data/hour%i/B_analytical.npy' %hour)
B_sc = np.load('example_data/hour%i/B_sc.npy' %hour)
r_field = np.load('example_data/hour%i/r_field.npy' %hour)
r0 = np.load('example_data/hour%i/r0.npy' %hour)

# %% find spacecraft positions
nx,ny,nz = 10,10,10     # resolution of 3d reconstruction grid

[positions, times] = get_positions()
r = positions[:,hour,:]
r = r - np.mean(r,axis=0)
r_sc = r + np.matlib.repmat(r0,9,1)  
n = nx*ny*nz

# %% find shape of all tetrahedral subsets of spacecraft
poly_indices = find_polys(r_sc,4)
m = np.shape(poly_indices)[0]
L_save = np.zeros([m])
E_save = np.zeros([m])
P_save = np.zeros([m])
for j in range(m):
        # j is row of indices array to use points of 
        # (a,b,c,d)-> use s/c numbered a,b,c,d to define tetrahedron
        index = list(poly_indices[j])
        [R, L_save[j], E_save[j], P_save[j]] = calc_RLEP(r_sc[index,:])
        
# %% compute the curlometer on all reconstruction points 
B_estimate = np.zeros([n,3])
Dist2Bary = np.zeros([m])
tetra_passed = np.zeros([n])
# loop over every reconstructed point k
for k in range(n):
    if np.mod(k,100) == 0:
        print('%.1f%% of points computed' %(100*k/n) )

    X0 = r_field[k,:]
    B_recon = np.zeros([m,3])
    for j in range(m):
        index = list(poly_indices[j])
        # find the distance from the barycenter of each tetrahedron to each reconstructed point
        Dist2Bary[j] = la.norm(X0 - np.mean(r_sc[index,:],axis=0))
        # find 1st order curlometer solution using this particular subset of s/c
        DivB1, B_recon[j,:], J_curl1 = CurlometerO1(r_sc[index,:],B_sc[index,:],X0)

    # determine which tetrahedra are well-shaped and nearby to the chosen point k
    ind_good = np.logical_and(np.sqrt(E_save**2 + P_save**2) < chi_thres, Dist2Bary < L_coeff*L_save)
    tetra_passed[k] = np.sum(ind_good)
    
    # use only these well-shaped and nearby tetrahedron in the estimate for B at this point k
    if tetra_passed[k] != 0:
        B_estimate[k,:] = np.mean(B_recon[ind_good],axis=0)
        
# %% compute error
Error = 100* np.divide(la.norm(B_estimate - B_analytical,axis=1), la.norm(B_estimate,axis=1), out=np.ones_like(la.norm(B_estimate,axis=1)), where = la.norm(B_estimate,axis=1)!=0)

endTime = datetime.now()
print('Time to execute:',endTime-startTime)

# %% output key
'''
B_analytical    - true value of B at all points on grid around spacecraft
B_estimate      - value of B at all points on grid if only averaging over tetrahedra which pass both criteria
B_sc            - synthetic B measurements taken by the 9 spacecraft

r_field         - absolule position of the locations of B reconstruction (which are centered around spacecraft)
r_sc            - absolute position of the nine spacecraft which take B measurements

Error           - percent error at all points in space if using criteria to select subsets to average over
tetra_passed    - reports the number of tetrahedra that passed both criteria at each location in space
''' 

# %% plot reconstructed field example
r_field_cube = np.reshape(r_field,(nx,ny,nz, 3))
B_field_cube = np.reshape(B_analytical,(nx,ny,nz, 3))
B_recon_cube = np.reshape(B_estimate,(nx,ny,nz, 3))
z_slice = int(nz/2)
SCALE = .8

fig = plt.figure(figsize=(8,6))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

L2 = plt.quiver(r_field_cube[:,:,z_slice,0], r_field_cube[:,:,z_slice,1], B_recon_cube[:,:,z_slice,0], B_recon_cube[:,:,z_slice,1],scale=SCALE, headaxislength=3, headwidth=3, headlength=3, facecolor='r',alpha=0.3)
L3 = plt.quiver(r_field_cube[:,:,z_slice,0], r_field_cube[:,:,z_slice,1], B_field_cube[:,:,z_slice,0], B_field_cube[:,:,z_slice,1],scale=SCALE, headaxislength=3, headwidth=3, headlength=3, facecolor='k',alpha=1)
plt.scatter(r_sc[:,0], r_sc[:,1], s=100, color='b',marker='+', alpha=1)

plt.xlabel(r'$x$ (km)' ,fontsize=20)
plt.ylabel(r'$y$ (km)',fontsize=20)
plt.grid(b=True, which='major', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.05)
plt.quiverkey(L2, 0.35, 1.02, SCALE/20, r'Reconstruction',labelpos ='E' )
plt.quiverkey(L3, 0.65, 1.02, SCALE/20, r'True Field',labelpos ='E' )
plt.title('Curlometer Reconstruction: Hour %i\n' %(hour),fontsize=20, wrap=True)
plt.tight_layout()
plt.savefig('figures/Brecon_hour%i.png' %hour, format='png',dpi = 600)
