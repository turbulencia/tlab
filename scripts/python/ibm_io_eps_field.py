#!/usr/bin/env python3
import numpy as np
import netCDF4  as nc 
import matplotlib.pyplot as plt
"""
#############################################################################
# Version 0.01 - Written by Jonathan Kostelecky
#---------------------------------------------------------------------------#
Script to create the needed geometry field (eps-/epsp-fields), 
from eps2d.nc field as int1 or float8.

Further explenations
eps2d.nc exists in same folder
   - with the specifictions:
         - format:   netcdf
         - variable: eps
         - size:     grid.nx*grid.nz
         - 0s for fluid and 1s for solid 
   - otherwise: create first eps2d.nc on your own!
   
ToDo: - implement bitwise representation of eps0.1

Created  on & by: 2023/08/01 - J. Kostelecky (j.kostelecky@posteo.de)

Modified on & by: ...        - ...
#############################################################################
"""
#---------------------------------------------------------------------------#
# choose options here
#---------------------------------------------------------------------------#
# choose
eps_format = 'int' # 'int' or 'real' (write in i1, f8)
write_eps  = True  # write eps0.1 
write_epsp = True  # write epsp0.1 for pressure grid staggering

# available roughness files
eps2d_present = True   # if True, eps2d.nc exists  
eps3d_present = False  # if True, eps0.1   exists 

# print roughness statistics from 2d
rstats_2d = False  # from 2d, gamma_f and gamma_s
rstats_3d = True   # gamma_s(y), to 'gamma_s.txt' from eps0.1

#---------------------------------------------------------------------------#
# supply informations for the header
#---------------------------------------------------------------------------#
# grid
nx = 1024
ny = 272
nz = nx

# simulation parameters
ni = 0      # iteration step
ti = 0      # physical time
re = 125000 # reynolds number 
nu = re**-1 # viscosity
fr = 0      # froude number
ro = 1      # rossby number
hmax = 40   # max height to compute gamma_s(y[:hmax])

#---------------------------------------------------------------------------#
# Don't change anything below here
#---------------------------------------------------------------------------#
# functions
def write_eps_from2d(eps2d, header, etype='real', epsp=False):

    # header
    hint   = header[0]
    hfloat = header[1]
    
    # output file
    if epsp: 
        o_fname = 'epsp0.1'
        e2d = stagger_eps2d(eps2d)
    else:
        o_fname = 'eps0.1'
        e2d = eps2d

    o_f = open(o_fname,'wb')
    
    # write header of ofile
    print('--------------------------------------------------')
    print('Write eps field   :', o_fname)
    if etype == 'int': hint[0] = 20
    hint.tofile(o_f)
    if etype == 'real': hfloat.tofile(o_f)
   
    # write eps from xz-eps-plane
    nxy = hint[1]*hint[2]
    
    # buffered io (save RAM)
    for k in range(hint[3]):
        i_dat = np.zeros((hint[1], hint[2]))
        for i in range(hint[1]):  
            height = int(e2d[i,k])
            i_dat[i,:height] = 1. 
        i_dat = i_dat.reshape((nxy), order='F')
        if etype == 'real':
            i_dat[:].astype(type_f8).tofile(o_f)
        elif etype == 'int':
            i_dat[:].astype(type_i1).tofile(o_f)        
    o_f.close()
    print('Write  to file    : DONE')
    return
#---------------------------------------------------------------------------#
def stagger_eps2d(eps2d):
    # horizontal staggering - move right interfaces one step to the left,
    # don't touch vertical direction, since no staggering is applied here
    
    nx = eps2d.shape[0]
    nz = eps2d.shape[1]
    
    e2d = np.zeros((nx,nz))
    e2d[:,:] = eps2d[:,:]
    for i in range(nx):
        for k in range(nz):
            if (eps2d[i,k] == 0.) and (eps2d[i,k-1] != 0.):
                e2d[i,k-1] = 0
    for k in range(nz):
        for i in range(nx):
            if (eps2d[i,k] == 0) and (eps2d[i-1,k] != 0):
                e2d[i-1,k] = 0
    return e2d
#---------------------------------------------------------------------------#    
def roughness_statistics_2d(eps2d):
    
    nx = eps2d.shape[0]
    nz = eps2d.shape[1]
    
    e = np.zeros((nx, nz))
    
    # compute fractions
    for i in range(nx):
        for k in range(nz):
            if eps2d[i,k] > 0: e[i,k] = 1
    for i in range(nx-1):
        for k in range(nz-1):
            if (e[i,k] == 0.) and (e[i,k+1] == 1.): e[i,k+1] *= 0.5
            if (e[i,k] == 1.) and (e[i,k+1] == 0.): e[i,k  ] *= 0.5
    for k in range(nz-1):
        for i in range(nx-1):
            if (e[i,k] == 0.) and (e[i+1,k] > 0.): e[i+1,k] *= 0.5
            if (e[i,k] > 0.) and (e[i+1,k] == 0.): e[i  ,k] *= 0.5
    gamma_s = e.sum()/(nx * nz)*100
    
    print('--------------------------------------------------')
    print('Solid fraction [%] : ', np.round(       gamma_s,5))
    print('Fluid fraction [%] : ', np.round( 100 - gamma_s,5))

    return
#---------------------------------------------------------------------------#    
def roughness_statistics_3d(eps):
    # volume approach
    for k in range(eps.shape[2]):     # in x
        for j in range(eps.shape[1]):
            for i in range(eps.shape[0]-1):
                if (eps[i,j,k] == 0) and (eps[i+1,j,k] != 0):
                    eps[i+1,j,k] *= 0.5
                if (eps[i,j,k] != 0) and (eps[i+1,j,k] == 0):
                    eps[i,j,k] *= 0.5
    for i in range(eps.shape[0]):     # in z
        for j in range(eps.shape[1]):
            for k in range(eps.shape[2]-1):
                if (eps[i,j,k] == 0) and (eps[i,j,k+1] != 0):
                    eps[i,j,k+1] *= 0.5
                if (eps[i,j,k] != 0) and (eps[i,j,k+1] == 0):
                    eps[i,j,k] *= 0.5
        for k in range(eps.shape[2]): # in y
            for j in range(eps.shape[1]-1):
                if (eps[i,j,k] == 0) and (eps[i,j+1,k] != 0):
                    eps[i,j+1,k] *= 0.5
                if (eps[i,j,k] != 0) and (eps[i,j+1,k] == 0):
                    eps[i,j,k] *= 0.5
    # gamma_s
    gamma_s = eps.mean(axis=(0,2))
    
    # save
    name = 'gamma_s.txt'
    np.savetxt(name, gamma_s)
    
    # clear
    del eps
    
    return gamma_s
#---------------------------------------------------------------------------#
def read_binary_field(grid=[1,1,1], hmax=50): # buffered
    fname = 'eps0.1'
    a     = np.zeros((grid[0],hmax,grid[2]))
    fxy   = np.zeros((grid[0], grid[1]))

    # load
    print('--------------------------------------------------')
    print('Reading binary file: ', grid)
    if eps_format == 'int': 
        seek = header_size_int
        dtype = type_i1
    elif eps_format == 'real': 
        seek = header_size
        dtype = type_f8
    sk = seek
    for i in range(grid[2]):
        fxy[:,:]  = read_binary(fname, dtype, np.prod(grid[0]*grid[1]), sk).reshape((grid[0],grid[1]), order='F')
        a[:,:,i] = fxy[:,:hmax]
        sk += grid[0]*grid[1]
    # a[:,:,:] = read_binary(fname, dtype, np.prod(grid), seek).reshape((grid[0],grid[1],grid[2]), order='F') # unbuffered
    print('file: ', fname)
    return a[:,:,:]
def read_binary(fname, t, count, seek):
    f = open(fname,'rb')
    f.seek(seek,0)
    rec = np.fromfile(f, t, count)
    f.close()
    return rec
#---------------------------------------------------------------------------#
# data types
type_i1 = np.dtype('<i1') 
type_i4 = np.dtype('<i4') 
type_f8 = np.dtype('<f8')

# build header
head_params_int   = 5 
head_params_float = 4
header_size       = 52
header_size_int   = 20
# 
head_int   = np.array([header_size, nx, ny, nz, ni], type_i4)
head_float = np.array([ti, nu, fr, ro], type_f8)

# read eps2d from nc file
if eps2d_present:
    key   = 'eps2d'
    enc   = nc.Dataset(key+'.nc', 'r', format='netCDF-4')
    eps2d = enc[key][:,:]
if eps3d_present:
    write_eps = False
    eps3d = read_binary_field(grid=[nx,ny,nz], hmax=hmax)
    eps2d = eps3d.sum(axis=1)

# write eps3d based on eps2d (buffered IO)
if write_eps:
    write_eps_from2d(eps2d, [head_int, head_float], etype=eps_format, epsp=False)
if write_epsp:
    write_eps_from2d(eps2d, [head_int, head_float], etype=eps_format, epsp=True)  

# statistics from 2d
if rstats_2d:
    roughness_statistics_2d(eps2d)
if rstats_3d:
    eps3d = read_binary_field(grid=[nx,ny,nz], hmax=hmax)
    gamma_s = roughness_statistics_3d(eps3d)
    
    # plot
    plt.figure()
    plt.title('gamma_s (solid fraction -- volume approach)')
    plt.xlabel('grid nodes')
    plt.xlim(0,hmax+10)
    plt.ylabel('gamma_s')
    plt.plot(gamma_s)
    plt.grid(True)
    plt.show()