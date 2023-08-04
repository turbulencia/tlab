#!/usr/bin/env python3
import numpy as np
import netCDF4  as nc 
import matplotlib.pyplot as plt
import struct
"""
#############################################################################
# Version 0.01 - Written by Jonathan Kostelecky
#---------------------------------------------------------------------------#
IO-script for geometry fields (eps/epsp) from eps2d.nc/eps0.1 (2d/3d)
Convert/Create/Read/Write...

Further explenations

- eps2d.nc exists in same folder
   - with the specifictions:
         - format:   netcdf
         - variable: eps
         - size:     grid.nx*grid.nz
         - values:   0s if only liquid in the vertical, otherwise an integer,
                     the value corresponds to the extent of the object 
                     in grid points from the ground
   - otherwise: create first eps2d.nc/eps0.1 on your own! 
- epsp0.1 can be created from eps2d.nc/eps0.1
- 3d-fields are stored in bitwise (converted to int1)/int1/float8 format
- IO is buffered to handle very large data with low RAM usage
   
ToDo:     - ...

Cautious: - script only works properly if objects are 
            located on the lower domain boundary. 
            Do not use for closed channel flow with objects
            on upper and lower domain boundary

Created  on & by: 2023/08/03 - J. Kostelecky (j.kostelecky@posteo.de)

Modified on & by: ...        - ...
#############################################################################
"""
#---------------------------------------------------------------------------#
# Choose here
#---------------------------------------------------------------------------#
# IO - data format of eps0.1/epsp0.1
# ('bit','int', 'real' {IO in bitwise (i1), i1, f8})
i_format = 'int' # [only considered of eps3d_present = True]
o_format = 'int' 

# input options (available roughness files)
eps2d_present = True   # if True, eps2d.nc exists  
eps3d_present = False  # if True, eps0.1   exists

# output options 
write_eps  = True  # write eps0.1 
write_epsp = True  # write epsp0.1 for pressure grid staggering

# statistics options (solid/fluid fractions)
rstats_2d = True   # from 2d, gamma_f and gamma_s
rstats_3d = True   # gamma_s(y), to 'gamma_s.txt' from eps0.1

#---------------------------------------------------------------------------#
# Supply parameters (needed for the header)
#---------------------------------------------------------------------------#
# grid
nx = 1024
ny = 272
nz = 1024

# simulation parameters
ni = 0      # iteration step
ti = 0      # physical time
re = 125000 # reynolds number 
nu = re**-1 # viscosity
fr = 0      # froude number
ro = 1      # rossby number

# max height to compute gamma_s(y[:hmax])
hmax = 45   

#---------------------------------------------------------------------------#
# Don't change anything below here!
#---------------------------------------------------------------------------#
# functions
def write_eps_from2d(eps2d, header, etype='real', epsp=False):

    # header
    hint   = header[0]
    hfloat = header[1]
    
    nx  = hint[1]
    ny  = hint[2]
    nz  = hint[3]
    nxy = nx*ny
    
    # data size
    bsize = nxy//8
    # rsize = bsize*8
    
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
    if etype == 'bit': 
        if np.mod(nx,8) != 0:
            raise ValueError('ERROR. Restriction mod(imax/8) = 0!')
        else:
            hint[0]  = 20
            hint[1] /= 8
    elif etype == 'int': 
        hint[0] = 20
    hint.tofile(o_f)
    if etype == 'real': hfloat.tofile(o_f)
   
    # buffered io (save RAM)
    for k in range(nz):
        print(f'\r{k+1}/{nz} write nxy planes to disk', end='')
        i_dat = np.zeros((nx,ny))
        for i in range(nx):
            height = int(e2d[i,k])
            i_dat[i,:height] = 1. 
        i_dat = i_dat.reshape((nxy), order='F')
        if etype == 'bit':
            # convert to bitwise-int1
            i_dat_bit = np.zeros(bsize)
            i_dat_bit = bit2int(i_dat_bit,i_dat)
        if etype == 'real':
            i_dat[:].astype(type_f8).tofile(o_f)
        elif etype == 'int':
            i_dat[:].astype(type_i1).tofile(o_f)    
        elif etype == 'bit':
            i_dat_bit[:].astype(type_i1).tofile(o_f)    
    o_f.close()
    print()
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
def read_binary_field(grid=[1,1,1], hmax=50, etype='real'): # buffered
    a     = np.zeros((grid[0],hmax,grid[2]))
    fxy   = np.zeros((grid[0], grid[1]))

    # load
    print('--------------------------------------------------')
    print('Reading binary file: ', grid)
    if etype == 'int': 
        seek = header_size_int
        dtype = type_i1
    elif etype == 'real': 
        seek = header_size
        dtype = type_f8
    sk = seek
    for i in range(grid[2]):
        print(f'\r{i+1}/{nz} read nxy planes from disk', end='')
        fxy[:,:]  = read_binary(fname, dtype, np.prod(grid[0]*grid[1]), sk).reshape((grid[0],grid[1]), order='F')
        a[:,:,i] = fxy[:,:hmax]
        sk += grid[0]*grid[1]
    print()
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
# bitwise functions (cf. ibm_eps_bit.py script)
def read_binary_field_bit(grid=[8,1,1], hmax=-1):
    print('--------------------------------------------------')
    print('Reading binary file: ', grid)
    # data size
    bsize = grid[0]*grid[1]//8
    rsize = bsize*8
    # header
    sk    = header_size_int
    # buffered
    fxy_int = np.zeros(bsize)
    a       = np.zeros((grid[0],hmax,grid[2]))
    # read
    for i in range(grid[2]):
        print(f'\r{i+1}/{nz} read nxy planes from disk', end='')
        f = open(fname,'rb')
        f.seek(sk,0)
        fxy_int = np.fromfile(f, type_i1, bsize)
        f.close()
        # convert
        fxy_bit = np.zeros(rsize)
        fxy_bit = int2bit(fxy_bit, fxy_int) 
        fxy_bit = fxy_bit.reshape((grid[0],grid[1]), order='F')
        # crop
        a[:,:,i] = fxy_bit[:,:hmax]
        sk += bsize
    print()
    return a[:,:,:]
#---------------------------------------------------------------------------#
def int2bit(out,data): # option 2 (bit faster then option 1)
    bsize = data.size
    for i in range(bsize):
        ip = i * 8
        by   = struct.pack('b',data[i])
        by2b = ''.join(format(ord(by), '08b') for byte in by)
        j = 0
        for k in range(-1,-9,-1):
            out[j+ip] = int(str(by2b)[k])
            j += 1
    return out
#---------------------------------------------------------------------------#
def bit2int(out,data):
    bsize = out.size
    for i in range(bsize):
        by = []
        ip = i * 8
        j  = 0
        for k in range(7,-1,-1):
            by.append(int(data[k+ip]))
            j += 1
        bit = "".join(str(i) for i in by)
        if bit[0] == '1':
            out[i] = int(bit,2)-256
        else:    
            out[i] = int(bit,2)
    return out
#---------------------------------------------------------------------------#
# specifications
#---------------------------------------------------------------------------#
fname = 'eps0.1'

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

#---------------------------------------------------------------------------#
# input of geometry field
#---------------------------------------------------------------------------#
# read eps2d from nc file
if eps2d_present:
    key   = 'eps2d'
    enc   = nc.Dataset(key+'.nc', 'r', format='netCDF-4')
    eps2d = enc[key][:,:]

# read eps3d from eps0.1
if eps3d_present:
    if i_format == 'bit':
        # restriction
        if np.mod(nx,8) != 0:
            raise ValueError('ERROR. Restriction mod(imax/8) = 0!')
        else:
            eps3d = read_binary_field_bit(grid=[nx,ny,nz], hmax=hmax)
    else:
        eps3d = read_binary_field(grid=[nx,ny,nz], hmax=hmax, etype=i_format)
    # generate eps2d-field
    eps2d = eps3d.sum(axis=1)

#---------------------------------------------------------------------------#
# output of geometry field
#---------------------------------------------------------------------------#  
# write eps3d based on eps2d (buffered IO)
if write_eps:
    write_eps_from2d(eps2d, [head_int, head_float], etype=o_format, epsp=False)
if write_epsp:
    write_eps_from2d(eps2d, [head_int, head_float], etype=o_format, epsp=True)  

#---------------------------------------------------------------------------#
# generate statistics of geometry field
#---------------------------------------------------------------------------#  
if rstats_2d:
    roughness_statistics_2d(eps2d)
if rstats_3d:
    if write_eps: # read written eps0.1 field
        eps3d = read_binary_field(grid=[nx,ny,nz], hmax=hmax, etype=i_format)
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

#---------------------------------------------------------------------------#    
# Debug - Plot sclices of eps3d
#---------------------------------------------------------------------------#    
# for i in range(0,40,5):
#     plt.figure()
#     plt.imshow(eps3d[:,i,:])
#     plt.show()