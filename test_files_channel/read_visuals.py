import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
from scipy import integrate
# import matplotlib.colors as mcolors
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 

#-----------------------------------------------------------------------------#
# path to flow fields
path    = str(os.path.dirname(__file__) + '/../test_little_channel/' )

# index
index_flow = 1

# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')

#-----------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read pressure field (usual scalar field)
f = open(path +'Pressure000000.1','rb')
f.seek(52,0)
pre = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
pre = pre.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# # read pressure field (single precission, no header)
# f = open(path +'Pressure000000','rb')
# f.seek(0,0)
# pre = np.fromfile(f, np.dtype('<f4'), grid.nx*grid.ny*grid.nz)
# pre = pre.reshape((grid.nx,grid.ny,grid.nz),order='F')
# f.close()

# simulation propteries
re_cl      = 4226                 # centerline Re-number
re_tau_lam = re_cl**0.88 * 0.116  # friction   Re-number (lam.)
nu         = 1 / re_cl            # kinematic viscosity
# ro    = 1                    # Rossby number
# pr    = 1                    # Prandtl number
# rho   = 1                    # density

#---------------------------------------------------------------------------#
# pressure gradient
pre_xy = pre.mean(axis=2)
dpdx   = np.zeros([grid.nx,grid.ny]) 
for i in range(0,grid.ny):
    dpdx[:,i] = mp.derivative1(grid.x, pre_xy[:,i], grid.nx)

#---------------------------------------------------------------------------#
# plots - 2d fields

# --------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('pressure field',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('xy - plane')
plt.pcolormesh(grid.x,grid.y,pre[:,:,0].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('xy - plane')
plt.pcolormesh(grid.x,grid.y,dpdx[:,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#

# mean pressure
plt.figure(figsize=size)
plt.xlabel('$y$')
plt.ylabel("$p$")
plt.plot(grid.y, pre.mean(axis=(0,2)), label='pre')
plt.legend()
plt.grid(True)
plt.show()
# 
plt.figure(figsize=size)
plt.xlabel('$x$')
plt.ylabel("$p$")
plt.plot(grid.x, pre.mean(axis=(1,2)), label='pre')
plt.legend()
plt.grid(True)
plt.show()
# 
plt.figure(figsize=size)
plt.xlabel('$x$')
plt.ylabel("$p$")
plt.plot(grid.x, dpdx[:,0], label='dpdx')
# plt.plot(grid.x, dpdx.mean(axis=(0,1)), label='dpdx_mean')
plt.legend()
plt.grid(True)
plt.show()

