# import numpy as np
# import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# from scipy import integrate
# import matplotlib.colors as mcolors
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 
import sys
#---------------------------------------------------------------------------#
# path to flow fields
path    = str(os.path.dirname(__file__) + '/../test_little_channel/' )
# path    = str(os.path.dirname(__file__) + '/../test_parallel_channel/' )
# path    = str(os.path.dirname(__file__) + '/../test_yamo_180/' )
 
# index
index_flow   = 100
index_scal   = index_flow
# forcing_flow = 'cfr'
forcing_flow = 'cpg'

# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud' #'nearest'
figs    = 'figs' 
plt.close('all')

#---------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read flow fields
flow = mp.Field(path,var='flow',index=index_flow, forcing=forcing_flow)
flow.read_3d_field()

# read scalar fields
scal = mp.Field(path,var='phi1',index=index_scal)
scal.read_3d_field()


# bulk velocity of the flow
ub, upar = mp.ubulk(flow.u, grid.y)
#---------------------------------------------------------------------------#
# plots - 2d flow fields
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
'''
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u')
plt.xlabel("z")
plt.ylabel("y")
plt.pcolormesh(grid.z,grid.y,flow.u[10,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity v')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.v[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity w')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.w[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
'''
#---------------------------------------------------------------------------#
# u-mean
plt.figure(figsize=size)
plt.xlabel("u_mean-velocity")
plt.ylabel("y")
plt.xlim(0,flow.u.mean(axis=(0,2)).max())
plt.ylim(0,grid.y.max())
plt.grid('True')
plt.plot(flow.u.mean(axis=(0,2)), grid.y, marker='.',label='u_mean')
plt.legend(loc=1)
plt.show()

#---------------------------------------------------------------------------#
# plots - 2d scal fields
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- scalar')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,scal.phi1[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
# phi-mean
plt.figure(figsize=size)
plt.xlabel("phi1_mean")
plt.ylabel("y")
plt.xlim(scal.phi1.mean(axis=(0,2)).min(),scal.phi1.mean(axis=(0,2)).max())
plt.ylim(0,grid.y.max())
plt.grid('True')
plt.plot(scal.phi1.mean(axis=(0,2)), grid.y, marker='.',label='phi_mean')
plt.legend(loc=1)
plt.show()

