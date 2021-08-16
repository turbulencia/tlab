# %%
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# import matplotlib. as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_parallel/' ) # name = 'flow.20.1' # file = str(path+name)
index = 0
#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u')
plt.xlabel("z")
plt.ylabel("y")
plt.pcolormesh(grid.z,grid.y,flow.u[grid.nx//2,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
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
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("z")
plt.pcolormesh(grid.x,grid.z,flow.u[:,1,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# u-mean
plt.figure(figsize=size)
plt.xlabel("u-velocity")
plt.ylabel("y")
plt.plot(flow.u.mean(axis=(0,2)), grid.y, marker='.',label='u_mean')
plt.legend(loc=1)
plt.show()