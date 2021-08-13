# %%
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_little_channel/' ) # name = 'flow.20.1' # file = str(path+name)
index = 10

#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

# orignial field
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()

# # u_mod field 
# f = open(path +'dil.1','rb')
# f.seek(52,0)
# dil = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
# dil = dil.reshape((grid.nx,grid.ny,grid.nz),order='F')
# f.close()

# %%
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
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
# plt.figure(figsize=size)
# plt.title('2d-plot -- xy-plane -- velocity w')
# plt.xlabel("x")
# plt.ylabel("y")
# plt.pcolormesh(grid.x,grid.y,flow.w[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()

# %%
#---------------------------------------------------------------------------#
# u-mean
plt.figure(figsize=size)
plt.xlabel("u-velocity")
plt.ylabel("y")
plt.plot(flow.u.mean(axis=0), grid.y, marker='.',label='u_mean')
plt.legend(loc=1)
plt.show()

# %%
#---------------------------------------------------------------------------#
# vertical grid spacing
dy = grid.y[1:] - grid.y[:-1]
plt.figure(figsize=size)
plt.xlabel("nodes")
plt.ylabel("delta_y -- grid spacing")
plt.plot(np.arange(1,grid.ny), dy, marker='.',label='dy')
plt.legend(loc=1)
plt.show()


sys.exit()