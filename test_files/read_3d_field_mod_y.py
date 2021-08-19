# %%
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_parallel/' ) # name = 'flow.20.1' # file = str(path+name)
index = 0

#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

# orignial field
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()

# u_mod field 
f = open(path +'fld_mod.1','rb')
f.seek(52,0)
u_mod = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
u_mod = u_mod.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# read eps field
f = open(path +'eps0.1','rb')
f.seek(52,0)
eps = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
eps = eps.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# sys.exit()


#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 120 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u_mod')
plt.xlabel("z")
plt.ylabel("y")
#
plt.pcolormesh(grid.z,grid.y,u_mod[0,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
# plt.figure(figsize=size)
# plt.title('2d-plot -- xy-plane -- velocity u_mod')
# plt.xlabel("x")
# plt.ylabel("y")
# #
# plt.pcolormesh(grid.x,grid.y,u_mod[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
# plt.colorbar()
# plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u_mod')
plt.xlabel("x")
plt.ylabel("z")
#
plt.pcolormesh(grid.x,grid.z,u_mod[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#

# sys.exit()

# u_mod = u_mod * eps
# %%
plt.figure(figsize=size)
plt.xlabel("y")
plt.ylabel("v-velocity")
plt.vlines(grid.y[9], ymin=-0.002, ymax=0.002)
for i in range(107,117):
    plt.plot(grid.y[:20],u_mod[0,:20,i], marker='.',label='z-node='+str(i))
plt.legend(loc=1)
plt.grid(True)
plt.show()

# small check 
v     = flow.v * (1. - eps)
v_mod = u_mod  * (1. - eps)
res   = v - v_mod
print(str(res.sum()))


sys.exit()

# %%

plt.figure(figsize=size)
plt.xlabel("y")
plt.ylabel("v-velocity")
for i in range(7,11):
    plt.plot(grid.z[:],u_mod[0,i,:], marker='.',label='y-node='+str(i))
plt.legend(loc=1)
plt.grid(True)
plt.show()

# %%
sys.exit()


#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u')
plt.xlabel("z")
plt.ylabel("y")
#
plt.pcolormesh(grid.z,grid.y,flow.w[0,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
#
plt.pcolormesh(grid.x,grid.y,flow.w[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("z")
#
plt.pcolormesh(grid.x,grid.z,flow.w[:,1,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()





# %%
