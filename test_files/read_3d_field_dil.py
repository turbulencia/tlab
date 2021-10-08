import numpy as np
import sys
import matplotlib.pyplot as plt
from   matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import os
import my_pylib as mp
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 
# import os
#-----------------------------------------------------------------------------#
# file
index        = 0
# forcing_flow = 'cpg'
# path  = '/home/jonathan/simulations/data/channel_flow/cpg_scal/'
path = str(os.path.dirname(__file__) + '/../test_little/' )
# path  = '/home/jonathan/simulations/data/channel_flow/cpg_test/'

# plot settings 
font = {'family':'monospace', 'weight':'bold', 'size':14}; rc('font', **font) # font
plt.rcParams['figure.dpi'] = 210
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read flow fields
flow = mp.Field(path,var='flow',index=index)#, forcing=forcing_flow)
flow.read_3d_field()

# # read scalar fields
# scal1 = mp.Field(path,var='phi1',index=index)#, forcing=forcing_flow)
# scal1.read_3d_field()

# read dudx, dvdy, dwdz, dil, eps
f = open(path +'dudx.1','rb')
f.seek(52,0)
dudx = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
dudx = dudx.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'dvdy.1','rb')
f.seek(52,0)
dvdy = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
dvdy = dudx.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'dwdz.1','rb')
f.seek(52,0)
dwdz = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
dwdz = dudx.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'eps0.1','rb')
f.seek(52,0)
eps = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
eps = eps.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'dil.1','rb')
f.seek(52,0)
dil = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
dil = dil.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'result1.1','rb')
f.seek(52,0)
result1 = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
result1 = result1.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'result2.1','rb')
f.seek(52,0)
result2 = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
result2 = result2.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#---------------------------------------------------------------------------#

dil_py = - (dudx + dvdy + dwdz)

err = dil - dil_py

sys.exit()












#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u_mod')
plt.xlabel("z")
plt.ylabel("y")
#
plt.pcolormesh(grid.z,grid.y,dil[0,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u_mod')
plt.xlabel("x")
plt.ylabel("y")
#
plt.pcolormesh(grid.x,grid.y,dil[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u_mod')
plt.xlabel("x")
plt.ylabel("z")
#
plt.pcolormesh(grid.x,grid.z,dil[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#

sys.exit()



plt.figure(figsize=size)
plt.xlabel("z")
plt.ylabel("w-velocity")
for i in range(0,10):
    plt.plot(grid.z,dil[-1,i,:], marker='.',label='y-node='+str(i))
plt.legend(loc=1)
plt.show()

# small check 
w     = flow.w * (1. - eps)
w_mod = u_mod  * (1. - eps)
res   = w - w_mod
print(str(res.sum()))

# plt.figure(figsize=size)
# plt.xlabel("z")
# plt.ylabel("w-velocity")
# for i in range(0,5):
#     plt.plot(grid.z,u_mod_glob[-1,i,:], marker='.',label='y-node='+str(i))
#     plt.plot(grid.z,u_mod_loc[-1,i,:], marker='x',label='y-node='+str(i))
# plt.legend(loc=1)
# plt.show()


# plt.figure(figsize=size)
# plt.xlabel("z")
# plt.ylabel("w-velocity")
# for i in range(0,5):
#     plt.plot(grid.z,u_mod_glob[-1,i,:] - u_mod_loc[-1,i,:], marker='.',label='y-node='+str(i))
#     # plt.plot(grid.z,u_mod_loc[-1,i,:], marker='x',label='y-node='+str(i))
# plt.legend(loc=1)
# plt.show()

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