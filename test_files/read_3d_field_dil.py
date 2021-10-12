import numpy as np
from   scipy import interpolate
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
index        = 20000
# forcing_flow = 'cpg'
# path  = '/home/jonathan/simulations/data/channel_flow/cpg_scal/'
path = str(os.path.dirname(__file__) + '/../test_parallel/' )
# path  = '/home/jonathan/simulations/data/channel_flow/cpg_test/'

# plot settings 
font = {'family':'monospace', 'weight':'bold', 'size':14}; rc('font', **font) # font
plt.rcParams['figure.dpi'] = 210
size    = (8,6)
shading = 'nearest'#'gouraud'
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
dvdy = dvdy.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#
f = open(path +'dwdz.1','rb')
f.seek(52,0)
dwdz = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
dwdz = dwdz.reshape((grid.nx,grid.ny,grid.nz),order='F')
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
f = open(path +'fld_mod.1','rb')
f.seek(52,0)
fld_mod = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
fld_mod = fld_mod.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#---------------------------------------------------------------------------#
sys.exit()
#---------------------------------------------------------------------------#

# 2d plot - yz
plt.figure(figsize=size)
plt.title(r'2d-plot -- yz-plane -- velocity u-mod')
plt.xlabel(r'$z$')
plt.ylabel(r'$y$')
#
plt.pcolormesh(grid.z,grid.y,fld_mod[0,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()


# =============================================================================
# spline functions
# =============================================================================
x_spl = np.zeros(8)
y_spl = np.zeros(8)

# left half
for i in range(8,11):
    x_spl[i-8] = grid.z[i]
    y_spl[i-8] = flow.w[0,6,i]
    
x_spl[3] = grid.z[11]
x_spl[4] = grid.z[20]
y_spl[3] = 0
y_spl[4] = 0

for i in range(21,24):
    x_spl[i-16] = grid.z[i]
    y_spl[i-16] = flow.w[0,6,i]

# spline interpolation
tck = interpolate.splrep(x_spl, y_spl, s=0, k=3)

# new
x_new = grid.z[11:21]
y_new = interpolate.splev(x_new, tck, der=0)

# plot 
alpha=0.5
plt.figure(figsize=size)
plt.grid(True)
plt.plot(grid.z[8:24],fld_mod[0,6,8:24], marker='.', alpha=alpha)
# plt.plot(grid.z[8:24],flow.w[0,6,8:24],  marker='x', alpha=alpha)
# plt.plot(x_spl, y_spl, marker='o',alpha=alpha)
plt.plot(x_new, y_new, marker='o')
# plt.legend(loc=1)
plt.show()



# plot 
plt.figure(figsize=size)
plt.grid(True)
plt.plot(grid.z[8:24],fld_mod[0,6,8:24],marker='.')
plt.plot(grid.z[11],eps[0,6,11]-1,marker='x')
plt.plot(grid.z[20],eps[0,6,20]-1,marker='x')
plt.legend(loc=1)
plt.show()




# =============================================================================
# spline functions
# =============================================================================
x_spl = np.zeros(8)
y_spl = np.zeros(8)

# left half
for i in range(0,3):
    x_spl[i] = - grid.y[3-i]
    
x_spl[4] = grid.y[9]

for i in range(0,3):
    x_spl[5+i] = grid.y[10+i]
    y_spl[5+i] = flow.v[0,10+i,11]

# spline interpolation
tck = interpolate.splrep(x_spl, y_spl, s=0, k=3)

# new
x_new = grid.y[:10]
y_new = interpolate.splev(x_new, tck, der=0)

# plot 
alpha=0.5
plt.figure(figsize=size)
plt.grid(True)
plt.plot(grid.y[:20],fld_mod[0,:20,11],marker='.',alpha=alpha)
plt.plot(grid.y[:20],flow.v[0,:20,11],marker='x',alpha=alpha)
plt.plot(x_new, y_new, marker='o',alpha=alpha)
# plt.plot(x_spl, y_spl, marker='o')
# plt.legend(loc=1)
plt.show()












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