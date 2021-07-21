import numpy as np
import matplotlib.pyplot as plt
import os
import my_pylib as mp
import sys
#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_parallel/' ) 
# path  = str(os.path.dirname(__file__) + '/../test_little/' ) 
index   = 1
nob_max = 4
#---------------------------------------------------------------------------#
# read grid field
grid = mp.DnsGrid(path+'grid')

# read eps field
f = open(path +'eps0.1','rb')
f.seek(52,0)
eps = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
eps = eps.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# read nobi field
f = open(path +'nobi0.1','rb')
f.seek(52,0)
nobi = np.fromfile(f, np.dtype('<f8'), grid.ny*grid.nz)
nobi = nobi.reshape((grid.ny,grid.nz),order='F')
f.close()

# read nobj field
f = open(path +'nobj0.1','rb')
f.seek(52,0)
nobj = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.nz)
nobj = nobj.reshape((grid.nx,grid.nz),order='F')
f.close()

# read nobk field
f = open(path +'nobk0.1','rb')
f.seek(52,0)
nobk = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny)
nobk = nobk.reshape((grid.nx,grid.ny),order='F')
f.close()
'''
# 3d nobi field
f = open(path +'nobi3d.1','rb')
f.seek(52,0)
nobi3d = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobi3d = nobi3d.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# 3d nobj field
f = open(path +'nobj3d.1','rb')
f.seek(52,0)
nobj3d = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobj3d = nobj3d.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# 3d nobk field
f = open(path +'nobk3d.1','rb')
f.seek(52,0)
nobk3d = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobk3d = nobk3d.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
'''
sys.exit()
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 180 
size    = (8,6)
shading = 'nearest'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nob fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('nobi')
plt.pcolormesh(grid.z,grid.y,nobi[:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj')
plt.pcolormesh(grid.x,grid.z,nobj[:,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk')
plt.pcolormesh(grid.x,grid.y,nobk[:,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()

plt.show()
#---------------------------------------------------------------------------#
# sys.exit()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nob_3d fields (for debugging)',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('nobi_3d')
plt.pcolormesh(grid.z,grid.y,nobi3d[0,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_3d')
plt.pcolormesh(grid.x,grid.z,nobj3d[:,0,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_3d')
plt.pcolormesh(grid.x,grid.y,nobk3d[:,:,0].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()