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
#---------------------------------------------------------------------------#
# # read nobi field
# f = open(path +'nobi0.1','rb')
# f.seek(52,0)
# nobi = np.fromfile(f, np.dtype('<f8'), grid.ny*grid.nz)
# nobi = nobi.reshape((grid.ny,grid.nz),order='F')
# f.close()

f = open(path +'nobi_b0.1','rb')
f.seek(52,0)
nobi_b = np.fromfile(f, np.dtype('<f8'), nob_max*grid.ny*grid.nz)
nobi_b = nobi_b.reshape((nob_max,grid.ny,grid.nz),order='F')
f.close()

f = open(path +'nobi_e0.1','rb')
f.seek(52,0)
nobi_e = np.fromfile(f, np.dtype('<f8'), nob_max*grid.ny*grid.nz)
nobi_e = nobi_e.reshape((nob_max,grid.ny,grid.nz),order='F')
f.close()
#---------------------------------------------------------------------------#
# # read nobj field
# f = open(path +'nobj0.1','rb')
# f.seek(52,0)
# nobj = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.nz)
# nobj = nobj.reshape((grid.nx,grid.nz),order='F')
# f.close()

f = open(path +'nobj_b0.1','rb')
f.seek(52,0)
nobj_b = np.fromfile(f, np.dtype('<f8'), nob_max*grid.nx*grid.nz)
nobj_b = nobj_b.reshape((grid.nx,nob_max,grid.nz),order='F')
f.close()

f = open(path +'nobj_e0.1','rb')
f.seek(52,0)
nobj_e = np.fromfile(f, np.dtype('<f8'), nob_max*grid.nx*grid.nz)
nobj_e = nobj_e.reshape((grid.nx,nob_max,grid.nz),order='F')
f.close()
#---------------------------------------------------------------------------#
# # read nobk field
# f = open(path +'nobk0.1','rb')
# f.seek(52,0)
# nobk = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny)
# nobk = nobk.reshape((grid.nx,grid.ny),order='F')
# f.close()

f = open(path +'nobk_b0.1','rb')
f.seek(52,0)
nobk_b = np.fromfile(f, np.dtype('<f8'), nob_max*grid.nx*grid.ny)
nobk_b = nobk_b.reshape((grid.nx,grid.ny,nob_max),order='F')
f.close()

f = open(path +'nobk_e0.1','rb')
f.seek(52,0)
nobk_e = np.fromfile(f, np.dtype('<f8'), nob_max*grid.nx*grid.ny)
nobk_e = nobk_e.reshape((grid.nx,grid.ny,nob_max),order='F')
f.close()

# sys.exit(0)
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
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
plt.suptitle('nobi_b fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi_b layer 1')
plt.pcolormesh(grid.z,grid.y,nobi_b[0,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi_b layer 2')
plt.pcolormesh(grid.z,grid.y,nobi_b[1,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi_b layer 3')
plt.pcolormesh(grid.z,grid.y,nobi_b[2,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobi_e fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi_e layer 1')
plt.pcolormesh(grid.z,grid.y,nobi_e[0,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi_e layer 2')
plt.pcolormesh(grid.z,grid.y,nobi_e[1,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi_e layer 3')
plt.pcolormesh(grid.z,grid.y,nobi_e[2,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobj_b fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_b  layer 1')
plt.pcolormesh(grid.x,grid.z,nobj_b[:,0,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_b layer 2')
plt.pcolormesh(grid.x,grid.z,nobj_b[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_b layer 3')
plt.pcolormesh(grid.x,grid.z,nobj_b[:,3,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobj_e fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_e  layer 1')
plt.pcolormesh(grid.x,grid.z,nobj_e[:,0,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_e layer 2')
plt.pcolormesh(grid.x,grid.z,nobj_e[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj_e layer 3')
plt.pcolormesh(grid.x,grid.z,nobj_e[:,2,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobk_b fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,3,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_b  layer 1')
plt.pcolormesh(grid.x,grid.y,nobk_b[:,:,0].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_b  layer 2')
plt.pcolormesh(grid.x,grid.y,nobk_b[:,:,1].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_b  layer 3')
plt.pcolormesh(grid.x,grid.y,nobk_b[:,:,2].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,4)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_b layer 4')
plt.pcolormesh(grid.x,grid.y,nobk_b[:,:,3].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
# plt.subplot(2,3,5)
# plt.xlabel('$x$')
# plt.ylabel("$y$")
# plt.title('nobk_b layer 5')
# plt.pcolormesh(grid.x,grid.y,nobk_b[:,:,4].T, shading=shading, cmap='RdBu_r')
# plt.colorbar()
# #
# plt.subplot(2,3,6)
# plt.xlabel('$x$')
# plt.ylabel("$y$")
# plt.title('nobk_b  layer 6')
# plt.pcolormesh(grid.x,grid.y,nobk_b[:,:,5].T, shading=shading, cmap='RdBu_r')
# plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobk_e fields',fontsize='xx-large')
#--------------------------#
plt.subplot(2,3,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_e  layer 1')
plt.pcolormesh(grid.x,grid.y,nobk_e[:,:,0].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_e  layer 2')
plt.pcolormesh(grid.x,grid.y,nobk_e[:,:,1].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_e  layer 3')
plt.pcolormesh(grid.x,grid.y,nobk_e[:,:,2].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,4)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk_e  layer 4')
plt.pcolormesh(grid.x,grid.y,nobk_e[:,:,3].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
# plt.subplot(2,3,5)
# plt.xlabel('$x$')
# plt.ylabel("$y$")
# plt.title('nobk_e  layer 5')
# plt.pcolormesh(grid.x,grid.y,nobk_e[:,:,4].T, shading=shading, cmap='RdBu_r')
# plt.colorbar()
#
# plt.subplot(2,3,6)
# plt.xlabel('$x$')
# plt.ylabel("$y$")
# plt.title('nobk_e  layer 6')
# plt.pcolormesh(grid.x,grid.y,nobk_e[:,:,5].T, shading=shading, cmap='RdBu_r')
# plt.colorbar()
#
plt.show()
#############################################################################
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
sys.exit(0)
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#############################################################################
# 3d nobi_b and nobi_e fields
f = open(path +'nobi3d_b.1','rb')
f.seek(52,0)
nobi3d_b = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobi3d_b = nobi3d_b.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

f = open(path +'nobi3d_e.1','rb')
f.seek(52,0)
nobi3d_e = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobi3d_e = nobi3d_e.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# 3d nobj_b and nobj_e fields
f = open(path +'nobj3d_b.1','rb')
f.seek(52,0)
nobj3d_b = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobj3d_b = nobj3d_b.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

f = open(path +'nobj3d_e.1','rb')
f.seek(52,0)
nobj3d_e = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobj3d_e = nobj3d_e.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

# 3d nobk_b and nobk_e fields
f = open(path +'nobk3d_b.1','rb')
f.seek(52,0)
nobk3d_b = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobk3d_b = nobk3d_b.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()

f = open(path +'nobk3d_e.1','rb')
f.seek(52,0)
nobk3d_e = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
nobk3d_e = nobk3d_e.reshape((grid.nx,grid.ny,grid.nz),order='F')
f.close()
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobi_b 3d fields for debugging',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi3d_b layer 1')
plt.pcolormesh(grid.z,grid.y,nobi3d_b[0,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi3d_b layer 2')
plt.pcolormesh(grid.z,grid.y,nobi3d_b[1,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi3d_b layer 3')
plt.pcolormesh(grid.z,grid.y,nobi3d_b[2,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobi_e 3d fields for debugging',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi3d_e layer 1')
plt.pcolormesh(grid.z,grid.y,nobi3d_e[0,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi3d_e layer 2')
plt.pcolormesh(grid.z,grid.y,nobi3d_e[1,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobi3d_e layer 3')
plt.pcolormesh(grid.z,grid.y,nobi3d_e[2,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobj_b 3d fields for debugging',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj3d_b  layer 1')
plt.pcolormesh(grid.x,grid.z,nobj3d_b[:,0,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj3d_b layer 2')
plt.pcolormesh(grid.x,grid.z,nobj3d_b[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj3d_b layer 3')
plt.pcolormesh(grid.x,grid.z,nobj3d_b[:,2,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobj_e 3d fields for debugging',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj3d_e  layer 1')
plt.pcolormesh(grid.x,grid.z,nobj3d_e[:,0,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj3d_e layer 2')
plt.pcolormesh(grid.x,grid.z,nobj3d_e[:,1,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('nobj3d_e layer 3')
plt.pcolormesh(grid.x,grid.z,nobj3d_e[:,2,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobk_b 3d fields for debugging',fontsize='xx-large')
#--------------------------#
plt.subplot(2,3,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_b  layer 1')
plt.pcolormesh(grid.x,grid.y,nobk3d_b[:,:,0].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_b  layer 2')
plt.pcolormesh(grid.x,grid.y,nobk3d_b[:,:,1].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_b  layer 3')
plt.pcolormesh(grid.x,grid.y,nobk3d_b[:,:,2].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,4)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_b layer 4')
plt.pcolormesh(grid.x,grid.y,nobk3d_b[:,:,3].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.subplot(2,3,5)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_b layer 5')
plt.pcolormesh(grid.x,grid.y,nobk3d_b[:,:,4].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
# plt.subplot(2,3,6)
# plt.xlabel('$x$')
# plt.ylabel("$y$")
# plt.title('nobk3d_b  layer 6')
# plt.pcolormesh(grid.x,grid.y,nobk3d_b[:,:,5].T, shading=shading, cmap='RdBu_r')
# plt.colorbar()
#
plt.show()
#---------------------------------------------------------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('nobk_e 3d fields for debugging',fontsize='xx-large')
#--------------------------#
plt.subplot(2,3,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_e  layer 1')
plt.pcolormesh(grid.x,grid.y,nobk3d_e[:,:,0].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,2)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_e  layer 2')
plt.pcolormesh(grid.x,grid.y,nobk3d_e[:,:,1].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,3)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_e  layer 3')
plt.pcolormesh(grid.x,grid.y,nobk3d_e[:,:,2].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,3,4)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_e  layer 4')
plt.pcolormesh(grid.x,grid.y,nobk3d_e[:,:,3].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.subplot(2,3,5)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('nobk3d_e  layer 5')
plt.pcolormesh(grid.x,grid.y,nobk3d_e[:,:,4].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
# plt.subplot(2,3,6)
# plt.xlabel('$x$')
# plt.ylabel("$y$")
# plt.title('nobk3d_e  layer 6')
# plt.pcolormesh(grid.x,grid.y,nobk3d_e[:,:,5].T, shading=shading, cmap='RdBu_r')
# plt.colorbar()
#
plt.show()


