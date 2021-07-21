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
# plot settings 
plt.rcParams['figure.dpi'] = 180 
size    = (8,6)
shading = 'nearest'
figs    = 'figs' 
plt.close('all')
# --------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('eps field - immersed objects in flow domain',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$x$')
plt.ylabel("$y$")
plt.title('xy - plane')
plt.pcolormesh(grid.x,grid.y,eps[:,:,15].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('yz - plane')
plt.pcolormesh(grid.z,grid.y,eps[0,:,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$x$')
plt.ylabel("$z$")
plt.title('xz - plane')
plt.pcolormesh(grid.x,grid.z,eps[:,0,:].T, shading=shading, cmap='RdBu_r')
plt.colorbar()
#
plt.show()