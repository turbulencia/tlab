import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import my_pylib as mp


#---------------------------------------------------------------------------#
# read data txt
path     = str(os.path.dirname(__file__) + '/../test_parallel/' ) 

with open(path+'nobj_pro0.txt', 'r') as f:  
    data = f.readlines()
    nobj0 = np.zeros(len(data))
    
    for i in range(0, len(data)):    
        nobj0[i] = data[i]
    del data
with open(path+'nobj_pro1.txt', 'r') as f:  
    data = f.readlines()
    nobj1 = np.zeros(len(data))
    
    for i in range(0, len(data)):    
        nobj1[i] = data[i]
    del data
with open(path+'nobj_pro2.txt', 'r') as f:  
    data = f.readlines()
    nobj2 = np.zeros(len(data))
    
    for i in range(0, len(data)):    
        nobj2[i] = data[i]
    del data
with open(path+'nobj_pro3.txt', 'r') as f:  
    data = f.readlines()
    nobj3 = np.zeros(len(data))
    
    for i in range(0, len(data)):    
        nobj3[i] = data[i]
    del data


nobj0 = nobj0.reshape((96,32),order='F')
nobj1 = nobj1.reshape((96,32),order='F')
nobj2 = nobj2.reshape((96,32),order='F')
nobj3 = nobj3.reshape((96,32),order='F')

# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 200 
size    = (8,6)
shading = 'nearest'
figs    = 'figs' 
plt.close('all')
#--------------------------#
plt.figure(figsize=size)
plt.rcParams.update({'font.size':11})
plt.suptitle('Immersed objects in flow domain',fontsize='xx-large')
#--------------------------#
plt.subplot(2,2,1)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('nobj - pro0')
plt.pcolormesh(grid.z[:32],grid.y,nobj0, shading=shading, cmap='RdBu_r')
plt.colorbar()
#--------------------------#
plt.subplot(2,2,2)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('nobj - pro1')
plt.pcolormesh(grid.z[:32],grid.y,nobj1, shading=shading, cmap='RdBu_r')
plt.colorbar()
# #--------------------------#
plt.subplot(2,2,3)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('nobj - pro2')
plt.pcolormesh(grid.z[:32],grid.y,nobj2, shading=shading, cmap='RdBu_r')
plt.colorbar()
# #--------------------------#
plt.subplot(2,2,4)
plt.xlabel('$z$')
plt.ylabel("$y$")
plt.title('nobj - pro3')
plt.pcolormesh(grid.z[:32],grid.y,nobj3, shading=shading, cmap='RdBu_r')
plt.colorbar()
# #
plt.show()