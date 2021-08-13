# %%
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
import netCDF4  as nc 

#---------------------------------------------------------------------------#
# read grid and flow fields
path = str(os.path.dirname(__file__) + '/../test_little_channel/' )
grid = mp.DnsGrid(path+'grid')

# %%
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#-----------------------------------------------------------------------------#
# plot vertical grid spacing
dy = grid.y[1:] - grid.y[:-1]
plt.figure(figsize=size)
plt.grid(True)
plt.xlim(0,grid.ny)
plt.ylim(0,dy.max())
plt.xlabel("nodes")
plt.ylabel("delta_y -- grid spacing")
plt.plot(np.arange(1,grid.ny), dy, marker='.',label='dy')
plt.legend(loc=1)
plt.show()
#-----------------------------------------------------------------------------#
# plot vertical grid distribution of nodes
plt.figure(figsize=size)
plt.grid(True)
plt.xlim(0,grid.ny)
plt.ylim(0,grid.y.max())
plt.xlabel("nodes")
plt.ylabel("y-position")
plt.plot(np.arange(0,grid.ny), grid.y, marker='.',label='y')
plt.legend(loc=1)
plt.show()
# %%
