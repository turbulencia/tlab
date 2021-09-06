import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# import netCDF4  as nc 

#---------------------------------------------------------------------------#
# read grid and flow fields
# path = str(os.path.dirname(__file__) + '/../test_little_channel/' )
path = str(os.path.dirname(__file__) + '/../test_remesh/' )
grid = mp.DnsGrid(path+'grid')

# grid spacing
dy = np.diff(grid.y)
dx = np.diff(grid.x)
# dy = grid.y[1:] - grid.y[:-1]
# dx = grid.x[1:] - grid.x[:-1]

# positions of mid nodes (dy is plotted here)
ym = (grid.y[:-1] + grid.y[1:]) / 2 

# print information
print('--------------------------------------------------')
print('stretched grid information in vertical direction')
print('origin   :', grid.y[0])
print('end      :', grid.y[-1])
print('min step :', dy.min())
print('max step :', dy.max())
print('min/max  :', dy.min()/dy.max())
print('--------------------------------------------------')
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 210 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#-----------------------------------------------------------------------------#
# plot vertical grid spacing
plt.figure(figsize=size)
plt.grid(True)
plt.xlim(0,ym.max())
plt.ylim(0,1.2*dy.max())
plt.xlabel("y mid-node postions")
plt.ylabel("delta_y")
plt.plot(grid.y[0],dy[0],   'go', color='red',   label='origin')
plt.plot(grid.y[-1],dy[-1], 'go', color='black', label='end')
plt.plot(ym, dy, marker='.',label='dy')
plt.legend(loc=1)
plt.show()
#-----------------------------------------------------------------------------#
'''# plot vertical grid spacing in visc. units (channel flow)
ymax_new = 2
ynew     = grid.y/grid.y[-1]*ymax_new

# in case of channel flow
re_tau = 180

# grid spacing
dxnew = np.diff(grid.x) * re_tau
# dynew = np.diff(grid.y) * re_tau
dynew = np.diff(ynew) * re_tau
dznew = np.diff(grid.z) * re_tau

# positions of mid nodes (dy is plotted here)
ymnew = (ynew[:-1] + ynew[1:]) / 2 


print('stretched grid information in viscous units with remeshing')
print('with remeshing y_max = 2*delta (delta=1)')
print('min step :', dynew.min())
print('max step :', dynew.max())
print('--------------------------------------------------')

plt.figure(figsize=size)
plt.grid(True)
plt.xlim(0,2)
plt.ylim(0,1.2*dynew.max())
plt.xlabel("y mid-node postions / delta")
plt.ylabel("delta_y^+")
plt.plot(ynew[0],dynew[0],   'go', color='red',   label='origin')
plt.plot(ynew[-1],dynew[-1], 'go', color='black', label='end')
plt.plot(ymnew, dynew, marker='.',label='dy')
plt.legend(loc=1)
plt.show()
'''
