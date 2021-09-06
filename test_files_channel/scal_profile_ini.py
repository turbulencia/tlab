# import numpy as np
# import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 

#-----------------------------------------------------------------------------#
# path 
path    = str(os.path.dirname(__file__) + '/../test_little_channel/' )

# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')

#-----------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

#---------------------------------------------------------------------------#
# dns.ini settings 

# ProfileScalar1=Linear
ThickScalar1 = 2
DeltaScalar1 = 1.0
YCoorScalar1 = 0.0
MeanScalar1  = 1.0

# compute phi_profile
y_start   = grid.y[0]
y_end     = grid.y[-1]
ycenter   = y_start + y_end * YCoorScalar1
yrel      = grid.y - ycenter
xi        = yrel / ThickScalar1

amplify   = - xi

s_profile = MeanScalar1 + DeltaScalar1 * amplify 

#---------------------------------------------------------------------------#
# plot phi(y) -- ini profile
plt.figure(figsize=size)
plt.xlabel("phi_mean ini-velocity")
plt.ylabel("y")
plt.xlim(s_profile.min(),s_profile.max())
plt.ylim(0,grid.y.max())
plt.grid('True')
plt.plot(s_profile, grid.y, marker='.',label='phi_mean')
plt.legend(loc=1)
plt.show()