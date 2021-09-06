import numpy as np
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
# ProfileVelocity    = parabolic
VelocityX            = 0.0
YCoorVelocity        = 0.5#1.0	# y-coordinate of profile reference point, relative to the total scale, equation (2.1).
DeltaVelocity        = 1	# Reference profile difference, equation (2.1).

# compute ThickVelocity to ensure u(y=0)=0
ThickVelocity        = grid.y.max() * YCoorVelocity / (2 * np.sqrt(1 - VelocityX / DeltaVelocity))
# in this case: ThickVelocity = 0.5 * grid.y.max() * YCoorVelocity
print('--------------------------------------------------')
print('ThickVelocity [u(0)=0]: ', str(round(ThickVelocity,5))) # Reference profile thickness, equation (2.1).

# compute u_profile
y_start   = grid.y[0]
y_end     = grid.y[-1]
ycenter   = y_start + y_end * YCoorVelocity
yrel      = grid.y - ycenter
xi        = yrel / ThickVelocity

amplify   = ( 1 - xi/2 ) * ( 1 + xi/2 )

u_profile = VelocityX + DeltaVelocity * amplify 

#---------------------------------------------------------------------------#
# plot u(y) -- ini profile
plt.figure(figsize=size)
plt.xlabel("u_mean ini-velocity")
plt.ylabel("y")
plt.xlim(0,u_profile.max())
plt.ylim(0,grid.y.max())
plt.grid('True')
plt.plot(u_profile, grid.y, marker='.',label='u_mean')
plt.legend(loc=1)
plt.show()