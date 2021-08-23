import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
from scipy import integrate
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_little_channel/' ) # name = 'flow.20.1' # file = str(path+name)
index = 100

#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

# original field
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()

#  read field 
# f = open(path +'field.1','rb')
# f.seek(52,0)
# field = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
# field = field.reshape((grid.nx,grid.ny,grid.nz),order='F')
# f.close()

#---------------------------------------------------------------------------#
# bulk velocity of the flow
um = flow.u.mean(axis=(0,2)) # average in x
ub = integrate.simpson(um, grid.y)
print('=========================================')
print('bulk velocity:        ', ub)

# analytical solution
ucl      = 1                                            # centerline velocity 
ycl      = grid.y.max() / 2                             # centerline position
u_par    = - (ucl / ycl**2 ) * (grid.y - ycl)**2 + ucl # parabolic ini velocity profile
ub_exact = (2/3) * grid.y.max() * ucl                   # exact bulk velocity
print('bulk velocity (exact):', ub_exact)
print('error [%]:            ', (ub - ub_exact)/ub_exact * 100)
print('=========================================')
# sys.exit()
#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# 2d plot
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity v')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.v[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
# plt.figure(figsize=size)
# plt.title('2d-plot -- xy-plane -- velocity w')
# plt.xlabel("x")
# plt.ylabel("y")
# plt.pcolormesh(grid.x,grid.y,flow.w[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()

# #---------------------------------------------------------------------------#
# # u-mean
# plt.figure(figsize=size)
# plt.xlabel("u_mean-velocity")
# plt.ylabel("y")
# plt.xlim(0,flow.u.mean(axis=(0,2)).max())
# plt.ylim(0,grid.y.max())
# plt.grid('True')
# plt.plot(flow.u.mean(axis=(0,2)), grid.y, marker='.',label='u_mean')
# plt.legend(loc=1)
# plt.show()

# # v-mean
# plt.figure(figsize=size)
# plt.xlabel("v_mean-velocity")
# plt.ylabel("y")
# # plt.xlim(0,flow.v.mean(axis=(0,2)).max())
# plt.ylim(0,grid.y.max())
# plt.grid('True')
# plt.plot(flow.v.mean(axis=(0,2)), grid.y, marker='.',label='v_mean')
# plt.legend(loc=1)
# plt.show()