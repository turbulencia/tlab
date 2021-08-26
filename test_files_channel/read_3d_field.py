import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
from scipy import integrate
# import matplotlib.colors as mcolors
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 

#-----------------------------------------------------------------------------#
# path to flow fields
path    = str(os.path.dirname(__file__) + '/../test_little_channel/' )

# index
index_flow = 1

# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')

#-----------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read flow fields
flow = mp.Field(path,var='flow',index=index_flow)
flow.read_3d_field()

#  read field 
# f = open(path +'flow.rand.1','rb')
# f.seek(52,0)
# field = np.fromfile(f, np.dtype('<f8'), grid.nx*grid.ny*grid.nz)
# field = field.reshape((grid.nx,grid.ny,grid.nz),order='F')
# f.close()

# simulation propteries
re_cl      = 4226                 # centerline Re-number
re_tau_lam = re_cl**0.88 * 0.116  # friction   Re-number (lam.)
nu         = 1 / re_cl            # kinematic viscosity
# ro    = 1                    # Rossby number
# pr    = 1                    # Prandtl number
# rho   = 1                    # density

#---------------------------------------------------------------------------#
# bulk velocity of the flow
um = flow.u.mean(axis=(0,2)) # average in x
ub = integrate.simpson(um, grid.y)
print('--------------------------------------------------')
print('Computing bulk velocity ...')
print('bulk velocity:         ', ub)

# analytical solution
ucl      = 1                                            # centerline velocity 
ycl      = grid.y.max() / 2                             # centerline position
u_par    = - (ucl / ycl**2 ) * (grid.y - ycl)**2 + ucl  # parabolic ini velocity profile
ub_exact = (2/3) * grid.y.max() * ucl                   # exact bulk velocity
print('bulk velocity (exact): ', ub_exact)
print('error [%]:             ', (ub - ub_exact)/ub_exact * 100)
print('--------------------------------------------------')

#---------------------------------------------------------------------------#
# plots - 2d flow fields

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

plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity w')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.w[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()

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