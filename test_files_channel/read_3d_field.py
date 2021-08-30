import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
from scipy import integrate
# import matplotlib.colors as mcolors
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 

#---------------------------------------------------------------------------#
# path to flow fields
# path    = str(os.path.dirname(__file__) + '/../test_little_channel/' )
# path    = str(os.path.dirname(__file__) + '/../test_parallel_channel/' )
path    = str(os.path.dirname(__file__) + '/../test_yamo_180/' )
 
# index
index_flow = 200001

# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')

#---------------------------------------------------------------------------#
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

# forcing type
forcing_cpg = False  # constant pressure gradient
forcing_cfr = True # constant flow rate

re = 5000           # reynoldsnumber in dns.ini file
if forcing_cpg:
    re_tau = re
    re_cl  = (re_tau/0.116)**(1/0.88)
    fcpg   = ( re_tau / re_cl )**2    # constant streamwise pressure gradient
if forcing_cfr:
    re_cl  = re 
    re_tau = re_cl**0.88 * 0.116  
    
nu    = 1 / re_cl  # kinematic viscosity, always!
# ro    = 1          # Rossby number
# pr    = 1          # Prandtl number
rho   = 1          # density
delta = 1          # channel half height
#---------------------------------------------------------------------------#
# bulk velocity of the flow
um = flow.u.mean(axis=(0,2)) # average in x
ub = (1 / grid.y.max()) * integrate.simpson(um, grid.y)
print('--------------------------------------------------')
print('bulk velocity:         ', ub)

# analytical solution
ucl      = 1 # 15.63                                       # centerline velocity 
ycl      = grid.y.max() / 2                             # centerline position
u_par    = - (ucl / ycl**2 ) * (grid.y - ycl)**2 + ucl # parabolic ini velocity profile
ub_exact = (2/3) * ucl                   # exact bulk velocity
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
plt.title('2d-plot -- yz-plane -- velocity u')
plt.xlabel("z")
plt.ylabel("y")
plt.pcolormesh(grid.z,grid.y,flow.u[10,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
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
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity w')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.w[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()

#---------------------------------------------------------------------------#
# u-mean
plt.figure(figsize=size)
plt.xlabel("u_mean-velocity")
plt.ylabel("y")
plt.xlim(0,flow.u.mean(axis=(0,2)).max())
plt.ylim(0,grid.y.max())
plt.grid('True')
plt.plot(flow.u.mean(axis=(0,2)), grid.y, marker='.',label='u_mean')
plt.legend(loc=1)
plt.show()

    
# plt.figure(figsize=size)
# plt.xlabel("u_mean-velocity")
# plt.ylabel("y")
# plt.xlim(0,1)#flow.u.mean(axis=(0,2)).max())
# plt.ylim(0,grid.y.max())
# plt.grid('True')
# for i in range(0,100,10):
#     plt.plot(flow.u[i,:,grid.nz//2], grid.y, marker='.',label=str(i))
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