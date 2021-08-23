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
index = 0
#---------------------------------------------------------------------------#
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
# shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

# original field
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()

# read ubulk from dns.out
header_dns_out = 3 
with open(path+'dns.out', 'r') as f:
    dns_out = f.readlines()
    ubulk = np.zeros(len(dns_out)-header_dns_out)
    i = 0
    for line in dns_out[header_dns_out:]: 
        data     = line.split()
        ubulk[i] = float(data[-1])
        i += 1
        
# plot ubulk
plt.figure(figsize=size)
plt.title('streamwise bulk velocity')
plt.xlabel("iterations")
plt.ylabel("ubulk")
plt.plot(np.arange(ubulk.size), ubulk, 'o',label="delta1")
plt.legend(loc=1)
plt.grid(True)
plt.show()
        
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
#---------------------------------------------------------------------------#
sys.exit()
#---------------------------------------------------------------------------#
# some contant flow rate calculation
ucl2   = 0.3
end    = 100 
delta1 = np.zeros(end) 
delta2 = np.zeros(end) 

# for i in range(end):
#     print('interative step:  ', i)
#     u_par2    = - (ucl2 / ycl**2 ) * (grid.y - ycl)**2 + ucl2
#     ub2       = (2/3) * grid.y.max() * ucl2
#     print('bulk velocity current:', ub2)
#     diff = (ub_exact - ub2)
    
#     u_par3 = np.zeros(grid.ny)
#     u_par3[:] = u_par2[:] + diff
    
#     ub3 = integrate.simpson(u_par3, grid.y)
#     print('bulk velocity new:    ', ub3)
    
    
#     ucl2 = ucl2 + diff
    
#     delta1[i] = diff
u_par2    = - (ucl2 / ycl**2 ) * (grid.y - ycl)**2 + ucl2
for i in range(end):
    print('interative step:  ', i)
    ub2       = (2/3) * grid.y.max() * ucl2
    print('bulk velocity current:', ub2)
    diff = (ub_exact - ub2)
    
    u_par3 = np.zeros(grid.ny)
    u_par3[:] = u_par2[:] + diff
    
    ub3 = integrate.simpson(u_par3, grid.y)
    print('bulk velocity new:    ', ub3)
    
    
    ucl2 = ucl2 + diff
    u_par2 = u_par3

    
    delta1[i] = diff

ucl2      = 0.3
u_par2    = - (ucl2 / ycl**2 ) * (grid.y - ycl)**2 + ucl2
for i in range(end):
    print('interative step:  ', i)
    ub2       = (2/3) * grid.y.max() * ucl2
    print('bulk velocity:', ub2)
    diff = (ub_exact - ub2)
    
    ucl_diff = (3/2) * diff / grid.y.max()
    
    u_par_diff = - (ucl_diff / ycl**2 ) * (grid.y - ycl)**2 + ucl_diff
    
    ucl2   = ucl2 + diff
    u_par2 = u_par2 + u_par_diff
    
    delta2[i] = diff
    
#---------------------------------------------------------------------------#
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# 2d plot
plt.figure(figsize=size)
plt.title('...')
plt.xlabel("x")
plt.ylabel("y")
plt.plot(np.arange(end), delta1, 'o',label="delta1")
plt.plot(np.arange(end), delta2, 'x',label="delta2")
plt.legend(loc=1)
plt.show()


