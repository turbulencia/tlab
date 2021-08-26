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
path_out  = str(os.path.dirname(__file__) + '/' ) # name = 'flow.20.1' # file = str(path+name)

index = 1000
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
    itera = np.zeros(len(dns_out)-header_dns_out)
    i = 0
    dns_out = np.delete(dns_out,range(header_dns_out), axis=0)
    for line in dns_out[:]:
        data     = line.split()
        ubulk[i] = float(data[-1])
        itera[i] = int(data[1])
        i += 1
        
# plot ubulk
plt.figure(figsize=size)
plt.title('streamwise bulk velocity')
plt.xlabel("iterations")
plt.ylabel("ubulk")
plt.plot(itera, ubulk, 'o',label="ubulk")
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
ub = (1 / grid.y.max()) * integrate.simpson(um, grid.y)
print('--------------------------------------------------')
print('bulk velocity:         ', ub)

# analytical solution
ucl      = 1                                            # centerline velocity 
ycl      = grid.y.max() / 2                             # centerline position
u_par    = - (ucl / ycl**2 ) * (grid.y - ycl)**2 + ucl # parabolic ini velocity profile
ub_exact = (2/3) * ucl                   # exact bulk velocity
print('bulk velocity (exact): ', ub_exact)
print('error [%]:             ', (ub - ub_exact)/ub_exact * 100)
print('--------------------------------------------------')
#---------------------------------------------------------------------------#