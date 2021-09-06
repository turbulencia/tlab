import numpy as np
import sys
import matplotlib.pyplot as plt
from   matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import os
import my_pylib as mp
import netCDF4  as nc 
import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 
import os
#-----------------------------------------------------------------------------#
# file
index = 0 # all from 200000 - 500000
# path = str(os.path.dirname(__file__) + '/' )
name  = 'avg'
path  = '/home/jonathan/simulations/data/channel_flow/cpg/'
file  = str(path+name+str(index)+'.nc')

# plot settings 
font = {'family':'monospace', 'weight':'bold', 'size':14}; rc('font', **font) # font
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')

#-----------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read avg.nc 
avg  = nc.Dataset(file, 'r')
# avg.variables.keys() # see all variables stored

# flow params (output from reading flow field)
re_tau_ini = 209.99999101975882
re_cl      = 5035.0230425920145
visc       = 0.00019860882294695577
fr         = 1.0
ro         = 1.0
rho        = 1.0
# delta      = y[-1]/2 # delta = 1
#-----------------------------------------------------------------------------#
# load variables - the first index is time, the second is vertical node

# time and vertical distance
t  = avg.variables['t'][:]
dt = np.diff(t)
y  = avg.variables['y'][:]
 
# velocities
u = avg.variables['rU'][:,:].mean(axis=0)
v = avg.variables['rV'][:,:].mean(axis=0)
w = avg.variables['rW'][:,:].mean(axis=0)

# Re-stresses
uu = avg.variables['Rxx'][:,:].mean(axis=0)
vv = avg.variables['Ryy'][:,:].mean(axis=0)
ww = avg.variables['Rzz'][:,:].mean(axis=0)
uv = avg.variables['Rxy'][:,:].mean(axis=0)
uw = avg.variables['Rxz'][:,:].mean(axis=0)
vw = avg.variables['Ryz'][:,:].mean(axis=0)

# TKE
tke = avg.variables['Tke'][:,:]

# pressure
p = avg.variables['rP'][:,:]

#-----------------------------------------------------------------------------#
# plots - averages

# u - velocity
plt.figure(figsize=size)
plt.title(r'Streamwise velocity')
plt.xlabel(r'$y$')
plt.ylabel(r'$u_i$')
plt.ylim(0,1.1)
plt.xlim(0,2)
plt.plot(y, upar, color='orange', label=r'$\overline{u}_{ini}$')
plt.plot(y, u, color='blue', label=r'$\overline{u}$')
plt.axhline(y=ub, xmin=0, xmax=2,color='orange', linestyle='--', label=r'$u_{b}$')
plt.axhline(y=2./3., xmin=0, xmax=2,color='blue', linestyle='--', label=r'$u_{b,ini}$')
plt.legend(loc=1)
# plt.xscale('log')
plt.grid(True)
plt.show()
# plt.savefig(".svg")

# Reynolds stresses
plt.figure(figsize=size)
plt.title(r'Reynolds stresses')
plt.xlabel(r'$y$')
plt.ylabel(r'$\overline{u_i^\prime u_j^\prime}$')
# plt.ylim(0,1.1)
plt.xlim(0,2)
plt.plot(y, uu, label=r'$\overline{u^\prime u^\prime}$')
plt.plot(y, vv, label=r'$\overline{v^\prime v^\prime}$')
plt.plot(y, ww, label=r'$\overline{w^\prime w^\prime}$')
plt.plot(y, uv, label=r'$\overline{u^\prime v^\prime}$')
plt.plot(y, uw, label=r'$\overline{u^\prime w^\prime}$')
plt.plot(y, vw, label=r'$\overline{v^\prime w^\prime}$')
plt.legend(loc=1)
# plt.xscale('log')
plt.grid(True)
plt.show()
# plt.savefig(".svg")