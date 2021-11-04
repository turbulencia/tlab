# import numpy as np
import sys
import matplotlib.pyplot as plt
from   matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import os
import my_pylib as mp
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 
# import os
#-----------------------------------------------------------------------------#
# file
index        = 38001
forcing_flow = 'cpg'
# path  = '/home/jonathan/simulations/data/channel_flow/cpg_scal/'
path = str(os.path.dirname(__file__) + '/../test_parallel/' )
# path  = '/home/jonathan/simulations/data/channel_flow/cpg_test/'

# plot settings 
font = {'family':'monospace', 'weight':'bold', 'size':14}; rc('font', **font) # font
plt.rcParams['figure.dpi'] = 210
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read flow fields
flow = mp.Field(path,var='flow',index=index, forcing=forcing_flow)
flow.read_3d_field()

# # read scalar fields
# scal1 = mp.Field(path,var='phi1',index=index, forcing=forcing_flow)
# scal1.read_3d_field()

# # read scalar fields
# scal2 = mp.Field(path,var='phi2',index=index, forcing=forcing_flow)
# scal2.read_3d_field()

# bulk velocity of the flow
ub = mp.ubulk(flow.u, grid.y)
#---------------------------------------------------------------------------#
sys.exit()
#---------------------------------------------------------------------------#
# plots - 2d flow fields
plt.figure(figsize=size)
plt.title(r'velocity u')
plt.xlabel(r'$x/\delta$')
plt.ylabel(r'$y/\delta$')
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,17].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
plt.figure(figsize=size)
plt.title(r'velocity u')
plt.xlabel(r'$z/\delta$')
plt.ylabel(r'$y/\delta$')
plt.pcolormesh(grid.z,grid.y,flow.u[10,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
plt.figure(figsize=size)
plt.title(r'velocity u')
plt.xlabel(r'$z/\delta$')
plt.ylabel(r'$y/\delta$')
plt.pcolormesh(grid.z,grid.y,flow.u.mean(axis=0), shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# # plots - 2d scalar fields
# plt.figure(figsize=size)
# plt.title(r'scalar 1')
# plt.xlabel(r'$x/\delta$')
# plt.ylabel(r'$y/\delta$')
# plt.pcolormesh(grid.x,grid.y,scal1.phi1[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()
# #
# plt.figure(figsize=size)
# plt.title(r'scalar 1')
# plt.xlabel(r'$z/\delta$')
# plt.ylabel(r'$y/\delta$')
# plt.pcolormesh(grid.z,grid.y,scal1.phi1[10,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()
# #
# plt.figure(figsize=size)
# plt.title(r'scalar 2')
# plt.xlabel(r'$x/\delta$')
# plt.ylabel(r'$y/\delta$')
# plt.pcolormesh(grid.x,grid.y,scal2.phi2[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()
# #
# plt.figure(figsize=size)
# plt.title(r'scalar 2')
# plt.xlabel(r'$z/\delta$')
# plt.ylabel(r'$y/\delta$')
# plt.pcolormesh(grid.z,grid.y,scal2.phi2[10,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()
#---------------------------------------------------------------------------#
# plot means
# u-mean
plt.figure(figsize=size)
plt.title(r'mean velcocity u')
plt.xlabel(r'$\overline{u}$')
plt.ylabel(r'$y/\delta$')
plt.xlim(0,1)
plt.ylim(0,2)
plt.grid('True')
plt.plot(flow.u.mean(axis=(0,2)), grid.y, marker='.',label=r'$\overline{u}$')
plt.legend(loc=1)
plt.show()
# # phi-mean
# plt.figure(figsize=size)
# plt.title(r'mean scalar')
# plt.xlabel(r'$\overline{\Phi}_i$')
# plt.ylabel(r'$y/\delta$')
# plt.xlim(0,1)
# plt.ylim(0,2)
# plt.grid('True')
# plt.plot(scal1.phi1.mean(axis=(0,2)), grid.y, marker='.',label=r'$\overline{\Phi}_1$')
# plt.plot(scal2.phi2.mean(axis=(0,2)), grid.y, marker='.',label=r'$\overline{\Phi}_2$')
# plt.legend(loc=1)
# plt.show()