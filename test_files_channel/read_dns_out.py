import numpy as np
import sys
import matplotlib.pyplot as plt
from   matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# import os
import my_pylib as mp
# import netCDF4  as nc 
# import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 
# import os
# from scipy import integrate
# import matplotlib.colors as mcolors
#---------------------------------------------------------------------------#
# path to dns.out file
# path = str(os.path.dirname(__file__) + '/../test/' )
path  = '/home/jonathan/simulations/data/channel_flow/cpg_test/'

# plot settings 
font = {'family':'monospace', 'weight':'bold', 'size':14}; rc('font', **font) # font
plt.rcParams['figure.dpi'] = 210 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid 
grid = mp.DnsGrid(path+'grid')

# read dns_out
out = mp.DnsOut(path+'dns.out')
# sys.exit()
#---------------------------------------------------------------------------#
# plots

# plot ubulk - iteration steps
plt.figure(figsize=size)
plt.title(r'streamwise bulk velocity')
plt.xlabel(r"iterations")
plt.ylabel(r'$u_b$')
plt.ylim(0.64,0.69)
plt.plot(out.data[:,0], out.data[:,8],color='blue',label=r'$u_b$')
plt.hlines(y=out.data[:,8].mean(), xmin=out.data[0,0], xmax=out.data[-1,0], color='blue',  linestyle='--', label=r'$\overline{u}_{b}$')
plt.hlines(y=2./3.,                xmin=out.data[0,0], xmax=out.data[-1,0], color='green', linestyle='--', label=r'$u_{b,ini}$')
plt.legend(loc=1)
plt.grid(True)
plt.show()

# # plot ubulk - time
# plt.figure(figsize=size)
# plt.title(r'streamwise bulk velocity')
# plt.xlabel(r'time')
# plt.ylabel(r'$u_b$')
# plt.ylim(0.6,0.7)
# plt.plot(out.data[:,1], out.data[:,8],label=r'$u_b$')
# plt.legend(loc=1)
# plt.grid(True)
# plt.show()

# # plot time_dt
# plt.figure(figsize=size)
# plt.title(r'time steps')
# plt.xlabel(r'iterations')
# plt.ylabel(r'$\Delta t$')
# plt.ylim(0.01,0.03)
# plt.plot(out.data[:,0], out.data[:,2],label=r'$\Delta t$')
# plt.legend(loc=1)
# plt.grid(True)
# plt.show()

# plot ubulk - iteration steps
plt.figure(figsize=size)
plt.title(r'Forcing')
plt.xlabel(r"iterations")
plt.ylabel(r'$F_1$')
plt.ylim(0.0005,0.003)
plt.plot(out.data[:,0], out.data[:,9], color='blue',label=r'$F_1$')
plt.hlines(y=out.data[:,9].mean(), xmin=out.data[0,0], xmax=out.data[-1,0], color='blue', linestyle='--', label=r'$\overline{F}_{1}$')
plt.legend(loc=1)
plt.grid(True)
plt.show()