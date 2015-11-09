import netCDF4 as nc   # Reading NetCDF4 files.
import numpy   as np   # For array operations.
import pylab   as pl   # For plotting.
from matplotlib import rc

pl.ion()

datafile_nc = nc.Dataset('./avg80100-128000.nc', 'r')

t=datafile_nc.variables['t'][:]
y=datafile_nc.variables['y'][:]
f=datafile_nc.variables['Tke'][:,:] # the first index is time, the second is vertical node

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Times']})
rc('font', size=18)
rc('xtick', direction='out')
rc('ytick', direction='out')
rc('xtick.major', size=5)
rc('xtick.minor', size=2)
rc('ytick.major', size=5)
rc('ytick.minor', size=2)

pl.figure(figsize=(12,10))
ax = pl.subplot(111)

# plot 10 profiles in the given range of times
for it in range(0,np.size(t),np.size(t)/10):
    pl.plot(f[it,:],y)

pl.xlabel(r'Tke')
pl.ylabel(r'y')
ax.spines['right'].set_visible(False)
ax.get_yaxis().tick_left()
ax.spines['top'].set_visible(False)
ax.get_xaxis().tick_bottom()

pl.show()
