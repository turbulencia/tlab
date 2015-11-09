import netCDF4 as nc   # Reading NetCDF4 files.
import numpy   as np   # For array operations.
import pylab   as pl   # For plotting.

datafile_nc = nc.Dataset('./avg80100-128000.nc', 'r')

t=datafile_nc.variables['t'][:]
y=datafile_nc.variables['y'][:]
f=datafile_nc.variables['Tke'][:,:] # the first index is time, the second is vertical node

# plot 10 profiles in the given range of times
for it in range(0,np.size(t),np.size(t)/10):
    pl.plot(f[it,:],y)

pl.show()
