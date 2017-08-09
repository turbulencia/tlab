import netCDF4 as nc   # Reading NetCDF4 files.
import numpy   as np   # For array operations.
import pylab   as pl   # For plotting.
from matplotlib import rc

#pl.ion()

datafile_nc = nc.Dataset('./tower000001-000010.nc', 'r')

t=datafile_nc.variables['t'][:]
x=datafile_nc.variables['x'][:]
y=datafile_nc.variables['y'][:]
z=datafile_nc.variables['z'][:]
f=datafile_nc.variables['v'][:,:,:,:] # indeces: 1. time
                                      #          2. z (horizontal)
                                      #          3. x (horizontal)
                                      #          4. y (vertical)
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Times']})
rc('font', size=16)
rc('xtick', direction='out')
rc('ytick', direction='out')
rc('xtick.major', size=5)
rc('xtick.minor', size=2)
rc('ytick.major', size=5)
rc('ytick.minor', size=2)

pl.figure(figsize=(6,5))
ax = pl.subplot(111)

stridex=4
stridez=4
# plot all columns
for iz in range(0,np.size(z),stridez):
    for ix in range(0,np.size(x),stridex):
        print(iz,ix)
        pl.plot(f[0,iz,ix,:],y)

pl.xlabel(r'u')
pl.ylabel(r'y')
ax.spines['right'].set_visible(False)
ax.get_yaxis().tick_left()
ax.spines['top'].set_visible(False)
ax.get_xaxis().tick_bottom()

pl.tight_layout(pad=0.1, w_pad=0.1)
pl.show()
