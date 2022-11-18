import netCDF4 as nc            # Reading NetCDF4 files.
import numpy   as np            # Array operations.
import matplotlib.pyplot as plt # Plot data

from matplotlib import rc       # Globally setup figure options
rc('text',       usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font',       family='serif', size=12)
rc('grid',       linestyle='dotted')
rc('axes',       grid=True)

# At home, screen 27''
rc('figure',     dpi=200)
rc('savefig',    dpi=100)

# pl.ion()

data = nc.Dataset('./avg0.nc', 'r')
data.set_auto_mask(False)

t=data.variables['t'][:]
y=data.variables['y'][:]
f=data.variables['Tke'][:,:] # the first index is time, the second is vertical node

plt.figure(figsize=(5,4))
axes = plt.gca()

its = [ 0 ]                                         # Plot just the first time
#its = list(range(0,np.size(t),int(np.size(t)/10)))  # plot 10 profiles in the given range of times
for it in its:
    plt.plot(f[it,:],y)

plt.xlabel(r'function $f(y)$')
plt.ylabel(r'height $y$')
axes.spines['right'].set_visible(False)
axes.get_yaxis().tick_left()
axes.spines['top'].set_visible(False)
axes.get_xaxis().tick_bottom()

plt.show()
