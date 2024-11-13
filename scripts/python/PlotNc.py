import netCDF4 as nc            # Reading NetCDF4 files.
import numpy   as np            # Array operations.
import sys
import matplotlib.pyplot as plt # Plot data

from matplotlib import rc       # Globally setup figure options
rc('text',       usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font',       family='serif', size=12)
# rc('grid',       linestyle='dotted')
# rc('axes',       grid=True)

# At home, screen 27''
rc('figure',     dpi=200)
rc('savefig',    dpi=100)

# pl.ion()

scale = 0.840   # km

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 list-of-files.")
    quit()

setoffiles = sorted(sys.argv[1:])

file = setoffiles[0]
name = file.removesuffix('.nc')

data = nc.Dataset(file, 'r')
data.set_auto_mask(False)

x=data.variables['x'][:]
y=data.variables['y'][:]
z=data.variables['z'][:]
a=data.variables[name][:,:,:]

if ( len(z) == 1 ):
    x1 = x*scale
    x2 = y*scale
    f = a[0,:,:]

if ( len(y) == 1 ):
    x1 = x*scale
    x2 = z*scale
    f = a[:,0,:]

if ( len(x) == 1 ):
    x1 = y*scale
    x2 = z*scale
    f = a[:,:,0]

# mean, std = np.mean( a ), np.std( a )
print(np.min(f),np.max(f))
fig = plt.figure( figsize=(8,8) )
ax = plt.gca()
plt.pcolormesh(x1,x2,f,shading='auto',cmap='Blues_r',vmin=0.0,vmax=0.0005)
# plt.contourf(x1,x2,a)
# plt.axis('equal')
ax.set_aspect('equal', adjustable='box')
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_position(('data',-0.05))
ax.spines['right'].set_visible(False)
ax.spines['left'].set_position(('data',-0.05))
ax.set_xlim([x1[0],x1[-1]])
ax.set_ylim([x2[0],x2[-1]]) #*0.65])
plt.xlabel(r'distance (km)')
plt.ylabel(r'distance (km)')
# plt.colorbar()
# plt.tight_layout(pad=2)
# plt.title(name)
plt.savefig("{}.jpg".format(name),dpi=300,bbox_inches='tight')
plt.show()
