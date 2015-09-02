import numpy
import struct
from pylab import *
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import colors
#import matplotlib as mpl

ion()

# data to be edited
#
#nx = 5120/2
nx = 3619
ny = 840

filenames = ["rcr93500.C11"]#, "rcr93500.Cww"]

scalel = (6.6666/5120.) *(0.005*70000**3)**(1./4.) 
scaley = (0.005*70000**3)**(1./4.)

# read files and accumulate the data
#
a = numpy.zeros((ny*nx))
for file in filenames:
    print "Adding data from file", file
    fin = open(file,"rb")
    raw = fin.read(nx*ny*4)
    a = a + numpy.array(struct.unpack('>{}f'.format(nx*ny), raw))
    fin.close()
    
a = a.reshape((ny,nx))
a = a /len(filenames)

# axis information
#
x = numpy.zeros((nx))
for index in range (len(x)):
    x[index]=(index-1)*scalel

y = numpy.loadtxt("y.dat")
y = y *scaley

# obtaining pseudo- integral length scales
#
l = numpy.zeros(ny)
for j in range (len(y)):
    l[j] = scalel *( a[j,:].sum() - 0.5 )

# plotting
#
rc('text', usetex=True)
#rc('font', family='sans-serif')
matplotlib.rcParams.update({'font.size': 30})

#colors = ['White', 'DarkOrange', 'DarkRed', 'LightGreen', 'DarkGreen', 'LightBlue', 'DarkBlue']
#colors = ['1.0','0.75','0.50','0.25','0.0']
#colors = ['1.0', '#fff5f0', '#fcbca2', '#fb6b4b', '#cb181d', '#67000d']
#nc = len(colors)
nc = 5
levels = MaxNLocator(nbins=nc).tick_values(0., 1.)
cmap = get_cmap('Blues')
#cmap = matplotlib.colors.ListedColormap(colors)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#define plot size in inches (width, height) & resolution(DPI)
fig=figure(figsize=(8, 10), dpi=100)
ax=fig.add_subplot(111)
pcolormesh(x[1:], y[1:], a[1:,1:], cmap=cmap, norm=norm)
plot(l[1:],y[1:],'k',linewidth=3)
xlim(1., 3000.) #max(z[:nx-2]))
ylim(1., 3000.) #max(y[1:]))

ax.set_xscale("log")
ax.set_yscale("log")

xlabel("$\lambda^+$")
ylabel("$z^+$")

# cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
# cb = colorbar(ax, cax = cbaxes)
colorbar(orientation='horizontal')

grid(linewidth=2)
minorticks_off()
tick_params(bottom='off',top='off',left='off',right='off')

show()
