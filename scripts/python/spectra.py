import numpy
import struct
from matplotlib.pyplot import *
from matplotlib.colors import *
from matplotlib.ticker import MaxNLocator

ion()

# data to be edited
#
#nx = 5120/2
nx = 3619
ny = 1024

hp = 679

filenames = ["rsp123300.Ev1"]#, "rsp123300.Eww"]
tag       = '(f)'
fname     = 'E1wStable.png'

# colors and levels
#
#colors = ['0.75', 'White', '#f1eef6', '#d0d1e6', '#a6bddb', '#74a9cf', '#2b8cbe', '#045a8d'] #PuBu
colors = ['0.75','White', 'DarkSalmon', '#aa0000', 'LightBlue', 'SteelBlue', 'LightGreen', '#008c00']
#colors =        ['White', 'DarkSalmon', '#aa0000', 'LightBlue', 'SteelBlue', 'LightGreen', '#008c00']#, '#fec44f', '#d95f0e']
cmap = matplotlib.colors.ListedColormap(colors)
#cmap = get_cmap('PuBu')

#levels =         [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] # auto
#levels =         [ -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5] # auto
#levels =         [ 0.0, 1., 4., 7., 10., 13., 16., 19., 22., 25.] # auto
levels = [ -0.1, -0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] # cross
nc = len(levels)-1

# scales to normalize data
#
length = 1./(0.005*70000**3)**(1./4.)     # length scales
scalek = length
scaley = 1. /length
scalel = 6.6666 /length

#scalee = 1. /( (0.005*length) **(2./3.) ) # psd scale; can be over-written below to change normalization
#scalee = 1. /0.005

# read files and accumulate the data
#
a = numpy.zeros((ny*nx))
for file in filenames:
    print "Reading data from file", file,"..."
    fin = open(file,"rb")
    raw = fin.read(nx*ny*4)
    a = a + numpy.array(struct.unpack('>{}f'.format(nx*ny), raw))
    fin.close()
    
a = a.reshape((ny,nx))
a = a /len(filenames)

# Files contain only half of the spectra
#
a = a *2.

# obtaining the total energy
#
print "Calculating energy..."
e = numpy.zeros(ny)
for j in range (len(e)):
    e[j] = a[j,0] + a[j,1:].sum()

scalee = 1. /abs(e).max() # normalizing with maximum energy
    
# filter
#
print "Smoothing with Daniell window of size 3..."
scratch = numpy.zeros((ny,nx))
for j in range (ny):
    for i in range (1,nx-1):
        scratch[j,i] = ( a[j,i-1] + 2. * a[j,i] + a[j,i+1] ) / 4.
# print "Smoothing with Daniell window of size 11..."
# scratch = numpy.zeros((ny,nx))
# for j in range (ny):
#     for i in range (5,nx-5):
#         scratch[j,i] = ( a[j,i-5] + 2.*a[j,i-4] + 2.*a[j,i-3] + 2.*a[j,i-2] + 2.*a[j,i-1] + 2.*a[j,i] + 2.*a[j,i+1] + 2.*a[j,i+2] + 2.*a[j,i+3] + 2.*a[j,i+4] + a[j,i+5] ) /20.

a = scratch

# check transformation
# 
aux = numpy.zeros(ny)
for j in range (len(aux)):
    aux[j] = a[j,0] + a[j,1:].sum()

error = numpy.zeros(ny)
for j in range (len(error)):
    if e[j] == 0.:
        error[j] = 0.
    else:
        error[j] = abs(aux[j]-e[j])/e[j]
print 'Relative error after smoothing:', error.max()

e = aux

scalee = 1./abs(e).max() # normalizing with maximum energy
print 'Maximum normalized value:', e.max() *scalee

# axis information
#
x = numpy.linspace(0, nx-1, num=nx)

y = numpy.loadtxt("y.dat")

z = numpy.zeros(nx)
for i in range (len(z)-1): # z[nx] would be infinity; we skip it
    z[i]=1./x[nx-i-1]

# creating premultiplied spectrum
#
b = numpy.zeros((ny,nx))
for j in range (len(y)):
    for i in range (len(x)):
        b[j,i] = a[j,nx-i-1]*x[nx-i-1]

# check transformation
# 
aux = numpy.zeros(ny)
for j in range (len(aux)):
    aux[j] += a[j,0]
    for i in range (len(z)-1):
        aux[j] += b[j,i]*z[i]

error = numpy.zeros(ny)
for j in range (len(error)):
    if e[j] == 0.:
        error[j] = 0.
    else:
        error[j] = abs(aux[j]-e[j])/e[j]
print "Relative error in premultiplied spectra:", error.max()

# obtaining integral length scales
#
l = numpy.zeros(ny)
for j in range (len(y)):
#    l[j] = a[j,0] /e[j]  *scalel / 2. 
#    l[j] = scalel * a[j,:].sum() / b[j,:].sum() 
    # for i in range (1,len(z)):
    #     l[j] += a[j,i]/x[i] /e[j]  *scalel / 2.
    for i in range (1,len(z)):
        l[j] += a[j,i] /x[i] / e[j]  *scalel *2.

# scaling
#
a = a *scalee 
b = b *scalee

x = x *scalek
y = y *scaley
z = z *scalel

# obtaining the maximum for a given height
#
mz = numpy.zeros(ny)
for j in range (len(y)):
    i = b[j,:].argmax() 
    mz[j] = z[i]

# obtaining the maximum for a given wavelength
#
ml = numpy.zeros(nx)
for i in range (len(z)):
    j = b[:,i].argmax() 
    ml[i] = y[j]

# z_ml = numpy.zeros((2,nx)); ml = numpy.zeros((2,nx)); i_ml = numpy.zeros(2, dtype=numpy.int)
# for i in range (len(z)-4):
#     jref = 0
#     while ( y[jref] < 0.2 *z[i] ):
#         jref += 1
#     jmin = 0; jmax = jref                           # lower branch
#     j = jmin + b[jmin:jmax+1,i].argmax()
#     if ( b[j,i] > b[j-1,i] and b[j,i] > b[j+1,i] ): # local maxima
#         z_ml[0,i_ml[0]] = z[i]; ml[0,i_ml[0]] = y[j]; i_ml[0] +=1
#     jmin = jref; jmax = ny                          # upper branch
#     j = jmin + b[jmin:jmax+1,i].argmax()
#     if ( b[j,i] > b[j-1,i] and b[j,i] > b[j+1,i] ): # local maxima
#         z_ml[1,i_ml[1]] = z[i]; ml[1,i_ml[1]] = y[j]; i_ml[1] +=1
        
# plotting array b
#
print "Plotting..."

rc('text', usetex=True)
rc('font', family='serif')
matplotlib.rcParams.update({'font.size': 35})

#levels = MaxNLocator(nbins=nc).tick_values(b.min(), b.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig=figure(figsize=(9,9), dpi=100)       #define plot size in inches (width, height) & resolution(DPI)
fig.subplots_adjust(bottom=0.1,left=0.1) #space for labels
ax=fig.add_subplot(111)
pcolormesh(z[:nx-2], y[1:], b[1:,:nx-2], cmap=cmap, norm=norm)

#colorbar(ticks=None, orientation='horizontal', aspect=25, format='%1.g')

xlim(10., 3000.)
ylim(1.,  3000.)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_aspect(1)

xlabel("$\lambda / z_\kappa$",labelpad=-4)
ax.set_xticks([10,20,30,40,50,60,70,80,90,100,
               200,300,400,500,600,700,800,900,1000,
               2000])
ax.tick_params(axis='x', pad=6, width=2)
#matplotlib.pyplot.gca().axes.get_xaxis().set_ticklabels([]) #remove labels

#ylabel("$z / z_\kappa$",      labelpad=-4)
ax.set_yticks([1,2,3,4,5,6,7,8,9,10,
               20,30,40,50,60,70,80,90,100,
               200,300,400,500,600,700,800,900,1000,
               2000])
ax.tick_params(axis='y',         width=2)
matplotlib.pyplot.gca().axes.get_yaxis().set_ticklabels([]) #remove labels

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

plot([hp,hp],   [1,3000], 'k--',linewidth=2)
plot([10,3000], [hp,hp],  'k--',linewidth=2)
plot([1,3000],  [0.2,600],'k--',linewidth=2)
#plot(l[1:],     y[1:],    'k-', linewidth=2)
#plot(mz[1:],     y[1:],    'k-', linewidth=2)
plot(z[:nx-4],   ml[:nx-4],'k-', linewidth=2)
# for n in range (len(i_ml)):
#     plot(z_ml[n,:i_ml[n]],   ml[n,:i_ml[n]],'k-', linewidth=2)

text(0.05, 0.9, tag, transform = ax.transAxes, fontweight='bold', fontsize='large')

show()

fig.savefig(fname, bbox_inches='tight', dpi=200)
#y[141]
#numpy.savetxt('/home/zmaw/m300041/develop/papers/14blm.cbl/figs/aux0.dat',numpy.transpose([z[:nx-1],b[141,:nx-1]]))
