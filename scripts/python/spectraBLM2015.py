import numpy
import struct
from matplotlib.pyplot import *
from matplotlib.colors import *
from matplotlib.ticker import MaxNLocator

ion()

# variables
directory = '/scratch/local1/m300041/CONVECTION/CBL3D/5120x1024x5120-70000-20141002/postprocess.timeaveraged/'
filenames = ["rsp123300.Evv"]#, "rsp123300.Eww"]
#directory = '/scratch/local1/m300041/CONVECTION/PLATE3D/5120x0840x5120-70000-20131021/postprocess.timeaveraged/'
#filenames = ["rsp7000.Euu", "rsp7000.Eww"]
flag_max  = 1 # number of local maxima to look for

nx = 3619
hp = 679

#mode = "frame"
#mode = "lines"
mode = "data"
data_row = 1; data_col = 2 # counter start at 0

# constants
#
tags = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']
rows = 2
cols = (len(tags)-1) /rows +1

#colors = ['0.75','White', 'DarkSalmon', '#aa0000', 'LightBlue', 'SteelBlue', 'LightGreen', '#008c00']
colors = ['0.75','White', '#c49494', '#944a4a', '#b0c4de', '#4478ad', '#a3ca83', '#3f7a0d']
cmap = matplotlib.colors.ListedColormap(colors)

levels = [ -0.1, -0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] # cross spectra
nc = len(levels)-1

# axes
#
y  = numpy.loadtxt(directory+"y.dat")
ny = len(y)

x  = numpy.linspace(0, nx-1, num=nx)

# scales to normalize data
#
length = 1./(0.005*70000**3)**(1./4.)     # length scales
scalek = length
scaley = 1. /length
scalel = 6.6666 /length

# read files and accumulate the data
#
a = numpy.zeros((ny*nx))
for file in filenames:
    print "Reading data from file", file,"..."
    fin = open(directory+file,"rb")
    raw = fin.read(nx*ny*4)
    a = a + numpy.array(struct.unpack('>{}f'.format(nx*ny), raw))
    fin.close()
    
a = a.reshape((ny,nx))
a = a /len(filenames)

# Files contain only half of the spectra
#
a = a *2.

# filter
#
print "Smoothing with Daniell window of size 3..."
scratch = numpy.zeros((ny,nx))
for j in range (ny):
    for i in range (1,nx-1):
        scratch[j,i] = ( a[j,i-1] + 2. * a[j,i] + a[j,i+1] ) / 4.

a = scratch

# obtaining the total energy to normilze the spectra
#
print "Calculating energy and normalizing spectra..."
e = numpy.zeros(ny)
for j in range (len(e)):
    e[j] = a[j,0] + a[j,1:].sum()

scalee = 1. /abs(e).max() # normalizing with maximum energy

# creating premultiplied spectrum
#
b = numpy.zeros((ny,nx))
for j in range (len(y)):
    for i in range (len(x)):
        b[j,i] = a[j,nx-i-1]*x[nx-i-1]

z = numpy.zeros(nx)
for i in range (len(z)-1): # z[nx] would be infinity; we skip it
    z[i]=1./x[nx-i-1]

# scaling
#
a = a *scalee 
b = b *scalee

x = x *scalek
y = y *scaley
z = z *scalel

# obtaining the maximum for a given wavelength
#
if flag_max == 1:
    ml = numpy.zeros(nx)
    for i in range (len(z)):
        j = b[:,i].argmax() 
        ml[i] = y[j]
else:
    ml = numpy.zeros((2,nx)); z_ml = numpy.zeros((2,nx)); i_ml = numpy.zeros(2, dtype=numpy.int)
    for i in range (len(z)-4):
        jref = 0
        while ( y[jref] < 0.2 *z[i] ):
            jref += 1
        jmin = 0; jmax = jref                           # lower branch
        j = jmin + b[jmin:jmax+1,i].argmax()
        if ( b[j,i] > b[j-1,i] and b[j,i] > b[j+1,i] ): # local maxima
            z_ml[0,i_ml[0]] = z[i]; ml[0,i_ml[0]] = y[j]; i_ml[0] +=1
        jmin = jref; jmax = ny                          # upper branch
        j = jmin + b[jmin:jmax+1,i].argmax()
        if ( b[j,i] > b[j-1,i] and b[j,i] > b[j+1,i] ): # local maxima
            z_ml[1,i_ml[1]] = z[i]; ml[1,i_ml[1]] = y[j]; i_ml[1] +=1
                    
# plotting array b
#
print "Plotting "+mode+"..."

rc('text', usetex=True)
rc('font', family='serif')
matplotlib.rcParams.update({'font.size': 28})

#levels = MaxNLocator(nbins=nc).tick_values(b.min(), b.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig, ax = subplots(rows, cols, figsize=(15,10.3))
fig.subplots_adjust(left  =0.055, bottom=0.08,
                    right =1.,    top   =1.,
                    hspace=0.06 , wspace=0.04)

for itag in range (len(tags)): # 0, 1, 2...
    irow = itag /cols          # 0, 1, 2...
    icol = itag %cols          # 0, 1, 2...
    ax[irow,icol].set_xlim(10., 3000.)
    ax[irow,icol].set_ylim(1.,  3000.)
    #
    ax[irow,icol].set_adjustable("box")
    ax[irow,icol].set_xscale("log")
    ax[irow,icol].set_yscale("log")
    ax[irow,icol].set_aspect(1)
    if   mode == "data":
        if irow == data_row and icol == data_col:
            ax[irow,icol].pcolormesh(z[:nx-2], y[1:], b[1:,:nx-2], cmap=cmap, norm=norm)
        ax[irow,icol].set_axis_off()
    elif mode == "lines":
        if irow == data_row and icol == data_col:
            if numpy.shape(numpy.shape(ml))[0] == 2: # check is array ml has 1 or 2 dimensions
                for n in range (len(i_ml)):
                    ax[irow,icol].plot(z_ml[n,:i_ml[n]],   ml[n,:i_ml[n]],'k-')
            else:
                ax[irow,icol].plot(z[:nx-4],   ml[:nx-4],'k-')
        ax[irow,icol].set_axis_off()
    else:
        # X axes
        ax[irow,icol].set_xticks([10,20,30,40,50,60,70,80,90,100,
                             200,300,400,500,600,700,800,900,1000,
                             2000])
        ax[irow,icol].tick_params(axis='x', pad=6)#, width=2)
        if irow == rows-1:
            ax[irow,icol].set_xlabel("$\lambda / z_\kappa$",labelpad=-0.5)
        else:
            ax[irow,icol].axes.get_xaxis().set_ticklabels([]) #remove labels
        # Y axes
        ax[irow,icol].set_yticks([1,2,3,4,5,6,7,8,9,10,
                             20,30,40,50,60,70,80,90,100,
                             200,300,400,500,600,700,800,900,1000,
                             2000])
        ax[irow,icol].tick_params(axis='y', pad=2)
        if icol == 0:
            ax[irow,icol].set_ylabel("$z / z_\kappa$", labelpad=-4, position=(0,0.45))
        else:
            ax[irow,icol].axes.get_yaxis().set_ticklabels([]) #remove labels
        # lines
        ax[irow,icol].plot([hp,hp],   [1,3000], 'k--')
        ax[irow,icol].plot([10,3000], [hp,hp],  'k--')
        ax[irow,icol].plot([1,3000],  [0.2,600],'k--')
                    
        text(0.05, 0.9, tags[itag], transform = ax[irow,icol].transAxes, fontweight='bold',fontsize=30)#'large')
        
if mode == "frame":
    text(0.4 , 0.58, '$\lambda=5z$', rotation=45, transform = ax[1,0].transAxes)
    text(0.3 , 0.82, '$z=z_*$', transform = ax[1,0].transAxes)
    text(0.77, 0.40, '$\lambda=z_*$', rotation=90, transform = ax[1,0].transAxes)

#colorbar(ticks=None, orientation='horizontal', aspect=25, format='%1.g')

if   mode == "data":
    fig.savefig(mode+str(data_row)+str(data_col)+'.png', transparent=True, dpi=100)
elif mode == "frame":
    fig.savefig(mode+'.pdf', transparent=True, dpi=100)
else:
    fig.savefig(mode+str(data_row)+str(data_col)+'.pdf', transparent=True, dpi=100)
