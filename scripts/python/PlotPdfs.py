#!/usr/bin/python3

import numpy as np
import struct
import sys
import matplotlib.pyplot as plt
import os

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

# etype = ">" # big-endian
etype = "<" # little-endian

dtype = "f" # floating number
# dtype = 'B' # unsigned character, for gate files

# do not edit below this line

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 dimensions level list-of-files.")
    quit()

ndim = int(sys.argv[1])                 # 1 for univariate, 2 for bivariate...
level = int(sys.argv[2])                # -1 to plot whole map instead of individual levels
setoffiles = sorted(sys.argv[3:])

# obtain size from first file
# the first data is size and vertical coordinates
# the last level contains the global pdfs
# the last two entries per level are min/max values
fin = open(setoffiles[0], 'rb')
raw = fin.read( (1+ndim)*4 )
ny  = struct.unpack((etype+'{}i').format(1+ndim), raw)[0]
nb  = struct.unpack((etype+'{}i').format(1+ndim), raw)[1:]
raw = fin.read( ny*sizeofdata )
y   = np.array(struct.unpack((etype+'{}'+dtype).format(ny), raw))
fin.close()

print("Files with {} bins and {} levels.".format(nb,ny))

# reading data
nb_size = np.prod(list(nb)) + 2 *ndim
a = np.zeros((nb_size*(ny+1)),dtype=float)
for file in setoffiles:
    print("Processing file {} ...".format(file))
    fin = open(file, 'rb')
    fin.seek( (1+ndim)*4 + ny*sizeofdata ) # Skip the y-coordinates
    raw = fin.read()
    a   = a +np.array(struct.unpack((etype+'{}'+dtype).format(int(nb_size*(ny+1))), raw))
    fin.close()

a = a /len(setoffiles)
a = a.reshape((ny+1,nb_size))

# processing data
if ndim == 1:
    nb = nb[0]

    # normalinzing histograms to obtain pdf (s.t. it integrates to 1 using midpoint rule)
    for j in range(ny+1):
        samplesize = np.sum(a[j,:nb])
        samplestep = (a[j,nb+1]-a[j,nb]) /( nb -1 )
        if samplestep > 0: # otherwise the pdf is zero, by construction in Tlab
            a[j,:nb] = a[j,:nb] /( samplesize *samplestep )

    if level > 0:
        var1 = np.linspace(a[level-1,nb],a[level-1,nb+1],num=nb)
        plt.plot( var1, a[level-1,:nb])
        plt.xlabel("var1")
        plt.ylabel("pdf")
        if level <= ny:
            plt.title("height {:4.2f}".format(y[level-1]))
        else:
            plt.title("Global pdf")
        plt.tight_layout(pad=0.1)
        plt.savefig("{}.pdf".format(setoffiles[0]))
        plt.show()

    else:
        # axis information
        var1 = np.zeros((ny,nb),dtype=float)
        y_ex = np.zeros((ny,nb),dtype=float)
        for j in range(ny):
            var1[j,:] = np.linspace(a[j,nb],a[j,nb+1],num=nb)
            y_ex[j,:] = y[j]
        #   var1[j,:] = var1[j,:]-0.5*nbstep # colormesh uses coordinates for the corners of the region

        # choose an interval to define the color range
        nymin = int(ny/10)
        nymax = nymin *3
        levels=np.linspace(0.,np.amax(a[nymin:nymax,:nb]),num=20)

        #plt.pcolormesh(xy[0,:,:],xy[1,:,:],a[:ny,:nb],levels)
        plt.contourf(var1,y_ex,a[:ny,:nb],levels)    # this is faster
        #plt.contourf(xy[0,:,:],xy[1,:,:],np.log(a[:ny,:nb]),20)
        plt.xlabel("var1")
        plt.ylabel("height")
        plt.colorbar(label='pdf',format="%.2g")
        plt.tight_layout(pad=0.1)
        plt.savefig("{}.pdf".format(setoffiles[0]))
        plt.show()

if ndim == 2:
    if level > 0:
        var1 = np.linspace( a[level-1,nb[0]*nb[1]],   a[level-1,nb[0]*nb[1]+1], num=nb[0] )
        var2 = np.linspace( a[level-1,nb[0]*nb[1]+2], a[level-1,nb[0]*nb[1]+3], num=nb[1] )

        plt.contourf( var1, var2, a[level-1,:nb[0]*nb[1]].reshape(nb[1],nb[0]) )
        plt.xlabel("var1")
        plt.ylabel("var2")
        if level <= ny:
            plt.title("height {:4.2f}".format(y[level-1]))
        else:
            plt.title("Global pdf")
        plt.colorbar(label='pdf',format="%.2g")
        plt.tight_layout(pad=0.1)
        plt.savefig("{}.pdf".format(setoffiles[0]))
        plt.show()

    else:
        print("Undveloped option.")
        quit()
