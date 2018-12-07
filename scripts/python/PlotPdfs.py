#!/usr/bin/python3

import numpy as np
import struct
import sys
import matplotlib.pyplot as plt
import os

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

etype = ">" # big-endian
# etype = "<" # little-endian

dtype = "f" # floating number 
# dtype = 'B' # unsigned character, for gate files

# do not edit below this line

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 dimension list-of-files.")
    quit()

dimension  = sys.argv[1] # 1 for univariate, 2 for bivariate... So far, only univariate
setoffiles = sorted(sys.argv[2:])

# axis information
#
y  = np.loadtxt("y.dat")
ny = len(y)

# obtain size from first file
# the last level contains the global Pdf
# the last two entries per level are min/max values
nb = int( os.stat(setoffiles[0]).st_size /sizeofdata /(ny+1) ) -2

print("Files with {} bins and {} levels.".format(nb,ny))

# processing data
a = np.zeros(((ny+1)*(nb+2)),dtype=float)
for file in setoffiles:
    print("Processing file {} ...".format(file))
    fin = open(file, 'rb')
    raw = fin.read()
    a = a +np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
    fin.close()

a = a.reshape((ny+1,nb+2))
a = a /len(setoffiles)

# normalinzing histograms to obtain pdf (s.t. it integrates to 1 using midpoint rule)
for j in range(ny+1):
    samplesize = np.sum(a[j,:nb])
    samplestep = (a[j,nb+1]-a[j,nb]) /( nb -1 )
    a[j,:nb] = a[j,:nb] /( samplesize *samplestep )
    
# axis information
#
xy = np.zeros((2,ny,nb),dtype=float)
for j in range(ny):
    xy[0,j,:] = np.linspace(a[j,nb],a[j,nb+1],num=nb)
    xy[1,j,:] = y[j]
#   xy[0,j,:] = xy[0,j,:]-0.5*nbstep # colormesh uses coordinates for the corners of the region

# Choose an interval to define the color range
nymin = int(ny/10)
nymax = nymin *3
levels=np.linspace(0.,np.amax(a[nymin:nymax,:nb]),num=20)
#plt.pcolormesh(xy[0,:,:],xy[1,:,:],a[:ny,:nb],levels)
plt.contourf(xy[0,:,:],xy[1,:,:],a[:ny,:nb],levels)    # this is faster
#plt.contourf(xy[0,:,:],xy[1,:,:],np.log(a[:ny,:nb]),20)
plt.xlabel("variable")
plt.title("Pdf")

plt.ylabel("height")

axes = plt.gca()

plt.colorbar()

plt.tight_layout(pad=0.1)
plt.savefig("{}.pdf".format(setoffiles[0]))
plt.show()
