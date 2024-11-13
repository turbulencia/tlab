#!/usr/bin/python3

import numpy as np
import struct
import sys
import os.path
import matplotlib.pyplot as plt

sizeofdata = 4  # in bytes
# sizeofdata = 1  # for gate files

# etype = ">"     # big-endian
etype = "<"     # little-endian

dtype = "f"     # floating number
# dtype = 'B'     # unsigned character, for gate files

nx = 0          # number of points in Ox; if 0, then search tlab.ini
ny = 0          # number of points in Oy; if 0, then search tlab.ini
nz = 0          # number of points in Oz; if 0, then search tlab.ini

# do not edit below this line

WhatToPlot = [ 'single', 'y-distribution']

# getting grid size from tlab.ini, if necessary
if ( nx == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("imax="):
            nx = int(line.split("=",1)[1])
            break

if ( ny == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("jmax="):
            ny = int(line.split("=",1)[1])
            break

if ( nz == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("kmax="):
            nz = int(line.split("=",1)[1])
            break

print("Grid size is {}x{}x{}.".format(nx,ny,nz))

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,zy,ry] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

# processing data
if   ( planetype == 'xy' ):     # Number of modes
    nk = int(nx /2.)
elif ( planetype == 'zy' ):
    nk = int(nz /2.)
elif ( planetype == 'ry' ):
    nk = int( min(nx /2.,nz /2.) )
else:
    print("Usage: python $0 [xy,zy,ry] list-of-files.")
    quit()

a = np.zeros((ny*nk))
for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    raw = fin.read()
    a = a +np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
    fin.close()

a = a.reshape((ny,nk))
a = a /len(setoffiles)
a = a *2.                       # Files contain only half of the spectra

# obtaining the total energy
#
print("Calculating variance...")
e = np.zeros(ny)
for j in range (len(e)):
    e[j] = a[j,0] + a[j,1:].sum()
# scalee = 1. /abs(e).max()     # normalizing with maximum energy

np.savetxt('variance.dat',e)

###############################################################################
# plot single spectra
#
tag = 'single'
if tag in WhatToPlot:
    j = int(ny/2)
    plt.figure( figsize=(4.,3.) )
    plt.plot(a[j,1:])
    plt.xlabel('mode number')
    plt.ylabel('mode energy')
    plt.yscale("log")
    plt.tight_layout(pad=0.1)
    plt.savefig("{}.pdf".format(setoffiles[0]))
    plt.show()

###############################################################################
# plot vertical distribution of spectra
#
tag = 'y-distribution'
if tag in WhatToPlot:
    if os.path.exists('y.dat'):
        y = np.loadtxt("y.dat")
    else:
        print("No grid file found, using grid node index instead.")
        y = np.linspace(1, ny, num=ny)

    k = np.linspace(0, nk-1, num=nk)    # mode number
    l = np.zeros(nk-1)                  # wavelength
    for i in range (len(k)-1):          # l[nk] would be infinity; we skip it
        l[i]=1./k[nk-i-1]

    b = np.zeros((ny,nk-1))             # creating premultiplied spectrum
    for j in range (len(y)):
        for i in range (len(k)-1):
            b[j,i] = a[j,nk-i-1]*k[nk-i-1]

    # k = k -0.5 # colormesh uses coordinates for the corners of the region
    # plt.pcolormesh(k,y,a)
    # plt.xlabel("mode number")
    # plt.title("spectrum")
    #plt.pcolormesh(l,y[1:],b[1:,:])    # skip first plane for log plot if it has height 0
    plt.contourf(l,y[1:],b[1:,:],20)    # this is faster
    plt.xlabel("normalized wavelength")
    plt.ylabel("height")
    plt.title("premultiplied spectrum")
    plt.xscale("log")
    # plt.yscale("log")
    # plt.colorbar()
    plt.tight_layout(pad=0.1)
    plt.savefig("{}-map.pdf".format(setoffiles[0]))
    plt.show()
