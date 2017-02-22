#!/usr/bin/python3

import numpy as np
import struct
import matplotlib.pyplot as plt
import sys

nx = 128 # number of points in Ox
ny =  96 # number of points in Oy
nz = 128 # number of points in Oz

# do not edit
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 [xy,xz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

sizeofdata = 4 # Single precision

fin = open('grid.'+planetype, 'rb')
raw = fin.read(nx*sizeofdata)
x1 = np.array(struct.unpack('>{}f'.format(nx), raw))
nx1= nx
if   ( planetype == 'xy' ):
    raw = fin.read(ny*sizeofdata)
    x2 = np.array(struct.unpack('>{}f'.format(ny), raw))
    nx2= ny
elif ( planetype == 'xz' ):
    raw = fin.read(nz*sizeofdata)
    x2 = np.array(struct.unpack('>{}f'.format(nz), raw))
    nx2= nz
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    raw = fin.read()
    a = np.array(struct.unpack('>{}f'.format(fin.tell()/sizeofdata), raw))
    a = a.reshape((nx2,nx1))
    fin.close()

    plt.pcolormesh(x1,x2,a)
    # plt.contourf(x1,x2,a)
    plt.axis('equal')
    axes = plt.gca()
    axes.set_xlim([x1[0],x1[nx1-1]])
    axes.set_ylim([x2[0],x2[nx2-1]])
    plt.colorbar()
    plt.title(file)

    # plt.tight_layout(pad=0.1, w_pad=0.4, h_pad=0.1)
    # plt.savefig("{}.pdf".format(file))
    plt.show()
