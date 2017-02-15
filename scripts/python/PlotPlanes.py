#!/usr/bin/python3

import numpy as np
import struct
import matplotlib.pyplot as plt
import sys

nx = 128 # number of points in Ox
nz = 128 # number of points in Oz

# do not edit
if ( len(sys.argv) == 1 ):
    print("Usage: python $0 list-of-files.")
    quit()

setoffiles = sorted(sys.argv[1:])

sizeofdata = 4 # Single precision

fin = open('grid.plane', 'rb')
raw = fin.read(nx*sizeofdata)
x = np.array(struct.unpack('>{}f'.format(nx), raw))
raw = fin.read(nz*sizeofdata)
z = np.array(struct.unpack('>{}f'.format(nz), raw))
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    raw = fin.read()
    a = np.array(struct.unpack('>{}f'.format(fin.tell()/sizeofdata), raw))
    a = a.reshape((nz,nx))
    fin.close()

    plt.pcolormesh(x,z,a)
    # plt.contourf(x,z,a)
    plt.axis('equal')
    axes = plt.gca()
    axes.set_xlim([x[0],x[nx-1]])
    axes.set_ylim([z[0],z[nz-1]])
    plt.colorbar()
    plt.title(file)

    # plt.tight_layout(pad=0.1, w_pad=0.4, h_pad=0.1)
    # plt.savefig("{}.pdf".format(file))
    plt.show()
