#!/usr/bin/python3

import numpy as np
import struct
import sys
import matplotlib.pyplot as plt

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

etype = ">" # big-endian
# etype = "<" # little-endian

dtype = "f" # floating number 
# dtype = 'B' # unsigned character, for gate files

nx = 0 # number of points in Ox; if 0, then search dns.ini
ny = 0 # number of points in Oy; if 0, then search dns.ini
nz = 0 # number of points in Oz; if 0, then search dns.ini

# do not edit below this line

# getting grid size from dns.ini, if necessary
if ( nx == 0 ):
    for line in open('dns.ini'):
        if line.lower().replace(" ","").startswith("imax="):
            nx = int(line.split("=",1)[1])
            break

if ( ny == 0 ):
    for line in open('dns.ini'):
        if line.lower().replace(" ","").startswith("jmax="):
            ny = int(line.split("=",1)[1])
            break
        
if ( nz == 0 ):
    for line in open('dns.ini'):
        if line.lower().replace(" ","").startswith("kmax="):
            nz = int(line.split("=",1)[1])
            break
        
print("Grid size is {}x{}x{}.".format(nx,ny,nz))

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 [xy,xz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

# processing data
fin = open('grid.'+planetype, 'rb')
raw = fin.read(nx*4)
x1 = np.array(struct.unpack(etype+'{}f'.format(nx), raw))
nx1= nx
if   ( planetype == 'xy' ):
    raw = fin.read(ny*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(ny), raw))
    nx2= ny
elif ( planetype == 'xz' ):
    raw = fin.read(nz*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(nz), raw))
    nx2= nz
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    raw = fin.read()
    a = np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
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
