#!/usr/bin/python3

import numpy as np   # For array operations.
import struct
import sys

nx = 128 # number of points in Ox
ny =  96 # number of points in Oy
nz = 128 # number of points in Oz

setofplanes = [ 1, 16 ]

# do not edit
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,xz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

sizeofdata = 4
sizeofheader = 0
if   ( planetype == 'xy' ):
    sizeofmask = len(str(nz)) 
elif ( planetype == 'xz' ):
    sizeofmask = len(str(ny)) 
else:
    print("Usage: python $0 [xy,xz] list-of-files.")
    quit()
        
def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

print("Processing file grid ...")
fin = open('grid', 'rb')
fout = open('grid.'+planetype,'wb')

fin.seek(56,0)
raw = fin.read(nx*8)
a = struct.unpack('>{}d'.format(nx), raw)
raw = struct.pack('>{}f'.format(nx),*a)
fout.write(raw)

if   ( planetype == 'xy' ):
    fin.seek(8,1)
    raw = fin.read(ny*8)
    a = struct.unpack('>{}d'.format(ny), raw)
    raw = struct.pack('>{}f'.format(ny),*a)
    fout.write(raw)
elif ( planetype == 'xz' ):
    fin.seek(8+ny*8+8,1)
    raw = fin.read(nz*8)
    a = struct.unpack('>{}d'.format(nz), raw)
    raw = struct.pack('>{}f'.format(nz),*a)
    fout.write(raw)

fout.close()
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    for plane in setofplanes:
        fout = open(file+'.'+planetype+tag(sizeofmask,plane),'wb')
        if   ( planetype == 'xy' ):
            fin.seek(sizeofheader +(plane-1)*nx*ny*sizeofdata,0)
            raw = fin.read(nx*ny*sizeofdata)
            fout.write(raw)
        elif ( planetype == 'xz' ):
            for k in range(nz):
                fin.seek(sizeofheader +k*nx*ny*sizeofdata +(plane-1)*nx*sizeofdata,0)
                raw = fin.read(nx*sizeofdata)
                fout.write(raw)
        fout.close()
    fin.close()
