#!/usr/bin/python3

import numpy as np   # For array operations.
import struct
import sys

setofplanes = [ 1, 16 ]

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

sizeofheader = 0
# sizeofheader = 36 # for gate files

etype = ">" # big-endian
# etype = "<" # little-endian

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
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,xz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

if   ( planetype == 'xy' ):
    sizeofmask = len(str(nz)) 
elif ( planetype == 'xz' ):
    sizeofmask = len(str(ny)) 
else:
    print("Usage: python $0 [xy,xz] list-of-files.")
    quit()

# further initialization
def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

# processing data
print("Processing file grid ...")
fin = open('grid', 'rb')
fout = open('grid.'+planetype,'wb')

fin.seek(56,0)
raw = fin.read(nx*8)
a = struct.unpack(etype+'{}d'.format(nx), raw)
raw = struct.pack(etype+'{}f'.format(nx),*a)
fout.write(raw)

if   ( planetype == 'xy' ):
    fin.seek(8,1)
    raw = fin.read(ny*8)
    a = struct.unpack(etype+'{}d'.format(ny), raw)
    raw = struct.pack(etype+'{}f'.format(ny),*a)
    fout.write(raw)
elif ( planetype == 'xz' ):
    fin.seek(8+ny*8+8,1)
    raw = fin.read(nz*8)
    a = struct.unpack(etype+'{}d'.format(nz), raw)
    raw = struct.pack(etype+'{}f'.format(nz),*a)
    fout.write(raw)

fout.close()
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    for plane in setofplanes:
        fout = open(file+'.'+planetype+tag(sizeofmask,plane),'wb')
        fin = open(file, 'rb')
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
