#!/usr/bin/python3

import numpy as np   # For array operations.
import struct
import sys

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
    print("Usage: python $0 list-of-files.")
    quit()

setoffiles = sorted(sys.argv[1:])

sizeofdata = 8 # doble precission

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    fout= open(file+'.little','wb')

    entries = 5                 # reading meta-data in integer form
    raw = fin.read(entries*4)
    a = struct.unpack('>{}i'.format(entries), raw)
    print(a)
    raw = struct.pack('<{}i'.format(entries),*a)
    fout.write(raw)

    entries = (a[0]-entries*4) //sizeofdata  # reading parameters
    print("Converting %i parameters..." % entries)
    raw = fin.read(entries*sizeofdata)
    a = struct.unpack('>{}d'.format(entries), raw)
    print(a)
    raw = struct.pack('<{}d'.format(entries),*a)
    fout.write(raw)

    print("Converting fields...")
    for k in range(nz):
        raw = fin.read(nx*ny*sizeofdata)
        a = struct.unpack('>{}d'.format(nx*ny), raw)
        raw = struct.pack('<{}d'.format(nx*ny),*a)
        fout.write(raw)

    fout.close()
    fin.close()
