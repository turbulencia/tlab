#!/usr/bin/python3

import numpy as np   # For array operations.
import sys

nx = 5120 # number of points in Ox
ny = 2048 # number of points in Oy
nz = 5120 # number of points in Oz

setoflines = [ [1740,   1],   # Give [j,k]-duples to process
               [1740,1280],
               [1740,2560],
               [1740,3840] ]

# do not edit
if ( len(sys.argv) == 1 ):
    print("Usage: python $0 list-of-files.")
    quit()

setoffiles = sorted(sys.argv[1:])

sizeofdata = 8
sizeofheader = 52
sizeofmaskJ = len(str(ny)) 
sizeofmaskK = len(str(nz))

def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

print("Processing file grid ...")
fin = open('grid', 'rb')
fin.seek(56,0)
raw = fin.read(nx*sizeofdata)
fout = open('grid.I','wb')
fout.write(raw)
fout.close()
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    for line in setoflines:
        fin.seek(sizeofheader+(line[1]-1)*nx*ny*sizeofdata+(line[0]-1)*nx*sizeofdata,0)
        raw = fin.read(nx*sizeofdata)
        fout = open(file+'.I.'+tag(sizeofmaskJ,line[0])+'.'+tag(sizeofmaskK,line[1]),'wb')
        fout.write(raw)
        fout.close()
    fin.close()
