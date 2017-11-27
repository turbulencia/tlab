#!/usr/bin/python3

import numpy as np   # For array operations.
import sys

np = 10  # number of particles

sizeofmask = 6

# do not edit
sizeofdata   = 4
sizeofheader = 0

def itnumber(filename):
    main = filename.split(".",1)[1]
    return int(main)

def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

if ( len(sys.argv) == 1 ):
    print("Usage: python $0 list-of-files.")
    quit()

filenames = sorted(sys.argv[1:],key=itnumber)
        
for file in filenames:
    print("Processing %s..." % file )
    fin = open(file, 'rb')
    fin.seek(0,2)      # go to the file end.
    eof = fin.tell()   # get the end of file location
    fin.seek(0,0)      # go back to file beginning
    name =     file.split(".",1)[0]
    time = int(file.split(".",1)[1]) - int( eof /( 3*np*sizeofdata) )
    while(fin.tell() != eof):
        time = time +1
        raw = fin.read(3*np*sizeofdata)
        fout = open(name+'.'+tag(sizeofmask,time)+'.vtk','wb')
        fout.write(raw)
        fout.close()
    fin.close()
    
