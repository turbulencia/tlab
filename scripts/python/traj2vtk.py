#!/usr/bin/python3

import numpy
import struct
import sys

np = 0   # number of particles; if 0, then search dns.ini

# do not edit
sizeofdata   = 4
sizeofheader = 0
sizeofmask   = 6

if ( len(sys.argv) == 1 ):
    print("Usage: python $0 number-of-variables list-of-files.")
    quit()

if ( np == 0 ):
    for line in open('dns.ini'):
        if "trajectorynumber" in line.lower():
            np = int(line.split("=",1)[1])
print("{} trajectories.".format(np))

def itnumber(filename):
    main = filename.split(".")[1]
    return int(main)

def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

nv = int(sys.argv[1])
filenames = sorted(sys.argv[2:],key=itnumber)
        
filetypes = []
for name in filenames:
    type = name.split(".")[0]
    if not (any(type in s for s in filetypes)):
        filetypes.append(type)
        
filetimes = []
for name in filenames:
    time = name.split(".")[1]
    if not (any(time in s for s in filetimes)):
        filetimes.append(time)

a = numpy.zeros((nv,np),dtype=float)

for time in filetimes:
    for type in filetypes:
        file = type+'.'+time
        print("Processing %s ..." % ( file ) )
    
        fin = [ open(file+'.'+str(iv+1), 'rb') for iv in range(nv) ]
        [ fin[iv].seek(0,2) for iv in range(nv) ]
        eof = [ fin[iv].tell() for iv in range(nv) ]
        [ fin[iv].seek(0,0) for iv in range(nv) ]
        for iv in range(nv):
            if eof[iv] != eof[0]:
                sys.exit("File sizes mismatch")

        itime = int(time) - int( eof[0] /( (np+1)*sizeofdata) )
        while(fin[0].tell() != eof[0]):
            itime = itime +1
            rtime = [ fin[iv].read(sizeofdata) for iv in range(nv) ] # Read physical time
            # rtime1= [ struct.unpack('>{}f'.format(1), rtime[iv]) for iv in range(nv) ]
            # print(rtime1)
            raw = [ fin[iv].read(np*sizeofdata) for iv in range(nv) ] # Read np particles
            for iv in range(nv):
                a[iv] = numpy.array(struct.unpack('>{}f'.format(np), raw[iv]))
            rawout = struct.pack('>{}f'.format(np*nv),*numpy.transpose(a).reshape(np*nv).tolist())

            fout = open(type+'.'+tag(sizeofmask,itime)+'.vtk','wb')
            fout.write(rawout)
            fout.close()
        
        [ fin[iv].close() for iv in range(nv) ]
    
