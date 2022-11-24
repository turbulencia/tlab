#!/usr/bin/python3

import numpy
import struct
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc

np = 0   # total number of particles in the files; if 0, then search dns.ini
varnames=['x-coordinate','y-coordinate','z-coordinate',
          'scalar']
#          'x-vorticity','y-vorticity','z-vorticity']
    
# do not edit
nv = len(varnames)

sizeofdata   = 4
sizeofheader = 0
sizeofmask   = 6

rc('text',       usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font',       family='serif', size=12)

matplotlib.rcParams.update({'font.size': 12})

if ( len(sys.argv) == 1 ):
    print("Usage: python $0 number-of-particles-to-plot list-of-files.")
    quit()

if ( np == 0 ):
    for line in open('dns.ini'):
        if "trajnumber" in line.lower():
            np = int(line.split("=",1)[1])

def itnumber(filename):
    main = filename.split(".")[1]
    return int(main)

def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

np2plot = int(sys.argv[1])
print("Processing {} out of {} trajectories.".format(np2plot,np))

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
fin = open(filetypes[0]+'.'+filetimes[0]+'.1')
fin.seek(0,2)           
nt = len(filetimes) *int( fin.tell() /( (np+1)*sizeofdata ) )
print('Processing {:d} times ...'.format(nt))
    
a = numpy.zeros((nv,nt,np2plot+1),dtype=float) # adding space for the time

it=0
for time in filetimes:
    for type in filetypes:
        file = type+'.'+time
        print('Processing {} ...'.format( file ) )
    
        fin = [ open(file+'.'+str(iv+1), 'rb') for iv in range(nv) ]
        [ fin[iv].seek(0,2) for iv in range(nv) ]
        eof = [ fin[iv].tell() for iv in range(nv) ]
        [ fin[iv].seek(0,0) for iv in range(nv) ]
        for iv in range(nv):
            if eof[iv] != eof[0]:
                sys.exit("File sizes mismatch")

        offset = 0
        while( offset < eof[0] ):
            [ fin[iv].seek(offset) for iv in range(nv) ]
            raw   = [ fin[iv].read((np2plot+1)*sizeofdata) for iv in range(nv) ] # Read time & np particles
            for iv in range(nv):
                a[iv,it] = numpy.array(struct.unpack('<{}f'.format(np2plot+1), raw[iv]))
            offset= offset +(np+1)*sizeofdata
            it = it +1
        [ fin[iv].close() for iv in range(nv) ]

a = numpy.transpose(a,(0,2,1)) # transposing the last two dimensions to speedup plotting
for iv in range(nv):
    plt.figure()
    for ip in range(np2plot):
        plt.plot(a[iv,0,:],a[iv,ip+1,:])
    plt.title(varnames[iv])
    plt.xlabel(r"Time")
    plt.ylabel(r"Variable")
    plt.tight_layout(pad=0.1, w_pad=0.1)
    plt.savefig("{}.pdf".format(varnames[iv]))

plt.show()
