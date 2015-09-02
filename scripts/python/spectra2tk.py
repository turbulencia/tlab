#!/usr/bin/python2

import numpy
import struct

nx = 3619
ny = 1280

time      = 12000  # iteration number to be process
filetag   = "rsp"  
variables = ["uu", "vv", "ww", "pp", "11"]

# no need to modify below
if filetag == "xsp" or filetag == "zsp" or filetag == "rsp":
    nametag = "E"

if filetag == "xcr" or filetag == "zcr" or filetag == "rcr":
    nametag = "C"

# printing ascii data to screen
print 'RTIME = %i' % ( time )
print 'IMAX = 1'
print 'JMAX = %i' % ( ny )
print 'IVAR = 1'

line = 'I J K'
for var in variables:
    line = line + ' ' + nametag + var

print line

i = 0
a = numpy.empty(len(variables)*ny*nx)
for var in variables:
    file = filetag + str(time) + '.' + nametag + var

#    print "Reading data from file", file, "..."
    fin = open(file,"rb")
    raw = fin.read(nx*ny*4)
    a[i*nx*ny:(i+1)*nx*ny] = numpy.array(struct.unpack('>{}f'.format(nx*ny), raw))
    fin.close()
    i += 1

a = a.reshape((len(variables),ny,nx))

for j in range (ny):
    for i in range (nx):
        print '1', j+1, i+1, ' '.join(['{:0.4e}'.format(i) for i in a[0:len(variables),j,i]])
