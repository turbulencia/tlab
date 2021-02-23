#!/usr/bin/python3

import numpy as np
import struct
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.interpolate as interpolate

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

# etype = ">" # big-endian
etype = "<" # little-endian

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
if   ( planetype == 'xy' ):
    raw = fin.read(nx*4)
    x1 = np.array(struct.unpack(etype+'{}f'.format(nx), raw))
    nx1= nx
    raw = fin.read(ny*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(ny), raw))
    nx2= ny
elif ( planetype == 'xz' ):
    raw = fin.read(nx*4)
    x1 = np.array(struct.unpack(etype+'{}f'.format(nx), raw))
    nx1= nx
    raw = fin.read(nz*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(nz), raw))
    nx2= nz
elif ( planetype == 'yz' ):
    raw = fin.read(ny*4)
    x1 = np.array(struct.unpack(etype+'{}f'.format(ny), raw))
    nx1= ny
    raw = fin.read(nz*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(nz), raw))
    nx2= nz
fin.close()

# translation
x1e = np.empty(nx1+1)
y1e = np.empty(nx1+1)
x1e[0:nx1] = x1[0:nx1]
x1e[nx1] = x1[nx1-1] + x1[1]-x1[0]
scale = x1e[nx1] -x1e[0]

# load the time data for the translation
n_loc, t_loc = np.load("./times.npy") # Iterations 80000, 80010, 80020,... from tower data
n_loc = n_loc.astype(int).tolist() # get the time

# define translation velocity
Utrans = 0.537284965

for file in setoffiles:
# get iteration number and correspoding time
    it = file.split(".")[1]
    t = t_loc[n_loc.index(int(it))]

    print("Processing file {:s} for time {:f}.".format(file,t))
    fin = open(file, 'rb')
    raw = fin.read()
    a = np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
    a = a.reshape((nx2,nx1))
    fin.close()

# do translation
    x1t = np.remainder(x1e -Utrans*t,scale)
    for j in range(nx2):
        y1e[0:nx1] = a[j,:]
        y1e[nx1] = a[j,0]
        f = interpolate.interp1d(x1e,y1e)
        a[j,:] = f(x1t)[0:nx1]

    fig = plt.figure(figsize=[x1[-1],x2[-1]])
#    plt.pcolormesh(x1,x2,a,vmin=-0.537284965,vmax=0.537284965)
    plt.pcolormesh(x1,x2,a,cmap='magma_r',norm=colors.PowerNorm(0.5,np.amin(a),6.0),shading='auto')
    # plt.imshow(a,origin='lower',cmap='magma_r',norm=colors.PowerNorm(0.5,np.amin(a),6.0))
    plt.axis('equal')
    plt.gca().set_axis_off()

    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.margins(0,0)
    plt.savefig("{}.jpg".format(file),dpi=1920./x1[-1])
    plt.close(fig)
