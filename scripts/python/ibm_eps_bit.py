#!/usr/bin/env python3
import numpy as np 
import os
import struct
"""
#############################################################################
# Version 0.01 - Written by Jonathan Kostelecky
#---------------------------------------------------------------------------#
Script to read/modify/create/write IBM geometry field (eps-field), 
which is stored in a bitwise manner as int1.

Restriction: mod(imax/8) = 0

Created  on & by: 2022/04/12 - J. Kostelecky (j.kostelecky@posteo.de)

Modified on & by: ...        - ...
#############################################################################
"""
#-----------------------------------------------------------------------------#
# functions to convert/modify eps-field (IBM geometry) [bitwise <--> int1]
#-----------------------------------------------------------------------------#
def int2bit_1(out,data): # option 1
    bsize = data.size
    for i in range(bsize):
        ip = i * 8
        binary = f'{data[i]:08b}'
        if binary[0] == '-':
            binary = f'{data[i]+256:08b}'
        j = 0
        for k in range(-1,-9,-1):
            out[j+ip] = int(binary[k])
            j += 1
    return out

def int2bit_2(out,data): # option 2 (bit faster then option 1)
    bsize = data.size
    for i in range(bsize):
        ip = i * 8
        by   = struct.pack('b',data[i])
        by2b = ''.join(format(ord(by), '08b') for byte in by)
        j = 0
        for k in range(-1,-9,-1):
            out[j+ip] = int(str(by2b)[k])
            j += 1
    return out

def bit2int_2(out,data):
    bsize = out.size
    for i in range(bsize):
        by = []
        ip = i * 8
        j  = 0
        for k in range(7,-1,-1):
            by.append(int(data[k+ip]))
            j += 1
        bit = "".join(str(i) for i in by)
        if bit[0] == '1':
            out[i] = int(bit,2)-256
        else:    
            out[i] = int(bit,2)
    return out

#-----------------------------------------------------------------------------#
# data specification of eps field
#-----------------------------------------------------------------------------#
# path to data
current_path = os.getcwd() + '/'
path         = current_path
fname        = 'eps0.1'

# data types (little endian)
type_i1 = np.dtype('<i1'); type_i4 = np.dtype('<i4'); type_f8 = np.dtype('<f8')
sizeofdata_int1 = 1; sizeofdata_int4 = 4; sizeofdata_float = 8

# header
head_params = 5 
head_size   = head_params * sizeofdata_int4

#-----------------------------------------------------------------------------#
# read
#-----------------------------------------------------------------------------#
# header
f = open(path + fname,'rb')
f.seek(0,0)
header = np.fromfile(f, type_i4, head_params)
f.close()
print('Header size           :', header[0])
print('Grid   size (nx*ny*nz):', header[1]*8,'x',header[2],'x',header[3])

# data size (attention: h[1] = grid.nx*8!)
bsize = np.prod(header[1:4])
rsize = bsize * 8

# read eps field as int1
f = open(path + fname,'rb')
f.seek(header[0],0)
data = np.fromfile(f, np.dtype('<i1'), bsize)
f.close()

#-----------------------------------------------------------------------------#
# convert to bitwise 
#-----------------------------------------------------------------------------#
eps = np.zeros(rsize)
eps = int2bit_1(eps,data) # eps = int2bit_2(eps,data) # faster
eps = eps.reshape((header[1]*8,header[2],header[3]),order='F') # (attention: h[1] = grid.nx*8!)

#-----------------------------------------------------------------------------#
# modify eps field here
#-----------------------------------------------------------------------------#
eps_new = eps
# ... ...

#-----------------------------------------------------------------------------#
# convert to int1
#-----------------------------------------------------------------------------#
eps_new_1d = np.zeros(rsize)
eps_new_1d = eps_new.reshape((rsize),order='F')

out = np.zeros(bsize)
out = bit2int_2(out,eps_new_1d)

#-----------------------------------------------------------------------------#
# write modified/new eps field
#-----------------------------------------------------------------------------#
file_out = fname[:-1] + str(2)
with open(file_out, "wb") as fout:
    print("Write out  file %s ..." % file_out)       
    header.astype('<i4').tofile(fout)
    out.astype('<i1').tofile(fout)
    fout.close()