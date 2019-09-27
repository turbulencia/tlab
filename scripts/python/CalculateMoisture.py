#!/usr/bin/python3

# Script to calculate the moisture field from the buoyancy and topdown fields
# according to Eq. 37 in Mellado et al., QJRMS 143:2403-2419, 2017.

import numpy as np
import struct
import sys

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

# etype = ">" # big-endian
etype = "<" # little-endian

dtype = "f" # floating number 
# dtype = 'B' # unsigned character, for gate files

# Physical parameters
Fb = 0.005 # Surface buoyancy flux
N2 = 3.0   # Buoyancy gradient in free atmosphere, positive when increasing with height
gammaq = -3.0 # Moisture lapse rate in free atmosphere, positive when decreasing with height

# do not edit below this line

# getting data from stdin
if ( len(sys.argv) <= 3 ):
    print("Usage: python $0 phi list-of-buyancy-files list-of-topdown-files.")
    quit()

phi = float(sys.argv[1])
numberofpairfiles = int(len(sys.argv[2:]) /2)
setoffiles1 = sorted(sys.argv[2:2+numberofpairfiles])
setoffiles2 = sorted(sys.argv[2+numberofpairfiles:])

# Calculating coefficients
L0 = (Fb/N2**1.5)**0.5

c1 = phi /( L0 *N2 )
c2 = 2.0 /( L0 *gammaq )

# Loop over the list of files
for i in range(numberofpairfiles):
    output = 'Moisture' +str(i)
    print("Using {} and {} to create {}".format(setoffiles1[i],setoffiles2[i],output))
    
    # Read buoyancy field
    fin = open(setoffiles1[i], 'rb')
    raw = fin.read()
    size = int(fin.tell()/sizeofdata)
    a = np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
    fin.close()
    
    # Read top-down field
    fin = open(setoffiles2[i], 'rb')
    raw = fin.read()
    if int(fin.tell()/sizeofdata) != size:
        print("Error: Scalar fields have different size")
        quit()
    b = np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
    fin.close()
    
    # Computing
    a = c1 *a + c2 *b
    
    # Writing
    fout = open(output, 'wb')
    raw = struct.pack((etype+'{}'+dtype).format(size),*a.tolist())
    fout.write(raw)
    fout.close()
