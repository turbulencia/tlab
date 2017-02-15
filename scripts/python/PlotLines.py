#!/usr/bin/python3

import numpy as np
import struct
import matplotlib.pyplot as plt
import sys

# do not edit
if ( len(sys.argv) == 1 ):
    print("Usage: python $0 list-of-files.")
    quit()

setoffiles = sorted(sys.argv[1:])

fin = open('grid.I', 'rb')
raw = fin.read()
x = np.array(struct.unpack('>{}d'.format(fin.tell()/8), raw))
fin.close()

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    raw = fin.read()
    a = np.array(struct.unpack('>{}d'.format(fin.tell()/8), raw))
    plt.plot(x,a)
    plt.title(file)
    plt.show()
    fin.close()
