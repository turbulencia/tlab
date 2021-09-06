# import numpy as np
# import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# from scipy import integrate
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
# path = str(os.path.dirname(__file__) + '/../test_little_channel/' ) # name = 'flow.20.1' # file = str(path+name)
path  = str(os.path.dirname(__file__) + '/../test_yamo_180/' )
index = 0
#---------------------------------------------------------------------------#
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
# shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')

# flow field
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()

# bulk velocity of the flow
ub = mp.ubulk(flow.u, grid.y)