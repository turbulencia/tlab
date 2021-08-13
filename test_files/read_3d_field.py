import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# import matplotlib.colors as mcolors

#---------------------------------------------------------------------------#
# path to 3d-fields
path  = str(os.path.dirname(__file__) + '/../test_parallel/' ) # name = 'flow.20.1' # file = str(path+name)
index = 10

#---------------------------------------------------------------------------#
# read grid and flow fields
grid = mp.DnsGrid(path+'grid')
flow = mp.Field(path,var='flow',index=index)
flow.read_3d_field()
# sys.exit()
#---------------------------------------------------------------------------#
# validate fortran routine
# ip_b =  grid.nx*grid.ny*(grid.nz//2-10) #+1
# for k in range(1,20):
#     # print("k = ", k)   
#     for j in range(1,20):
#         # print("     ip_b = ", ip_b)
#         field_3d_tlab[ip_b:ip_b+grid.nx-1] = 0
#         ip_b =  ip_b + grid.nx        
#     ip_b = grid.nx*grid.ny*(grid.nz//2-10+k)
# field_3d = field_3d_tlab.reshape((grid.nx,grid.ny,grid.nz),order='F') 

#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# 2d plot - yz
plt.figure(figsize=size)
plt.title('2d-plot -- yz-plane -- velocity u')
plt.xlabel("z")
plt.ylabel("y")
#
plt.pcolormesh(grid.z,grid.y,flow.u[grid.nx//2,:,:], shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xy
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
#
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,grid.nz//3].T, shading=shading, cmap='RdBu_r') #shading='gouraud',
plt.colorbar()
plt.show()
#---------------------------------------------------------------------------#
# 2d plot - xz
plt.figure(figsize=size)
plt.title('2d-plot -- xz-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("z")
#
plt.pcolormesh(grid.x,grid.z,flow.u[:,1,:], shading=shading, cmap='RdBu_r')
plt.colorbar()
plt.show()




#-----------------------------------------------------------------------------#
# plot lines
#-----------------------------------------------------------------------------#

# plt.figure(figsize=size)
# for i in range(0,10):
#     plt.plot(grid.z,field_3d[:,i,0])
# plt.show()


# plt.figure(figsize=size)
# for i in range(60,70):
#     plt.plot(grid.y,field_3d[0,:,i])
# plt.show()


#-----------------------------------------------------------------------------#
# plot with axis
#-----------------------------------------------------------------------------#

# plt.close("all")




# f = open(path +'epsi','rb')
# f.seek(0,0)
# rec = np.fromfile(f, np.dtype('<f8'), 128**2*96)
# rec = rec.reshape((128,96,128),order='F')
# f.close()