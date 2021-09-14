import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
#---------------------------------------------------------------------------#
# path to flow fields
path = str(os.path.dirname(__file__) + '/../test/' )

# indices and files
index_flow   = 0
index_fields = [1,2] #[1,2,3] # 3d case, here just 2d case
files        = []
for i in index_fields: 
    files.append(path + 'flow' + '.' + str(index_flow) + '.' + str(i))
            
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud' #'nearest'
figs    = 'figs' 
plt.close('all')

#---------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# data specification for flow data
sizeofdata   = 8  # in bytes
sizeofheader = 52
etype        = "<" # > big-endian; < little-endian
dtype        = "f" # float

planetype    = 'xy'
plane        = [1]#,100]
#---------------------------------------------------------------------------#
# rw specified planes from 3d-file with header
for file in files:
    print('--------------------------------------------------')
    print("Processing file %s ..." % file)
    for planenumber in plane:
        with open(file, "rb") as fin:
            header_int   = []
            header_float = []
            for i in range(5):
                header_int.append(np.fromfile(fin,etype+'i'+str(4),count=1))
            for i in range(4):
                header_float.append(np.fromfile(fin,etype+'f'+str(8),count=1))
            # new headers
            header_int_new = [52,128,128,1,0]
            header_float_new = [header_float[i] for i in range(len(header_float))]
            # read plane
            fin.seek(sizeofheader +(planenumber-1)*grid.nx*grid.ny*sizeofdata,0)
            raw_plane = np.fromfile(fin, dtype='<f8',count=grid.nx*grid.ny)
            # write out plane with header
            file_out = file + '.' + str(planenumber) + '.' + planetype
            with open(file_out, "wb") as fout:
                print("Write out  file %s ..." % file_out)       
                np.asarray(header_int_new  ).astype('<i4').tofile(fout)
                np.asarray(header_float_new).astype('<f8').tofile(fout)
                raw_plane.astype('<f8').tofile(fout)
#---------------------------------------------------------------------------#
# header definition:
    # 5 x type '<i4' 
        # headerlength      52
        # nx*ny*nz          128*128*128
        # iteration step    250000
    # 4 x type '<f8'
        # rtime             4351.52202238
        # visc (1/Re)       0.00019861
        # froude            1.
        # rossby            1.
#---------------------------------------------------------------------------#
# # rw specified planes from 3d-file
# for file in files:
#     print('--------------------------------------------------')
#     print("Processing file %s ..." % file)
#     for planenumber in plane:
#         fin      = open(file, 'rb')
#         file_out = file + '.' + str(planenumber) + '.' + planetype
#         print("Write out  file %s ..." % file_out)       
#         fout     = open(file_out,'wb')
#         if (planetype == 'xy'):
#             fin.seek(sizeofheader +(planenumber-1)*grid.nx*grid.ny*sizeofdata,0)
#             raw = fin.read(grid.nx*grid.ny*sizeofdata)
#             fout.write(raw)
#         fin.close()
#         fout.close()
#---------------------------------------------------------------------------#               
# # rw specified planes from 3d-file
# print("Processing file %s ..." % file)
# for planenumber in plane:
#     fin  = open(file, 'rb')
#     fout = open(file+'.'+str(planenumber)+'.'+planetype,'wb')
#     if ( planetype == 'xy' ):
#         fin.seek(sizeofheader +(planenumber-1)*grid.nx*grid.ny*sizeofdata,0)
#         raw = fin.read(grid.nx*grid.ny*sizeofdata)
#         fout.write(raw)
#     fin.close()
#     fout.close()
#---------------------------------------------------------------------------#