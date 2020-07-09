import netCDF4 as nc
import numpy as np
import sys

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 [flow,scal,tower] list-of-files.")
    quit()

filetype   = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

# processing data
if   ( filetype == 'tower' ):
    # setoflevels = np.arange(0,100,10,dtype=np.int16)
    setoflevels = np.array([0, 147, 276, 356, 419, 479, 537, 595, 652],dtype=np.int16)
    setofvars   = 'u v w p s'.split(' ')
# elif ( filetype == 'scal' ):
#     setoflevels = np.arange(0,100,10,dtype=np.int16)
# elif ( filetype == 'flow' ):
#     setoflevels = np.arange(0,100,10,dtype=np.int16)
else:
    print("Error: file type not supported.")
    quit()

for file in setoffiles:
    print("Processing file %s ..." % file)
    file_org = nc.Dataset(file,      'r')
    file_dst = nc.Dataset(file.replace('.nc','_r.nc'), 'w')

    # read independent variables from origin nc-file
    x_org = file_org.variables['x'][:]
    y_org = file_org.variables['y'][setoflevels]
    z_org = file_org.variables['z'][:]
    t_org = file_org.variables['t'][:]

    print(y_org)

    # create dimensions for destiny nc-file
    file_dst.createDimension('x',len(x_org))
    file_dst.createDimension('y',len(y_org))
    file_dst.createDimension('z',len(z_org))
    file_dst.createDimension('t',len(t_org))

    # create and write independent variables in destiny nc-file
    x_dst = file_dst.createVariable('x', 'f4', ('x',))
    y_dst = file_dst.createVariable('y', 'f4', ('y',))
    z_dst = file_dst.createVariable('z', 'f4', ('z',))
    t_dst = file_dst.createVariable('t', 'f4', ('t',))
    x_dst[:] = x_org[:]
    y_dst[:] = y_org[:]
    z_dst[:] = z_org[:]
    t_dst[:] = t_org[:]

    # read and write iteration numbers
    it_org = file_org.variables['it'][:]
    it_dst = file_dst.createVariable('it', 'i4', ('t',))
    it_dst[:] = it_org[:]

    # read and write other dependent variables
    for var in setofvars:
        var_org = file_org.variables[var][:,:,:,setoflevels]
        var_dst = file_dst.createVariable(var, 'f4', ('t','z','x','y',))
        var_dst[:] = var_org[:]

    file_dst.close()
