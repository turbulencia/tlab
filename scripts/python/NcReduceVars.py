import netCDF4 as nc
import sys

t_stride = 2    # Reducing also the number of times
y_stride = 2    # Reducing also the number of levels

# getting data from stdin
if ( len(sys.argv) <= 1 ):
    print("Usage: python $0 [flow,scal] list-of-files.")
    quit()

filetype   = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

# processing data
if   ( filetype == 'flow' ):
    setofvars = 'rU rV rW rR Rxx Ryy Rzz Rxy Rxz Ryz Wx Wy Wz Wx2 Wy2 Wz2 Prd Eps Trp'.split(' ')
elif ( filetype == 'scal' ):
    setofvars = 'rS rS2 Rsu Rsv Rsw rS_y'.split(' ')
else:
    print("Error: file type not supported.")
    quit()

for file in setoffiles:
    print("Processing file %s ..." % file)
    file_org = nc.Dataset(file,      'r')
    file_dst = nc.Dataset(file.replace('.nc','_r.nc'), 'w')

    # read independent variables from origin nc-file
    t_org = file_org.variables['t'][::t_stride]
    y_org = file_org.variables['y'][::y_stride]

    # create dimensions for destiny nc-file
    file_dst.createDimension('y',len(y_org))
    file_dst.createDimension('t',len(t_org))

    # create and write independent variables in destiny nc-file
    t_dst = file_dst.createVariable('t', 'f8', ('t',))
    y_dst = file_dst.createVariable('y', 'f8', ('y',))
    t_dst[:] = t_org[:]
    y_dst[:] = y_org[:]

    # read and write iteration numbers
    it_org = file_org.variables['it'][::t_stride]
    it_dst = file_dst.createVariable('it', 'i4', ('t',))
    it_dst[:] = it_org[:]

    # read and write other dependent variables
    for var in setofvars:
        var_org = file_org.variables[var][::t_stride,::y_stride]
        var_dst = file_dst.createVariable(var, 'f8', ('t','y',))
        var_dst[:,:] = var_org[:,:]

    file_dst.close()
