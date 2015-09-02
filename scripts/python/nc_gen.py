import gzip
import struct
import numpy as np
from netCDF4 import Dataset

#SHOULD BE UNIVERSAL FOR ALL OUR ENSIGHT FILES
ENSIGHT_GRID_OFFSET=644 # put 56 for normal grid file
ENSIGHT_DATA_OFFSET=244 # put 52 for normal 3d binary

#GENERATE THE DATASET (new)
ncf=Dataset('2048_004.nc','w',format='NETCDF4')
#PUT GLOBAL ATTRIBUTES 
ncf.Conventions="CF-1.6" #netcdf-cf convention for metadata in netcdf files which avizo claimes they are supporting
ncf.reynolds_number=125000.
ncf.schmidt_number=1.0
ncf.froude_numer=0.4
ncf.source='dns-5.6.1'
ncf.time_integration='RK-explicit4'
ncf.spatial_scheme=''

#HARDCODE DOMAIN SIZE 
 
ncf.createDimension('lon',  2048) 
ncf.createDimension('level', 150)
ncf.createDimension('lat',  2048)
ncf.createDimension('time',    7)

################################################################################
# CREATE VARIABLES
################################################################################  

#DIMENSION VARIABLES
####################
vlon=ncf.createVariable('lon',   'f', 'lon') 
vlon.long_name='streamwise distance normalized by Rossby Radius'
vlon.standard_name='x'
vlon.units='meter'

vlev=ncf.createVariable('level', 'f', 'level') 
vlev.long_name='vertical distance to the wall normalized by Rossby Radius'
vlev.standard_name='y'
vlev.units='meter'
vlev.positive='up'

vlat=ncf.createVariable('lat',   'f', 'lat') 
vlat.long_name='spanwise distance normalized by Rossby Radius'
vlat.standard_name='z'
vlat.units='meter'

vtim=ncf.createVariable('time',  'f', 'time')
vtim.long_name='time normalized with the inverse coriolis parameter'
vtim.standard_name='t'
vtim.units="seconds since 2001-1-1 0:0:0"
vtim.calendar='none'

#DATA VARIABLES
####################
vens=ncf.createVariable('LnEnstrophy', 'f', ('time','lon','level','lat'))
vens.long_name='Natural Logarithm of Flow Enstrophy normalized by the square of the inverse coriolis parameter'
vens.standard_name='LnEnstrophy'
vens.units='1'

vbuo=ncf.createVariable('Buoyancy',    'f', ('time','lon','level','lat'))
vbuo.long_name='Buoyancy'
vbuo.standard_name='Buoyancy'
vbuo.units='1'

vgbu=ncf.createVariable('GradBuoyancy','f', ('time','lon','level','lat'))
vgbu.long_name='Natural logarithm of gradient of Buoyancy'
vbuo.standard_name='LnBuoyancyGradient'
vbuo.units='1'

vu=ncf.createVariable('xwind','f',('time','lon','level','lat'))
vu.long_name='streamwise wind velocity normalized by the Gesotrophic Wind'
vu.standard_name='u'
vu.units='m/s'

vw=ncf.createVariable('zwind','f',('time','lon','level','lat'))
vw.long_name='streamwise wind velocity normalized by the Gesotrophic Wind'
vw.standard_name='u'
vw.units='m/s'

vv=ncf.createVariable('ywind','f',('time','lon','level','lat'))
vv.long_name='streamwise wind velocity normalized by the Gesotrophic Wind'
vv.standard_name='u'
vv.units='m/s'

# READ THE GRID
gfile=open("grid.ensight","rb")
gfile.seek(ENSIGHT_GRID_OFFSET)

in_type=np.dtype('>i4')
in_type.newbyteorder('S')
[nx,ny,nz] = np.fromfile(gfile,in_type,3)
isize_field=nx*ny*nz

in_type=np.dtype('>f4')
gfile.seek(ENSIGHT_GRID_OFFSET+12)
vlon[:]=np.fromfile(gfile,in_type,nx) 
vlev[:]=np.fromfile(gfile,in_type,ny)
vlat[:]=np.fromfile(gfile,in_type,nz)
vtim[:]=[0.,1.,2.,3.,4.,5.,6.]

times=("112000","112050","112100","112150", "112200","112500","113000")
vars=("LnEnstrophy","LnScalar1Gradient","Scalar1","VelocityX","VelocityY","VelocityZ")
#READ THE DATA
itime=0
for v in vars:
    for t in times:
        name=v+t
        dfile=open(name,'rb')
        dfile.seek(ENSIGHT_DATA_OFFSET)
        print np.fromfile(dfile,in_type,10)
        dfile.seek(ENSIGHT_DATA_OFFSET)  

        #THIS HAS TO BE FIXED 
        # if ( int(t) < 112300 ):
        #     itime=(int(t)-112000)/50
        # else:    
        #     itime=(int(t)-112500)/500 + 5
        # print v,t,itime 

        # PUT A LOOP HERE OVER THE Z DIMENSION 
        # isize_plane = nx*ny 
        #for iz in range(nz): 
        #    vens[itime,:,:,iz]=np.fromfile(dfile,in_type,isize_plane) 

        if ( v == "LnEnstrophy" ):
            vens[itime,:,:,:]=np.fromfile(dfile,in_type,isize_field)
        if ( v == "Scalar1" ):
            vbuo[itime,:,:,:]=np.fromfile(dfile,in_type,isize_field)
        if ( v == "LnScalar1Gradient"):
            vgbu[itime,:,:,:]=np.fromfile(dfile,in_type,isize_field)
        if ( v == "VelocityX" ):
            vu[itime,:,:,:]=np.fromfile(dfile,in_type,isize_field) 
        if ( v == "VelocityY" ):
            vw[itime,:,:,:]=np.fromfile(dfile,in_type,isize_field) 
        if ( v == "VelocityZ" ):
            vv[itime,:,:,:]=-np.fromfile(dfile,in_type,isize_field) 
        dfile.close
        
ncf.close()

