import struct
import array
import string
import numpy
import datetime

from netCDF4 import Dataset
from os import listdir
from os import system 

import os.path 
import os

user = os.environ['USER'] 
lock_fname = 'lock_tower_merge_for_{}'.format(user) 
if os.path.isfile(lock_fname):  
    print 'ERROR: cannot run tower_merge.py, locked by {}'.format(lock_fname)  
else :
    command = 'touch {}'.format(lock_fname) 
    system('touch {}'.format(lock_fname) )

quit()

def read_record(file, mode):
    reclen,=struct.unpack(mode[0]+'i',f.read(4))
    if ( mode[1] == 'i' or mode[1] == 'f'):
        size = reclen/4
    elif ( mode[1] == 'd' ):
        size = reclen/8
    else:
        print 'ERROR: unsupported data type for FORTRAN RECORD READ'
        exit
    read_string = mode[0]+mode[1]*size
    data=struct.unpack(read_string,f.read(reclen))
    recend,=struct.unpack(mode[0]+'i',f.read(4))
    if ( reclen != recend ) :
        print 'ERROR: FORTRAN RECORD READ - end tag does not match record size'
    return data

# SETUP DIRECTORY WHERE TO PUT PROCESSED DATA
#
timestamp=datetime.datetime.now().strftime("%Y%m%d-%I%M%S")
dump_dir='towerdump_'+timestamp
sys_command='mkdir ' + dump_dir
system(sys_command) 

# ##########################################################################################
# GET GRID AND TOWER STRIDE INFORMATION
# ##########################################################################################

f = open('dns.ini','r')
stride = -1
for line in f:
    if 'Stride' in line:
        stride = [int(i) for i in (line.split('=', 2)[1].split(','))]

if ( stride==-1) :
    print 'ERROR: Keyword stride not found in dns.ini'
    exit(1)
f.close()

f = open('grid','r') 
str=f.read(4)
i1,= struct.unpack('>i',str) # READ AS BIG ENDIAN
i2,= struct.unpack('<i',str) # READ AS LITTLE ENDIAN

if ( i1 == 12 ):
    endian='>'
    print 'BIG ENDIAN DATA'
    reclen = i1
elif ( i2 == 12):
    endian='<'
    print 'LITTLE ENDIAN DATA'
    reclen = i2
else:
    print 'ERROR: Cannot determine Endianness from first Integer in grid file'
    print '       Assuming BIG     gives:', i1
    print '       Assuming LITTLE  gives:', i2
    quit 

iread=endian+'i'
fread=endian+'f'
dread=endian+'d'

f.seek(0)
griddims = [i for i in read_record(f,iread)]
gridsize = [d for d in read_record(f,dread)]
xgrid    = [x for x in read_record(f,dread)]
ygrid    = [y for y in read_record(f,dread)]
zgrid    = [z for z in read_record(f,dread)]

print 'GRID:  ', griddims
print '       ', gridsize
print 'STRIDE:', stride

# ##########################################################################################
# CONSTRUCT TOWER INFORMATION
# ##########################################################################################

xtower = xgrid[0:griddims[0]:stride[0]]
ytower = ygrid[0:griddims[1]:stride[1]]
ztower = zgrid[0:griddims[2]:stride[2]]
itowerx = range(1,griddims[0]+1,stride[0])
itowery = range(1,griddims[1]+1,stride[1])
itowerz = range(1,griddims[2]+1,stride[2])
ntowerx = len(itowerx)
ntowery = len(itowery)
ntowerz = len(itowerz)

# ##########################################################################################
# INITIALIZE PROCESSING OF TOWERS
# ##########################################################################################

# Find out which time slices to process
#
slices = []
varcount = -1
from glob import glob
for f in glob('tower.mean.*.?'): 
    dummy = f.split('.',4)   
    print dummy 
    [start,end] =  dummy[2].split('-')   
    ivar = int(dummy[3])
    if ( ivar > varcount ): 
        varcount = ivar 
    if  [start,end] not in slices:
        slices.append([start,end])

end = sorted([int(i[1]) for i in slices], key=int)
start=sorted([int(i[0]) for i in slices], key=int)
slices=[ [start[i],end[i] ] for i in range(len(slices))]
start_iteration = min(start)
end_iteration   = max(end) 
ntimes = max(end) - min(start) + 1


print 'NTIMES:', ntimes, 'TOWER SIZE:', ntowery, 'SLICES:', len(end)

# ##########################################################################################
# BUILD NETCDF FILE 
# ##########################################################################################

for slc in slices:
    ncfile = Dataset('tower_new.nc','w',format='NETCDF4') 
    nc_xdim = ncfile.createDimension('x',ntowerx)
    nc_ydim = ncfile.createDimension('y',ntowery)
    nc_zdim = ncfile.createDimension('z',ntowerz)
    nc_tdim = ncfile.createDimension('t',slc[1]-slc[0]+1) 

    nc_xvar = ncfile.createVariable('x','f4',('x',) )
    nc_xvar.longname='streamwise distance normalized by Rossby Radius'
    nc_xvar.units='m'
    nc_xvar.axis='x'

    nc_yvar = ncfile.createVariable('y','f4',('y',) )
    nc_yvar.long_name='height normalized by Rossby Radius'
    nc_yvar.positive='up'
    nc_yvar.units='m'

    nc_zvar = ncfile.createVariable('z','f4',('z',) )
    nc_zvar.long_name = 'spanwise distance normalized by Rossby Radius'
    nc_zvar.units='m'
    nc_zvar.axis='z'

    nc_tvar = ncfile.createVariable('t','f4',('t',) )
    nc_tvar.long_name='time'
    nc_tvar.units='days since 0001-01-01 00:00'
    nc_tvar.calendar='none'
    nc_tvar.axis='T' 

    nc_itvar= ncfile.createVariable('it','i4',('t',) )
    nc_zvar.long_name='Iteration'
    nc_zvar.units='1'

    nc_uvar = ncfile.createVariable('u','f4',('t','z','x','y',),zlib=False,least_significant_digit=4)
    nc_uvar.description = 'Streamwise Velocity U, instantaneous'
    nc_uvar.units       = 'm s^-1; normalized by reference velocity'
    nc_vvar = ncfile.createVariable('v','f4',('t','z','x','y',),zlib=False,least_significant_digit=4)
    nc_vvar.description = 'VerticalVelocity W, instantaneous'
    nc_vvar.units       = 'm s^-1; normalized by reference velocity'
    nc_wvar = ncfile.createVariable('w','f4',('t','z','x','y',),zlib=False,least_significant_digit=4)
    nc_wvar.description = 'Spanwise Velocity V, instantaneous'
    nc_wvar.units       = 'm s^-1; normalized by reference velocity'
    nc_pvar = ncfile.createVariable('p','f4',('t','z','x','y',),zlib=False,least_significant_digit=4) 
    nc_pvar.description = 'Local Ageostrophic pressure, instantaneous'
    nc_pvar.units       = 'kg s^-2 m^-1; normalized by reference pressure'
    nc_svar = ncfile.createVariable('s','f4',('t','z','x','y',),zlib=False,least_significant_digit=4)
    nc_svar.description = 'Scalar concentration S, instantaneous'
    nc_svar.units       = '1; normalized by reference value at the boundary'
    
    nc_umvar = ncfile.createVariable('uM','f4',('t','y',),zlib=False,least_significant_digit=4)
    nc_umvar.description = 'Streamwise Velocity U, averaged'
    nc_umvar.units       = 'm s^-1; normalized by reference velocity'
    nc_vmvar = ncfile.createVariable('vM','f4',('t','y',),zlib=False,least_significant_digit=4)
    nc_vmvar.description = 'VerticalVelocity W, averaged'
    nc_vmvar.units       = 'm s^-1; normalized by reference velocity'
    nc_wmvar = ncfile.createVariable('wM','f4',('t','y',),zlib=False,least_significant_digit=4)
    nc_wmvar.description = 'Spanwise Velocity V, averaged'
    nc_wmvar.units       = 'm s^-1; normalized by reference velocity'
    nc_pmvar = ncfile.createVariable('pM','f4',('t','y',),zlib=False,least_significant_digit=4) 
    nc_pmvar.description = 'Local Ageostrophic pressure, averaged'
    nc_pmvar.units       = 'kg s^-2 m^-1; normalized by reference pressure'
    nc_smvar = ncfile.createVariable('sM','f4',('t','y',),zlib=False,least_significant_digit=4)
    nc_smvar.description = 'Scalar concentration S, averaged'
    nc_smvar.units       = '1; normalized by reference value at the boundary'
    
    
    # ##########################################################################################
    # PUT GRID DATA 
    # ##########################################################################################
    nc_xvar[:] = xtower
    nc_yvar[:] = ytower
    nc_zvar[:] = ztower 
    
    # ##########################################################################################
    # PROCESS TOWERS 
    # ##########################################################################################

    slc_string = '.'+string.zfill(slc[0],6) + '-' + string.zfill(slc[1],6)
    files = []
    vname = ['','u','v','w','p','s']
    print datetime.datetime.now().strftime("%Y%m%d %I:%M:%S%p") ,\
        ': processing time slice', slc[0],'-',slc[1],'/',end_iteration
    for ivar in range(1,int(varcount)+1):
        name = 'tower.mean' + slc_string + '.' + string.zfill(ivar,1)
        try:
            with open(name):
                files.append(name)
        except IOError:
            print 'ERROR problem opening file ', name
        for i in itowerx:
            for k in itowerz:
                name = 'tower.'+string.zfill(i,6)+'x'+string.zfill(k,6) \
                       + slc_string + '.' + string.zfill(ivar,1)
                try:
                    with open(name):
                        files.append(name)
                except IOError:
                    print 'ERROR problem opening file ', name

    tower_jmax = len(itowery) 
    ntimes_loc = slc[1] - slc[0] + 1
    iteration_rt=array.array('d',(-1. for i in range(0,ntimes_loc)))
    iteration_it=array.array('i',(-1  for i in range(0,ntimes_loc)))

    for f in files:
        datf = open(f,'r')
        ivar = int(f.split('.',4)[3])
        position = f.split('.',2)[1]
        if ( position == 'mean'):
            [xpos,zpos] = [-1,-1]
            idx_i = -1
            idx_k = -1
        else: 
            [xpos,zpos] = position.split('x',2)
            idx_i = int(xpos)/stride[0] 
            idx_k = int(zpos)/stride[2] 

        times_loc = slc[1] - slc[0] + 1
        data_raw = datf.read((tower_jmax+2)*8*times_loc)
        read_string = endian+ 'd'*(tower_jmax+2)*times_loc
        data_flt = struct.unpack(read_string,data_raw)
        data = numpy.reshape(data_flt,[times_loc,tower_jmax+2])
        
        # MOVE FILE TO dump directory for this batch 
        datf.close()
        sys_command='mv ' + f + ' ' + dump_dir 
        system(sys_command) 

        time_it = data[:,1]
        time_rt = data[:,0]
        t_index = [int(it-slc[0]+1) for it in time_it] 
        its=min(t_index); ite=max(t_index)+1;
        data = data[:,2:tower_jmax+2]

        for i in range(its,ite):
            if ( iteration_it[i] > 0 and int(time_it[i-its]) != iteration_it[i] ):
                print 'ERROR: Timestamp in file', f, '(', int(time_it[i-its]),')'
                print '       Does not agree with expected time:', iteration_it[i]
                exit(1)
            iteration_rt[i] = time_rt[i-its]
            iteration_it[i] = int(time_it[i-its])
            
        if ( idx_i >= 0 and idx_k >= 0 ):  
            if ( ivar == 1 ) : 
                nc_uvar[its:ite,idx_k,idx_i,:] = data
            if ( ivar == 2 ) :
                nc_vvar[its:ite,idx_k,idx_i,:] = data
            if ( ivar == 3 ) :
                nc_wvar[its:ite,idx_k,idx_i,:] = data
            if ( ivar == 4 ) :
                nc_pvar[its:ite,idx_k,idx_i,:] = data
            if ( ivar == 5 ) :
                nc_svar[its:ite,idx_k,idx_i,:] = data

        else:
            if ( ivar == 1 ) :
                nc_umvar[its:ite,0:tower_jmax] = data
            if ( ivar == 2 ) :
                nc_vmvar[its:ite,0:tower_jmax] = data
            if ( ivar == 3 ) :
                nc_wmvar[its:ite,0:tower_jmax] = data
            if ( ivar == 4 ) :
                nc_pmvar[its:ite,0:tower_jmax] = data
            if ( ivar == 5 ) :
                nc_smvar[its:ite,0:tower_jmax] = data
            
    # Write real time and iteration data
    nc_tvar[:] = iteration_rt
    nc_itvar[:]= iteration_it
                
    ncfile.close()
    sys_command='mv tower_new.nc tower'+string.zfill(slc[0],6)+'-'+string.zfill(slc[1],6)+'.nc' 
    system(sys_command) 


system('rm {}'.format(lock_fname))

exit 
