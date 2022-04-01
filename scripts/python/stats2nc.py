# read statistics from ASCII files and put them in netCDF
# Chiel van Heerwaarden,     2012 -- created 
# Cedrick.Ansorge@gmail.com, 2013 -- Some generalizations
#
import gzip as gz
import netCDF4
import subprocess
import array as arr
import os 
 
from pylab import *

def addGroup(gDict,gLine,ig):
    garr=gLine.split()
    if  size(garr) > 3:
        gName=garr[2]
        gDict[gName] = garr[3:]

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def avg2dict(avgtype,avgpath,jmax,ftype='',tstart=-1, tend=-1,tstep=-1,reynolds=1):

    print('TYPE: ', avgtype)
    print('PATH: ', avgpath)
    print('jmax: ', jmax)
    print('fTYPE:', ftype) 
    ftype_count = len(ftype) 

    ###########################################################
    # get times automatically from files in directory 
    ########################################################### 
    if ( tstart == -1 ) : 
        command = "find " + avgpath + ' -name \"' + avgtype + "[0-9]*\""+ftype
        # print(command)
        p=subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE) 
        flist = []
        ordered_list = []
        for file in p.stdout.readlines():
            dummy = file.strip()
            try:
                with open(dummy):
                    cstart=len('{}/{}'.format(avgpath,avgtype))
                    cend=len(dummy)-ftype_count             #strip off .gz
                    if ( is_number(dummy[cstart:cend]) ):
                        iter = int(dummy[cstart:cend])
                        ordered_list.append([iter,file.strip()])
            except IOError:
                print('ERROR - File', file, 'does not exist' )

        flist = [f[1] for f in sorted(ordered_list,key=lambda file:file[0]) ]

        retval = p.wait()
        ntimes = len(flist)
    else :
        ntimes = (tend - tstart) / tstep + 1


    print('FILES for', avgtype,':', ntimes )

    if ftype != '.nc' and ntimes > 0  : 
        print('  -- OLD (ASCII) VERSION --\n      first: {}\n      last:  {}'.format(flist[0],flist[-1]))
        d,vgd=avg2dict_asc(flist,ntimes,avgtype,avgpath,jmax,ftype,tstart, tend,tstep,reynolds)
        return d,vgd 
    elif ftype == 'nc' and ntimes > 0 :
        print(' -- NEW (NETCDF) VERSION --\n      first: {}\n      last:  {}'.format(flist[0],flist[-1]))
        mergeNC(avgtype,avgpath,'nc_all.nc') 
        return

def mergeNC(avgtype,avgpath,ofile,reynolds=1):

    command = 'ncrcat --4 -L 9 {}/{}*.nc {}.nc'.format(avgpath,avgtype,ofile)
    print(command) 
    os.system(command)

    # recalculate friction values for backward compatibility
    f_h=netCDF4.Dataset('{}.nc'.format(ofile),'r+')

    t=f_h['t']
    du=f_h['U_y1'][:,1]
    dw=f_h['W_y1'][:,1]

    vF=f_h.createVariable('FrictionVelocity','f8',('t',))
    vF.setncattr('Group','Friction')
    vA=f_h.createVariable('FrictionAngle','f8',('t',))
    vA.setncattr('Group','Friction') 

    nu=2./(reynolds*reynolds) 
    
    vF[:] = np.sqrt(nu*np.sqrt(du**2+dw**2))
    vA[:] = np.arctan2(dw,du)

    
    return

def avg2dict_asc(file_list,ntimes,avgtype,avgpath,jmax,ftype,tstart=-1, tend=-1,tstep=-1,reynolds=1):
    
  ftype_count = len(ftype) 
    
  ############################################################ 
  if ( ntimes == 0 ) : 
      return -1

  ig=0
  hdrStr_old=''
  update_vList=True
  updateGroups=True
  t_start=tstart
  
  for t in range(ntimes):
    if t_start == -1 :
        filename = file_list[t]
        #number starts after <path>/<file type>
        cstart=len('{}/{}'.format(avgpath,avgtype))
        cend=len(filename)-ftype_count       #strip off .gz
        filenum = filename[cstart:cend]
        tend = filenum
        if (t == 0):
            tstart = filenum
    else:
        filenum = tstart+t*tstep
        filename = '{}/{}{}{}'.format(avgpath,avgtype,filenum,gzip_str)

    # process the file
    if ftype == '.gz':
        f = gz.open(filename,'r')
    else:
        f= open(filename,'r')
   
    # retrieve the time
    # print(filename)
    datastring = f.readline().split()
    if datastring[0] == 'RTIME':
        time = datastring[2]
    elif datastring[0].decode('UTF-8') == 'RTIME' :
        time= float(datastring[2])
    else:
        print("ERROR getting time from file for type {}, it {}".format(avgtype,filenum))
        print(datastring)
        quit()

    # process the group items in the header
    datastring = f.readline()
    d0=datastring.split()[0]
    if type(d0) is bytes :
        d0=str(d0.decode('UTF-8') )

    if filenum == tstart:
        varGroups={}
        groupVars={}
        
    while d0 == 'GROUP' or d0 == 'IMAX' or d0 == 'JMAX':
        if d0 == 'GROUP'and updateGroups:
            addGroup(varGroups,datastring,ig)
            ig=ig+1
        datastring = f.readline()
        d0=datastring.split()[0]
        if type(d0) is bytes:
            d0=d0.decode('UTF-8')

    hdrStr=datastring
    hdr=datastring.split()
    data = f.readline().split()
    
    if hdrStr != hdrStr_old:   # if there is a first or new header, we need to update 
        hdrTotal=size(hdr)
        hdrProf=size(data)
        hdrBulk=hdrTotal-hdrProf
        hdrStr_old = hdrStr
        update_vList = True
        print('New HEADER for file {}\n  ...  {} VARIABLES total ( {} profiles +  {} bulk)'.format(filename,hdrTotal,hdrProf,hdrBulk))
    else:
        update_vList = False

    for n in range(hdrTotal * ( 1 if updateGroups else 0 )):
        v=hdr[n]
        for g in varGroups.keys():
            if v in varGroups[g]:
                groupVars[v]=g
                break
        updateGroups=False

    if(filenum == tstart):
    #initialize netCDF file
        avg = {}
        avg['Time'] = zeros(ntimes)
        avg['Iteration'] = zeros(ntimes)

        for n in range(hdrTotal):
            v=hdr[n]
            if n < hdrProf:
                avg[v] = zeros((ntimes, jmax))
            else:
                avg[v] = zeros(ntimes)

    elif update_vList:
    # variables are initialized, but the header has changed => need to check whether we have to define something new
        for n in range(hdrProf):
            v=hdr[n]
            if not v in avg:
                print('   ... new variable {} at it {} ({})'.format(v,filenum,t))
                if n < hdrProf:
                    avg[v] = zeros((ntimes,jmax))
                else:
                    avg[v] = zeros(ntimes)
                # if there are new variables, we will parse the group header next time to retrieve group informations
                updateGroups=True

    # process the data
    # first store the time
    avg['Time'][t] = time
    avg['Iteration'][t] = filenum

    for i in range(jmax):
        if i>0:
            data = f.readline().split() 

        if ( len(data) != hdrProf and len(data) != hdrTotal ) :
            print("ERROR: length of line ({}) neither matches expected number of profiles ({})".format(len(data),hdrProf))
            print("       nor expected number of total variables ({})".format(hdrTotal)) 
        # process the vertical profiles
        for n in range(hdrProf):
            avg[hdr[n]][t,i] = data[n]
  
        # process the time series
        if(len(data) == hdrTotal):
            for n in range(hdrProf, hdrTotal):
                avg[hdr[n]][t] = data[n]

    # calculate friction velocity and angle for backward compatibility
    if (  ( 'FrictionVelocity' in avg and not 'FrictionVelocity' in hdr ) or
          ( 'FrictionAngle'    in avg and not 'FrictionAngle'    in hdr ) ):
        dudy = avg['U_y1'][t,0]
        dwdy = avg['W_y1'][t,0]
        us = np.sqrt(np.sqrt(dudy**2+dwdy**2)/(reynolds*reynolds/2))
        al = np.arctan2(dwdy,dudy)*180/np.pi if us>1e-12 else 0. 

        print('Calculating u* and alpha:',t,us,al,np.average(avg['FrictionVelocity'][:t]),np.average(avg['FrictionAngle'][:t] ) )
        avg['FrictionVelocity'][t] = us
        avg['FrictionAngle'][t] = al 

    f.close()
    
  return avg,groupVars

#############################################################

def dict2nc(dict, ncname, flag=0,groups=-1):
  # process the dictionaries to netcdf files
  avgnc  = netCDF4.Dataset("{}.nc".format(ncname), "w")

  # retrieve dimensions
  time   = dict["Time"]
  ntimes = time.shape[0]

  y    = dict["Y"]
  jmax = y.shape[1]

  print("Creating netCDF file with ntimes = {} and jmax = {}".format(ntimes, jmax))
  
  # create dimensions in netCDF file
  dim_y = avgnc.createDimension('y', jmax)
  dim_t = avgnc.createDimension('t', ntimes)

  # create variables
  var_t = avgnc.createVariable('t', 'f8',('t',)) 
  var_t.units='seconds'
  var_y = avgnc.createVariable('y', 'f8',('y',)) 
  var_y.long_name='Height above Surface' 
  var_y.positive='up'
  var_y.standard_name='height' 
  var_y.units='level'

  var_it= avgnc.createVariable('it','i4',('t',))
  # store the data
  # first, handle the dimensions
  var_t[:]  = dict['Time'][:]
  var_y[:]  = dict['Y']   [0,:] 
  var_it[:] = [int(f) for f in dict['Iteration'][:]]
  
  # now make a loop through all vars in dict.
  for varname in dict.keys(): 
      if(not(  (varname == "Iteration") or (varname == "Y") or \
                   (varname == "Time") or (varname == "I") or (varname == "J") ) ):
          vardata = dict[varname]
          if(len(vardata.shape) == 2):
              if( (vardata.shape[0] == ntimes) and (vardata.shape[1] == jmax) ):
          #print("Storing {} in 2D (t,y) array".format(varname))
                  var_name = avgnc.createVariable(varname,'f8',('t','y',))
                  var_name[:,:] = vardata
          if(len(vardata.shape) == 1):
              if(vardata.shape[0] == ntimes):
          #print("Storing {} in 1D (t) array".format(varname))
                  var_name = avgnc.createVariable(varname,'f8',('t',))
                  var_name[:] = vardata

  if groups != -1:
      for v in groups.keys():
          if  v != 'dummy' :
              avgnc[v].setncattr('Group',groups[v])

  # close the file
  avgnc.close()

