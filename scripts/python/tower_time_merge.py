from glob import glob 
from netCDF4 import Dataset 
files=glob("tower*-*.nc")

f_list=[] 
for f in files:  
    fname=f.split('/')[-1]    
    fbase=fname.split('.')[0]
    it_start=int(fbase.split('-')[0][5:] )
    it_end = int(fbase.split('-')[1])
    f_list.append([fname,it_start,it_end])

f_list=sorted(f_list,key=lambda x:x[1])  

f_last=f_list[0] 
it_first=f_last[1]
end_last=f_last[2] 
ntime=f_last[2]-f_last[1] +1 
for f in f_list[1:]:  
    
    if int(f[1]) != f_last[2] +1:  
        print('') 
        print('WARNING: Gap in time detected: last end: {} ; this begin: {}'.format(f_last[2],f[1]))  
        print('')
    f_last=f 
    ntime += f_last[2]-f_last[1]+1

it_last=f_last[2] 

print('Processing {} time steps'.format(ntime))


out_name='allTowers.{}-{}.nc'.format(it_first,it_last)
print(out_name)
nc_out=Dataset(out_name,'w',clobber=True)

f_count=0 
it_loc=0 
for f in f_list:  
    nc_in=Dataset(f[0],'r',format='NETCDF4') 
    # New file -- initialization 
    if f_count == 0: 
        for d in ['x','y','z']: 
            nc_out.createDimension(d,size=len(nc_in.dimensions[d]))  
        nc_out.createDimension('t',size=ntime) 

        for name,variable in nc_in.variables.items(): 
            nc_out.createVariable(name,variable.datatype,variable.dimensions) 
            nc_out[name].setncatts(nc_in[name].__dict__) 

            if name in [ 'x','y','z' ]: 
                nc_out[name][:]=nc_in[name][:]
    
    nt_loc=len(nc_in.dimensions['t']) 
    it_end=it_loc+nt_loc 
    print('Writing data from file {} to local step {}-{}'.format(f,it_loc,it_end))


    for name,variable in nc_in.variables.items():  
        if name in ['u','v','w','p','s']:
            nc_out[name][it_loc:it_end,:,:,:]=nc_in[name][:,:,:,:]
        elif name in ['uM','vM','wM','pM','sM']: 
            nc_out[name][it_loc:it_end,:]=nc_in[name][:,:] 
        elif name in ['t','it']: 
            nc_out[name][it_loc:it_end]=nc_in[name][:]
    it_loc+=nt_loc 
    f_count+=1
