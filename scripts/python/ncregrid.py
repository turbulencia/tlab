'''
Interpolate and concatenate up to five netCDF4 files with different vertical grids. A cubic spline interpolation is performed. All netCDF4 files need to be located in the same 
directory and should follow the naming convention (i.e. first set: avg2000-4000.nc, avg1s2000-4000.nc,... and second set: avg4000-6000.nc, avg1s4000-6000.nc,...).
(Note: To perform a visual inspection of the results uncomment line 70-83 and adapt it to your needs.)

Bernhard Schulz, May 2019
'''

import warnings
import netCDF4 as nc
import numpy   as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

''' define iterations to be merged '''
it_1 = '3500'								# first index 
it_2 = '6500'								# second index
it_3 = '10100'								# third index
number_scalars = 4							# maximum 5

'''-----------------------------------------------------------------------------------------------------------------------'''
'''-----------------------------------------------------------------------------------------------------------------------'''
'''-----------------------------------------------------------------------------------------------------------------------'''
''' load data; _1 for SMALL GRID; _2 for LARGE GRID; _12 CONCATENATED FILES '''

path_avg_1 = './avg' +  it_1 + '-' + it_2 + '.nc';	avg_1 = nc.Dataset(path_avg_1, 'r')
path_avg_2 = './avg' +  it_2 + '-' + it_3 + '.nc';	avg_2 = nc.Dataset(path_avg_2, 'r')
print( 'merge ', path_avg_1, 'with', path_avg_2 )

''' create and read scalars '''
list_scalars = ['avg1s', 'avg2s', 'avg3s', 'avg4s', 'avg5s']
avg1s_1 = 'nan'; avg2s_1 = 'nan'; avg3s_1 = 'nan'; avg4s_1 = 'nan'; avg5s_1 = 'nan'	# create dummy varibles 
avg1s_2 = 'nan'; avg2s_2 = 'nan'; avg3s_2 = 'nan'; avg4s_2 = 'nan'; avg5s_2 = 'nan'
list_1 = [avg1s_1, avg2s_1, avg3s_1, avg4s_1, avg5s_1]
list_2 = [avg1s_2, avg2s_2, avg3s_2, avg4s_2, avg5s_2]
for idx in range(number_scalars):
	path_1 = './' + list_scalars[idx] +  it_1 + '-' + it_2 + '.nc'
	list_1[idx] = nc.Dataset(path_avg_1, 'r')

	path_2 = './' + list_scalars[idx] +  it_2 + '-' + it_3 + '.nc'
	list_2[idx] = nc.Dataset(path_avg_2, 'r')
	print( 'interpolate and concatenate ', path_1, 'with ', path_2 )


''' names for concatenated files '''
avg_12 = 'avg' +  it_1 + '-' + it_3 + '.nc' 

avg1s_12 = 'nan'; avg2s_12 = 'nan'; avg3s_12 = 'nan'; avg4s_12 = 'nan'; avg5s_12 ='nan'
list_12 = [avg1s_12, avg2s_12, avg3s_12, avg4s_12, avg5s_12]
for idx in range(number_scalars):
	list_12[idx] = list_scalars[idx] +  it_1 + '-' + it_3 + '.nc' 

#print('list of variables', avg_1.variables.keys())

'''-----------------------------------------------------------------------------------------------------------------------'''
'''-----------------------------------------------------------------------------------------------------------------------'''
'''-----------------------------------------------------------------------------------------------------------------------'''
''' define functions '''

def grid_12(y_1,A_1,y_2,A_2,varname):
	''' perform  n-th order (k=3) spline interpolation '''

	dim_t = len(A_1[:,0])
	dim_y = len(y_2)
	A_12 = np.zeros( (dim_t,dim_y) )
	for idx in range(dim_t):
		tck = interp.splrep(y_1, A_1[idx,:], k=3)
		A_12[idx,:] = interp.splev(y_2, tck)

#	''' quick and dirty check '''
#	diff = abs(abs(A_2[0,:]) - abs(A_12[-1,:]))
#	relative_diff = 100.0*np.nanmax(diff)/A_2[0,np.nanargmax(diff)]
#	
#	if relative_diff > 0.1 and np.nanmax(diff) > 0.1 and y_2[np.nanargmax(diff)] > -19.0:	# adapt this line to your needs
#            print('WARNING: profiles deviate by ',np.nanmax(diff),' at height ',y_2[np.nanargmax(diff)])
#            plt.figure('Profiles for ' + str(varname) + ' deviate by ' + str( round(np.nanmax(diff),2) ) + ' at height ' + str(y_2[np.nanargmax(diff)]) )
#            plt.plot(A_1[idx,:],y_1,'b',label='initial data')
#            plt.plot(A_12[idx,:],y_2,'r--',label='interpolated data')
#            plt.xlabel(str(varname)+' at iteration'+str(idx))
#            plt.ylabel('height')
#            plt.legend()
#            plt.show()

	return A_12 

def interpolate(data_1,data_2,species):
	''' interpolate data_1 into vertical grid of data_2 '''

	dictonary = {}
	for varname in data_1.variables:	
		print('interpolate and concatenate',str(species),str(varname))
		y_1 = data_1.variables['y'] 
		y_2 = data_2.variables['y'] 
		var_1 = data_1.variables[str(varname)] 
		var_2 = data_2.variables[str(varname)] 
	
		if len(var_1.shape) == 1 and varname == 'y':
			dictonary[str(varname)] = y_2
	
		if len(var_1.shape) == 1 and varname != 'y':
			if len(var_1) != len(data_1.variables['it']): 
				warnings.warn('Warning: length of varname_1 is unequal number of iterations')
			dictonary[str(varname)] = np.concatenate( (var_1,var_2), axis=0 ) 
	
		if len(var_1.shape) == 2:
			var_12 = str(varname)+'_12'
			var_12 = grid_12(y_1,var_1, y_2,var_2,varname)
			
			''' concatenate data '''
			dictonary[str(varname)] = np.concatenate( (var_12,var_2), axis=0 )
			#print('shape of',str(varname)+'_1:',var_1.shape, str(varname)+'_12:',var_12.shape, str(varname)+'_2:',var_2.shape, str(varname)+':',dictonary[str(varname)].shape)

	return dictonary

def writeNC(dictonary,name_12):
	''' process the dictionary to netcdf files '''

	avgnc  = nc.Dataset(name_12, 'w')
	
	''' dimension '''
	ntimes = len(dictonary['t'])
	jmax = len(dictonary['y'])
	
	print("Creating netCDF file with ntimes = {} and jmax = {}".format(ntimes, jmax))
	
	''' create dimensions in netCDF file '''
	dim_y = avgnc.createDimension('y', jmax)
	dim_t = avgnc.createDimension('t', ntimes)
	
	''' create variables '''
	var_t = avgnc.createVariable('t', 'f8',('t',)) 
	var_t.units='Days since 0000-01-01 00:00'
	var_y = avgnc.createVariable('y', 'f8',('y',)) 
	var_y.long_name='Height above Surface' 
	var_y.positive='up'
	var_y.standard_name='height' 
	var_y.units='level'
	var_it= avgnc.createVariable('it','i4',('t',))
	
	''' first, handle the dimensions '''
	var_t[:]  = dictonary['t'][:]
	var_y[:]  = dictonary['y'][:] 
	var_it[:] = [int(f) for f in dictonary['it'][:]]
	
	''' now make a loop through all vars '''
	for varname in list(dictonary):
	  if( not( (varname == "it") or (varname == "y") or (varname == "t") ) ):
	    vardata = dictonary[varname]
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
	
	avgnc.close()
	return

'''-----------------------------------------------------------------------------------------------------------------------'''
'''-----------------------------------------------------------------------------------------------------------------------'''
'''-----------------------------------------------------------------------------------------------------------------------'''
''' run code '''

''' interpolate and concatenate data '''
dict_avg = interpolate(avg_1,avg_2,'avg')
writeNC(dict_avg,avg_12)

for idx in range(number_scalars):
	dictonary = interpolate(list_1[idx], list_2[idx], list_scalars[idx] )
	writeNC(dictonary, list_12[idx])

