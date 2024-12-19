'''
Interpolate and concatenate up to five netCDF4 files with different vertical grids. A cubic spline interpolation is performed. 

Bernhard Schulz, May 2019
Raphael Pistor, 2023
JP Mellado, 2024, simplifying to just do it for one file. 
'''

import warnings
import netCDF4 as nc
import numpy   as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import sys

if ( len(sys.argv) == 1 ):
    print("Usage: python $0 nc-file-1 nc-file-2 nc-file-out.")
    quit()

file1 = sys.argv[1]
file2 = sys.argv[2]
fileOut = sys.argv[3]

print( 'Interpolate and concatenate '+file1+' with '+file2+'...' )

avg_1 = nc.Dataset(file1, 'r')
avg_2 = nc.Dataset(file2, 'r')

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

dict_avg = interpolate(avg_1,avg_2)
writeNC(dict_avg,fileOut)

