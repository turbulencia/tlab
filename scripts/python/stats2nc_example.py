#!/opt/local/bin/python2.7
# script to process DNS statistics
# Chiel van Heerwaarden, 2011
# Cedrick Ansorge, 2015
#
import stats2nc

pathes=['.']
types= ['avg','avg1s','avg2s','int']
jmax = 512

for p in pathes:
    for t in types:
        d = stats2nc.avg2dict(t,p, jmax,0)
#        d = stats2nc.avg2dict(t,p, jmax,1) if gzipped files
        if ( d != -1 ):
            tstart=d['Iteration'][0]
            tend  =d['Iteration'][len(d['Iteration'][:])-1]
            stats2nc.dict2nc(d, '{}/{}{}-{}'.format(p,t,tstart,tend))
        else:
            'ERROR - in directory:', p
            'ERROR   creating the dictionary with stats2nc seems to have failed'
            'ERROR   for files of type', t
