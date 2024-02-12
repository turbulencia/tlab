#!/bin/sh
#
# habu, 4 tasks_per_node and network.MPI set to css0
# marcellus, 8 tasks_per_node and network.MPI set to csss
#
############################################################################
#@ environment = COPY_ALL; MP_EUILIB=us; XLSMPOPTS=parthds=1:spins=0:yields=0:schedule=static;MP_CPU_USE=unique;AIXTHREAD_SCOPE=S;AIXTHREAD_MNRATIO=8:8;MP_SHARED_MEMORY=yes;MP_PULSE=0;RT_GRQ=ON;MP_CSS_INTERRUPT=yes;MP_INTRDELAY=100 
#@ output = LL_out.$(jobid)
#@ error = LL_err.$(jobid)
#@ initialdir = /scr/mellado
#@ job_type = parallel
#@ notify_user =  mellado@mae.ucsd.edu
#@ network.MPI = csss,not_shared,US
#@ node_usage = not_shared
#@ notification = error
#@ notification = complete
#@ class = batch
########################### CHANGE ONLY THIS ###############################
#@ account_no=AFOSRSA2
#@ tasks_per_node = 8
########################### END NQS DIRECTIVES #############################
#@ queue
########################### DO NOT MODIFY ##################################

. $HOME/bin/CONFIG
export MP_INTRDELAY=100

# Change Account 
#newacct $ACCOUNT

# Change this to reflect your home directory
MSFHOME=/u/home/$LOGNAME
MSTOR=jules-hip0
ABORT="no"

MASSDIR=$MSFHOME/$JOBDIR
LOCALDIR=/scr/$LOGNAME/$JOBDIR
TRASHDIR=/scr/$LOGNAME/$JOBDIR-transfered
   
# Create Working Directories
if [ ! -d /tmp/$LOGNAME  ]; then
   mkdir /tmp/$LOGNAME 
fi
if [ ! -d $LOCALDIR ]; then
   mkdir  $LOCALDIR
fi

cd $LOCALDIR

# If not pt* files, then set DNSFIRST to yes
tmp=`ls pt*[0-9]`
if [ $? = 0 ]; then
    DNSFIRST="no"
else
    DNSFIRST="yes"
fi

# Check if temp_rest is already in Job Directory
if [ ! -e temp_rest ]; then
    # Try to create links to temp_rest and scal_rest
    dns.pre $LOCALDIR $STEP
fi

touch *
    
# Run Simulation
simulation $TIMESTAMP 

# Abort on severall errors
if [ -e tlab.error ]; then
    ABORT="yes"
fi
	
if [ -e tlab.error.0 ]; then
    ABORT="yes"
fi

if [ -e core ]; then
    touch dns.core
    rm -f core
    ABORT="yes"
fi

if [ $DNSFIRST = "no" ]; then
    /bin/rm -f *_rest
else
    mv temp_rest pt0
    mv scal_rest sc0
fi

# Compress and tar statistics
DATE=`date +%m%d%y-%H%M%S`
dns.stats dns-stat-$DATE.tar 

# Move to mass storage 
dns.archive $TRASHDIR $MASSDIR

# Requeue simulation
if [ $ABORT = "no" ]; then
	if [ -e $LOCALDIR/dns.nqs.new-vars ]; then
	    . $LOCALDIR/dns.nqs.new-vars
	fi

	if [ -e $LOCALDIR/tlab.ini ]; then
	    ITIME=`awk -F"=" '{ 
				if ( $1 == "End" ) 
				    {
				    print $2 
				    }
			    }' $LOCALDIR/tlab.ini` 
	else
	    echo "Error getting max time"
	    exit 1
	fi
 
	if [ $ITIME -lt $MAXITER ]; then
	    # Submit Script
	    qsend -name $NAME -time $TIME -maxiter $MAXITER -mem $MEM \
	          -script $SCRIPT -jobdir $JOBDIR -step $STEP -np $NP

	fi       
fi

