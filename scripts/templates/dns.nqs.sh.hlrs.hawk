#!/bin/bash -vx
#
#PBS -m abe
#PBS -M cansorge@uni-koeln.de
########################### CHANGE ONLY THIS ###############################
########################### END LSF DIRECTIVES #############################
############################################################################
########################### DO NOT MODIFY ##################################
export TOOLS_HOME=$HOME/bin

export PATH=$PATH:$TOOLS_HOME
export CONFIG_FILE=$TOOLS_HOME/CONFIG.default

. $CONFIG_FILE

ABORT="no"

# Change this to reflect your directories
LOCALDIR=$JOBDIR
TRASHDIR=$JOBDIR-transfered

# Create Working Directories
mkdir -p $LOCALDIR

# Avoid catching lustre problems 
# (write logs [tlab.log, dns.out, dns.obs] to home, afterwards copy back)
export DNS_LOGGER_PATH=$HOME/logs/logs-$TIMESTAMP
mkdir -p $DNS_LOGGER_PATH

#####################################################################
# USE CUSTOM LFS STRIPING TO OVERRIDE THE UNFORTUNATE DEFAULT VALUES  
#    --stripe_count -1: use all available OSF's 
#    --stripe_size 32M: use striping across I/O devices every 32MB 
lfs setstripe --stripe-count -1 --stripe-size 32M $LOCALDIR
#
# FORCE USE OF THE UFS PROTOCOL FOR [PARALLEL] I/O 
export ROMIO_FSTYPE_FORCE="ufs:" 

cd $LOCALDIR

dns.pre $LOCALDIR $STEP

touch *

module load aocc/2.1.0 mpt/2.23 fftw/3.3.8 amd-libm

# Write list of used nodes
cat $PBS_NODEFILE > nodefile.txt

# Run Simulation
case $RUNMODE in
  "preprocess" )
    preprocess $TIMESTAMP
    ABORT="yes"
    ;;
  "simulation" )
    simulation $TIMESTAMP $NP $PINNING
    if [ $? -ne 0 ]; then
       ABORT="yes"
    fi
    ;;
  "postprocess" )
    postprocess $TIMESTAMP
    ABORT="yes"
  ;;
esac

# Abort on several errors
if [ -f tlab.err ]; then
    ABORT="yes"
fi

if [ -f tlab.err.0 ]; then
    ABORT="yes"
fi

stat -t core* >/dev/null 2>&1 && ABORT="yes"

# Clean
if [ -f tlab.ini ]; then
    cp tlab.ini tlab.ini-$TIMESTAMP
fi
LOGFILES="tlab.ini.bak tlab.log dns.out dns.obs dns.les partinfos.txt mapping.txt nodefile.txt"
for FILE in $LOGFILES; do
    if [ -f $FILE ]; then
        mv $FILE $FILE-$TIMESTAMP
    fi
done

# Copy log-files from home to localdir
cd $DNS_LOGGER_PATH
for FILE in $LOGFILES; do
    if [ -f $FILE ]; then
        mv $FILE $FILE-$TIMESTAMP
        cp $FILE-$TIMESTAMP $LOCALDIR
    fi
done
cd $LOCALDIR

# Organize statistics
STATSDIR=stats-$TIMESTAMP
if [ ! -e $STATSDIR ]; then
    mkdir $STATSDIR

    if [ -f tlab.ini ]; then
	cp tlab.ini $STATSDIR
    fi

    LIST=`ls | egrep 'avg[a-zA-Z]*[0-9]|[xzr](sp|cr)[0-9]|pdf[0-9]|int[0-9]|kin[0-9]'`
    echo "Moving statistic files into $STATSDIR"
    if [ -n "$LIST" ]; then
       mv $LIST $STATSDIR
    fi

else
    echo "$STATSDIR exists. Aborting"
fi

# Organize planes
PLANESDIR=planes-$TIMESTAMP
if [ ! -e $PLANESDIR ]; then
    mkdir $PLANESDIR
    mv planes?.* $PLANESDIR

else
    echo "$PLANESDIR exists. Aborting"

fi

if [ $RUNMODE = "simulation" ];then

# Requeue simulation
    if [ $ABORT = "no" ]; then
	if [ -f $LOCALDIR/dns.nqs.new-vars ]; then
	    . $LOCALDIR/dns.nqs.new-vars
	fi

	if [ -f $LOCALDIR/tlab.ini ]; then
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
	    qsend -name $NAME -queue $QUEUE -time $TIME -maxiter $MAXITER -mem $MEM \
	          -script $SCRIPT -jobdir $JOBDIR -step $STEP -np $NP \
                  -rankspernode $RANKSPERNODE -commmode $COMMMODE \
		  -pinning $PINNING

	fi
    fi
fi
