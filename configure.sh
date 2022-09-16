#!/bin/bash

WDIR=`pwd` 
CMAKE=cmake
BUILDS="BIG LITTLE PARALLEL DEBUG"
kernel=`uname`
host=`hostname` 

case $kernel in 
    Darwin ) 
	SYST=macbook ;; 
    AIX )
	SYST=blizzard;; 
    Linux )
	case $host in 
	    thunder* ) 
		SYST=thunder;;
            mpipc* )
                SYST=mpipc;; 
            juqueen* )
                SYST=juqueen;; 
            laptop*|home* )
                SYST=archlinux;; 
            eetac* )
                SYST=eetac;; 
            mi* )
                SYST=uhhpc;; 
            *.cen.uni-hamburg.de )
                SYST=uhhpc-c;; 
	esac 
esac

for B in $BUILDS; do
    dir=build_$B 
    echo -e "\nCreating $dir..." 
    if [ -d $dir ] ; then 
       rm -rf $dir
    fi
    mkdir $dir; cd $dir;  
    $CMAKE ../src -DSYST=$SYST -DBUILD_TYPE=$B #| grep 'flags' | cut -d':' -f2
    cd .. 
done

exit 0 
