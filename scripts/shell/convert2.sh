#!/bin/sh

if [ $# -eq 0 ]; then
   echo "Usage: $0 file"
   exit 1
fi

TMPDIR=./convert

if [ ! -d $TMPDIR ]; then 
   mkdir $TMPDIR
fi

for file in $*; do

# correct comment symbol
    sed 's/'isord'/'imode_fdm'/g' $file > $file.1

#move to TMPDIR
    mv $file $TMPDIR
    mv $file.1 $file
done 

