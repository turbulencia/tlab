#!/bin/sh

if [ $# -eq 0 ]; then
   echo "Usage: $0 file"
   exit 1
fi

TMPDIR=./convert

for file in $*; do

    name=`echo $file | awk -F. '{print $1}'`
    if [ ! -d $TMPDIR ]; then 
       mkdir $TMPDIR
    fi

# correct comment symbol
    sed 's/^c/!/g'     $file   > $file.1
    sed 's/^!    /!/g' $file.1 > $file.2

# correct the continuation of lines
    awk 'BEGIN { FS="" ; RS="\n     \\$" } { if (NR>1) {printf "      "}; print $0 "&" }' $file.2 > $file.3

# eliminate the last ampersand
    awk '{ if (NR<LINES) print $0 }' LINES=`cat $file.3 | wc -l` $file.3 > $name.f90

#move to TMPDIR
    mv $file.1 $file.2 $file.3 $file $TMPDIR

done 

