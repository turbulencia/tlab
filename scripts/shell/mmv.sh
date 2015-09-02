#!/bin/sh

if [ $# -eq 0 ]; then
   echo "Usage: $0 string1 string2"
   exit 1
fi

LIST=`ls -d *$1*`
for ORG in $LIST; do
  DES=`echo $ORG | sed s/$1/$2/`
  mv $ORG $DES
done

