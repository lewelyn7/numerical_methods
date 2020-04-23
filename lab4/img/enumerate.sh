#!/bin/bash
IFS=$'\n'
names=`find -name "*.png" | sort`
if [ $1 -e 0 ] 
then
        it=0
else
        it=$1
fi
for i in $names
do
        mv "$i" "$it.jpg"
        it=$(($it + 1))
done
