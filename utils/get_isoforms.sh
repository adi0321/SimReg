#!/bin/bash

if [ $# -lt 1 ]
then
	echo -e "Usage: \n"
	echo "argv[1] - GTF file"
	
	exit 1
fi

#-F is field separator

#awk -F";" '{ per[$2] += 1 } END { for (i in per) print i}' $1

awk -F";" '{ per[$2] += 1 } END { for (i in per) print i}' $1 | sed -e 's/transcript_id//g' | sed -e 's/"//g' | sed -e 's/ //g'
