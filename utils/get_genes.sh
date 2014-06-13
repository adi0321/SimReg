#!/bin/bash

if [ $# -lt 1 ]
then
	echo -e "Usage: \n"
	echo 'argv[1] - known Gene GTF (e.g., /data3/serghei/trip_test/knownGeneGnfAtlas2_chr21_only.gtf)'
	exit 1
fi

awk -F"\"" '{ per[$2] += 1 } END { for (i in per) print i}' $1
