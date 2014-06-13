#!/bin/bash

if [ $# -lt 3 ]
then
        echo -e "Usage: \n"
        echo 'argv[1] - list with gene names'
        echo 'argv[2] - annotation GTF (e.g., knownGene....gtf)'
		echo 'argv[3] - output file'
	exit 1
fi

gene_names=`cat $1`
gtf=$2
out=$3

rm -rf $out
touch $out

for gene in $gene_names
do
	grep -w $gene $gtf >> $out
done

