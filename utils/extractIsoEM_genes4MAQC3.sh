#!/bin/bash
#10/11/2014
#This script is used for IsoEM to parse the gene file estimates

if [ $# -lt 1 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - Genes Frequency Estimates (two column file)"
	echo ""
	exit 1
fi

inFile=$1
outFile=${1}.MAQC3.txt

if [ -f ${outFile} ];then
	rm -rf $outFile
fi

probes=`awk '{print $1}' /data1/adrian/MAQC3_Estimations/hs.Ensembl.70_ENCODE.v15.Homo_sapiens.GRCh37.70.gtf.MAQC3_genes.txt`

for g in $probes
do
	
		#echo "Extract gene: $g"
		gene=`grep $g $inFile`
		
		if [[ -z "${gene}" ]]; then
			# -z string - True if the length of string is zero.
			geneFreq=0;
		else
		
			geneFreq=`grep $g $inFile | awk '{print $2}'`
		fi
		
		echo -e "$g\t$geneFreq" >> $outFile

done
