#!/bin/bash

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

while read probes
do
	#echo "$probes"
	geneFreq=0
	
	for t in $probes
	do
		#echo "Extract: $t"
		
		genePrefix=${t:0:4}
		#echo "$genePrefix"
		
		if [ "$genePrefix" = "ENSG" ]; then
			g=$t
			continue
		fi
				
		#for each transcript
		#echo "Transcript $t"
		
		trFreq=`grep $t $inFile`
		#1st check if transcript exists in IsoEM results
		
		if [[ -z "${trFreq}" ]]; then
			# -z string - True if the length of string is zero.
			trFreq=0;
		else
		
			trFreq=`grep $t $inFile | awk '{print $2}'`
		fi
		
			#echo -e "$t\t$trFreq"
		
		geneFreq=$(echo "$geneFreq + $trFreq" | bc -l)

	done
	
	echo -e "$g\t$geneFreq" >> $outFile
	#exit 7
	
done < /data1/adrian/MAQC3_Estimations/hs.Ensembl.70_ENCODE.v15.Homo_sapiens.GRCh37.70.gtf.MAQC3_genes.txt




