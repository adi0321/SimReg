#!/bin/bash

if [ $# -lt 1 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - Genes Frequency Estimates (two column file)"
	echo ""
	exit 1
fi

annotation=/import1/GENCODE_v19/gencode.v19.annotation.gtf
inFile=$1
outFile=${1}.MAQC.txt

if [ -f ${outFile} ];then
	rm -rf $outFile
fi

probes=`awk '{print $1}' /import3/MAQC/TaqMan/HBRR_ENS_expression.1`

for g in $probes
do
	
	#g="ENSG00000131095"
	echo "Extract gene: $g"
	geneFreq=0
	#Find transcript names that belong to this gene:
		
	transcripts=`grep $g $annotation | awk '($3=="exon") || ($2=="exon")  {ID=$12;L[ID]+=$5-$4+1} END{for(i in L){print i}}'`

	#for each transcript
	for t in $transcripts
	do

		t=`echo $t| cut -d';' -f 1`
		#-d is the delimiter and -f is the filed(token) number
		
		t=`echo ${t:1:-1}`
		#remove quotes (first and last character)
		
			#echo "Transcript $t"
		
		trFreq=`grep $t $inFile | awk '{print $2}'`
			#echo "Tr freq: $trFreq"
		
		#Solve scientific notations since bc cannot do arithmetic wiht e numbers
		trFreq=`awk '{ printf("%.20f", $1)}' <<< $trFreq`
			#echo "Tr freq: $trFreq"
		
			#echo -e "$t\t$trFreq"	
		
		geneFreq=$(echo "$geneFreq + $trFreq" | bc -l)
			
	done
	
	echo -e "$g\t$geneFreq" >> $outFile
	#exit 7
done




