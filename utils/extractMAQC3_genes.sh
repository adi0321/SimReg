#!/bin/bash
#Extract genes and transcripts from MAQC3 in this format
#GeneName Transcript1 Transcript2
echo $0

if [ $# -lt 1 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - Annotation GTF File"
	echo ""
	exit 1
fi

#annotation=/import1/GENCODE_v19/gencode.v19.annotation.gtf
annotation=$1
outFile=${1}.MAQC3_genes.txt

if [ -f ${outFile} ];then
	rm -rf $outFile
	rm -rf genesNotFound.txt
fi

#probes=`awk '{print $1}' /import3/MAQC/TaqMan/UHRR_ENS_expression.1`

#probes=`head /import3/MAQC3/qPCR/ID.txt | awk '{print $2}'` 
probes=`awk '{print $2}' /import3/MAQC3/qPCR/ID.txt`

for g in $probes
do
	
	#g="ENSG00000131095"
	#echo "Extract gene: $g"

	if [ "$g" = "ORF" ]; then
		continue
	fi	
	
	#Find transcript names that belong to this gene:
	unset transcripts
	transcripts=`grep $g $annotation | awk '($3=="exon") || ($2=="exon")  {ID=$12;L[ID]+=$5-$4+1} END{for(i in L){print i}}'`

	if [[ -z "${transcripts}" ]];then
		echo "$g" >> genesNotFound.txt
		echo "$g" >> $outFile
		continue
	else
		echo -ne "$g" >> $outFile
	fi
	#for each transcript
	for t in $transcripts
	do

		t=`echo $t| cut -d';' -f 1`
		#-d is the delimiter and -f is the filed(token) number
		
		t=`echo ${t:1:-1}`
		#remove quotes (first and last character)
		
		#echo "Transcript $t"
		echo -ne "\t$t" >> $outFile
			
	done
	
	echo "" >> $outFile
	#exit 7
done




