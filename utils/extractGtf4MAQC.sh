#!/bin/bash

annotation=/import1/GENCODE_v19/gencode.v19.annotation.gtf

outFile=gencode.v19.annotation.MAQC.gtf

if [ -f ${outFile} ];then
	rm -rf $outFile
fi

probes=`awk '{print $1}' /import3/MAQC/TaqMan/UHRR_ENS_expression.1`

for g in $probes
do
	
	echo "Extract gene: $g"
		
	grep $g $annotation >> $outFile

done




