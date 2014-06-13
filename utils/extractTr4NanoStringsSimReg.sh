#!/bin/bash

if [ $# -lt 1 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - Transcript Frequency Estimates (two column file)"
	echo ""
	exit 1
fi

inFile=$1
outFile=${1}.NanoStrings.txt

if [ -f ${outFile} ];then
	rm -rf $outFile
fi

probes=`ls /import1/GENCODE_v19/nanoString/`


if [ -d nanoStrings ];then
	rm -rf nanoStrings
fi

mkdir nanoStrings

for p in $probes
do
	
	tr=`cat /import1/GENCODE_v19/nanoString/$p`
	#echo "For each probe p = $p"
	#echo "tr=$tr"

	for t in $tr
	do
		#echo "Extract transcript: $t"
		grep $t $inFile >> ./nanoStrings/$p
	done
done

#Sum the counts from each probe and report the results
for p in $probes
do

	#probename=$(basename "$p")
	#extension="${filename##*.}"
	probeName="${p%.*}"
	#echo "probeName=$probeName"
	
	awk -v outFile=$outFile -v probeName=$probeName '
	BEGIN{
		FREQ=0;
	}
	{
			FREQ+=$2;
	}
	END{
		printf("%s %s\n", probeName, FREQ)>>outFile
	}
	' ./nanoStrings/$p
	
done

echo ""
echo "Done!"

