#!/bin/bash
#This script is used only to extract subFA from subGTF

if [ $# -lt 2 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - GTF: Input gtf file with transcripts"
	echo "argv[2] - FA: Transcripts sequences"
	echo ""
	exit 1
fi

gtf=$1
fa=$2

#First, get file without path:
filename=$(basename "$gtf")
extension="${filename##*.}"
filename="${filename%.*}"

out=$PWD/${filename}.fa

if [ -s $out ];then
	rm -rf $out
fi

tr_names=`awk -F";" '{ per[$2] += 1 } END { for (i in per) print i}' $gtf | sed -e 's/transcript_id//g' | sed -e 's/"//g' | sed -e 's/ //g'`

for tr in $tr_names
do
	# the $? variable that is set to the return status of the command. So you'd have:
	line=$(grep $tr $fa)
	if [ $? -eq 1 ]
	then
		echo "$tr not found" >> notFoundTranscripts.txt
	else
		grep -A 1 $tr $fa >> $out
	fi
done
