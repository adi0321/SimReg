#!/bin/bash
hostname=`hostname`
echo -e "\nRun: $0 on $hostname ... "

readonly DEBUG=0
readonly SCRIPT_DIR=$(dirname $0)

if [ $# -lt 1 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - GTF: Input gtf file with transcripts"
	echo ""
	exit 1
fi

gtf=$1

trLen="${gtf%.*}.trLen.txt"

#Compute transcripts lengths:
awk '($3=="exon") || ($2=="exon")  {ID=$12;L[ID]+=$5-$4+1} END{for(i in L){print i, L[i]}}' $gtf | sort | awk '{print $1, $2}' > $trLen