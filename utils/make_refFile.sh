#!/bin/bash

#============================================================================
# Project     : Life Technologies Collaborative Research Grant on "Software for Robust Transcript Discovery and Quantification"
# Author      : Adrian Caciula
# Version     : 
# Copyright   : This program is property of Georgia State University and Life Technologies
# Description : Create reference file for Grinder
# Created     : 
# Updated     : 04/29/2014 
#===========================================================================

readonly SCRIPT_DIR=$(dirname $0)

hostname=`hostname`
echo -e "\n\tRun: $0 on $hostname ... "

##For string comparison, use if [ "$s1" == "$s2" ] 
##For the a contains b, use if [[ $s1 == *"$s2"* ]]
if [ "$hostname" == "alan" ]
then
	mRna=/data/malta/hg18/myKnownGeneMrna.txt
elif [ "$hostname" == "rna1" -o "$hostname" == "cnv1" ]
then
	mRna=/import1/UCSC/hg18/myKnownGeneMrna.txt
	#mRna=/import1/UCSC/hg19/hg19KnownGeneMrna.txt
else
	echo -e "\ntError: in $0: Unknown host!"
	echo -e "\t\tPlease provide the correct path"
	exit 1
fi


if [ $# -lt 1 ]
then
	echo ""
	echo "Too few parameters:"
	echo 'argv[1] - gtf file (full path - or use $PWD in front)'
	echo -e "\nOther Default Inputs:"
	echo "   mRna = $mRna"
	echo ""
	exit 1
fi

gtf=$1

echo -e "\n\t Get transcripts names ..."
tr_names=`${SCRIPT_DIR}/get_isoforms.sh $gtf`

output_dir=$(dirname ${gtf})
ref_file=${output_dir}/`basename "$gtf" .gtf`.fa
	
#remove previous reference file
rm -rf $ref_file
	
echo -e "\n\t Extract fa sequences from $mRna \n\t   to $ref_file ... "	

for transcr in $tr_names
do
	if [ "$(grep $transcr $mRna)" ]
	then
		#The return value is as in C: return 1 if true and 0 if false.	
	
		grep -w -A1 $transcr $mRna >> $ref_file
	else
		echo "Error: Sequence for transcript $transcr was not found in $mRna"
		exit 0
	fi
done

echo -e "\n\tDone!"
