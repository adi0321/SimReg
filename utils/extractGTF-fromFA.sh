#!/bin/bash

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
out=$PWD/$(basename "$gtf")

if [ -s $out ];then
	rm -rf $out
fi

#Get only transcript names for which we have a sequence
tr_names=`awk '{if(substr($1,0,1)==">") print substr($1,2,17)}' $fa`

temp_arr=( $tr_names )
number_of_tr=${#temp_arr[@]}
echo "# of transcripts: $number_of_tr"

for tr in $tr_names
do

	grep $tr $gtf >> $out

done
