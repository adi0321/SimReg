#!/bin/bash
hostname=`hostname`
echo -e "\nRun: $0 on $hostname ... "

if [ $# -lt 2 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - GTF: Input gtf file with transcripts"
	echo "argv[2] - 0-Random, 1-Uniform"
	echo ""
	exit 1
fi

readonly DEBUG=0
readonly code=~/SimRegBest/

gtf=$1
option_freq=$2

tr_names=`awk -F";" '{ per[$2] += 1 } END { for (i in per) print i}' $gtf | sed -e 's/transcript_id//g' | sed -e 's/"//g' | sed -e 's/ //g'`
temp_arr=( $tr_names )
number_of_tr=${#temp_arr[@]}

#echo "tr_names = $tr_names"
#echo "number_of_tr = $number_of_tr"

#First, get file without path:
filename=$(basename "$gtf")
extension="${filename##*.}"
filename="${filename%.*}"

out=$PWD/$filename.tr.freq

	if [ -s $out ];then
		rm -rf $out*
	fi

touch $out


	
	##Randomly distribute frequecy to transcripts 
	#tr_freq=`~/code/scientific_notation_calculator $current_gene_freq $number_of_tr $option_freq`
	tr_freq=`${code}/utils/scientific_notation_calculator 1 $number_of_tr $option_freq`
	
			if [ $DEBUG -ne 0 ]; then
				echo "number of transcripts: $number_of_tr"
			fi


	arr=($tr_freq)
	i=0
	for t in $tr_names
	do
			if [ $DEBUG -ne 0 ]; then
				echo -e "$t\t${arr[$i]}"
			fi
			
		echo -e "$t\t${arr[$i]}" >> $out
		i=$((i+1))
	done

echo "Done Random"


##Normalize: (they are already normalized)
#compute sum of all
#sum_of_all=`awk '{ sum += $2 } END { print sum }' $out`

#echo "Sum of all tr. freq. is $sum_of_all"
#echo "out=$out"
#awk -v var=$sum_of_all '{print $1"\t"$2/var}' $out
#awk -v var=$sum_of_all '{print $1"\t"$2/var}' $out > ${out}.norm



