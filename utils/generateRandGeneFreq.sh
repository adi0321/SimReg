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

genes_names=`${code}/utils/get_genes.sh $gtf`
temp_arr=( $genes_names )
number_of_genes=${#temp_arr[@]}

#number_of_genes=`echo $genes_names | wc -l`
#echo "genes_names = $genes_names"
#echo "number_of_genes = $number_of_genes"

#genes_freq=/import1/UCSC/hg18/gnfHumanAtlas2_column1.cluster_concentrations

#First, get file without path:
filename=$(basename "$gtf")
extension="${filename##*.}"
filename="${filename%.*}"

out=$PWD/$filename.genes.freq
rm -rf $out*
touch $out


	
	##Randomly distribute frequecy to genes 
	#tr_freq=`~/code/scientific_notation_calculator $current_gene_freq $number_of_tr $option_freq`
	genes_freq=`~/code/scientific_notation_calculator 1 $number_of_genes $option_freq`
	
			if [ $DEBUG -ne 0 ]; then
				echo "number of genes: $number_of_genes"
			fi


	arr=($genes_freq)
	i=0
	for t in $genes_names
	do
		#echo -e "$t\t${arr[$i]}"
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



