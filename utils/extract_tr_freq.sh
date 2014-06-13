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

genes_names=`${SCRIPT_DIR}/get_genes.sh $gtf`
#genes_names="219823_at"

genes_freq=/import1/UCSC/hg18/gnfHumanAtlas2_column1.cluster_concentrations
#genes_freq=


#First, get file without path:
filename=$(basename "$gtf")
extension="${filename##*.}"
filename="${filename%.*}"

out=$PWD/$filename.tr.freq
rm -rf $out*
touch $out
touch $out.norm


for gene in $genes_names
do

	current_gene_freq=`grep -w $gene $genes_freq | awk '{print $2}'` 
	
			if [ $DEBUG -ne 0 ]; then
				echo -e "\n\nfor gene $gene"
				echo "current_gene_freq = $current_gene_freq"
			fi
	#get gene and extract transcript names
	tr=`grep -w $gene $gtf | awk -F";" '{ per[$2] += 1 } END { for (i in per) print i}' | sed -e 's/transcript_id//g' | sed -e 's/"//g' | sed -e 's/ //g'`
	
	number_of_tr=0
	for t in $tr
	do
		number_of_tr=$((number_of_tr+1))

	done

		
	
	##Randomly distribute frequecy to transcripts 
	tr_freq=`${SCRIPT_DIR}/scientific_notation_calculator $current_gene_freq $number_of_tr 0`
			
			if [ $DEBUG -ne 0 ]; then
				echo "transcripts are $tr"
				echo "number of transcripts: $number_of_tr"
				grep -w $gene $genes_freq
				echo "tr_freq = $tr_freq"
				
				if [[ $gene == "219823_at" ]]; then
					echo "gene = $gene"
					exit 7
				fi
				
				if [[ $gene == "823_at" ]];then
					echo "gene = "$gene
					echo "number_of_tr = $number_of_tr"
					exit 7
				fi
			fi


	arr=($tr_freq)
	i=0
	for t in $tr
	do
		#echo -e "$t\t${arr[$i]}"
		echo -e "$t\t${arr[$i]}" >> $out
		i=$((i+1))
	done
	
done

##Normalize:
#compute sum of all
sum_of_all=`awk '{ sum += $2 } END { print sum }' $out`

echo "Sum of all tr. freq. is $sum_of_all"
#echo "out=$out"
#awk -v var=$sum_of_all '{print $1"\t"$2/var}' $out
awk -v var=$sum_of_all '{print $1"\t"$2/var}' $out > ${out}.norm



