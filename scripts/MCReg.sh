#!/bin/bash

hostname=`hostname`
echo -e "\nRun: $0 on $hostname ... "
if [ $# -lt 6 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - gtf file"
	echo "argv[2] - d values file"
	echo "argv[3] - obs values file"
	echo "argv[4] - total number of observed reads"
	echo "argv[5] - output file"
	echo "argv[6] - Read Classes Size File"
	echo ""
	exit 1
fi

readonly DEBUG=0
readonly SCRIPT_DIR=$(dirname $0)

gtf=$1
d_values=$2
obsVFile=$3
total_no_reads=$4 #total number of reads

#out_file="../0.mcreg.iso.estimates"
out_file=$5
rcSize_file=$6

if [ -s $out_file ];then
	rm -rf $out_file
fi

echo "Input parameters are:"
echo "1. gtf = $gtf"
echo "2. d_values = $d_values"
echo "3. obsVFile = $obsVFile"
echo "4. total_no_reads= $total_no_reads"
echo "5. out_file = $out_file"


rm -rf ./temp/

mkdir temp

#Compute transcripts lengths:
#This step was done in compute_obsRC_CC.cpp
#rm -rf ./tr_length.txt
#awk '($3=="exon") || ($2=="exon")  {ID=$12;L[ID]+=$5-$4+1} END{for(i in L){print i, L[i]}}' $gtf | sort | awk '{print $1, $2}' > tr_length.txt



if [ -s "$d_values" ]
then
awk '{if($1=="Component:"){ 

		prefix=$2
		sum_obs_freq=$11
		nr_tr=$14
		
		#for i=nr_tr+1
		for (i = 15; i <= NF; i++){
			printf("%s\n", $i)>>"./temp/"prefix"_tr_names.txt"
		}
		
		#printf("%s\t%s\n", $2, $5)>>"./temp/no_reads.txt"
		printf("%s\n", $5)>>"./temp/"prefix"_no_reads.txt"
		printf("%s\n", sum_obs_freq)>>"./temp/"prefix"_sum_obs_freq.txt"
		
		
	}
	else if($1=="---"){
		#print "---"
	}
	else{
			#printf("Number of transcripts is: %d \n",nr_tr);
		#i=1 is observed read class frequency
		#i=2 is read class name
		#i=3 is the first d value
		for (i = 3; i <= nr_tr+2; i++){
			printf("%.10f ",$i)>>"./temp/"prefix"_d_values.txt";
		}
		
		printf("%f\n",$1)>>"./temp/"prefix"_o_values.txt";
		printf("\n")>>"./temp/"prefix"_d_values.txt";
	}
}' $d_values
else
	echo -e "\nERROR: $d_values is empty.\n"
	exit 7
fi


cd ./temp/

components=`ls *_tr_names.txt`
total_sum_reads_portion=0; #Sum of all (L*f/n) from each component

##Here we should actually have only one component
##so we may ignore this for loop for now since it should iterate only once

for c in $components
do
	transcripts=`cat $c`
	SUBSTRING=`echo $c| cut -d'_' -f 1`
	#-d is the delimiter and -f is the filed(token) number
				if [ $DEBUG -ne 0 ]; then
					echo "Transcripts = $transcripts"
					echo "SUBSTRING = $SUBSTRING"
				fi
	for t in $transcripts
	do
		grep $t ../tr_length.txt | awk '{print $2}' >> ${SUBSTRING}_tr_lengths.txt 
	done
	
	#array with transcripts lengths
	tr_lengths_array=(`cat ${SUBSTRING}_tr_lengths.txt`)	
	
	nt=`wc -l ${c} | awk '{print $1}'`
	s2=`echo $c| cut -d'_' -f 1` #component id
	nrc=`wc -l ${s2}_d_values.txt | awk '{print $1}'`
	sum_obs_freq=`cat ${s2}_sum_obs_freq.txt`
	
			if [ $DEBUG -ne 0 ]; then
				
				echo "Tr lengths : ${tr_lengths_array[@]}"
				echo "nt=$nt"
				echo "s2 = $s2"
				echo "nrc (Number of reads classes) = $nrc" ##number of reads classes
				echo "sum_obs_freq = $sum_obs_freq"
				
					#exit 7
			fi
	
	echo -e "\nRun QP Python solver for component id $s2\n"
	#python ~/code/wt/regression/qp_v3.py $nt $nrc $PWD/${s2}_d_values.txt $PWD/${s2}_o_values.txt $PWD/${s2}_tr_lengths.txt > ${s2}.qp.log
	python ${SCRIPT_DIR}/qp_v3.py $nt $nrc $PWD/${s2}_d_values.txt $obsVFile $PWD/${s2}_tr_lengths.txt $rcSize_file
	wait

	#Here we need to compute all those previous files again and run solver again
	#Why we cannot run the solver in the previous step? (I don't understand this question, what step?)
	
	##Now we need to simulate again observed reads
	##(The simulation is done in run-mcreg.sh)
	
	
	#Afer we run the solver we need to compute n divided by the average length for this component
	#n / sum(L_i*f_i)
	
	tr_estim=`cat ${s2}_d_values.iso.estimates`
	
	#if [[ "$s2" == "1018" ]]; then
	#	exit 7
	#fi
	
	##Extract f from solver estimates and multiply with the corresponding transcript length
	i=0 #used for array index
	local_avg_tr_len=0 #Average transcripts lengths
	
	for t in $tr_estim
	do
		if [[ "$t" == "[" ]]
		then
			continue
		else
			tr_i=`echo $t | cut -d ']' -f 1`	
			
			first_ch=${tr_i:0:1}
			#${string:position:length}
				
				if [ $DEBUG -ne 0 ]; then
					echo "Transcript i is $tr_i"
					echo "First Character = $first_ch"
				fi
			
			if [[ "$first_ch" == "[" ]]
			then
				tr_i=`echo $tr_i | cut -d '[' -f2`
			fi
			
			#Solve scientific notations since bc cannot do arithmetic with e numbers
			negNum=$(echo "$tr_i < 0" | bc -l)
			#negativeNumber should be 1 if true and 0 if false
			
			if [[ $negNum == "1" ]]; then
				tr_i=0
			else
				tr_i=`awk '{ printf("%.10f", $1)}' <<< $tr_i`
			fi
			
			
			
				if [ $DEBUG -ne 0 ]; then
					echo "New2 tr_i is $tr_i"
				fi
				
			#TO DO: for running time improvement we can collect here the estimates from the solver
			#solver_estimates[$comp_id][$tr_name]=$tr_freq
			#we could use here a 2d array and store all estimations so we can remove the next for loop  
			
			#Multiply each f with the corresponding transcript length
			local_avg_tr_len=$(echo "$local_avg_tr_len + ${tr_lengths_array[$i]} * $tr_i" | bc -l)
			
			i=$[$i + 1]
		fi
	done
	
	#Number of reads divided by average local length
	no_reads_c=`cat ${s2}_no_reads.txt` #Number of reads in the current component
	
	total_sum_reads_portion=$(echo "$total_sum_reads_portion + $no_reads_c / $local_avg_tr_len" | bc -l)
			
											if [ $DEBUG -ne 0 ]; then
												echo "no_reads_c = $no_reads_c"
												echo "local_avg_tr_len = $local_avg_tr_len"
												echo "$no_reads_c / $local_avg_tr_len" | bc -l
												echo "Sum for all reads portions total_sum_reads_portion = $total_sum_reads_portion"
													#exit 7
											fi	
	
	#delete array
	unset tr_lengths_array
	#reset variable
	no_reads_c=0

done

##For each component collect the local results and compute the overall values

components=`ls *_d_values.iso.estimates`
for c in $components
do
	comp_id=`echo $c| cut -d'_' -f 1`
	#array with transcript names
	arr_tr_names=(`cat ${comp_id}_tr_names.txt`)
	no_reads_c=`cat ${comp_id}_no_reads.txt`
	tr_estim=`cat $c`
	i=0
		
	declare -A solver_estimates
	local_avg_tr_len=0 #Average transcripts lengths
	
	#array with transcripts lengths
	tr_lengths_array=(`cat ${comp_id}_tr_lengths.txt`)	
		
		if [ $DEBUG -ne 0 ]; then
			echo -e "\n\nComponent=$c"
			echo "comp_id = $comp_id"
			echo "arr_tr_names=${arr_tr_names[@]}"	
			echo "no_reads_c = $no_reads_c"
			echo "tr_estim = $tr_estim"
			echo "i=$i"
			echo "Tr lengths : ${tr_lengths_array[@]}"
		fi
	

	for t in $tr_estim
	do
		if [[ "$t" == "[" ]]
		then
			continue
		else
			tr_i=`echo $t | cut -d ']' -f 1`	
			
			first_ch=${tr_i:0:1}
			#${string:position:length}
				
				if [ $DEBUG -ne 0 ]; then
					echo "Transcript i is $tr_i"
					echo "First Character = $first_ch"
				fi
			
			if [[ "$first_ch" == "[" ]]
			then
				tr_i=`echo $tr_i | cut -d '[' -f2`
			fi
			
			
			negNum=$(echo "$tr_i < 0" | bc -l)
			#negativeNumber should be 1 if true and 0 if false
			if [[ $negNum == "1" ]]; then
				tr_i=0
			else
				#Solve scientific notations since bc cannot do arithmetic wiht e numbers
				tr_i=`awk '{ printf("%.10f", $1)}' <<< $tr_i`
			fi

			
				if [ $DEBUG -ne 0 ]; then
					echo "New2 tr_i is $tr_i"
				fi
				
			solver_estimates[$i]=$tr_i
			
			local_avg_tr_len=$(echo "$local_avg_tr_len + ${tr_lengths_array[$i]} * $tr_i" | bc -l)
			
						if [ $DEBUG -ne 0 ]; then
							echo "${tr_lengths_array[$i]} * $tr_i" | bc -l
							echo "tr_i = $tr_i"
							echo "tr_lengths_array[$i] = ${tr_lengths_array[$i]}"
							echo "Local average transcripts lengths (L_i * F_i) = $local_avg_tr_len"
						fi
						
			i=$[$i + 1]
		fi
	done
	
							if [ $DEBUG -ne 0 ]; then
								echo "Solver estimates frequencies: ${solver_estimates[@]}"
								echo "Local average transcripts lengths (L_i * F_i) = $local_avg_tr_len"
									#exit 7
							fi
							
	i=0 ##Reset i
	for t in "${arr_tr_names[@]}"
	do
	
		#New Formula (11/7/2013)
		if [[ $total_sum_reads_portion != 0 ]]; then
			frequency=$(echo "${solver_estimates[$i]} * ( ${no_reads_c} / $local_avg_tr_len ) / $total_sum_reads_portion" | bc -l)
		else
			frequency=0
		fi
			if [ $DEBUG -ne 0 ]; then
				echo "Where:"
				echo -e "\tTranscript t=$t and index i=$i"
				echo -e "\tsolver_estimates[$i] = ${solver_estimates[$i]}"
				echo -e "\tno_reads_c = ${no_reads_c}"
				echo -e "\t* local_avg_tr_len = $local_avg_tr_len"
				echo -e "\t / total_sum_reads_portion = $total_sum_reads_portion"
				echo -e "\t$t = $frequency"
				
				#exit 7
			fi
			
		i=$[$i + 1]
		
		echo -e "$t\t${frequency}">>$out_file
		#echo "out_file = $out_file"
		
	done
	
	#delete arrays
	unset arr_tr_names
	unset solver_estimates
	unset tr_lengths_array
	
	#if [[ "$comp_id" == "53" ]]
	#then
	#	exit 7
	#fi
	
done


echo -e "\nMCReg.sh Done!\n"

