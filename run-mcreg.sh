#!/bin/bash

#****************************************************************************
#*	Credits:	Adrian Caciula - Georgia State University, Atlanta, GA		*
#*												*****************************
#*	Date Created: 11/10/2013					*
#*	Last Update: 06/12/2014						*
#************************************************

#This is an extension of best MCReg version 21

##ToDO: Add option for D values (give the file if you have it and do not compute it anymore)

readonly DEBUG=0
readonly SCRIPT_DIR=$(dirname $0)

startTotalTime=`date +%s`

hostname=`hostname`

echo -e "\nRunning: $0 v.2.1 on $hostname ... "

usage() { echo "
Usage: 
   $0 [options]* [-G <GTF_File>] [-F <FA_File>] [-m <int>] [-d <int>] {-1 <m1> -2 <m2> -R <Ref_File> | -S <SAM_File>}

Mandatory options:   
   -G <GTF_File> File with known genes and isoforms in GTF format (Full path to file)                                                       
   -m <int> Fragment length mean                   
   -d <int> Fragment length standard deviation
   -F <FA_File> Reference FASTA File with transcript sequences
   -1 <m1>  FA or FASTQ File containing sequences paired with mates in <m2>
   -2 <m2>  FA or FASTQ File containing sequences paired with mates in <m1>
   -R <Ref_File> Reference base name (path with prefix to indexes)
   -S <SAM_File> File with aligned reads to trnascripts - Bowtie output - (Full path to file)
   (The following Bowtie command is used to align observed reads if no alignment is given:
   bowtie -k 60 -p 12 --chunkmbs 128 [REF_File] {-f | -g} -1 [obs_pair1_file] -2 [obs_pair2_file] -I [fragInsLenMin] -X [fragInsLenMax] -S [outFile])
	-C <Path> - Path to connected components (from observed reads)
	
Optional arguments:
   -l <int> Read length [ default: 100 ]
   -r <int> Number of reads per fragment [ default: 2 ]
   -g	Use Grinder Simulator for generation Monte Carlo reads [ default: generate-reads]
   " 1>&2; exit 7; }

if (($# == 0)); then
        usage
fi

#TO DO: What about if the sequence of candidate transcripts. This must be extracted from the gtf (See cufflinks it has a script that do such extractions)
#Maybe use that script if -F option is not provided???
##(in other words we will try to move -F option as optional arguments)

   
#If your var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":G:m:d:R:S:F:C:l:r:1:2:t:g" o; do
    case "${o}" in
        G)
            GTF_File=${OPTARG}
            ;;
        m)
            mean=${OPTARG}
            ;;
        d)
			deviation=${OPTARG}
			;;
		R)
			REF_File=${OPTARG}
			;;
		S)
			SAM_File=${OPTARG}
			;;
		l)
			read_length=${OPTARG}
			;;
		r)
			#reads per fragment
			rpf=${OPTARG}
			;;
		g)
			simulator="grinder"
			;;
		1)
			obs_pair1_file=${OPTARG}
			;;
		2)
			obs_pair2_file=${OPTARG}
			;;
		F)
			FA_File=${OPTARG}
			;;
		C)
			CC_Path=${OPTARG}
			;;
		t)
			precision=${OPTARG}
			#tuning
			;;
		*)
            echo -e "\nERROR: Invalid option: -$OPTARG" >&2
			usage
            ;;
    esac
done
shift $((OPTIND-1))

# -z string
#		True if the length of string is zero.

if [[ -z "${read_length}" ]];then
	read_length=100;
fi

if [[ -z "${rpf}" ]];then
	rpf=2;
fi

if [[ -z "${simulator}" ]];then
	#simulator="generate-reads";
	simulator="sim-reads";
fi

if [[ -z "${precision}" ]];then
	precision=0;
fi

##Only check for several flags here
errs=0
declare -A option=(
    [GTF_File]="-G"
	[FA_File]="-F"
    [mean]="-m"
    [deviation]="-d"
)
for var in "${!option[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo -e "ERROR: specify a value for $var with ${option[$var]}"
        ((errs++))
    fi
done
((errs > 0)) && usage

if [ ! -f ${GTF_File} ];then
	echo "Error: File $GTF_File not found!"
	exit 1
fi

if [ ! -f ${FA_File} ];then
	echo "Error: File $FA_File not found!"
	exit 1
fi

					echo ""
					echo "Input Parameters:"
					echo "   GTF_File: $GTF_File"
					echo "   FA_File: $FA_File"
					echo "   Reference: $REF_File"
					echo "   Mean fragment length: $mean"
					echo "   Fragment standard deviation: $deviation"
					echo "   Read length: $read_length"
					echo "   Number of reads per fragment: $rpf"
					echo "   Reads Simulator: $simulator"
					echo "   Precision Level = $precision"


if [ -d "./tmp" ]; then
  # Control will enter here if $DIRECTORY exists.
  echo -e "\nRemoving previous temporary (./tmp) directory . . . \n"
  rm -rf ./tmp/
fi

mkdir ./tmp/

cd ./tmp/

#if SAM File with transcriptome alignment DOES NOT Exists then
if [[ -z "${SAM_File}" ]]; then

	echo "No SAM file provided -- Reads will be mapped to given reference using Bowtie"
	errs=0
	declare -A option=(
		[REF_File]="-R"
		[obs_pair1_file]="-1"
		[obs_pair2_file]="-2"
	)
	for var in "${!option[@]}"; do
		if [[ -z "${!var}" ]]; then
			echo -e "ERROR: specify a value for $var with ${option[$var]}"
			((errs++))
		fi
	done
	
	((errs > 0)) && echo "or use -S option for transcriptome alignment" &&usage

	if [ ! -f $obs_pair1_file ];then
		echo "Error: File $obs_pair1_file not found!"
		exit 7
	fi
	
	if [ ! -f $obs_pair2_file ];then
		echo "Error: File $obs_pair2_file not found!"
		exit 7
	fi
	
					echo ""
					echo "Input Parameters for Bowtie Mapping:"
					echo "   Reference: $REF_File"
					echo "   Pair 1 sequences file: $obs_pair1_file"
					echo "   Pair 2 sequences file: $obs_pair2_file"

 
	##Step1 - Create components from observed reads 
		##Stept 1.A - Map Observed reads using bowtie
	
	echo -e "\nMapping all observed reads, using Bowtie, to $GTF_File ..."
	
	#detect file type
	filename=$(basename "$obs_pair1_file")
	extension="${filename##*.}"
	echo -e "\tInput reads file extension = $extension"

	if [[ $extension == "fa" ]]; then
		input_option="-f"
	else
		input_option="-q"
	fi
	
	#Fragment insert length range
	fil_min=$(($mean-(4*$deviation)))
	fil_max=$(($mean+(4*$deviation)))
	
	echo -e "\tInput Option is: $input_option"
	echo -e "\tFragment Insert Length Min: $fil_min"
	echo -e "\tFragment Insert Length Max: $fil_max"
	

	start_time=`date +%s`
	echo "bowtie --best -v 3 -k 60 -p 16 --chunkmbs 128 $REF_File ${input_option} -1 $obs_pair1_file -2 $obs_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_k60v3best.sam"
	bowtie --best -v 3 -k 60 -p 16 --chunkmbs 128 $REF_File ${input_option} -1 $obs_pair1_file -2 $obs_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_k60v3best.sam
	
	#-v <int> Report alignments with at most <int> mismatches.
	#bowtie -k 60 -p 12 -v 0 --chunkmbs 128 $REF_File ${input_option} -1 $obs_pair1_file -2 $obs_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_60multiAligns.sam
	wait
	end_time=`date +%s`
	echo Done Mapping observed reads! Execution Time: `expr $end_time - $start_time`s.
	echo ""
	
	#Does this SAM file requires sorting?
	SAM_File=$PWD/bowtie_OBS_k60v3best.sam
fi

if [ ! -f ${SAM_File} ];then
	echo "Error: File $SAM_File not found!"
	exit 1
fi

echo "SAM_File = $SAM_File"


mkdir ./components
cd ./components
	
##################################################
#	Results File HEADER	################################
##################################################	
touch ../results.txt
echo -e "C.ID #Tr. #ObsRCmp\tTr.Name\tLclFreq.\tAvgTrLen">>../results.txt
#Where: 
#		1. C.ID - component ID
#		2. #Tr. - number of transcripts in current component
#		3. #ObsRCmp - number of observed reads per component
#		6. AvgTrLen - average transcript length in current component
#####################################################################################################
 
 
 #Correlation for each component
 touch ../comp_correl.txt
 echo -e "CmpID Correlation\t#Transcripts" >> ../comp_correl.txt
 
if [[ -z "${CC_Path}" ]]; then
	#If path to connected components is not given then

	#create the file for isolated genes with only one transcript
	touch ../singleTrGenes.txt
	
	##Step 1.B - Compute Connected Components
	echo -e "\nComputing Connected Components..."
	start_time=`date +%s`
	${SCRIPT_DIR}/lib/compute_obsRC_CC $SAM_File $GTF_File $FA_File
	
	if [[ $? -eq "1" ]];
	then
		echo "Error in compute_obsRC_CC"
		exit 1
	fi
	
	end_time=`date +%s`
	tsec=`expr $end_time - $start_time`
	echo "Done Computing Connected Components! Execution Time:" 
	${SCRIPT_DIR}/utils/seconds.sh $tsec
	echo ""

	CC_Path=$PWD
else
	#Else if path to connected components is given we'll need to copy the singleTrGenes.txt file
	cp ${CC_Path}/../singleTrGenes.txt ../
fi


##############################
#echo -e "\nStep 2 - Run MCReg for each component (that has more than 1 transcript)"
##############################

	if [ $DEBUG -ne 0 ]; then
		components=`ls $CC_Path`
		echo "$components" | head
	fi

	
#split the components into procs, where procs is the number of processors used
#procs=10
#First get the total number of components
#tNComp=`ls $CC_Path | wc -l`
#echo "Total number of components: $tNComp"

echo "Path to Connected Components: $CC_Path"

#for p in $(seq 1 $procs)
#do
	#Simple implementation
	#echo "p=$p"
	#components=`ls -d $CC_Path/${p}*`
	#components=`ls -d $CC_Path/*`
	#this step is now done in simReg.sh
	
	#echo "components: "
	#echo $components | head
	
		if [ $DEBUG -ne 0 ]; then
			echo "GTF_File = $GTF_File"
			echo "FA_File = $FA_File"
			echo "read_length = $read_length"
			echo "rpf = $rpf"
			echo "REF File = $REF_File"
			echo "precision = $precision"
		fi
		
	#${SCRIPT_DIR}/scripts/simReg.sh -m $mean -d $deviation -l $read_length -r $rpf -G $GTF_File -F $FA_File -s $SCRIPT_DIR -C $CC_Path
	${SCRIPT_DIR}/simreg -G $GTF_File -F $FA_File -m $mean -d $deviation -l $read_length -t $precision -s $SCRIPT_DIR -C $CC_Path

#done

echo "Move out from components directory"
echo "Current directory before moving out: $PWD"
cd ..

out_file=../MCReg.iso.estimates

if [[ -s ${out_file} ]]; then
	rm -rf ${out_file}
fi

pwd

echo -e "\nProcess the final results file"
#1st Create the reasults.txt file inside tmp by concatenating all the results file from each connected component
components=`ls $CC_Path`
for c in $components
do
	if [ -f $PWD/components/${c}/results.txt ];then
		cat $PWD/components/${c}/results.txt >> results.txt
	fi
done


#2nd
#Add the results from singleTrGenes.txt --- here we need to extract the transcript length for each gene
#And concatenate the results at the end to results.txt

#Compute total_sum_reads_portion		
	
# Load files in C++
${SCRIPT_DIR}/lib/combineFiles "singleTrGenes.txt" "trLen.txt" "results.txt"

echo "Compute Total Sum for Reads Portion of each Component"

total_sum_reads_portion=`awk '
BEGIN{
ID=-1;
flag=0;
}
{
	if(flag==0)
	{
		flag=1;
		#skip the header line
		next;
	}
	
	if(ID != $1)
	{
		ID=$1;
		#print "dolar 3 = ", $3;
		#print "dolar 6 = ", $6;
		
		if( ($6 != 0) && (length($6) != 0))
			total+=$3/$6;
	}
}
END{print total}' results.txt`

#echo -e "\nCompute total_sum_reads_portion=$total_sum_reads_portion + $no_reads_c / $local_avg_tr_len"

echo "total_sum_reads_portion=$total_sum_reads_portion"

##Compute the estimates and report the frequency
#What about the transcripts that are not covered by observed reads
#We also need to report gene frequency 
#maybe load all tr names in a structure and then only add the freq. and the rest are zero

#compute and print tr frequency
${SCRIPT_DIR}/lib/computeTrFreq $GTF_File "results.txt" $total_sum_reads_portion

endTotalTime=`date +%s`
totalTimeSec=`expr $endTotalTime - $startTotalTime`
echo "Finished! Total Execution Time:"
${SCRIPT_DIR}/utils/seconds.sh $totalTimeSec

