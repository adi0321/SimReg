#!/bin/bash

#****************************************************************************
#*	Credits:	Adrian Caciula - Georgia State University, Atlanta, GA		*
#*			*
#*								*****************************
#*	Date Created: 11/10/2013					*
#*	Last Update: 11/10/2013
#* Map observed reads and create the connected components						*
#************************************************

##ToDO: Add option for D values (give the file if you have it and do not compute it anymore)

readonly DEBUG=0
readonly SCRIPT_DIR=$(dirname $0)

hostname=`hostname`

echo -e "\nRunning: $0 v.2.1 on $hostname ... "

usage() { echo "
Usage: 
   $0 [options]* [-m <int>] [-d <int>] {-1 <m1> -2 <m2> -R <Ref_File> | -S <SAM_File>}

Mandatory options:                                                       
   -m <int> Fragment length mean                   
   -d <int> Fragment length standard deviation
   -1 <m1>  FA or FASTQ File containing sequences paired with mates in <m2>
   -2 <m2>  FA or FASTQ File containing sequences paired with mates in <m1>
   -R <Ref_File> Reference base name (transcript sequencies and indexes)
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
while getopts ":m:d:R:S:C:l:r:1:2:" o; do
    case "${o}" in
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
		1)
			obs_pair1_file=${OPTARG}
			;;
		2)
			obs_pair2_file=${OPTARG}
			;;
		C)
			CC_Path=${OPTARG}
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


if [ -d "./tmp" ]; then
  # Control will enter here if $DIRECTORY exists.
  echo -e "\nRemove previous temporary (./tmp) directory\n"
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
		[mean]="-m"
		[deviation]="-d"
	)
	for var in "${!option[@]}"; do
		if [[ -z "${!var}" ]]; then
			echo -e "ERROR: specify a value for $var with ${option[$var]}"
			((errs++))
		fi
	done
	echo "or use -S option for transcriptome alignment"
	((errs > 0)) && usage

	if [ ! -f $obs_pair1_file ];then
		echo "Error: File $obs_pair1_file not found!"
		exit 7
	fi
	
	if [ ! -f $obs_pair2_file ];then
		echo "Error: File $obs_pair2_file not found!"
		exit 7
	fi
	
					echo ""
					echo "Input Parameters:"
					echo "   GTF_File: $GTF_File"
					echo "   Mean fragment length: $mean"
					echo "   Fragment standard deviation: $deviation"
					echo "   Read length: $read_length"
					echo "   Number of reads per fragment: $rpf"
					echo "   Reads Simulator: $simulator"
					echo ""
					echo "Other Input Parameters:"
					echo "   REF_File: $REF_File"
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
	#bowtie -k 60 -p 12 --chunkmbs 128 $REF_File ${input_option} -1 $obs_pair1_file -2 $obs_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_60multiAligns.sam
	
	#-v <int> Report alignments with at most <int> mismatches.
	bowtie -k 60 -p 12 -v 0 --chunkmbs 128 $REF_File ${input_option} -1 $obs_pair1_file -2 $obs_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_60multiAligns.sam
	wait
	end_time=`date +%s`
	echo Done Mapping observed reads! Execution Time: `expr $end_time - $start_time`s.
	echo ""
	
	#Does the SAM file requires sorting?
	SAM_File=$PWD/bowtie_OBS_60multiAligns.sam
fi

if [ ! -f ${SAM_File} ];then
	echo "Error: File $SAM_File not found!"
	exit 7
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
	${SCRIPT_DIR}/../lib/compute_obsRC_CC $SAM_File
	end_time=`date +%s`
	echo Done Computing Connected Components! Execution Time: `expr $end_time - $start_time`s.
	echo ""
	CC_Path=$PWD
else
	#Else if path to connected components is given we'll need to copy the singleTrGenes.txt file
	cp ${CC_Path}/../singleTrGenes.txt ../
fi

