#!/bin/bash

#readonly SCRIPT_DIR=$(dirname $0)
readonly SCRIPT_DIR=~/SimRegBest/

#Variables that have not been inserted in the arg option
read_length=100
#Reads Per Fragment
rpf=2

hostname=`hostname`

echo -e "\nRunning: $0 on $hostname ... "

usage() { echo "
Usage: 
   $0 [options]* [-G <GTF_File>]

Mandatory options:   
	-G <GTF_File> File with known genes and isoforms in GTF format (Full path to file) 
	-F <FA_File> Reference FASTA File with transcript sequences   

Optional arguments:
   -s Simulator for generating Monte Carlo reads: grinder or generate-reads [ default: generate-reads (IsoEM-Simulator)]
   -c <int> Coverage [ default: 100 ]
   " 1>&2; exit 1; }

if (($# == 0)); then
        usage
fi


#If your var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":G:F:m:d:s:c:" o; do
    case "${o}" in
        G)
            GTF_File=${OPTARG}
            ;;
		F)
			FA_File=${OPTARG}
			;;
		m)
            mean=${OPTARG}
            ;;
        d)
			deviation=${OPTARG}
			;;
		s)
			simulator=${OPTARG}
			;;
		c)
			coverage=${OPTARG}
			;;
		*)
            echo -e "\nERROR: Invalid option: -$OPTARG" >&2
			usage
            ;;
    esac
done
shift $((OPTIND-1))

if [[ -z "${simulator}" ]];then
	simulator="generate-reads";
fi

if [[ -z "${coverage}" ]];then
	coverage=100;
fi

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


if [ -d "./results" ]; then
  # Control will enter here if $DIRECTORY exists.
  echo -e "\nRemove previous temporary (./results) directory\n"
  rm -rf ./results/
fi

mkdir ./results

echo -e "Run\t\tSimReg\t\tRSEM\t\tIsoEM" > Correlations.txt
#This will delete the previous correlations file

echo -e "\nCompute Total Exon Length (required for computing the total number of reads that need to be simulated): " 
#get total_exon_length:
total_exon_length=`awk '{s+=($5-$4)+1} END {print s}' $GTF_File` 

echo "total_exon_length=$total_exon_length"
#Calculate number of reads
numberOfReads=`echo "($coverage * $total_exon_length) / ($rpf * $read_length)" | bc`

echo "Total number of reads to be simulated: $numberOfReads"

for i in {1..10} 
do
	echo "i=$i"
	echo "Generate Random True Frequency"
	#Extract and normalize true transcript frequencies:
	${SCRIPT_DIR}/utils/extract_tr_freq_uniform.sh $GTF_File 1
	
	#Save true frequencies used 
		#First, get file without path:
		filename=$(basename "$GTF_File")
		extension="${filename##*.}"
		filename="${filename%.*}"
		cp $PWD/$filename.tr.freq.norm $PWD/results/$filename.${i}.tr.freq.norm

	#Simulate reads:	
		if [[ $simulator == "generate-reads" ]]; then
	
				if [ -s $PWD/simulation.properties ];then
					rm -rf $PWD/simulation.properties
				fi
	
				
			#Step 7: Create simulation.properties				
			echo "gtfFile=$GTF_File">>$PWD/simulation.properties
			echo "isoformSequencesFile=$FA_File">>$PWD/simulation.properties
			echo "clusterDistribution = customWeights,$PWD/$filename.tr.freq.norm">>$PWD/simulation.properties
			echo "fragmentLengthDistribution=normal,${mean},${deviation}">>$PWD/simulation.properties
			echo "fragmentStartingPositionDistribution=uniform">>$PWD/simulation.properties
			echo "isoformDistribution=uniform">>$PWD/simulation.properties
			echo "numberOfReads=$numberOfReads">>$PWD/simulation.properties
			echo "readLength=$read_length">>$PWD/simulation.properties
			echo "randomNumberGeneratorSeed=123">>$PWD/simulation.properties
			echo "readsPerFragment=$rpf">>$PWD/simulation.properties
			echo "firstReadOrigin=random">>$PWD/simulation.properties
				
			rm -rf ./nr=*
		
			end_time=`date +%s`
			echo Done: `expr $end_time - $start_time` s.

		
			#Step 8
			echo -e "\nGenerating Reads using \"generate-reads\"..."
			start_time=`date +%s`
			generate-reads 
			wait
			end_time=`date +%s`
			echo Done Generating Simulated Reads! Execution time: `expr $end_time - $start_time` s.
			echo ""
			
			mc_pair1_file=`ls *paired_1.fastq`
			mc_pair2_file=`ls *paired_2.fastq`
			
			#fastq files
			input_option="-q"
			
			refFile=`ls *paired.sam`
			
		else
			echo "Simulating reads using Grinder..."
			
			command -v grinder >/dev/null 2>&1 || { echo >&2 "
			ERROR: grinder NOT Found. 
			Install Grinder Simulator: http://sourceforge.net/projects/biogrinder/
			Aborting..."; exit 7; }
	
			if [ -f "$PWD/grinder.log" ];then
				rm -rf $PWD/grinder*
			fi
			
			grinder -reference_file ${cID}.fa -coverage_fold 1000 -insert_dist ${mean} uniform ${deviation}
			#grinder -reference_file $PWD/../2genes.fa -coverage_fold 100 -insert_dist 300 uniform 30 -abundance_file ./2genes.tr.freq.norm -mutation_dist uniform 0.1 >& $PWD/grinder.log

			wait
			
			#Step 2.2: Split grinder file into m1 and m2
			echo -e "\nSplitting Grinder files into <m1> and <m2> ... " 
			${SCRIPT_DIR}/utils/split_grinder $PWD/grinder-reads.fa
	
			mc_pair1_file=$PWD/grinder-reads_pair_1.fa
			mc_pair2_file=$PWD/grinder-reads_pair_2.fa
		
			#fasta files
			input_option="-f"
		
			refFile="grinder-reads.fa"
		fi
	
	

	
	
	

	



	#Phase II:
	#1. Map Observed Reads
	if [ -f "$PWD/bowtie_OBS_k60v3best.sam" ];then
		rm -rf $PWD/bowtie*
	fi
	
	bowtie --best -v 3 -k 60 -p 36 --chunkmbs 128 ../2genes -f -1 $PWD/grinder-reads_pair_1.fa -2 $PWD/grinder-reads_pair_2.fa -I 180 -X 420 -S $PWD/bowtie_OBS_k60v3best.sam
	wait

	#2. Run SimReg
	#1. SimReg
	~/MCReg_v23/run-mcreg.sh -G $PWD/../2genes.gtf -S $PWD/bowtie_OBS_k60v3best.sam -C $PWD/../SimReg/components/ -F $PWD/../2genes.fa -R $PWD/../2genes -1 $PWD/grinder-reads_pair_1.fa -2 $PWD/grinder-reads_pair_2.fa -m 300 -d 30 -l 100 > mcreg.log
	wait
	
	cp $PWD/MCReg.iso.estimates $PWD/results/MCReg.${i}.iso.estimates
	
	#Run RSEM
	/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --no-qualities $PWD/grinder-reads_pair_1.fa $PWD/grinder-reads_pair_2.fa $PWD/../RSEM/2genes $PWD/rsem
	
	cp $PWD/rsem.isoforms.results $PWD/results/rsem.${i}.isoforms.results
	
	#Run IsoEM
	/home/sahar/isoem-1.1.2/bin/convert-iso-to-genome-coords $PWD/bowtie_OBS_k60v3best.sam ../2genes.gtf $PWD/bowtie_OBS_k60v3best.sam.genome --ignore-pairing

	#run isoem given the frag len distr.
	isoem -G $PWD/../2genes.gtf -m 300 -d 30 $PWD/bowtie_OBS_k60v3best.sam.genome -o IsoEM
	
	
	#Compute Correlation
	awk '{print $2}' ./MCReg.iso.estimates > ./mcreg.iso.estimates_noNames
	sort $PWD/2genes.tr.freq.norm | awk '{print $2}' > ./true.tr.freq.noNames
	
	SimRegCorrel=`python ${SCRIPT_DIR}/scripts/correl.py ./true.tr.freq.noNames ./mcreg.iso.estimates_noNames`
	echo "SimRegCorrel = $SimRegCorrel"
	
	
	#tail -n +2 shows all lines of report from the second line onwards:
	tail -n +2 $PWD/rsem.isoforms.results | sort | awk '{print $6}' > ./rsem.isoforms.results.justFreq
	RSEMCorrel=`python ${SCRIPT_DIR}/scripts/correl.py ./true.tr.freq.noNames ./rsem.isoforms.results.justFreq`
	echo "RSEMCorrel = $RSEMCorrel"
		
	#3. IsoEM Correlation:
	sort ./IsoEM.iso_estimates | awk '{print $2}' > ./IsoEM.iso_estimates.justFreq
	IsoEMCorrel=`python ${SCRIPT_DIR}/scripts/correl.py ./true.tr.freq.noNames ./IsoEM.iso_estimates.justFreq`
	echo "IsoEMCorrel = $IsoEMCorrel"
	echo -e "${i}\t${SimRegCorrel}\t${RSEMCorrel}\t${IsoEMCorrel}" >> Correlations.txt
	
	#exit 7
	
done
