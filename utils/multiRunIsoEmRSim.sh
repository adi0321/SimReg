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
   $0 [options]* [-G <GTF_File>] [-F <Path to FA_File>] [-R <Path to Reference Files>] [-m <int>] [-d <int>]

Mandatory options:   
	-G <GTF_File> File with known genes and isoforms in GTF format (Full path to file) 
	-F <FA_File> Reference FASTA File with transcript sequences   
	-R <Ref_File> Reference base name (transcript sequencies and indexes)
	-I <Ref File> RSEM Indexes

Optional arguments:
   -s Simulator for generating Monte Carlo reads: grinder or generate-reads [ default: generate-reads (IsoEM-Simulator)]
   -c <int> Coverage [ default: 100 ]
   " 1>&2; exit 1; }

if (($# == 0)); then
        usage
fi


#If your var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":G:F:R:I:m:d:s:c:" o; do
    case "${o}" in
        G)
            GTF_File=${OPTARG}
            ;;
		F)
			FA_File=${OPTARG}
			;;
		R)
			REF_File=${OPTARG}
			;;
		I)
			RSEM_Index=${OPTARG}
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
	[REF_File]="-R"
	[RSEM_Index]="-I"
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
#numberOfReads=`echo "($coverage * $total_exon_length) / ($rpf * $read_length)" | bc`
numberOfReads=40000000
echo "Total number of reads to be simulated: $numberOfReads"

fil_min=$(($mean-(3*$deviation)))
fil_max=$(($mean+(3*$deviation)))

for i in {1..10} 
do
	echo "i=$i"
	
	#Simulate reads:	
		if [[ $simulator == "generate-reads" ]]; then
			
			echo "Generate Random GENE True Frequency"
			${SCRIPT_DIR}/utils/generateRandGeneFreq.sh $GTF_File 0
	
			#Save true frequencies used 
			#First, get file without path:
			filename=$(basename "$GTF_File")
			extension="${filename##*.}"
			filename="${filename%.*}"
			cp $PWD/$filename.genes.freq $PWD/results/$filename.${i}.genes.freq
	
	
				if [ -s $PWD/simulation.properties ];then
					rm -rf $PWD/simulation.properties
				fi
	
				
			#Step 7: Create simulation.properties				
			echo "gtfFile=$GTF_File">>$PWD/simulation.properties
			echo "isoformSequencesFile=$FA_File">>$PWD/simulation.properties
			echo "clusterDistribution = customWeights,$PWD/$filename.genes.freq">>$PWD/simulation.properties
			echo "fragmentLengthDistribution=normal,${mean},${deviation}">>$PWD/simulation.properties
			echo "fragmentStartingPositionDistribution=uniform">>$PWD/simulation.properties
			echo "isoformDistribution=uniform">>$PWD/simulation.properties
			echo "numberOfReads=$numberOfReads">>$PWD/simulation.properties
			echo "readLength=$read_length">>$PWD/simulation.properties
			echo "randomNumberGeneratorSeed=123">>$PWD/simulation.properties
			echo "readsPerFragment=$rpf">>$PWD/simulation.properties
			echo "firstReadOrigin=random">>$PWD/simulation.properties
				
			rm -rf ./nr=*
		
		
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
			
			trueTrFreq=`ls *iso_weights`
			#fastq files
			input_option="-q"
			
			refFile=`ls *paired.sam`

		else
			echo "Simulating reads using Grinder..."
			
			command -v grinder >/dev/null 2>&1 || { echo >&2 "
			ERROR: grinder NOT Found. 
			Install Grinder Simulator: http://sourceforge.net/projects/biogrinder/
			Aborting..."; exit 7; }
	
			echo "Generate Random Transcripts True Frequency"
			${SCRIPT_DIR}/utils/generateRandTrFreq.sh $GTF_File 0
	
			#Save true frequencies used 
			#First, get file without path:
			filename=$(basename "$GTF_File")
			extension="${filename##*.}"
			filename="${filename%.*}"
			cp $PWD/$filename.tr.freq $PWD/results/$filename.${i}.tr.freq
			trueTrFreq=$PWD/$filename.tr.freq
	
			if [ -f "$PWD/grinder.log" ];then
				rm -rf $PWD/grinder*
			fi
			
			#grinder -reference_file ${cID}.fa -coverage_fold 100 -insert_dist ${mean} uniform ${deviation}
			#grinder -reference_file $FA_File -coverage_fold 100 -insert_dist 300 uniform 30 -abundance_file $trueTrFreq -mutation_dist uniform 0.1 >& $PWD/grinder.log
			grinder -reference_file $FA_File -total_reads $numberOfReads -insert_dist 300 uniform 30 -abundance_file $trueTrFreq -mutation_dist uniform 0.1 >& $PWD/grinder.log

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
	echo "Map Observed Reads"
	if [ -f "$PWD/bowtie_OBS_k60v3best.sam" ];then
		rm -rf $PWD/bowtie*
	fi
	
	bowtie --best -v 3 -k 60 -p 36 --chunkmbs 128 $REF_File $input_option -1 $mc_pair1_file -2 $mc_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_k60v3best.sam
	wait

	echo "Run SimReg"
	#1. SimReg
	nohup nice ~/SimRegBest/run-mcreg.sh -G $GTF_File -S $PWD/bowtie_OBS_k60v3best.sam -F $FA_File -R /import1/UCSC/hg18/indexes/hg18_ref_genome -1 $mc_pair1_file -2 $mc_pair2_file -m $mean -d $deviation -l 100 >& mcreg.log
	wait
	
	cp $PWD/MCReg.iso.estimates $PWD/results/MCReg.${i}.iso.estimates
	
	echo "Run RSEM"
	#/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --no-qualities $mc_pair1_file $mc_pair2_file $PWD/../RSEM/2genes $PWD/rsem
	/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --sam $PWD/bowtie_OBS_k60v3best.sam $RSEM_Index $PWD/rsem
	
	cp $PWD/rsem.isoforms.results $PWD/results/rsem.${i}.isoforms.results
	
	echo "Run IsoEM"
	/home/sahar/isoem-1.1.2/bin/convert-iso-to-genome-coords $PWD/bowtie_OBS_k60v3best.sam $GTF_File $PWD/bowtie_OBS_k60v3best.sam.genome --ignore-pairing

	#run isoem given the frag len distr.
	isoem -G $GTF_File -m $mean -d $deviation $PWD/bowtie_OBS_k60v3best.sam.genome -o IsoEM
	
	
	#Compute Correlation
	awk '{print $2}' ./MCReg.iso.estimates > ./mcreg.iso.estimates_noNames
	sort $trueTrFreq | awk '{print $2}' > ./true.tr.freq.noNames
	
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
	
	exit 7

done
