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
	-S <SAM_File> File with aligned reads to trnascripts - Bowtie output - (Full path to file)
	-1 <m1>  FA or FASTQ File containing sequences paired with mates in <m2>
	-2 <m2>  FA or FASTQ File containing sequences paired with mates in <m1>

Optional arguments:
   -s Simulator for generating Monte Carlo reads: grinder or generate-reads [ default: generate-reads (IsoEM-Simulator)]
   -c <int> Coverage [ default: 100 ]
   " 1>&2; exit 1; }

if (($# == 0)); then
        usage
fi


#If your var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":G:F:R:I:S:1:2:m:d:c:" o; do
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
		S)
			SAM_File=${OPTARG}
			;;
		m)
            mean=${OPTARG}
            ;;
        d)
			deviation=${OPTARG}
			;;
		c)
			coverage=${OPTARG}
			;;
		1)
			pair1_file=${OPTARG}
			;;
		2)
			pair2_file=${OPTARG}
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
	[SAM_File]="-S"
	[RSEM_Index]="-I"
	[mean]="-m"
    [deviation]="-d"
	[pair1_file]="-1"
	[pair2_file]="-2"
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
#numberOfReads=100000000
echo "Total number of reads to be simulated: $numberOfReads"

fil_min=$(($mean-(3*$deviation)))
fil_max=$(($mean+(3*$deviation)))

for i in {1..10} 
do
	echo "i=$i"
	
	


	echo "Run SimReg"
	#1. SimReg
	#nohup nice ~/SimRegBest/run-mcreg.sh -G $GTF_File -S $PWD/bowtie_OBS_k60v3best.sam -F $FA_File -R /import1/UCSC/hg18/indexes/hg18_ref_genome -1 $pair1_file -2 $pair2_file -m $mean -d $deviation -l 100 >& mcreg.log
	wait
	
	#cp $PWD/MCReg.iso.estimates $PWD/results/MCReg.${i}.iso.estimates
	
	echo "Run RSEM"
	#/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --no-qualities $mc_pair1_file $mc_pair2_file $PWD/../RSEM/2genes $PWD/rsem
	/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --sam $PWD/bowtie_OBS_k60v3best.sam $RSEM_Index $PWD/rsem
exit 7	
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
