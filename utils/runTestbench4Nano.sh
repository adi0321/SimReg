#!/bin/bash
#This script is for only one testcase (for multiple testcases see runTestbench2.sh)
#Creates 200 subsamples (Bootstrapping)

readonly SCRIPT_DIR=$(dirname $0)

usage() { echo "
Usage: 
   $0 [-d Path2Dir] 

Mandatory options:   
   -d <Path2Dir> Full Path to testbench directory
   -S <SAM_File> File with aligned reads to trnascripts - Bowtie output - (Full path to file)
   -G <GTF_File> File with known genes and isoforms in GTF format (Full path to file)                                                       
   -F <FA_File> Reference FASTA File with transcript sequences
   -R <Ref_File> Reference base name (path with prefix to indexes)
   -E <RSEM_Ref> Reference base name for RSEM indexes (path with prefix to RSEM indexes)
   " 1>&2; exit 7; }

if (($# == 0)); then
        usage
fi

#If var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":d:S:G:F:R:E:T:" o; do
    case "${o}" in
        d)
            Path2Dir=${OPTARG}
            ;;
		S)
			SAM_File=${OPTARG}
			;;
		G)
            GTF_File=${OPTARG}
            ;;
		F)
			FA_File=${OPTARG}
			;;
		R)
			REF_File=${OPTARG}
			;;
		E)
			RSEM_Ref=${OPTARG}
			;;
		T)
			trueFreq_File=${OPTARG}
			;;
		*)
            echo -e "\nERROR: Invalid option: -$OPTARG" >&2
			usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ ! -d "$Path2Dir" ]; then
  # Control will enter here if $DIRECTORY exists.
  echo "Directory: $Path2Dir does not exist"
  exit 1
fi


##Check for several flags:
errs=0
declare -A option=(
    [SAM_File]="-S"
	[GTF_File]="-G"
	[FA_File]="-F"
	[REF_File]="-R"
	[RSEM_Ref]="-E"
	[trueFreq_File]="-T"
)
for var in "${!option[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo -e "ERROR: specify a value for $var with ${option[$var]}"
        ((errs++))
    fi
done
((errs > 0)) && usage

#Variable Declaration
fil_min=180
fil_max=420
mean=200
deviation=50

					echo ""
					echo "Input Parameters:"
					echo "   GTF_File: $GTF_File"
					echo "   FA_File: $FA_File"
					echo "   Reference: $REF_File"
					echo "   RSEM indexes path: $RSEM_Ref"
					echo "   Mean fragment length: $mean"
					echo "   Fragment standard deviation: $deviation"
					echo "   File with true transcript frequencies: $trueFreq_File"
					
cd $Path2Dir

correlFile=${Path2Dir}/Correlations.txt~

if [ -d "${Path2Dir}/results~" ];then
	echo -e "\nRemoving previous results (./results~) directory . . . \n"
	rm -rf ${Path2Dir}/results~
fi

mkdir ${Path2Dir}/results~
echo -e "Run\t\tDate-Time\tCorrel-SimRegBest2\tCorrel-RSEM\tCorrel-IsoEM\tMPE-SimRegBest2\tMPE-RSEM\tMPE-IsoEM" > ${Path2Dir}/Correlations.txt~

cd ${Path2Dir}/results~
	

	#run 200 times to compute 95% CI (Confidence Interval)
	for i in {0..19} 
	do
		mkdir "run$i"
		cd run$i
		
		#For run 0 use the initial SAM file
		if [ $i == 0 ]; then
			echo -e "\nRun $i: Run tools for initial sam file: $SAM_File..."
			SAM=$SAM_File
		else
			echo -e "\nRun $i: SubSampling from $SAM_File..."
			
			startTime=`date +%s`
			
			${SCRIPT_DIR}/subSampling $SAM_File 
			#It outputs in subSample.sam
			SAM=$PWD/subSample.sam
			
			endTime=`date +%s`
			totalTimeSec=`expr $endTotalTime - $startTotalTime`
			echo "Done SAM! Execution Time:"
			${SCRIPT_DIR}/utils/seconds.sh $totalTimeSec
		fi

			#############################	
			#	Run SimRegBest2			#
			#############################
			echo "Running Latest SimReg version"
			startTime=`date +%s`
			
			mkdir ./SimReg
			cd ./SimReg
	
			${SCRIPT_DIR}/../run-mcreg.sh -G $GTF_File -S $SAM -F $FA_File -R $REF_File -m $mean -d $deviation > simreg.log

			#Extract NanoCounts
			${SCRIPT_DIR}/extractTr4NanoStringsSimReg.sh $PWD/MCReg.iso.estimates
			wait
			#Compute Correlation
			sort ./MCReg.iso.estimates.NanoStrings.txt | awk '{print $2}' > ./mcreg.iso.estimates_noNames

			if [ ! -f "${trueFreq_File}.justFreq" ]; then
				#For Grinder True Freq
				#tail -n +2 $trueFreq_File | sort -k 2 | awk '{print $3}' > ${trueFreq_File}.justFreq
				#remove header for Grinder ranks and sort based on the 2nd col
				
				#For NanoCounts
				sort $trueFreq_File | awk '{print $2}' > ${trueFreq_File}.justFreq
			fi
	
			SimRegCorrel=`python ${SCRIPT_DIR}/../scripts/correl.py ${trueFreq_File}.justFreq ./mcreg.iso.estimates_noNames`
			echo "SimRegCorrel = $SimRegCorrel"
	
			#Compute MPE
			SimRegMPE=`python ${SCRIPT_DIR}/../scripts/mpe.py ${trueFreq_File}.justFreq ./mcreg.iso.estimates_noNames`
			echo "SimRegMPE = $SimRegMPE"
	
			
			endTime=`date +%s`
			totalTimeSec=`expr $endTime - $startTime`
			echo "Done SimReg! Execution Time:"
			${SCRIPT_DIR}/seconds.sh $totalTimeSec
	
			#################	
			#	Run RSEM	#
			#################
	
			echo -e "\n\n## Running RSEM ... "
			startTime=`date +%s`
			
			mkdir ../RSEM
			cd ../RSEM
		
			if [ ! -d $RSEM_Ref ];then
				echo "Prepare RSEM refernece in $PWD"
				mkdir rsemIndex
				RSEM_Ref=$PWD/rsemIndex
				/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-prepare-reference $FA_File $RSEM_Ref > rsem.index.log
			else
				echo "RSEM indexes: OK"
			fi
			
			#/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --no-qualities $PWD/../grinder-reads_pair_1.fa $PWD/../grinder-reads_pair_2.fa $PWD/${i} $PWD/rsem > rsem.log

			#use SimReg Bowtie alignment
			/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --sam --paired-end $SAM $RSEM_Ref $PWD/rsem > rsem.log

			#Extract NanoCounts
			${SCRIPT_DIR}/extractTr4NanoStringsRSEM.sh $PWD/rsem.isoforms.results
			wait
			
			#Compute Correlation for Entire Genome (usually reads from Grinder where we know true freq for all tr.)
			#tail -n +2 shows all lines of report from the second line onwards:
			#tail -n +2 $PWD/rsem.isoforms.results | sort | awk '{print $6}' > ./rsem.isoforms.results.justFreq
			sort $PWD/rsem.isoforms.results.NanoStrings.txt | awk '{print $2}' > ./rsem.isoforms.results.justFreq
			
			RSEMCorrel=`python ${SCRIPT_DIR}/../scripts/correl.py ${trueFreq_File}.justFreq ./rsem.isoforms.results.justFreq`
			echo "RSEMCorrel = $RSEMCorrel"
			
			#Compute MPE
			RSEMMPE=`python ${SCRIPT_DIR}/../scripts/mpe.py ${trueFreq_File}.justFreq ./rsem.isoforms.results.justFreq`
			echo "RSEMMPE = $RSEMMPE"

			
			endTime=`date +%s`
			totalTimeSec=`expr $endTime - $startTime`
			echo "Done RSEM! Execution Time:"
			${SCRIPT_DIR}/seconds.sh $totalTimeSec
			
			#################	
			#	Run IsoEM	#
			#################						
			echo -e "\n\n## Running IsoEM ... "
			startTime=`date +%s`
			mkdir ../IsoEM
			cd ../IsoEM
			
			if [ ! -f ${SAM}.genome ]; then
				echo "Convert to genome coordinates"
				/home/sahar/isoem-1.1.2/bin/convert-iso-to-genome-coords $SAM $GTF_File ${SAM}.genome --ignore-pairing
			else
				echo "IsoEM: Genome coordinates: OK"
			fi
			
			
			#run isoem given the frag len distr.
			isoem -G $GTF_File -m $mean -d $deviation ${SAM}.genome -o IsoEM
			wait
			
			#Extract NanoCounts
			${SCRIPT_DIR}/extractTr4NanoStringsSimReg.sh IsoEM.iso_estimates
			wait

			#IsoEM Correlation:
			#sort ./IsoEM.iso_estimates | awk '{print $2}' > ./IsoEM.iso_estimates.justFreq
			sort ./IsoEM.iso_estimates.NanoStrings.txt | awk '{print $2}' > ./IsoEM.iso_estimates.justFreq
			
			IsoEMCorrel=`python ${SCRIPT_DIR}/../scripts/correl.py ${trueFreq_File}.justFreq ./IsoEM.iso_estimates.justFreq`
			echo "IsoEMCorrel = $IsoEMCorrel"

			#Compute MPE
			IsoEMMPE=`python ${SCRIPT_DIR}/../scripts/mpe.py ${trueFreq_File}.justFreq ./IsoEM.iso_estimates.justFreq`
			echo "IsoEMMPE = $IsoEMMPE"
			
			
			endTime=`date +%s`
			totalTimeSec=`expr $endTime - $startTime`
			echo "Done IsoEM! Execution Time:"
			${SCRIPT_DIR}/seconds.sh $totalTimeSec
			#################
	
	
			currDate=`date '+%m.%d-%T'`
			echo -e "Run${i}\t${currDate}\t${SimRegCorrel}\t${RSEMCorrel}\t${IsoEMCorrel}\t${SimRegMPE}\t${RSEMMPE}\t${IsoEMMPE}" >> $correlFile

			echo -e "\nDone simulated reads: $s \n\n"
			##Move out of directory $s
			cd ../..
	
		echo -e "Done run: $i\n\n"
	done
	