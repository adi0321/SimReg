#!/bin/bash

readonly SCRIPT_DIR=$(dirname $0)

usage() { echo "
Usage: 
   $0 [-d Path2Dir] 

Mandatory options:   
   -d <Path2Dir> Full Path to testbench directory
   " 1>&2; exit 7; }

if (($# == 0)); then
        usage
fi

#If var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":d:" o; do
    case "${o}" in
        d)
            Path2Dir=${OPTARG}
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

cd $Path2Dir

correlFile=${Path2Dir}/Correlations.txt~

if [ ! -d "${Path2Dir}/results~" ];then
	mkdir ${Path2Dir}/results~
	echo -e "Date\tTestCaseName\tObservedReads\tSimReg_LatestVersion\tSimReg_Version21\tRSEM\tIsoEM" > ${Path2Dir}/Correlations.txt~
fi



testbench=`ls -B`
#	-B, --ignore-backups
#              do not list implied entries ending with ~

echo -e "testbenches: \n$testbench"

for i in $testbench 
do
	
	echo "i=$i"
	cd ${Path2Dir}/$i
	
	#Check FA File
	
	if [ ! -f $PWD/${i}.fa ]; then
		echo "FASTA File: Missing!"
		echo -e "\tExtracting FASTA File (reference sequences)..."
		#Extract the reference fa sequences
		${SCRIPT_DIR}/utils/make_refFile.sh $PWD/${i}.gtf
		
		echo -e "\tCreating Bowtie indexes..."
		bowtie-build $PWD/${i}.fa $i > bowtie-build.log
		echo "Done"
	else
		echo "FASTA File: Ok!"
	fi
	
	if [ ! -f $PWD/${i}.tr.freq ];then
		echo "True Frequencies File: Missing!"
		echo -e "\tExtract and normalize true transcript frequencies..."
		${SCRIPT_DIR}/utils/extract_tr_freq.sh $PWD/${i}.gtf
	else
		echo "True Frequencies File: Ok!"
	fi
	
	#!!This is a temp statement
	rm -rf ./SimReg
	##Delete it after debuging ends
	
	if [ ! -d "SimReg" ];then
	
		echo "Preprocessed SimReg directory: Missing!"
		echo -e "\tProcessing the annotation ..."
		mkdir SimReg
		cd SimReg
			echo -e "\t###	Phase I:"
			echo -e "\t# Simulate reads with same parameters used by observed reads ... "
			fragLen=300
			readLen=100
			echo -e "\tFragment Lenght: $fragLen"
			echo -e "\tRead Length: $readLen"
			echo -e "\tSimulating reads ... "
			${SCRIPT_DIR}/utils/sim-reads $PWD/../${i}.fa $fragLen $readLen > sim-reads.log
			wait
			
			echo "#Map the reads:"
			#bowtie --best -v 0 -k 60 -p 36 --chunkmbs 128 ../$i -f -1 $PWD/simreg-reads-pair1.fa -2 $PWD/simreg-reads-pair2.fa -I 180 -X 420 -S $PWD/bowtie_SimReg_k60v0best.sam > bowtieSimReg.log

			bowtie -v 0 -a -p 12 --chunkmbs 128 ../$i -f -1 $PWD/simreg-reads-pair1.fa -2 $PWD/simreg-reads-pair2.fa -I 1 -X 1000 -S $PWD/bowtie_SimReg_k60v0best.sam > bowtieSimReg.log
			
			echo "#Step 3: Parsing SimReg alignment and create the Read Classes" 
			${SCRIPT_DIR}/scripts/makeSimRegRC_CC.sh -S $PWD/bowtie_SimReg_k60v0best.sam -1 $PWD/simreg-reads-pair1.fa
			wait
			echo "Done Preprocessing Step!"

		cd ..
	else
		echo "Preprocessed SimReg directory: Ok!"
	fi
	
	#obsReads=`ls -d grinder*/`
	obsReads="grinderErrorFree	grinderMutUniform0.1"
	echo "obsReads = $obsReads"
exit 7	
	for g in $obsReads
	do
	
		echo "g = $g"
		
		if [ ! -d ${Path2Dir}/${i}/${g}/ ]; then
			echo "Directory ${Path2Dir}/${i}/${g}/ : Missing!"
			echo "Simulate Observed Reads using Grinder:"
			mkdir ${Path2Dir}/${i}/${g}/
			cd ${Path2Dir}/${i}/${g}/
			
			if [ $g == "grinderErrorFree" ]; then
				echo "\tSimulating error FREE reads ..."
				grinder -reference_file ../${i}.fa -coverage_fold 100 -insert_dist 300 uniform 30 -abundance_file ../${i}.tr.freq.norm > grinder.log
				
			else
			
				echo -e "\tSimulating reads using mutation distribution: uniform 0.1 ... "
				grinder -reference_file ../${i}.fa -coverage_fold 100 -insert_dist 300 uniform 30 -abundance_file ../${i}.tr.freq.norm -mutation_dist uniform 0.1 > grinder.log 
			
			fi
			
			#Split paried-end reads into 2 files (each pair in a different file)
			${SCRIPT_DIR}/utils/split_grinder $PWD/grinder-reads.fa

			echo "Sort and extract true freq (required for computing the correlation)"
			tail -n +2 $PWD/grinder-ranks.txt | sort -k 2 | awk '{print $3}' > $PWD/grinder-ranks.justFreq.txt
				
			echo -e "Done simulating observed reads\n\n"
		else
			echo "Directory ${Path2Dir}/${i}/${g}/ : Ok!"
			cd ${Path2Dir}/${i}/${g}/
		fi
		
		
		currDate=`date '+%m.%d-%T'`
		echo "Current Date: $currDate"

	
		
		#############################	
		#	Run SimReg Version 23	#
		#############################
		echo "Running Latest SimReg version"
		
		if [ ! -d "./SimReg" ]; then
			mkdir ./SimReg
		fi
		
		cd ./SimReg
		
		if [ ! -f $PWD/bowtie_OBS_k60v3best.sam ]; 
		then	
			rm -rf $PWD/bowtieOBS.log
			
			echo "SAM File NOT Found! Map observed reads..."
			bowtie --best -v 3 -k 60 -p 36 --chunkmbs 128 $PWD/../../$i -f -1 $PWD/../grinder-reads_pair_1.fa -2 $PWD/../grinder-reads_pair_2.fa -I 180 -X 420 -S $PWD/bowtie_OBS_k60v3best.sam > bowtieOBS.log
		fi
		
		~/MCReg_v23/run-mcreg.sh -G $PWD/../../${i}.gtf -S $PWD/bowtie_OBS_k60v3best.sam -C $PWD/../../SimReg/components/ -F $PWD/../../${i}.fa -R $PWD/../../$i -1 $PWD/../grinder-reads_pair_1.fa -2 $PWD/../grinder-reads_pair_2.fa -m 300 -d 30 -l 100 > mcreg.log
	
		cp $PWD/MCReg.iso.estimates ../../../results~/${currDate}-${i}-${g}-MCReg.v23.iso.estimates
	
		#Compute Correlation
		awk '{print $2}' ./MCReg.iso.estimates > ./mcreg.iso.estimates_noNames
	
		if [ ! -f "../grinder-ranks.justFreq.txt" ]; then
			tail -n +2 ../grinder-ranks.txt | sort | awk '{print $3}' > ../grinder-ranks.justFreq.txt
		fi
		
		SimRegCorrel=`python ${SCRIPT_DIR}/scripts/correl.py ../grinder-ranks.justFreq.txt ./mcreg.iso.estimates_noNames`
		echo "SimRegCorrel = $SimRegCorrel"
		
		
		
		#############################	
		#	Run SimReg Version 21	#
		#############################
		if [ ! -d "../SimReg_v21" ]; then
			mkdir ../SimReg_v21
		fi
		cd ../SimReg_v21
	
		~/MCReg_v21/run-mcreg.sh -G $PWD/../../${i}.gtf -S $PWD/../SimReg/bowtie_OBS_k60v3best.sam -F $PWD/../../${i}.fa -R $PWD/../../$i -1 $PWD/../grinder-reads_pair_1.fa -2 $PWD/../grinder-reads_pair_2.fa -m 300 -d 30 -l 100 > mcreg_v21.log
		wait
	
		cp $PWD/MCReg.iso.estimates ../../../results~/${currDate}-${i}-${g}-MCReg.v21.iso.estimates
	
		#Compute Correlation
		awk '{print $2}' ./MCReg.iso.estimates > ./mcreg.iso.estimates_noNames
	
		SimRegCorrel_v21=`python ${SCRIPT_DIR}/scripts/correl.py ../grinder-ranks.justFreq.txt ./mcreg.iso.estimates_noNames`
		echo "SimRegCorrel_v21 = $SimRegCorrel_v21"
		
		
		
		#################	
		#	Run RSEM	#
		#################
		
		echo "## Running RSEM ... "
		
		if [ ! -d "../RSEM" ]; then
			echo "Directory ../RSEM: Missing!"
			mkdir ../RSEM
			cd ../RSEM
			
			echo "Prepare RSEM refernece in $PWD"
			
			/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-prepare-reference ${Path2Dir}/${i}/${i}.fa $PWD/$i > rsem.index.log
		else
			cd ../RSEM
		fi
		
		#/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --no-qualities $PWD/../grinder-reads_pair_1.fa $PWD/../grinder-reads_pair_2.fa $PWD/${i} $PWD/rsem > rsem.log
	
		#use our alignment
		/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --sam $PWD/../SimReg/bowtie_OBS_k60v3best.sam --paired-end $PWD/${i} $PWD/rsem > rsem.log
		
		cp $PWD/rsem.isoforms.results ../../../results~/${currDate}-${i}-${g}-rsem.isoforms.results
	
		#Compute Correlation
		#tail -n +2 shows all lines of report from the second line onwards:
		tail -n +2 $PWD/rsem.isoforms.results | sort | awk '{print $6}' > ./rsem.isoforms.results.justFreq
		RSEMCorrel=`python ${SCRIPT_DIR}/scripts/correl.py ../grinder-ranks.justFreq.txt ./rsem.isoforms.results.justFreq`
		echo "RSEMCorrel = $RSEMCorrel"
	

		#################	
		#	Run IsoEM	#
		#################
		if [ ! -d "../IsoEM" ]; then
			echo "Directory ../IsoEM does not exist!"
			mkdir ../IsoEM
		fi
		
		cd ../IsoEM
		
		if [ ! -f $PWD/bowtie_OBS_k60v3best.sam.genome ]; then
			/home/sahar/isoem-1.1.2/bin/convert-iso-to-genome-coords $PWD/../SimReg/bowtie_OBS_k60v3best.sam $PWD/../../${i}.gtf $PWD/bowtie_OBS_k60v3best.sam.genome --ignore-pairing
		fi
		
		#run isoem given the frag len distr.
		isoem -G $PWD/../../${i}.gtf -m 300 -d 30 $PWD/bowtie_OBS_k60v3best.sam.genome -o IsoEM
	
		cp $PWD/IsoEM.iso_estimates ../../../results~/${currDate}-${i}-${g}-IsoEM.iso_estimates
	
		#IsoEM Correlation:
		sort ./IsoEM.iso_estimates | awk '{print $2}' > ./IsoEM.iso_estimates.justFreq
		IsoEMCorrel=`python ${SCRIPT_DIR}/scripts/correl.py ../grinder-ranks.justFreq.txt ./IsoEM.iso_estimates.justFreq`
		echo "IsoEMCorrel = $IsoEMCorrel"
	
		#################
		
		echo -e "${currDate}\t${i}\t${g}\t${SimRegCorrel}\t${SimRegCorrel_v21}\t${RSEMCorrel}\t${IsoEMCorrel}" >> $correlFile
	
		echo -e "\nDone for observed reads: $g \n\n"
		
	done
	
	echo -e "Done test case: $i\n\n"
	
done
