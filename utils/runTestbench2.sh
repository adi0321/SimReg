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

#Variable Declaration
fil_min=180
fil_max=420
mean=300
deviation=30

cd $Path2Dir

correlFile=${Path2Dir}/Correlations.txt~

if [ -d "${Path2Dir}/results~" ];then
	rm -rf ${Path2Dir}/results~
fi

mkdir ${Path2Dir}/results~
echo -e "Run\tDate-Time\tTestCaseName\tObservedReads\tSimRegBest2\tRSEM\tIsoEM" > ${Path2Dir}/Correlations.txt~

testbench=`ls -B`
#	-B, --ignore-backups
#              do not list implied entries ending with ~

echo -e "Testbenches: \n$testbench"

cd ${Path2Dir}/results~

for t in $testbench 
do
	
	echo -e "\nCurrent Testbench = $t"
	GTF_File=${Path2Dir}/${t}/${t}.gtf
	FA_File=${Path2Dir}/${t}/${t}.fa
	REF_File=${Path2Dir}/${t}/${t}
	
	filename=$(basename "$GTF_File")
	extension="${filename##*.}"
	filename="${filename%.*}"
	RSEM_Ref=${Path2Dir}/${t}/rsemIndex/$filename
	
	mkdir $t
	cd $t
	
	#Check FA File
	if [ ! -f ${Path2Dir}/${t}/${t}.fa ]; then
		echo "FASTA File: Missing!"
		echo -e "\tExtracting FASTA File (reference sequences)..."
		#Extract the reference fa sequences
		${SCRIPT_DIR}/utils/make_refFile.sh ${Path2Dir}/${t}/${t}.gtf
		
		echo -e "\tCreating Bowtie indexes..."
		bowtie-build ${Path2Dir}/${t}/${i}.fa $i > bowtie-build.log
		echo "Done"
	else
		echo "FASTA File: Ok!"
	fi

	#for each testbench run 10 times
	for i in {0..9} 
	do
		mkdir "run$i"
		cd run$i
		echo -e "\nRun $i: Simulating reads using Grinder..."
		
		command -v grinder >/dev/null 2>&1 || { echo >&2 "
		ERROR: grinder NOT Found. 
		Install Grinder Simulator: http://sourceforge.net/projects/biogrinder/
		Aborting..."; exit 7; }

		echo "Generate Random Transcripts True Frequency"
		${SCRIPT_DIR}/generateRandTrFreq.sh $GTF_File 0
		
		simReads="grinderErrorFree	grinderMutUniform0.1"
		#echo "simReads = $simReads"

		for s in $simReads
		do

			#echo "s = $s"
			
			if [ ! -d $s ]; then
				mkdir $s
			fi
			
			cd $s
			
			if [ $s == "grinderErrorFree" ]; then
				echo -e "\tSimulating error FREE reads using Grinder ..."
				grinder -reference_file $FA_File -coverage_fold 100 -insert_dist $mean uniform $deviation -abundance_file ../${t}.tr.freq > grinder.log
			
			else
				echo -e "\tSimulating reads using mutation distribution: uniform 0.1 ... "
				grinder -reference_file $FA_File -coverage_fold 100 -insert_dist $mean uniform $deviation -abundance_file ../${t}.tr.freq -mutation_dist uniform 0.1 > grinder.log 
			fi
		
			#Split grinder pair end reads file into m1 and m2
			echo -e "\nSplitting Grinder files into <m1> and <m2> ... " 
			${SCRIPT_DIR}/split_grinder $PWD/grinder-reads.fa

			echo "Sort and extract true freq (required for computing the correlation)"
			tail -n +2 $PWD/grinder-ranks.txt | sort -k 2 | awk '{print $3}' > $PWD/grinder-ranks.justFreq.txt
			
			pair1_file=$PWD/grinder-reads_pair_1.fa
			pair2_file=$PWD/grinder-reads_pair_2.fa
	
			#fasta files
			input_option="-f"
			
			currDate=`date '+%m.%d-%T'`
			echo "Current Date: $currDate"

			echo "Map Observed Reads"
			bowtie --best -v 3 -k 60 -p 36 --chunkmbs 128 $REF_File $input_option -1 $pair1_file -2 $pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_k60v3best.sam
	
			#############################	
			#	Run SimRegBest2			#
			#############################
			echo "Running Latest SimReg version"
	
			mkdir ./SimReg
			cd ./SimReg
	
			${SCRIPT_DIR}/../run-mcreg.sh -G $GTF_File -S $PWD/../bowtie_k60v3best.sam -F $FA_File -R $REF_File -1 $pair1_file -2 $pair2_file -m $mean -d $deviation > simreg.log

			#Compute Correlation
			awk '{print $2}' ./MCReg.iso.estimates > ./mcreg.iso.estimates_noNames

			if [ ! -f "../grinder-ranks.justFreq.txt" ]; then
				tail -n +2 ../grinder-ranks.txt | sort | awk '{print $3}' > ../grinder-ranks.justFreq.txt
			fi
	
			SimRegCorrel=`python ${SCRIPT_DIR}/../scripts/correl.py ../grinder-ranks.justFreq.txt ./mcreg.iso.estimates_noNames`
			echo "SimRegCorrel = $SimRegCorrel"
	
	
			#################	
			#	Run RSEM	#
			#################
	
			echo -e "\n\n## Running RSEM ... "
			mkdir ../RSEM
			cd ../RSEM
		
			if [ ! -d "${Path2Dir}/${t}/rsemIndex" ];then
				echo "Prepare RSEM refernece in $PWD"
				mkdir ${Path2Dir}/${t}/rsemIndex/
				/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-prepare-reference $FA_File ${Path2Dir}/${t}/rsemIndex/$filename > ${Path2Dir}/${t}/rsemIndex/rsem.index.log
			else
				echo "RSEM indexes: OK"
			fi
			
			#/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --paired-end --no-qualities $PWD/../grinder-reads_pair_1.fa $PWD/../grinder-reads_pair_2.fa $PWD/${i} $PWD/rsem > rsem.log

			#use SimReg Bowtie alignment
			/home/code/IMPORT/RSEM/rsem-1.2.8/rsem-calculate-expression --sam $PWD/../bowtie_k60v3best.sam --paired-end $RSEM_Ref $PWD/rsem > rsem.log

			#Compute Correlation
			#tail -n +2 shows all lines of report from the second line onwards:
			tail -n +2 $PWD/rsem.isoforms.results | sort | awk '{print $6}' > ./rsem.isoforms.results.justFreq
			RSEMCorrel=`python ${SCRIPT_DIR}/../scripts/correl.py ../grinder-ranks.justFreq.txt ./rsem.isoforms.results.justFreq`
			echo "RSEMCorrel = $RSEMCorrel"

			#################	
			#	Run IsoEM	#
			#################						
			echo -e "\n\n## Running IsoEM ... "
			mkdir ../IsoEM
			cd ../IsoEM
			
			if [ ! -f $PWD/bowtie_k60v3best.sam.genome ]; then
				echo "Convert to genome coordinates"
				/home/sahar/isoem-1.1.2/bin/convert-iso-to-genome-coords $PWD/../bowtie_k60v3best.sam $GTF_File $PWD/../bowtie_k60v3best.sam.genome --ignore-pairing
			else
				echo "IsoEM: Genome coordinates: OK"
			fi
			
			
			#run isoem given the frag len distr.
			isoem -G $GTF_File -m $mean -d $deviation $PWD/../bowtie_k60v3best.sam.genome -o IsoEM

			#IsoEM Correlation:
			sort ./IsoEM.iso_estimates | awk '{print $2}' > ./IsoEM.iso_estimates.justFreq
			IsoEMCorrel=`python ${SCRIPT_DIR}/../scripts/correl.py ../grinder-ranks.justFreq.txt ./IsoEM.iso_estimates.justFreq`
			echo "IsoEMCorrel = $IsoEMCorrel"

			#################
	
			echo -e "${i}\t${currDate}\t${t}\t${s}\t${SimRegCorrel}\t${RSEMCorrel}\t${IsoEMCorrel}" >> $correlFile

			echo -e "\nDone simulated reads: $s \n\n"
			##Move out of directory $s
			cd ../..
		done
	
		cd ..
		echo -e "Done run: $i\n\n"
	done
	
	cd ..
	echo -e "Done test bench: $t\n\n"
done
