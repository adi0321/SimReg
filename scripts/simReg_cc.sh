#!/bin/bash

startCmpTime=`date +%s`

#This is only for one component
readonly DEBUG=0
		
		if [ $DEBUG -ne 0 ]; then
			for var in "$@"
			do
				echo "$var"
			done
		fi

#echo -e "\nRunning: $0 ... "

usage() { echo "
Usage: 
   $0 [options]* [-s <SCRIPT_DIR>] [-C <Path to components>] [-c <components>] [-G <GTF_File>] [-F <FA_File>]

Mandatory options:   
	-s <SCRIPT_DIR>
	-c <int> component ID
	-G <GTF_File> File with known genes and isoforms in GTF format (Full path to file)
	-F <FA_File> Reference FASTA File with transcript sequences
	-C <Path> - Path to connected componentes (from observed reads)
	-m <int> Fragment length mean                   
	-d <int> Fragment length standard deviation
	
	Optional arguments:
   -l <int> Read length [ default: 100 ]
   -r <int> Number of reads per fragment [ default: 2 ]
   -g Simulator for generating Monte Carlo reads [ default: generate-reads (IsoEM-Simulator)]
   " 1>&2; exit 7; }

if (($# == 0)); then
        usage
fi

#If your var is just a flag, without any additional argument, just leave the var, without the ":" following
while getopts ":s:c:C:G:F:l:r:m:d:t:g" o; do
    case "${o}" in
		s)
            SCRIPT_DIR=${OPTARG}
            ;;
        C)
			CC_Path=${OPTARG}
			;;
		c)
			cID=${OPTARG}
			;;
		G)
            GTF_File=${OPTARG}
            ;;
        F)
            FA_File=${OPTARG}
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
		m)
            mean=${OPTARG}
            ;;
        d)
			deviation=${OPTARG}
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

	if [ $DEBUG -ne 0 ]; then
		echo "GTF_File = $GTF_File"
		echo "FA_File = $FA_File"
		echo "read_length = $read_length"
		echo "rpf = $rpf"
		echo "CC_Path = $CC_Path"
		echo "simulator=$simulator"
		echo "precision = $precision"
		#exit 7
	fi
	
errs=0
declare -A option=(
    [SCRIPT_DIR]="-s"
	[CC_Path]="-C"
	[cID]="-c"
	[GTF_File]="-G"
	[FA_File]="-F"
)
for var in "${!option[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo -e "ERROR: specify a value for $var with ${option[$var]}"
        ((errs++))
    fi
done
((errs > 0)) && usage

#old
#components=`ls -d $CC_Path/*`
#new 04.03.2014
#components=`ls -d ${CC_Path}*`

#cID stands for component id (remember that this iterates only through components with more than 1 transcript)
#echo "Component ID passed to simReg_cc.sh"
#echo "$cID"


#for c in $components
#do
	
	CC_Path=${CC_Path}/$cID
	
	#Extract directory name from its path
	#cID=`basename $c`
	
	#echo "cID = $cID"
	#echo "CC_Path = $CC_Path"
	
	if [ -d "./$cID" ]; then
		# Control will enter here if $DIRECTORY exists.
		cd ./$cID
	else
		mkdir ./$cID
		cd ./$cID
	fi

	obsTrNames="${CC_Path}/obsTr.txt"
	
	if [ ! -f $obsTrNames ];then
		echo "Error in simReg.sh: File $obsTrNames not found!"
		exit 7
	fi

	tr_names=`cat $obsTrNames`
	
		if [ $DEBUG -ne 0 ]; then
			echo "$cID"
			echo "CC_Path = $CC_Path"
			pwd
			echo "obsTrNames = $obsTrNames"
			echo "tr_names = $tr_names"
			#exit 7
		fi
	
	
	echo "Create the gtf and fa file for each transcript in current component ${cID}"
	start_time=`date +%s`
	for tran in $tr_names
	do
		#echo "Transcript name: $tran"
		
		grep $tran $GTF_File >> ${cID}.gtf
		
		#grep -A 1 $tran $FA_File >> ${cID}.fa
		
		#Remember that the shell will do shell interpolation before the program gets to AWK, 
		#so ${tran} will be replaced by the value. To be more safe, you might want to use
		#awk '/'$tran'/{print "Hello"}' $FA_File | head
		awk '/'$tran'/{p=1;print;next} p&&/^>/{p=0};p' $FA_File >> ${cID}.fa 
	done
	end_time=`date +%s`
	echo Done: `expr $end_time - $start_time` s.
	
	#Check to make sure that both gtf and fa have the same number of transcripts
	#number of transcripts in the fa file
	#nTFa=`awk {} ${cID}.fa`
	#echo "nTFa=$nTFa"
	
	#Get total number of transcripts from the gtf
	tnt=`${SCRIPT_DIR}/utils/get_isoforms.sh ${cID}.gtf | wc -l`
	echo "tnt=$tnt"

		#if [[ $tnt -gt 50 ]]; then
		#	echo "Total Number of transcripts in component $cID is $tnt" >> ../../largeComponents.txt
			#exit 7
		#fi
	
	
	#echo -e "\nPreparing simulation.properties file for IsoEM simulator ..."
	start_time=`date +%s`
	
	
	#Compute average transcript frequency 
	avg_tf=`echo "1 / $tnt" | bc -l`
		
	#Get number of transcripts in each gene
	genes_names=`${SCRIPT_DIR}/utils/get_genes.sh ${cID}.gtf`
		
	touch $PWD/nt_gene.txt
	
	for gene in $genes_names
	do
		nt=`grep -w $gene ${cID}.gtf | ${SCRIPT_DIR}/utils/get_isoforms.sh | wc -l`
		echo "$gene $nt">>$PWD/nt_gene.txt
	done
		
	G_Freq_File=${cID}.cluster_concentrations
	
		if [ -s $G_Freq_File ];then
			rm -rf $G_Freq_File
		fi
	
	write_flag="false"

	nt_gene=`cat $PWD/nt_gene.txt`
	#nt_gene is a two column file - 1st col is the name of the gene while the second one is the number of tr in the gene
	for gene in $nt_gene
	do	
			if [[ $write_flag == "false" ]]; then
				gene_name=$gene
				write_flag="true"
			else
				gene_freq=`echo "$gene * $avg_tf" | bc -l`	
				echo -e "${gene_name}\t${gene_freq}" >> $G_Freq_File
				write_flag="false"
			fi
	done
		
		#echo -e "\nCompute number of Monte Carlo reads: " 
		#get total_exon_length:
		total_exon_length=`awk '{s+=($5-$4)+1} END {print s}' ${cID}.gtf` 

		#echo "total_exon_length=$total_exon_length"
		
		
		
		if [[ $simulator == "generate-reads" ]]; then
	
			#Calculate number of reads for coverage 1000
			#ofreads =(coverage x total_exon_length) / (reads_per_fragment * read_length)     (paired-end) 
			#ofreads=(1000 x $total_exon_length) / 200 	#200 because they are paired-end reads

			mcreads=`echo "(1000 * $total_exon_length) / ($rpf * $read_length)" | bc`
			#just to be consistent with Grinder --- we'll divide by $read_length only ???
			#mcreads=`echo "(1000 * $total_exon_length) / $read_length" | bc`
	
				if [ $DEBUG -ne 0 ]; then
					echo "Number of transcripts in each gene: $nt_gene"
					echo "Total number of transcripts: $tnt"
					echo "Average transcript frequnecy: $avg_tf"
					echo "File with genes frequencies: $G_Freq_File"
					echo "FA File with isoform sequences: ${cID}.fa"
					echo "Total exon length: $total_exon_length"
					echo "Total number of Monte Carlo reads to be generated: $mcreads"
					echo "simulator = $simulator"
					#exit 7
				fi
	
			if [ -s $PWD/simulation.properties ];then
				rm -rf $PWD/simulation.properties
			fi
	
			#Step 7: Create simulation.properties				
			echo "gtfFile=${cID}.gtf">>$PWD/simulation.properties
			echo "isoformSequencesFile=${cID}.fa">>$PWD/simulation.properties
			echo "clusterDistribution = customWeights,$G_Freq_File">>$PWD/simulation.properties
			echo "fragmentLengthDistribution=normal,${mean},${deviation}">>$PWD/simulation.properties
			echo "fragmentStartingPositionDistribution=uniform">>$PWD/simulation.properties
			echo "isoformDistribution=uniform">>$PWD/simulation.properties
			echo "numberOfReads=$mcreads">>$PWD/simulation.properties
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
			
		elif [[ $simulator == "grinder" ]]; then

			echo -e "\nGenerating Monte Carlo Reads using Grinder..."
			
			command -v grinder >/dev/null 2>&1 || { echo >&2 "
			ERROR: grinder NOT Found. 
			Install Grinder Simulator: http://sourceforge.net/projects/biogrinder/
			Aborting..."; exit 7; }
			
			grinder -reference_file ${cID}.fa -coverage_fold 1000 -insert_dist ${mean} uniform ${deviation}
			wait
			
			#Step 2.2: Split grinder file into m1 and m2
			echo -e "\nSplitting Grinder files into <m1> and <m2> ... " 
			${SCRIPT_DIR}/utils/split_grinder $PWD/grinder-reads.fa
	
			mc_pair1_file=$PWD/grinder-reads_pair_1.fa
			mc_pair2_file=$PWD/grinder-reads_pair_2.fa
		
			#fasta files
			input_option="-f"
		
			refFile="grinder-reads.fa"
		else
			#else use SimReg Simulator
			#echo -e "\nGenerating Monte Carlo Reads using SimReg Reads Simulator..."
			${SCRIPT_DIR}/utils/sim-reads $PWD/${cID}.fa $mean $read_length #> sim-reads.log
			wait

			mc_pair1_file=$PWD/simreg-reads-pair1.fa
			mc_pair2_file=$PWD/simreg-reads-pair2.fa
			
			#fasta files
			input_option="-f"
			
			#both paires have same reference (so we can use any)
			refFile=$mc_pair1_file
		fi
		
		#Step 9: Map reads to transcriptome using bowtie
		
	
			
		#Step 10: Make bowtie index
		#echo -e "\nBuilding Bowtie indexes for ${cID}.fa ...\n"
			
		start_time=`date +%s`
		bowtie-build ${cID}.fa $cID > /dev/null 2>&1
		wait
		end_time=`date +%s`
		#echo Done: execution time was `expr $end_time - $start_time` s.
		
		#echo -e "\nMapping Monte Carlo reads using Bowtie ...\n"
		#Fragment insert length range
		fil_min=$(($mean-(3*$deviation)))
		fil_max=$(($mean+(3*$deviation)))
		
		start_time=`date +%s`
		#bowtie -v 0 -k 30 -p 12 $cID $input_option -1 $mc_pair1_file -2 $mc_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_MC_60multiAligns.sam
		bowtie -v 0 -k 60 -p 4 --chunkmbs 128 $cID $input_option -1 $mc_pair1_file -2 $mc_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_MC_60multiAligns.sam > /dev/null 2>&1
		wait
		mcSAM_File=$PWD/bowtie_MC_60multiAligns.sam
		
				if [ $DEBUG -ne 0 ]; then
					echo "Reference fasta file: ${cID}.fa"
					echo -e "\tInput Option is: $input_option"
					echo -e "\tFragment Insert Length Min: $fil_min"
					echo -e "\tFragment Insert Length Max: $fil_max"

					echo "in $PWD"
				fi
		
		#end_time=`date +%s`
		#echo Done Mapping Monte Carlo Reads! Execution time: `expr $end_time - $start_time` s.
		#echo ""
		
		
#if the number of transcripts is 2 then we need to divide by the ratio

#TO DO(4/24/2014): this if else should be between formula and solver
#the formula program should also output the next observed
if [[ $tnt -eq 0 ]]; then

			echo "This component has only $tnt transcripts"
			#compute frequency for two transcripts
			echo "CC_Path = ${CC_Path}"

			#divide the counts from each single class and report the the frequencies
			#This also computes the Total number of reads
			mkdir ./temp/
			
			${SCRIPT_DIR}/lib/compute-freq-42t ${cID}.gtf $refFile $PWD/bowtie_MC_60multiAligns.sam ${CC_Path}/obsRCcounts.txt
			
			#Compute transcripts lengths:
			awk '($3=="exon") || ($2=="exon")  {ID=$12;L[ID]+=$5-$4+1} END{for(i in L){print i, L[i]}}' ${cID}.gtf | sort | awk '{print $1, $2}' > tr_length.txt
			
				if [ $DEBUG -ne 0 ]; then
					cat tr_length.txt
					#exit 7
				fi
	
else
		#Else run the solver
		
		echo "Computing Simulated Read Classes and D Values(compute_sRC_d)"
		
		start_time=`date +%s`
		${SCRIPT_DIR}/lib/compute_sRC_d ${cID}.gtf $refFile $PWD/bowtie_MC_60multiAligns.sam ${CC_Path}/obsRCcounts.txt >> sim-reads.log
		wait

		end_time=`date +%s`
		#echo Done Computing Simulated Read Classes and D values! Execution Time: `expr $end_time - $start_time` s.
		#echo ""

		
		#To DO (2014.02.13 - Add also observed classes - right now they are discarded - 
			##but since the results are good we don't worry too much about this)
		#echo "Prepare file for solver"
		
				#May 1, 2014
				#File obsRCsize.txt must be recalculated 
				rm -rf $PWD/obsRCsize.txt
				#Based on the new size
				#(This will be done in MCReg_CC_v2)
		
		#start_time=`date +%s`
		${SCRIPT_DIR}/lib/MCReg_CC_v2 ${cID}_d_values.txt ${cID}_read_classes2.txt ${CC_Path}/obsRCcounts.txt 1 1 > mcreg_cc.log
		wait
		#end_time=`date +%s`
		#echo Done MCReg_CC_v2 execution time: `expr $end_time - $start_time` s.
		#echo ""
		
		
		#Step 12: Run qp solver and compute the transcripts frequency
		#Total number of mapped observed reads
		total_obs_reads=`awk '{ total+=$NF} END{print total}' ${CC_Path}/obsRCcounts.txt`
		d_values_file=`ls d_value*`
		estimFile=$PWD/mcreg.iso.estimates
		
			if [ $DEBUG -ne 0 ]; then
				echo "total_obs_reads = $total_obs_reads"
				echo "d_values_file is: $d_values_file"
				echo "Estimates File = $estimFile"
			fi

		#start_time=`date +%s`
		${SCRIPT_DIR}/scripts/MCReg.sh $PWD/${cID}.gtf $PWD/${d_values_file} $PWD/temp/0_o_values.txt $total_obs_reads $estimFile ${CC_Path}/obsRCsize.txt #> MCReg.sh.log

		#end_time=`date +%s`
		#User regular observed file (the one that will be created in temp) ?
		#echo Done MCReg.sh in `expr $end_time - $start_time` s.

fi
		##Compute the correlation for transcripts
		#awk '{print $2}' ./mcreg.iso.estimates > ./mcreg.iso.estimates_noNames.txt
		#python ${SCRIPT_DIR}/scripts/correl.py ../../true.txt ./mcreg.iso.estimates_noNames.txt >> ./tr_correl.txt
			
		##Compute Squared Deviation for transcripts
		#paste ../../true.txt ./mcreg.iso.estimates_noNames.txt | awk '{ print (($1-$2)**2); }' >> ./tr_sqDev.txt
		#awk '{ sum += $1 } END { print sum }' ./tr_sqDev.txt >> ./tr_sqDevSum.txt
		
##########################################################################################
		##Now take the estimates from $out_file and simulate Observed Reads again
		#option_precision=1
#precision is given as an argument and represents the number of iterations
	if [[ $precision != "0" ]]; then
		
		echo -e "\n\n##############################################################"
		echo "Start Precision Estimation (Tuning-up current estimations):"
		echo " "
		
		#ToDo: define number of iterations
		#Number of Iterations:
		#nOfIt=
		
		echo "Bash version ${BASH_VERSION}..."
		echo " "
		

		
		echo "Initialization Step"
		echo " "
		
		echo "a) initialize aimed=observed"
		cat $PWD/temp/0_o_values.txt > $PWD/temp/aimed_reads.txt
		
		deltaSqAll=$((0))
		echo "delta = $deltaSqAll"
		

		if [ -f $PWD/bowtie_OBS_multiAligns.sam ];then
			rm -rf $PWD/all.estimates.txt
		fi
		
		echo -e "Run I:" >> ./all.estimates.txt
		cat $PWD/mcreg.iso.estimates >> ./all.estimates.txt
		echo -e "\n" >> ./all.estimates.txt
		
			mkdir ./precision
			cd ./precision

		#for i in 1 2 3
		for i in $(seq 1 $precision)
		do
			
			if [ -f $PWD/bowtie_OBS_multiAligns.sam ]; then
				rm -rf $PWD/bowtie_OBS_multiAligns.sam
				rm -rf grinder*
				rm -rf ./0/
				rm -rf ./temp/
				rm -rf $PWD/mcreg.iso.estimates
			fi
		
			echo -e "\n Sim $i: Generate Simulated Reads using F' as abundance file"
			
			#simulator used for simulated reads
			simulatorSR="Formula"
			
		if [[ $simulatorSR != "Formula" ]]; then
			#if its not formula then it's either Grinder simulator or IsoEM simulator (generate-reads)
			if [[ $simulatorSR == "Grinder" ]]; then
				
				echo "Using Grinder for simulated reads"
				
				command -v grinder >/dev/null 2>&1 || { echo >&2 "
				ERROR: grinder NOT Found. 
				Install Grinder Simulator: http://sourceforge.net/projects/biogrinder/
				Aborting..."; exit 1; }
			
				grinder -reference_file ../${cID}.fa -coverage_fold 100 -insert_dist 300 uniform 30 -abundance_file $estimFile 
				#>& $PWD/grinder.${i}.log.txt
				wait
			
				echo "#Split paried-end reads into 2 files (each pair in a different file)"
				${SCRIPT_DIR}/utils/split_grinder $PWD/grinder-reads.fa
		
				bowtie -k 60 -p 12 ../${cID} -f -1 grinder-reads_pair_1.fa -2 grinder-reads_pair_2.fa -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_multiAligns.sam
			
			else
				
				echo "Using IsoEM-simulator for simulated reads"
				
				nreads=`echo "(100 * $total_exon_length) / ($rpf * $read_length)" | bc`
				
				if [ -s $PWD/simulation.properties ];then
					rm -rf $PWD/simulation.properties
				fi

				#Create simulation.properties				
				echo "gtfFile=../${cID}.gtf">>$PWD/simulation.properties
				echo "isoformSequencesFile=../${cID}.fa">>$PWD/simulation.properties
				echo "clusterDistribution = customWeights,$estimFile">>$PWD/simulation.properties
				echo "fragmentLengthDistribution=normal,${mean},${deviation}">>$PWD/simulation.properties
				echo "fragmentStartingPositionDistribution=uniform">>$PWD/simulation.properties
				echo "isoformDistribution=uniform">>$PWD/simulation.properties
				echo "numberOfReads=$nreads">>$PWD/simulation.properties
				echo "readLength=$read_length">>$PWD/simulation.properties
				echo "randomNumberGeneratorSeed=123">>$PWD/simulation.properties
				echo "readsPerFragment=$rpf">>$PWD/simulation.properties
				echo "firstReadOrigin=random">>$PWD/simulation.properties
			
				rm -rf ./nr=*
				
				generate-reads

				pair1_file=`ls *paired_1.fastq`
				pair2_file=`ls *paired_2.fastq`
				
				bowtie -k 60 -p 12 ../${cID} -q -1 $pair1_file -2 $pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_multiAligns.sam
				
			fi
			

			${SCRIPT_DIR}/lib/compute_obsRC_CC $PWD/bowtie_OBS_multiAligns.sam
			#(no need to compute again everything ... we just need the counts)
			#but it's fine for now
		else
			
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			echo "Else compute the counts using formula"										# Update: 4/29/2014					#@@@@
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			
			echo -e "\nComputing NEXT observed counts..."
			
				if [ $DEBUG -ne 0 ]; then
					echo "mcSAM_File = $mcSAM_File"
					echo "estimFile = $estimFile"
					echo "mc_pair1_file = $mc_pair1_file"
					echo "total_obs_reads = $total_obs_reads"
					#exit 7
				fi
					
			start_time=`date +%s`
			#7/5/2014 - In order to keep the same classes we'll need to load them in memory
			${SCRIPT_DIR}/lib/compute-nextCounts2 $mcSAM_File $estimFile $mc_pair1_file $total_obs_reads ${CC_Path}/obsRCcounts.txt
			#argv[1] Mapped reads (sam file)
			#argv[2] Transcripts Freq.
			#argv[3] File with pair 1 to extract references
			end_time=`date +%s`
			
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			echo Done Computing NEXT observed counts! Execution Time: `expr $end_time - $start_time`s.								#@@@@
			echo ""																													#@@@@
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		fi

			##Use same d values and read classes
			##Obs 2014.07.05: But this time we may have different read classes (so why not use current read classes and d values?)
				##But we need the same size for observed bc we have to do the subtraction few steps below
				
			${SCRIPT_DIR}/lib/MCReg_CC_v2 ../${cID}_d_values.txt ../${cID}_read_classes2.txt $PWD/0/obsRCcounts.txt 1 1
			#${SCRIPT_DIR}/lib/MCReg_CC_v2 $PWD/0/0_d_values.txt $PWD/0/0_read_classes2.txt $PWD/0/obsRCcounts.txt 1 1			
			#>& mcreg_cc.log.txt

			#Run qp solver and compute the transcripts frequency
			#Total number of mapped observed reads
			total_obs_reads2=`awk '{ total+=$NF} END{print total}' $PWD/0/obsRCcounts.txt`
			echo "total_obs_reads2 = $total_obs_reads2"
	
			d_values_file=`ls d_value*`
			
			echo "d_values_file is: $d_values_file"

			#This script computes o_values 
			#This d_value file also contain the new observed values (the new d values file was computed above by MCReg_CC_v2)

			#obsRCsize should be the same if we keep the same obs read classes since is just 1 / number of tr in the class
			${SCRIPT_DIR}/scripts/MCReg.sh $PWD/../${cID}.gtf $PWD/${d_values_file} $PWD/temp/0_o_values.txt $total_obs_reads2 $PWD/mcreg.iso.estimates ${CC_Path}/obsRCsize.txt
				#argv 5 - output file
			echo "Run ${i}: Done MCReg.sh"

			##Now compute delta and deltaSquare
			##Simulated Reads = ./temp/0_o_values.txt
			##Observed Reads = ../temp/0_o_values.txt
			
			#Delta is s-o (simulated - observed) -- DELTA is a vector
			
			paste ./temp/0_o_values.txt ../temp/0_o_values.txt | awk '{ print $1-$2; }' > DELTA_${i}.txt			
			paste ./temp/0_o_values.txt ../temp/0_o_values.txt | awk '{ print (($1-$2)**2); }' > deltaSq_${i}.txt
			
			cp ./temp/0_o_values.txt ./S_${i}.txt
			
			##Sum of current delta squares
			cur_deltaSq=`paste ./temp/0_o_values.txt ../temp/0_o_values.txt | awk '{ print (($1-$2)**2)}' | awk '{ sum += $1 } END { print sum }' `
			
			echo "Current sum of delta (squares) = $cur_deltaSq"

			
			#Unfortunately bc doesn't support scientific notation.
			#cur_deltaSq=${cur_deltaSq_temp/[eE]+*/*10^}
			cur_deltaSq=`awk '{ printf("%.10f", $1)}' <<< $cur_deltaSq`
			echo "Current sum of delta (squares) (without scientific notation) = $cur_deltaSq"		
			
			
			if [ $DEBUG -ne 0 ]; then
				python ${SCRIPT_DIR}/scripts/correl.py ./temp/0_o_values.txt ../temp/0_o_values.txt >> ../temp/reads_correl.txt
			fi
			
			
			#this file contains all delta squares
			awk '{ sum += $1 } END { print sum }' deltaSq_${i}.txt >> ../temp/deltaSq.sum.txt
			

			##########
			#2nd - Compute aimed reads which is prev_aimed - delta/2
			
			#Use this one here:
			#Old command: Aimed = Prev_Aimed - delta/2
			
			#Other Old Aimed = Prev_aimed - 1/delta Sq
			
			##New: 2014.02.22 - New Aimed = O - 1/sigma x sum (DELTA/delta)
			
			#Compute sigma ---there should be i and not i-1 in the sigma formula
			#sigma = `awk '{ sum += 1/($1) } END { print sum }' ../temp/deltaSq.sum.txt`
			#We do not need sigma anymore because we normalize anyway before we give it to the solver
			
			
			#Old command
			#There is a problem in processing the file and then overwriting the results in the same file
			paste ../temp/aimed_reads.txt ./DELTA_${i}.txt | awk '{ if (($1-($2/2)) < 0) { print 0 } else {print ($1-($2/2))} }' > ../temp/aimed_reads2.txt
			cp -rf ../temp/aimed_reads2.txt ../temp/aimed_reads.txt
			
			#New aimed: O - DELTA / delta
			
			
			#Compute A_{i+2}
				#Step 1: Compute the numerator: DELTA_{i}*A_{i+1} + DELTA_i+1 * A_i
				#echo "DELTA$(($i-1)).txt"
				
				#prod1
				#paste ./DELTA_${i}.txt ../temp/aimed_reads_$(($i+1)).txt | awk '{ print ($1*$2) }' > ./prod1.txt
				#prod2
				#paste ./DELTA_$(($i+1)).txt ../temp/aimed_reads_${i}.txt | awk '{ print ($1*$2) }' > ./prod2.txt
				#sum
				#paste ./prod1.txt ./prod2.txt | awk '{print ($1 + $2)}' > ./sum.txt
				
				#Step 2: Compute the denominator: S_i - S_{i+1}
				#paste ./S_${i}.txt ./S_$(($i+1)).txt | awk '{print ($1-$2)}' > ./denom.txt
				
				#Step 3: Compute A_i
				#paste ./sum.txt ./denom.txt | awk '{if($2==0) print "0"; else print ($1 / $2);}' > ../temp/aimed_reads_$((i+2)).txt
				
			
			#Normalize aimed reads before giving to solver
			#compute sum of all
			sum_of_all=`awk '{ sum += $1 } END { print sum }' ../temp/aimed_reads.txt`

			echo "Sum of all read classes freq. is $sum_of_all"
			#awk -v var=$sum_of_all '{print $1"\t"$2/var}' $out
			awk -v var=$sum_of_all '{print $1/var}' ../temp/aimed_reads.txt > ../temp/aimed_reads.norm.txt
			
			method=1
			if [[ $method -eq "1" ]];then
				echo "Method 1"
				${SCRIPT_DIR}/scripts/MCReg.sh $PWD/../${cID}.gtf $PWD/../${d_values_file} $PWD/../temp/aimed_reads.norm.txt $total_obs_reads2 $PWD/mcreg.iso.estimates ${CC_Path}/obsRCsize.txt
			else
				echo "Method 2"
				#this is the method from MCReg v5
				#Compute o2 which is o0 - delta/2
				paste ../temp/0_o_values.txt ./DELTA_${i}.txt | awk '{ if (($1-($2/2)) < 0) { print 0 } else {print ($1-($2/2))} }' > ./0_o_values2.txt
				${SCRIPT_DIR}/scripts/MCReg.sh $PWD/../${cID}.gtf $PWD/${d_values_file} $PWD/0_o_values2.txt $total_obs_reads2  $PWD/mcreg.iso.estimates ${CC_Path}/obsRCsize.txt
			fi

			#copy the estimates file one level up
			cat $PWD/mcreg.iso.estimates > ../mcreg.iso.estimates
			
			echo -e "Run ${i}:" >> ../all.estimates.txt
			cat $PWD/mcreg.iso.estimates >> ../all.estimates.txt
			echo -e "\n" >> ../all.estimates.txt


		done
		
		##Exit dir precision
		cd ..
	fi #end if precision
		




			##Compute the correlation for transcripts
			awk '{print $2}' ./mcreg.iso.estimates > ./mcreg.iso.estimates_noNames
			
			##Extract true transcript frequency for current component
			#for each transcript extract true freq from ref
			
			ref=${GTF_File%.*}
			#echo "ref=$ref"
			
			#Execute this only if true transcript frequencies ${ref}.tr.freq.norm exists
			#pwd: precision directory
			
			if [ -f ${ref}.tr.freq.norm ];then
		
				for tran in $tr_names
				do
					#echo "Transcript name: $tran"
		
					grep $tran ${ref}.tr.freq.norm >> ./true.tr.freq
				done

				#sort and extract only the values
				sort ./true.tr.freq | awk '{print $2}' > ./true.tr.freq.noNames

				current_correl=`python ${SCRIPT_DIR}/scripts/correl.py ./true.tr.freq.noNames ./mcreg.iso.estimates_noNames`
				echo "current_correl = $current_correl"
			
				echo -e "${cID}\t${current_correl}\t${tnt}" >> ../../comp_correl.txt
				#tnt holds total number of transcripts
	
			
				##Compute Squared Deviation for transcripts
				paste ./true.tr.freq.noNames ./mcreg.iso.estimates_noNames | awk '{ print (($1-$2)**2); }' >> ./tr_sqDev.txt
			
				#we need the sum of tr_sqDev
				awk '{ sum += $1 } END { print sum }' ./tr_sqDev.txt >> ./tr_sqDevSum.txt
		

				#cp $PWD/mcreg.0.iso.estimates $PWD/mcreg.iso.estimates
			fi
				##~~~~~~~2014.02.04~~~~~~~~~~~~~~~~~
		
		
		
	
		#~~~New code Dec.18,2013~~~~~~
		#array with transcripts lengths
		tr_lengths_array=(`cat ./tr_length.txt | awk '{print $2}' `)
			echo "tr_lengths_array=${tr_lengths_array[@]}"
		#tr_estim=`cat ${cID}.mcreg.iso.estimates`
		tr_estim=(`awk '{print $2}' mcreg.iso.estimates`)
			echo "tr_estim=${tr_estim[@]}"
	
	
		##Extract f from solver estimates and multiply with the corresponding transcript length
		i=0 #used for array index
		local_avg_tr_len=0 #Average transcripts lengths
		
		for tr_i in ${tr_estim[@]}
		do		
			#echo "tr_i=$tr_i"
			#echo "tr_lengths_array[$i] = ${tr_lengths_array[$i]}"
			#Multiply each f with the corresponding transcript length
			local_avg_tr_len=$(echo "$local_avg_tr_len + ${tr_lengths_array[$i]} * $tr_i" | bc -l)
			
			i=$[$i + 1]
		done
				
		#Number of reads divided by average local length
		no_reads_c=0 
		#Initialize Number of reads in the current component
		no_reads_files=`ls ./temp/*_no_reads.txt`
		
		for f in $no_reads_files
		do
			temp_reads_c=`cat $f`
			no_reads_c=$[$no_reads_c + $temp_reads_c]
		done
		
		echo "Number of reads in component ${cID} is: $no_reads_c"
		
		tr_names=`awk '{print $1}' mcreg.iso.estimates`
			#echo "tr_names=${tr_names}"
			#here there are two elements per line (that's why we'll use variable i in the loop below)
		
		i=0
		##Header Structure: "C.ID #Tr. #ObsRCmp\tTr.Name\tLclFreq.\tAvgTrLen"
		for t in $tr_names
		do
			#here we may need to report inside the component
			#so instead of ../../results.txt we say just results.txt
			echo -e "${cID}\t$tnt\t${no_reads_c}\t${t}\t${tr_estim[$i]}\t$local_avg_tr_len">>results.txt
			#tnt is the total number of transcripts in this component
			
			i=$[$i + 1]
		done
		
		#Remove Monte Carlo Reads and Mappings
		rm -rf nr*
		rm -rf bowtie_MC_60multiAligns.sam
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		endTotalTime=`date +%s`
		totalTimeSec=`expr $endTotalTime - $startCmpTime`
		echo "Done component ${cID} in $totalTimeSec seconds"
		
		cd ..
		
		#if [[ ${cID} == "37" ]]; then
		#	echo "cID = $cID"
		#fi
		
#done


cd ..