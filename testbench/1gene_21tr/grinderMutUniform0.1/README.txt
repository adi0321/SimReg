Commands used for testing:

tail -n +2 $PWD/grinder-ranks.txt | sort | awk '{print $3}' > $PWD/grinder-ranks.txt.justFreq

### 11/18/2014
nohup nice ~/SimRegBest2/utils/runTestbench3.sh -d $PWD -S $PWD/SimReg/bowtie_OBS_k60v3best.sam -G $PWD/../1gene_21tr.gtf -F $PWD/../1gene_21tr.fa -R $PWD/../1gene_21tr -E $PWD/../rsemIndex/1gene_21tr -T $PWD/grinder-ranks.txt >& run.log




