### 8/28/2014 ###

cat /data1/adrian/MAQC/ENSG00000131095/ENSG00000131095.gtf > 2MAQC_Genes.gtf
cat /data1/adrian/MAQC/ENSG00000129244/E >> 2MAQC_Genes.gtf

cat /data1/adrian/MAQC/ENSG00000131095/ENSG00000131095.fa > 2MAQC_Genes.fa
cat /data1/adrian/MAQC/ENSG00000129244/ENSG00000129244.fa >> 2MAQC_Genes.fa

bowtie-build $PWD/2MAQC_Genes.fa 2MAQC_Genes

### Combine the alignment file
awk '{if($5!="*") print $0}' /data1/adrian/MAQC/ENSG00000129244/bowtie_OBS_k60.sam > bowtie_filtered.sam
awk '{if($5!="*") print $0}' /data1/adrian/MAQC/ENSG00000131095/bowtie_OBS_k60.sam >> bowtie_filtered.sam

sort -k 3,3 -k 4,4n bowtie_filtered.sam > bowtie_filtered.sam.sorted

~/SimRegBest/run-mcreg.sh -S $PWD/bowtie_filtered.sam -G $PWD/2MAQC_Genes.gtf -F $PWD/2MAQC_Genes.fa -R $PWD/2MAQC_Genes -1 /import3/MAQC/RNA-Seq/Illumina/Brain/SRX018974/SRR039628.fastq -2 /import3/MAQC/RNA-Seq/Illumina/Brain/SRX018975/SRR039629.fastq -m 200 -d 50 -l 50
## Sorted SAM doesn't work with SimReg


