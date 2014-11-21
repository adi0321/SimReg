After you installed all 3 dependecies, then you may try the following test case(s) to test SimReg:
(The C++ source files have been compiled on Ubuntu server so the program may work as-is). Try it first and if it doesn't work please report the bugs/or suggestions on how ti fix the bug to adrian.caciula@gmail.com

1. First change directory to 1gene_21tr (this small test case consists of a single gene with 21 transcripts)

2. Then you may choose any of the simulated reads: either error free reads or simulated reads with errors: mutations 0.1 uniform

Let's say you choose "grinderMutUniform0.1":

Step 1:
cd ./1gene_21tr/grinderMutUniform0.1/SimReg/

Step 2: Run SimReg:
$PWD/../../../../run-mcreg.sh -S $PWD/bowtie_OBS_k60v3best.sam -G $PWD/../../1gene_21tr.gtf -F $PWD/../../1gene_21tr.fa -R $PWD/../../1gene_21tr -m 300 -d 30

(use absolute path names for files)

The results will be available in MCReg.iso.estimates


