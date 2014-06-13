#!/bin/bash
hostname=`hostname`

#This script computes the correlation for each component

echo -e "\nRun: $0 on $hostname ... "
if [ $# -lt 2 ]
then
	echo ""
	echo -e "Usage: \n"
	echo "argv[1] - Full path to components"
	echo "argv[2] - Full path to file with true frequencies"
	echo ""
	exit 1
fi

cmpPath=$1
trueFreq=$2

cd $cmpPath

components=`ls`

#cID stands for component id (remember that this iterates only through components with more than 1 transcript)
for cID in $components
do

	cd ./$cID
	
	#Extract true frequency for this component
	tr_names=`awk '{print $1}' mcreg.iso.estimates`
	echo "tr_names=${tr_names}"
	#here there are two elemnts per line (that's why we'll use variable i in the loop below)
		
	i=0

	for t in $tr_names
	do
			
			grep $t $$trueFreq >> ./trueFreq.txt
			
			i=$[$i + 1]
	done
	
	
done