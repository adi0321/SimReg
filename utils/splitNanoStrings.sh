#!/bin/bash

if [[ $PWD == "/import1/GENCODE_v19" ]]
then
	echo "Running Script from $PWD"
else
	echo "You have to run the script from /import1/GENCODE_v19"
	exit 1
fi

awk '{

	outFile=$1
	for (i=3; i<=NF; i++){
		printf("%s ", $i)>>"./nanoString/"outFile".txt"
	}
		
}' /import1/GENCODE_v19/nanoCounts.txt

echo ""
echo "Done!"

