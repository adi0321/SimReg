#!/bin/bash

if [ $# -lt 1 ]
then
	echo -e "Usage: \n"
	echo "argv[1] - file with transcript freq"
	exit 1
fi

trFreq=$1
trNames=$PWD/trNames4Duplicate.txt
dupTrFreq=$PWD/freq4DuplicateTr.txt
#duplicate transcript frequencies

awk '{print $1}' $trFreq | sort | uniq -d > $trNames

tr=`cat $trNames`

for t in $tr
do

	#pattern=`${t}	0`
	#echo $pattern
	#sed '/${t}"\t"/d' $trFreq
	grep $t $trFreq | awk '{if($2!=0)print $0}' >> $dupTrFreq
done

#Now remove those transcripts from the freq file
for t in $tr
do
	#echo "t=$t"
	sed -i -e '/'$t'/d' $trFreq
done

#Append files
cat $trFreq $dupTrFreq >> ${trFreq}.temp

sort ${trFreq}.temp > $trFreq

rm -rf $dupTrFreq
rm -rf $trNames
rm -rf ${trFreq}.temp

#while read line
#do
	
#	echo $line
	
#done < $PWD/lines2BRemoved.txt



