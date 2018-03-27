#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR="../fwsw_results/sweed/"
cd $DIR
cd $OUTDIR
#Running SweeD


#for each lg run SweeD-P in the background

for i in LG*.vcf
do
	filename=$(echo $i | sed -e 's/.vcf//g')
	echo "~/Programs/sweed-master/SweeD-P -name ${filename} -input ${i} -grid 50 -isfs FreqFiles.txt -checkpoint 3600" 
	~/Programs/sweed-master/SweeD-P -name ${filename} -input ${i} -grid 50 -isfs FreqFiles.txt -threads 4 &
done