#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR="../fwsw_results"
cd $DIR
cd $OUTDIR
#Running SweepFinder2

#first, convert vcf to sweepfinder

Rscript ../scripts/vcf2sf.R

#then for each lg run SweepFinder

for i in SF2/*AF.txt
do
	filename=$(echo $i | sed -e 's/txt/out/g')
	echo "~/Programs/SF2/SweepFinder2 -s 50 ${i} ${filename}"
	~/Programs/SF2/SweepFinder2 -s 50 ${i} ${filename}
done