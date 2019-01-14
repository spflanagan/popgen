#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR="../fwsw_results"
cd $DIR
cd $OUTDIR
#Running SweepFinder2

#first, convert vcf to sweepfinder

Rscript ../scripts/204_vcf2sf.R

#calculate genome-wide allele frequency spectrum
echo "Calculating genome-wide allele frequency spectrum"
~/Programs/SF2/SweepFinder2 -f SF2/GenomeWideSpectrum.txt SF2/SpectFile.txt

#then for each lg run SweepFinder

for i in SF2/*AF.txt
do
	filename=$(echo $i | sed -e 's/txt/out/g')
	echo "~/Programs/SF2/SweepFinder2 -l 50 ${i} SF2/SpectFile.txt ${filename}"
	~/Programs/SF2/SweepFinder2 -l 50 ${i} SF2/SpectFile.txt ${filename}
done