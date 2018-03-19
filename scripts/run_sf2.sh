#!/bin/bash

cd ~/sf_ubuntushare/popgen/fwsw_results
#Running SweepFinder2

#first, convert vcf to sweepfinder

Rscript ../scripts/vcf2sf.R

#then for each lg run SweepFinder
cd ~/Programs/SF2/
for i in fwsw_results/sf2/*AF.txt
do
	./SweepFinder2 -s 50 ${i} ${i}.out.txt

done