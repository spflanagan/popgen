#!/bin/bash

# This uses GNU parallel

cd "${0%/*}" # move to location of script
cd ../fwsw_results/dadi_results/ # move to dadi output location


############ SET THESE PARAMS ############
pops=('FLLG' 'FLCC' 'ALFW' 'ALST' 'LAFW' 'TXFW' 'TXCC')
projs=(70 60 72 70 72 46 61)

rangeX=1
rangeY=2

############ RUN THE ANALYSIS ############
# the python script takes these arguments in this order:
# snps_file, popi, prji, popj, prjj, x, y

for ((i=0; i<(${#pops[@]}-1); ++i)); do
	for ((j=(i+1); j<${#pops[@]}; ++j)); do
		sem -j -4 "python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps '${pops[$i]}' ${projs[$i]} '${pops[$j]}' ${projs[$j]} $rangeX $rangeY 2>&1 ${pops[$i]}-${pops[$j]}_${rangeX}.log"
	done
	sem --wait
done



