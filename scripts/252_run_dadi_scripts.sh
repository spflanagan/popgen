#!/bin/bash

# This uses GNU parallel
# Simply run this script - there is no need for nohup or &
# It may take a while to complete though
# Note: a number of messages from GNU parallel will pop up

cd "${0%/*}" # move to location of script


############ SET THESE PARAMS ############
pops=('FLLG' 'FLCC' 'ALFW' 'ALST' 'LAFW' 'TXFW' 'TXCC')
models=('IM2NG') # 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')
# 'SI' 'IM' 'AM' 'SC' 'SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 
rangeX=2
rangeY=10


############ CREATE POP COMBOS TO RUN ############
combos=()
for ((i=0; i<(${#pops[@]}-1); ++i)); do
	for ((j=(i+1); j<${#pops[@]}; ++j)); do
		combos+=("${pops[$i]}_${pops[$j]}")
	done
done


############ RUN THE ANALYSIS ############
# The scripts take rangeX and rangeY and pass them on to the python file
for ((x=${rangeX}; x<${rangeY}; x++)); do
	y=$(( x+1 ))
	for ((i=0; i<(${#combos[@]}); ++i)); do
		for ((mod=0; mod<(${#models[@]}); ++mod)); do
			outfile="../fwsw_results/dadi_results/${combos[$i]}/${combos[$i]}_${x}.${models[$mod]}.optimized.txt"
			if [[ ! -f $outfile ]]; then
				echo "sem -j -4 ./dadi_scripts/${combos[$i]}_${models[$mod]}.sh ${x} ${y}"
				#sem -j -4 --bg ./dadi_scripts/${combos[$i]}_${models[$mod]}.sh ${x} ${y}
			else
				echo "${outfile} already exists"
			fi 
		done
	done
done


