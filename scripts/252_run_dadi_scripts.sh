#!/bin/bash

# This uses GNU parallel
# Run this script in the background and save output to a log to avoid disruptions
# It may take a while to complete though
# Note: task IDs start from 0

cd "${0%/*}" # move to location of script

taskStart=${1} 
taskEnd=${2}

############ SET THESE PARAMS ############
pops=('FLLG' 'FLCC' 'ALFW' 'ALST' 'LAFW' 'TXFW' 'TXCC')
models=('SI' 'IM' 'AM' 'SC' 'SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 'IM2NG' 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')

rangeX=5
rangeY=10

outdir="../fwsw_results/dadi_results/"
outend=".optimized.txt"


############ CREATE POP COMBOS TO RUN ############
if [ "$ARG1" == "nofile" ]; then
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
				outfile="${outdir}${combos[$i]}/${combos[$i]}_${x}.${models[$mod]}${outend}"
				if [[ ! -f $outfile ]]; then
					echo "sem -j -4 ./dadi_scripts/${combos[$i]}_${models[$mod]}.sh ${x} ${y}"
					#sem -j -4 --bg ./dadi_scripts/${combos[$i]}_${models[$mod]}.sh ${x} ${y}
				else
					echo "${outfile} already exists"
				fi 
			done
		done
	done
else
	file="../fwsw_results/$ARG1"
	pops=()
	mod=()
	while read -r value1 value2
	do
		pops+=($value1)
		mod+=($value2)
	done <"$file"
	
	for ((x=${rangeX}; x<${rangeY}; x++)); do
		y=$(( x+1 ))
		for ((i=0; i<(${#pops[@]}); ++i)); do
			outfile="${outdir}${pops[$i]}/${pops[$i]}_${x}.${mod[$i]}${outend}"
			if [[ ! -f $outfile ]]; then
				echo "sem -j -4 ./dadi_scripts/${pops[$i]}_${mod[$i]}.sh ${x} ${y}"
				sem -j -4 --bg ./dadi_scripts/${pops[$i]}_${mod[$i]}.sh ${x} ${y}
			else
				echo "${outfile} already exists"
			fi 	
		done
	done


fi


