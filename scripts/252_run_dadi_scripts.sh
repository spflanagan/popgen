#!/bin/bash

# This uses GNU parallel
# Run this script in the background and save output to a log to avoid disruptions
# It may take a while to complete though
# Note: task IDs start from 0

cd "${0%/*}" # move to location of script

taskStart=${1} 
taskEnd=${2}

############ SET THESE PARAMS ############
fw_pops=('TXFW' 'LAFW' 'ALFW' 'FLLG')
sw_pops=('TXCC' 'ALST' 'ALST' 'FLCC')
models=('SI' 'IM' 'AM' 'SC' 'SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 'IM2NG' 'IM2mG' 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')
 
rangeX=1
rangeY=10

outdir="../fwsw_results/dadi_results/"
outend=".optimized.txt"


############ CREATE POP COMBOS TO RUN ############
tasks=()
for ((i=0; i<(${#fw_pops[@]}); ++i)); do
	for ((mod=0; mod<(${#models[@]}); ++mod)); do		
		tasks+=("${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}")
	done
done

len=${#tasks[@]}

############ RUN WITH SANITY CHECKS! ############
if [[ "$len" -gt "$taskEnd" ]]; then
	for ((x=${taskStart}; x<=${taskEnd}; x++)); do
		
		outfile="${outdir}${tasks[$x]}${outend}"
		if [[ ! -f $outfile ]]; then
			echo "sem -j -4 ./dadi_scripts/${tasks[$x]}.sh ${rangeX} ${rangeY}"
			sem -j -4 --bg ./dadi_scripts/${tasks[$x]}.sh ${rangeX} ${rangeY}
		else
			echo "${outfile} already exists"
		fi 	

	done
else
	echo "Ending taskID (${taskEnd}) greater than number of tasks (${len}): change the range and start qsub from 1" 1>&2
fi

	


