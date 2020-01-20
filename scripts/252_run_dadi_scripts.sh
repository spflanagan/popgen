#!/bin/bash

# This uses GNU parallel
# Simply run this script - there is no need for nohup or &
# Note: a number of messages from GNU parallel will pop up

cd "${0%/*}" # move to location of script


############ SET THESE PARAMS ############
pops=('FLLG' 'FLCC' 'ALFW' 'ALST' 'LAFW' 'TXFW' 'TXCC')  
models=('SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 'IM2NG' 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')
# 'SI' 'IM' 'AM' 'SC'
rangeX=2
rangeY=3


############ CREATE POP COMBOS TO RUN ############
combos=()
for ((i=0; i<(${#pops[@]}-1); ++i)); do
	for ((j=(i+1); j<${#pops[@]}; ++j)); do
		combos+=("${pops[$i]}_${pops[$j]}")
	done
done


############ RUN THE ANALYSIS ############
# The scripts take rangeX and rangeY and pass them on to the python file

for ((i=0; i<(${#combos[@]}); ++i)); do
	for ((mod=0; mod<(${#models[@]}); ++mod)); do
		sem -j -4 ./dadi_scripts/${combos[$i]}_${models[$mod]}.sh ${rangeX} ${rangeY}
	done
done


