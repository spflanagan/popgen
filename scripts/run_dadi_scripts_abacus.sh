#!/bin/bash

# This is adapted from 252_run_dadi_scripts.sh
# This script is for running on abacus using qsub
# Usage, where S is the starting ID and N is ending ID:
# qsub -t S-N -cwd -S /bin/bash run_dadi_scripts_abacus.sh
cd "${0%/*}" # move to location of script
taskID=$SGE_TASK_ID 

### --- SGE SETTINGS --- ###

#$ -M sarah.flanagan@canterbury.ac.nz
#$ -m abe
#$ -cwd
#$ -o /home/spf50/jobs/
#$ -e /home/spf50/jobs/

############ SET THESE PARAMS ############
pops=('FLLG' 'FLCC' 'ALFW' 'ALST' 'LAFW' 'TXFW' 'TXCC')   #
models=('SI' 'IM' 'AM' 'SC' 'SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 'IM2NG' 'IM2mG' 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')
 
rangeX=12
rangeY=13


############ CREATE POP COMBOS TO RUN ############
tasks=()
for ((i=0; i<(${#pops[@]}-1); ++i)); do
	for ((j=(i+1); j<${#pops[@]}; ++j)); do
		for ((mod=0; mod<(${#models[@]}); ++mod)); do		
			tasks+=("${pops[$i]}_${pops[$j]}_${models[$mod]}")
		done
	done
done

len=${#tasks[@]}

############ RUN WITH SANITY CHECKS! ############
if [[ "$len" -gt "$taskID" ]]; then
	./dadi_scripts_abacus/${tasks[$taskID]}.sh ${rangeX} ${rangeY}
else
	echo "taskID (${taskID}) greater than number of tasks (${len}): change the range and start qsub from 1" 1>&2
fi





