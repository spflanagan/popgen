#!/bin/bash

# SET PARAMS
# need plink-style inputs to run admixture, I've got some in the admixture directory - include path otherwise
PLINK_NAME="fwsw_all_filt.bed" # bed is working, ped isn't
ADMIX_DIR="/home/sarah/Programs/admixture_linux-1.3.0/"
OUT_DIR="../fwsw_results/admixture/"
TEST_K=true


cd "${0%/*}" # move to dir of script
cd "${OUT_DIR}" # move to output directory

if [ "$TEST_K" = true ]; then

	for K in {1..16}
	do
		"${ADMIX_DIR}"admixture --cv ${PLINK_NAME} $K | tee log${K}.out
	done
	grep -h CV log*out > K_CVs.txt # this we will load into R to plot

fi

