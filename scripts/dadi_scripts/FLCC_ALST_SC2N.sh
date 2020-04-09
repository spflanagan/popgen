#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLCC' 60 'ALST' 70 $rangeX $rangeY SC2N 2>&1 > FLCC_ALST/FLCC-ALST_SC2N_${rangeX}_${rangeY}.log
