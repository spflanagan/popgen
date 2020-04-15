#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'ALST' 70 $rangeX $rangeY SI2N 2>&1 > FLLG_ALST/FLLG-ALST_SI2N_${rangeX}_${rangeY}.log