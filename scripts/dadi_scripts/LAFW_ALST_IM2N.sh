#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py LAFW_ALST.dadi.snps 'LAFW' 30 'ALST' 30 $rangeX $rangeY IM2N 2>&1 > LAFW_ALST/LAFW_ALST_IM2N_${rangeX}_${rangeY}.log
