#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'ALFW' 72 'LAFW' 72 $rangeX $rangeY IM2N 2>&1 > ALFW_LAFW/ALFW-LAFW_IM2N_${rangeX}_${rangeY}.log