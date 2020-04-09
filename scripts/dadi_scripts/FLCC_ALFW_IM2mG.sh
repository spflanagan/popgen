#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLCC' 60 'ALFW' 72 $rangeX $rangeY IM2mG 2>&1 > FLCC_ALFW/FLCC-ALFW_IM2mG_${rangeX}_${rangeY}.log
