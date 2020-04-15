#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'LAFW' 72 $rangeX $rangeY IM2mG 2>&1 > FLLG_LAFW/FLLG-LAFW_IM2mG_${rangeX}_${rangeY}.log