#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'TXFW' 46 $rangeX $rangeY IM2NG 2>&1 > FLLG_TXFW/FLLG-TXFW_IM2NG_${rangeX}_${rangeY}.log
