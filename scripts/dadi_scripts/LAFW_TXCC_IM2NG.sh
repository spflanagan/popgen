#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'LAFW' 72 'TXCC' 61 $rangeX $rangeY IM2NG 2>&1 > LAFW_TXCC/LAFW-TXCC_IM2NG_${rangeX}_${rangeY}.log