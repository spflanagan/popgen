#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'TXCC' 61 $rangeX $rangeY AM2NG 2>&1 > FLLG_TXCC/FLLG-TXCC_AM2NG_${rangeX}_${rangeY}.log
