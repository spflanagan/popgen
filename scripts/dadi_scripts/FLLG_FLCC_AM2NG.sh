#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py FLLG_FLCC.dadi.snps 'FLLG' 30 'FLCC' 30 $rangeX $rangeY AM2NG 2>&1 > FLLG_FLCC/FLLG_FLCC_AM2NG_${rangeX}_${rangeY}.log
