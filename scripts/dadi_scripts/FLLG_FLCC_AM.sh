#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'FLCC' 60 $rangeX $rangeY AM 2>&1 > FLLG_FLCC/FLLG-FLCC_AM_${rangeX}_${rangeY}.log
