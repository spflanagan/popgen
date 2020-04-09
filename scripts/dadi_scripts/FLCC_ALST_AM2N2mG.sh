#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLCC' 60 'ALST' 70 $rangeX $rangeY AM2N2mG 2>&1 > FLCC_ALST/FLCC-ALST_AM2N2mG_${rangeX}_${rangeY}.log
