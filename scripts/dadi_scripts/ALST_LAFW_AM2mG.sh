#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'ALST' 70 'LAFW' 72 $rangeX $rangeY AM2mG 2>&1 > ALST_LAFW/ALST-LAFW_AM2mG_${rangeX}_${rangeY}.log
