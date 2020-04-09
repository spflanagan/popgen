#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'ALST' 70 'TXFW' 46 $rangeX $rangeY SCG 2>&1 > ALST_TXFW/ALST-TXFW_SCG_${rangeX}_${rangeY}.log
