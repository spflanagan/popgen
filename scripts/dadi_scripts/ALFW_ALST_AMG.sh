#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py ALFW_ALST.dadi.snps 'ALFW' 30 'ALST' 30 $rangeX $rangeY AMG 2>&1 > ALFW_ALST/ALFW_ALST_AMG_${rangeX}_${rangeY}.log
