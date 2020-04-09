#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'ALFW' 72 'ALST' 70 $rangeX $rangeY AM2m 2>&1 > ALFW_ALST/ALFW-ALST_AM2m_${rangeX}_${rangeY}.log
