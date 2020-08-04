#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py ALFW_ALST.dadi.snps 'ALFW' 30 'ALST' 30 $rangeX $rangeY IM2m 2>&1 > ALFW_ALST/ALFW_ALST_IM2m_${rangeX}_${rangeY}.log
