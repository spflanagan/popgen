#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'LAFW' 72 'TXFW' 46 $rangeX $rangeY AM2N2m 2>&1 > LAFW_TXFW/LAFW-TXFW_AM2N2m_${rangeX}_${rangeY}.log