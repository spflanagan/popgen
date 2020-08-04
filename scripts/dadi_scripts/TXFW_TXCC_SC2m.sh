#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py TXFW_TXCC.dadi.snps 'TXFW' 30 'TXCC' 30 $rangeX $rangeY SC2m 2>&1 > TXFW_TXCC/TXFW_TXCC_SC2m_${rangeX}_${rangeY}.log
