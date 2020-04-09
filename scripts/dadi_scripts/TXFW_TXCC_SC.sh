#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'TXFW' 46 'TXCC' 61 $rangeX $rangeY SC 2>&1 > TXFW_TXCC/TXFW-TXCC_SC_${rangeX}_${rangeY}.log
