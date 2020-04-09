#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLCC' 60 'TXCC' 61 $rangeX $rangeY IM2m 2>&1 > FLCC_TXCC/FLCC-TXCC_IM2m_${rangeX}_${rangeY}.log
