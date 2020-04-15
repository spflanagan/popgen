#!/bin/bash
cd "${0%/*}"
cd ../../fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
python ../../scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'ALFW' 72 $rangeX $rangeY SC2N2m 2>&1 > FLLG_ALFW/FLLG-ALFW_SC2N2m_${rangeX}_${rangeY}.log