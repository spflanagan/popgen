#!/bin/bash
cd "${0%/*}"
cd /share/data/people/spf50/fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'FLCC' 60 $rangeX $rangeY SC2mG 2>&1 > FLLG_FLCC/FLLG-FLCC_SC2mG_${rangeX}_${rangeY}.log
