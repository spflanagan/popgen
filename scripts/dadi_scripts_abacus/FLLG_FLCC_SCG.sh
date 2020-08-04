#!/bin/bash
cd "${0%/*}"
cd /share/data/people/spf50/fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py FLLG_FLCC.dadi.snps 'FLLG' 30 'FLCC' 30 $rangeX $rangeY SCG 2>&1 > FLLG_FLCC/FLLG_FLCC_SCG_${rangeX}_${rangeY}.log
