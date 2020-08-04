#!/bin/bash
cd "${0%/*}"
cd /share/data/people/spf50/fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py TXFW_TXCC.dadi.snps 'TXFW' 30 'TXCC' 30 $rangeX $rangeY SCG 2>&1 > TXFW_TXCC/TXFW_TXCC_SCG_${rangeX}_${rangeY}.log
