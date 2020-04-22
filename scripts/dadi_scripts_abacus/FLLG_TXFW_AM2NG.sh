#!/bin/bash
cd "${0%/*}"
cd /share/data/people/spf50/fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLLG' 70 'TXFW' 46 $rangeX $rangeY AM2NG 2>&1 > FLLG_TXFW/FLLG-TXFW_AM2NG_${rangeX}_${rangeY}.log
