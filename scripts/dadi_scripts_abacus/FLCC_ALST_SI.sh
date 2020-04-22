#!/bin/bash
cd "${0%/*}"
cd /share/data/people/spf50/fwsw_results/dadi_results/
rangeX=$1
rangeY=$2
/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py fwsw75.dadi.snps 'FLCC' 60 'ALST' 70 $rangeX $rangeY SI 2>&1 > FLCC_ALST/FLCC-ALST_SI_${rangeX}_${rangeY}.log
