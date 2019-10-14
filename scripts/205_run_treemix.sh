#!/bin/bash
# Date: 9 August 2017
# Author: Sarah P. Flanagan
# updated: 2019-15-10

ROOT="FLPB"
PREFIX="fwsw_"
VCF="../filter_rad_20191014@1654/14_filtered/radiator_data_20191014@1710.vcf"

cd "${0%/*}" # move to location of script
cd ../fwsw_results/treemix

Rscript ../../R/vcf2treemix.R ${VCF} "poplist" "fwsw_treemix" "../../../gwscaR"
gzip -c fwsw_treemix > fwsw_treemix.gz

# run treemix
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -se -o ${PREFIX}k100b #without the root
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -o ${PREFIX}k100b${ROOT}r
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 1 -o ${PREFIX}k100b${ROOT}rm1
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 2 -o ${PREFIX}k100b${ROOT}rm2
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 3 -o ${PREFIX}k100b${ROOT}rm3
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 4 -o ${PREFIX}k100b${ROOT}rm4
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 5 -o ${PREFIX}k100b${ROOT}rm5

# tests for treeness
threepop -i ${PREFIX}treemix.gz -k 100 >> ${PREFIX}threepop.txt
fourpop -i ${PREFIX}treemix.gz -k 100 >> ${PREFIX}fourpop.txt
