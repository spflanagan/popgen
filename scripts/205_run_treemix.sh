#!/bin/bash
# Date: 9 August 2017
# Author: Sarah P. Flanagan
# updated: 2017-12-18

ROOT="FLPB"
PREFIX="p4_"

cd ~/sf_ubuntushare/popgen/fwsw_results/treemix

treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -o ${PREFIX}k100b${ROOT}r
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 1 -o ${PREFIX}k100b${ROOT}rm1
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 2 -o ${PREFIX}k100b${ROOT}rm2
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 3 -o ${PREFIX}k100b${ROOT}rm3
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 4 -o ${PREFIX}k100b${ROOT}rm4
treemix -i ${PREFIX}treemix.gz -k 100 -bootstrap -root ${ROOT} -se -m 5 -o ${PREFIX}k100b${ROOT}rm5


threepop -i ${PREFIX}treemix.gz -k 100 >> ${PREFIX}threepop.txt
fourpop -i ${PREFIX}treemix.gz -k 100 >> ${PREFIX}fourpop.txt