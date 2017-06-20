#!/bin/bash
# Date: 19 June 2017

cd ~/sf_ubuntushare/popgen/fwsw_results/treemix

treemix -i fwsw.treemix.gz -k 100 -bootstrap -o fwsw.k100bFLLGr
treemix -i fwsw.treemix.gz -k 100 -bootstrap -m 1 -o fwsw.k100bFLLGrm1
treemix -i fwsw.treemix.gz -k 100 -bootstrap -m 2 -o fwsw.k100bFLLGrm2
treemix -i fwsw.treemix.gz -k 100 -bootstrap -m 3 -o fwsw.k100bFLLGrm3
treemix -i fwsw.treemix.gz -k 100 -bootstrap -m 4 -o fwsw.k100bFLLGrm4
treemix -i fwsw.treemix.gz -k 100 -bootstrap -m 5 -o fwsw.k100bFLLGrm5

threepop -i fwsw.treemix.gz -k 500
fourpop -i fwsw.treemix.gz