#!/bin/bash
# Date: 9 August 2017

cd ~/sf_ubuntushare/popgen/fwsw_results/treemix

# treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLLG -se -o fwsw.k100bFLLGr
# treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLLG -se -m 1 -o fwsw.k100bFLLGrm1
# treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLLG -se -m 2 -o fwsw.k100bFLLGrm2
# treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLLG -se -m 3 -o fwsw.k100bFLLGrm3
# treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLLG -se -m 4 -o fwsw.k100bFLLGrm4
# treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLLG -se -m 5 -o fwsw.k100bFLLGrm5

# threepop -i fwsw.treemix.gz -k 100 >> fwsw.threepop.txt
# fourpop -i fwsw.treemix.gz -k 100 >> fwsw.fourpop.txt

treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLPB -se -o fwsw.k100bFLPBr
treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLPB -se -m 1 -o fwsw.k100bFLPBrm1
treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLPB -se -m 2 -o fwsw.k100bFLPBrm2
treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLPB -se -m 3 -o fwsw.k100bFLPBrm3
treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLPB -se -m 4 -o fwsw.k100bFLPBrm4
treemix -i fwsw.treemix.gz -k 100 -bootstrap -root FLPB -se -m 5 -o fwsw.k100bFLPBrm5
