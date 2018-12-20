#!/bin/bash

#run from popgen/fwsw_results
populations -P ./stacks -O ./stacks/populations_subset100 -M ../fwsw_sub_map.txt -t 6 -p 7 -r 1.0 --fstats --smooth_fstats --bootstrap_fst --plink --vcf --genepop
wait 
populations -P ./stacks -O ./stacks/populations_subset75 -M ../fwsw_sub_map.txt -t 6 -p 7 -r 0.75 --fstats --smooth_fstats --bootstrap_fst --plink --vcf --genepop
wait
populations -P ./stacks -O ./stacks/populations_subset50 -M ../fwsw_sub_map.txt -t 6 -p 7 -r 0.5 --fstats --smooth_fstats --bootstrap_fst --plink --vcf --genepop
