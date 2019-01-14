#!/bin/bash

#run from popgen/fwsw_results
RUN100=false
RUN75=false
RUN50=false
RUN95=true
RUN90=true
RUN85=true

if [ "$RUN100" = true ]; then
	populations -b 2 -P ./stacks -M ../fwsw_sub_map.txt -t 5 -p 7 -r 1.0 --fstats --fstats --bootstrap_fst --plink --vcf --genepop
	mv ./stacks/batch_2.fst* ./stacks/populations_subset100
	mv ./stacks/batch_2.phistats* ./stacks/populations_subset100
	mv ./stacks/batch_2.plink* ./stacks/populations_subset100
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset100
	mv ./stacks/batch_2.vcf ./stacks/populations_subset100
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset100
	mv ./stacks/batch_2.genepop ./stacks/populations_subset100
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset100
fi
wait 

if [ "$RUN75" = true ]; then
	populations -b 2 -P ./stacks -M ../fwsw_sub_map.txt -t 5 -p 7 -r 0.75 --fstats --fstats --bootstrap_fst --plink --vcf --genepop
	mv ./stacks/batch_2.fst* ./stacks/populations_subset75
	mv ./stacks/batch_2.phistats* ./stacks/populations_subset75
	mv ./stacks/batch_2.plink* ./stacks/populations_subset75
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset75
	mv ./stacks/batch_2.vcf ./stacks/populations_subset75
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset75
	mv ./stacks/batch_2.genepop ./stacks/populations_subset75
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset75
fi
wait

if [ "$RUN50" = true ]; then
	populations -b 2 -P ./stacks -M ../fwsw_sub_map.txt -t 5 -p 7 -r 0.5 --fstats --fstats --bootstrap_fst --plink --vcf --genepop
	mv ./stacks/batch_2.fst* ./stacks/populations_subset50
	mv ./stacks/batch_2.phistats* ./stacks/populations_subset50
	mv ./stacks/batch_2.plink* ./stacks/populations_subset50
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset50
	mv ./stacks/batch_2.vcf ./stacks/populations_subset50
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset50
	mv ./stacks/batch_2.genepop ./stacks/populations_subset50
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset50
fi
wait

if [ "$RUN95" = true ]; then
	populations -b 2 -P ./stacks -M ../fwsw_sub_map.txt -t 5 -p 7 -r 0.95 --fstats --fstats --bootstrap_fst --plink --vcf --genepop
	mv ./stacks/batch_2.fst* ./stacks/populations_subset95
	mv ./stacks/batch_2.phistats* ./stacks/populations_subset95
	mv ./stacks/batch_2.plink* ./stacks/populations_subset95
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset95
	mv ./stacks/batch_2.vcf ./stacks/populations_subset95
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset95
	mv ./stacks/batch_2.genepop ./stacks/populations_subset95
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset95
fi
wait

if [ "$RUN90" = true ]; then 
	populations -b 2 -P ./stacks -M ../fwsw_sub_map.txt -t 5 -p 7 -r 0.90 --fstats --fstats --bootstrap_fst --plink --vcf --genepop
	mv ./stacks/batch_2.fst* ./stacks/populations_subset90
	mv ./stacks/batch_2.phistats* ./stacks/populations_subset90
	mv ./stacks/batch_2.plink* ./stacks/populations_subset90
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset90
	mv ./stacks/batch_2.vcf ./stacks/populations_subset90
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset90
	mv ./stacks/batch_2.genepop ./stacks/populations_subset90
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset90
fi
wait

if [ "$RUN85" = true ]; then
	populations -b 2 -P ./stacks -M ../fwsw_sub_map.txt -t 5 -p 7 -r 0.85 --fstats --fstats --bootstrap_fst --plink --vcf --genepop
	mv ./stacks/batch_2.fst* ./stacks/populations_subset85
	mv ./stacks/batch_2.phistats* ./stacks/populations_subset85
	mv ./stacks/batch_2.plink* ./stacks/populations_subset85
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset85
	mv ./stacks/batch_2.vcf ./stacks/populations_subset85
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset85
	mv ./stacks/batch_2.genepop ./stacks/populations_subset85
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset85
fi

