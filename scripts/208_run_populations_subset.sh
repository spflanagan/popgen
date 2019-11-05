#!/bin/bash

#run from popgen/fwsw_results
RUN100=false
RUN75=false
RUN50=false
RUN95=false
RUN90=false
RUN85=false
RUNALL=false
WHITELIST=false
WL75=true

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
wait

if [ "$RUNALL" = true ]; then
	populations -b 2 -P ./stacks -M ../fwsw_pops_map.txt -t 5 -p 16 -r 0.75 --fstats --fstats --plink --vcf --genepop --treemix --structure
	mv ./stacks/batch_2.fst* ./stacks/populations_all
	mv ./stacks/batch_2.phistats* ./stacks/populations_all
	mv ./stacks/batch_2.plink* ./stacks/populations_all
	mv ./stacks/batch_2.sumstats* ./stacks/populations_all
	mv ./stacks/batch_2.vcf ./stacks/populations_all
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_all
	mv ./stacks/batch_2.genepop ./stacks/populations_all
	#mv *batch_2*treemix* ./stacks/populations_all
	#mv *batch_2*structure* ./stacks/populations_all
	mv ./stacks/batch_2.populations.log ./stacks/populations_all
fi

if [ "$WHITELIST" = true ]; then
	populations -b 2 -P ./stacks -M ./stacks/strata.filtered.txt -t 5 -p 16 -r 0.75 -W ./stacks/whitelist.txt --fstats --plink --vcf --genepop --bootstrap_fst
	mv ./stacks/batch_2.fst* ./stacks/populations_whitelist
	mv ./stacks/batch_2.phistats* ./stacks/populations_whitelist
	mv ./stacks/batch_2.plink* ./stacks/populations_whitelist
	mv ./stacks/batch_2.sumstats* ./stacks/populations_whitelist
	mv ./stacks/batch_2.vcf ./stacks/populations_whitelist
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_whitelist
	mv ./stacks/batch_2.genepop ./stacks/populations_whitelist
	mv ./stacks/batch_2.populations.log ./stacks/populations_whitelist
fi

if [ "$WL75" = true ]; then
	# we're running the whitelisted sub 75 SNPs but with all populations
	populations -b 2 -P ./stacks -M ../fwsw_alt_map.txt -t 4 -W ./stacks/populations_subset75/pruned_snps.txt -p 1 -r 0.01 -R 0.01 --min-maf 0.001 --fstats --plink --vcf
	mv ./stacks/batch_2.fst* ./stacks/populations_subset75/all_pops
	mv ./stacks/batch_2.plink* ./stacks/populations_subset75/all_pops
	mv ./stacks/batch_2.sumstats* ./stacks/populations_subset75/all_pops
	mv ./stacks/batch_2.vcf ./stacks/populations_subset75/all_pops
	mv ./stacks/batch_2.haplotypes.tsv ./stacks/populations_subset75/all_pops
	mv ./stacks/batch_2.populations.log ./stacks/populations_subset75/all_pops
fi
