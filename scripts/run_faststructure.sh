#!/bin/bash

#run from popgen
for i in {1..12}
do
	echo "Starting K = $i"
	python ~/Programs/fastStructure/structure.py -K $i \
	--input=/home/sarah/sf_ubuntushare/popgen/sw_results/faststructure/subset.structure.recode \
		--output=/home/sarah/sf_ubuntushare/popgen/sw_results/faststructure/pruned_out_simple --full --format=str --seed=100

	
done
python ~/Programs/fastStructure/chooseK.py --input=sw_results/faststructure/pruned_out_simple 
  #  Model complexity that maximizes marginal likelihood = 2
 #   Model components used to explain structure in data = 5

