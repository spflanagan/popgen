#!/bin/bash

#run from scripts
cd ../sw_results/bayenv2
for filename in ./pruned_snps/*; do
	echo "Running bayenv2 on $filename"
	~/Programs/bayenv_2/bayenv2 -i $filename -m matrix.1.out \
		-e env_data_bayenv.txt -p 12 -k 100000 -n 5 -t -c -X -f \
		-r 628398 -o ./pruned_snps/$filename.freq
done
