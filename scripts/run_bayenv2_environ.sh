#!/bin/bash

#run from scripts
cd ../sw_results/environmental_assoc/new_bayenv/
for filename in ./pruned_snps/*; do
	echo "Running bayenv2 on $filename"
	~/Programs/bayenv2/bayenv2 -i $filename -m matrix.txt \
		-e standardized.env -p 12 -k 100000 -n 5 -t -c -X -f \
		-r 628398 -o $filename.freq
done
