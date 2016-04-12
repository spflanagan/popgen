#!/bin/bash

#run from scripts
cd ../sw_results/bayenv2/new_bayenv/
for filename in ./temp_var/*; do
	echo "Running bayenv2 on $filename"
	~/Programs/bayenv_2/bayenv2 -i $filename -m matrix.txt \
		-e wod_tempvar_data_std.txt -p 12 -k 100000 -n 1 -t -c -f \
		-r 628398 -o $filename.freq
		#	-e env_data_bayenv_std.txt -p 12 -k 100000 -n 5 -t -c -X -f \

done
