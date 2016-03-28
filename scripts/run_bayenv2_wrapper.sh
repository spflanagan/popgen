#!/bin/bash

cd ./sw_results/environmental_assoc/new_bayenv/
cp ~/Programs/bayenv2/calc_bf.sh .
cp ~/Programs/bayenv2/bayenv2 .

./calc_bf.sh all.SNPSFILE env_data_bayenv_std.txt matrix.txt 12 100000 5

rm calc_bf.sh
rm bayenv2
