#!/bin/bash


#~/Programs/PCAdaptPackage/ped2pcadapt ./sw_results/stacks/populations/subset.ped ./sw_results/pcadapt/pcadaptgen
#echo "Running PCAdapt with:"
#for i in {1..12}
#do
#	echo "K=$i"
#	~/Programs/PCAdaptPackage/PCAdapt -i ./sw_results/pcadapt/pcadaptgen -K $i -o ./sw_results/pcadapt/pcadapt.$i > ./sw_results/pcadapt/pcadapt.$i.log

#done

#cp sw_results/pcadapt/pcadapt.*.stats ~/Programs/PCAdaptPackage/Rscripts

#choose k using R 
R
#source("~/Programs/PCAdapt/Rscripts/get_errors.R")
#look for where the min value first occurs. K=4.
#then run multiple times 

for i in {1..10}
do	
	echo "Run $i"
	~/Programs/PCAdaptPackage/PCAdapt -i ./sw_results/pcadapt/pcadaptgen -K 4 -o ./sw_results/pcadapt/pcadapt.4_$i
done
