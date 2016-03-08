#!/bin/bash
DIRECTORY="/home/sarah/sf_ubuntushare/popgen/sw_results/faststructure/"
#echo "logistic k = 2"
#python ~/Programs/fastStructure/structure.py -K 2 \
#	--input=${DIRECTORY}subset.structure.recode \
#		--output=${DIRECTORY}res.k2.output_log --prior=logistic --full \
#		--format=str --seed=100
#echo "logistic k = 3"
#python ~/Programs/fastStructure/structure.py -K 3 \
#	--input=/${DIRECTORY}subset.structure.recode \
#		--output=${DIRECTORY}res.k3.output_log --prior=logistic \
#		--full --format=str --seed=100
#echo "logistic k = 4"
#python ~/Programs/fastStructure/structure.py -K 4 \
#	--input=${DIRECTORY}subset.structure.recode \
#		--output=${DIRECTORY}res.k4.output_log --prior=logistic --full \
#		--format=str --seed=100
#echo "logistic k = 5"
#python ~/Programs/fastStructure/structure.py -K 5 \
#	--input=${DIRECTORY}subset.structure.recode \
#		--output=${DIRECTORY}res.k5.output_log --prior=logistic \
#		--full --format=str --seed=100

echo "running distruct"
python ~/Programs/fastStructure/distruct.py -K 2 --input=${DIRECTORY}res.k2.output_log --output=${DIRECTORY}res.k2_distruct.svg

python ~/Programs/fastStructure/distruct.py -K 3 --input=${DIRECTORY}res.k3.output_log --output=${DIRECTORY}res.k3_distruct.svg

python ~/Programs/fastStructure/distruct.py -K 4 --input=${DIRECTORY}res.k4.output_log --output=${DIRECTORY}res.k4_distruct.svg

python ~/Programs/fastStructure/distruct.py -K 5 --input=${DIRECTORY}res.k5.output_log --output=${DIRECTORY}res.k2_distruct.svg


