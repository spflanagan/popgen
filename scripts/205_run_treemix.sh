#!/bin/bash
# Date: 9 August 2017
# Author: Sarah P. Flanagan
# updated: 2020-06-02

ROOT="FLAB"
PREFIX="fwsw"
VCF="../filter_rad_20191014@1654/14_filtered/radiator_data_20191014@1710.vcf"
METHOD="migrations" # one of convert, unrooted, rooted, or migrations
bootreps=100

cd "${0%/*}" # move to location of script
cd ../fwsw_results/treemix
# convert vcf to input format
if [ "$METHOD" == "convert" ]; then
	Rscript ../../R/vcf2treemix.R ${VCF} "poplist" "fwsw_treemix" "../../../gwscaR"
	gzip -c fwsw_treemix > fwsw_treemix.gz
fi
# run treemix without a root and generate a consensus tree
if [ "$METHOD" == "unrooted" ]; then
	echo "Running treemix without a root"
	[ ! -d "unrooted" ] && mkdir -p unrooted
	# run bootstraps
	for (( i=1; i<=$bootreps; i++ ))
	do
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -o unrooted/${PREFIX}_${i} 
	done
	# create consensus
	Rscript ../../R/202_combine_treemix.R unrooted/${PREFIX} unrooted/${PREFIX}_cat.tre 100
	sumtrees.py --unrooted -o unrooted/${PREFIX}_boottree.txt -F 'newick' unrooted/${PREFIX}_cat.tre
fi


# run treemix with a specified root and generate a consensus tree

if [ "$METHOD" == "rooted" ]; then
	echo "Running treemix with root ${ROOT}"
	[ ! -d "rooted" ] && mkdir -p rooted
	# run bootstraps
	for (( i=1; i<=$bootreps; i++ ))
	do
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -o rooted/${PREFIX}_${ROOT}_${i} 
	done
	# create consensus
	Rscript ../../R/202_combine_treemix.R rooted/${PREFIX}_${ROOT} rooted/${PREFIX}_${ROOT}_cat.tre 100 
	sumtrees.py --rooted -o rooted/${PREFIX}_${ROOT}_boottree.txt -F 'newick' rooted/${PREFIX}_${ROOT}_cat.tre
fi

# run treemix with migration edges
if [ "$METHOD" == "migrations" ]; then
	echo "Running treemix with migration edges"
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted/${PREFIX}_${ROOT}_consensus.tre -root ${ROOT} -se -m 1 -o ${PREFIX}_${ROOT}_m1
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted/${PREFIX}_${ROOT}_consensus.tre -root ${ROOT} -se -m 2 -o ${PREFIX}_${ROOT}_m2
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted/${PREFIX}_${ROOT}_consensus.tre -root ${ROOT} -se -m 3 -o ${PREFIX}_${ROOT}_m3
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted/${PREFIX}_${ROOT}_consensus.tre -root ${ROOT} -se -m 4 -o ${PREFIX}_${ROOT}_m4
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted/${PREFIX}_${ROOT}_consensus.tre -root ${ROOT} -se -m 5 -o ${PREFIX}_${ROOT}_m5
	# tests for treeness
	threepop -i ${PREFIX}_treemix.gz -k 100  >> ${PREFIX}_threepop.txt
	fourpop -i ${PREFIX}_treemix.gz -k 100 >> ${PREFIX}_fourpop.txt

fi

