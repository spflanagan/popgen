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
# copy the phylip consense program to this directory
cp ~/Programs/phylip-3.697/exe/consense ./

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
	rm outtree; rm outfile
	cp unrooted/${PREFIX}_cat.tre intree
	# NOW RUN CONSENSE from phylip & convert with R -- saved as unrooted_consensus.newick
	# create plot
	#R -e "png("../../figs/treemix_unrooted_consensus.png",height=8,width=8,units="in",res=300); plot(ape::read.tree("unrooted/fwsw_boottree.txt"));dev.off()"
	#
	# run threepop and fourpop
	threepop -i ${PREFIX}_treemix.gz -k 100  >> ${PREFIX}_threepop.txt
	fourpop -i ${PREFIX}_treemix.gz -k 100 >> ${PREFIX}_fourpop.txt
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
	# create consensus using consense
	Rscript ../../R/202_combine_treemix.R rooted/${PREFIX}_${ROOT} rooted/${PREFIX}_${ROOT}_cat.tre 100 
	rm outtree; rm outfile
	cp rooted/${PREFIX}_${ROOT}_cat.tre intree # Then use consense & convert with R -- saved as rooted_consensus.newick
	
	# convert it in R
	#R -e "ape::write.tree(ape::read.tree('rooted/${PREFIX}_${ROOT}_boottree.newick'),'rooted/${PREFIX}_${ROOT}_consensus.tre',digits=1)"
fi


# run treemix with migration edges
if [ "$METHOD" == "migrations" ]; then
	echo "Running treemix with migration edges"
	[ ! -d "migrations" ] && mkdir -p migrations
	# maximum likelihood trees
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf unrooted_consensus.newick -se -o unrooted/${PREFIX}_ML_consensus 
	# with migration edges
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted_consensus.newick -root ${ROOT} -se -m 0 -o migrations/${PREFIX}_${ROOT}_m0
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted_consensus.newick -root ${ROOT} -se -m 1 -o migrations/${PREFIX}_${ROOT}_m1
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted_consensus.newick -root ${ROOT} -se -m 2 -o migrations/${PREFIX}_${ROOT}_m2
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted_consensus.newick -root ${ROOT} -se -m 3 -o migrations/${PREFIX}_${ROOT}_m3
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted_consensus.newick -root ${ROOT} -se -m 4 -o migrations/${PREFIX}_${ROOT}_m4
	treemix -i ${PREFIX}_treemix.gz -k 100 -tf rooted_consensus.newick -root ${ROOT} -se -m 5 -o migrations/${PREFIX}_${ROOT}_m5
	# bootstraps
	for (( i=1; i<=$bootreps; i++ ))
	do
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -tf rooted_consensus.newick -root ${ROOT} -se -m 0 -o migrations/${PREFIX}_${ROOT}_m0_${i}
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -tf rooted_consensus.newick -root ${ROOT} -se -m 1 -o migrations/${PREFIX}_${ROOT}_m1_${i}
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -tf rooted_consensus.newick -root ${ROOT} -se -m 2 -o migrations/${PREFIX}_${ROOT}_m2_${i}
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -tf rooted_consensus.newick -root ${ROOT} -se -m 3 -o migrations/${PREFIX}_${ROOT}_m3_${i}
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -tf rooted_consensus.newick -root ${ROOT} -se -m 4 -o migrations/${PREFIX}_${ROOT}_m4_${i}
		treemix -i ${PREFIX}_treemix.gz -k 100 -bootstrap -tf rooted_consensus.newick -root ${ROOT} -se -m 5 -o migrations/${PREFIX}_${ROOT}_m5_${i}
	done
	

fi

