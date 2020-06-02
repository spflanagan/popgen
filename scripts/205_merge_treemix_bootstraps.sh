#!/bin/bash


# Author: Sarah Flanagan 
# Date: 2 June 2020
# Name: 205_merge_treemix_bootstraps.py
# Usage: 	205_merge_treemix_bootstraps.sh IN THE DIRECTORY WHERE YOU WANT IT TO BE RUN
# Make sure sumtrees (from Dendropy) is installed

# Inspired by script from: 
# https://github.com/mgharvey/misc_python/blob/master/bin/TreeMix/treemix_tree_with_bootstraps.py
# Original Name: treemix_tree_with_bootreps.py 
# Original Author: Michael G. Harvey
# Original Date: 13 May 2013
# Dependencies:
# treemix (https://code.google.com/p/treemix/)
# sumtrees package in dendropy (http://pythonhosted.org/DendroPy/scripts/sumtrees.html) 
# updated sumtrees link (2 June 2020): https://dendropy.org/programs/sumtrees.html



############################
# Define these parameters! #
############################

out_dir="treemix/"
stems=( "fwsw_k100bFLPBrm1" "fwsw_k100bFLPBrm2" "fwsw_k100bFLPBrm3" "fwsw_k100bFLPBrm4" "fwsw_k100bFLPBrm5") #"fwsw_k100bFLPBr") #
bootreps=10



for stem in "${stems[@]}"
do
	# remove old ones
	rm ${out_dir}${stem}_cat_trees.tre
	# Use R to parse, reformat, and combine trees
	Rscript ../R/202_combine_treemix.R ${out_dir}${stem} ${out_dir}${stem}_cat_trees.tre 
	# run sumtrees
	sumtrees.py --rooted -t ${out_dir}${stem}.treeout.tre -o ${out_dir}${stem}_boottree.txt -F 'newick' ${out_dir}${stem}_cat_trees.tre
done
