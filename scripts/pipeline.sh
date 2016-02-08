#!/bin/bash
#execute this file from /Documents/Sarahs/Popgen/


#bowtie alignment
#to build a bowtie index:
#./bowtie2-build ./allpaths_cms1.scaff.fa scovelli_allpaths
#cd samples/
#sh pop2align_newgenome_12Aug2014.sh 2> bowtie.out.txt


#ref_map.pl (stacks program). Database interaction turned off for now.

#cd ../
echo "running ref_map.sh"
sh ref_map.sh


#now run rxstacks for population-based corrections

echo "running rxstacks"
rxstacks -b 1 -P ./stacks/ -o ./cor_stacks/ --conf_lim 0.25 --prune_haplo \ 
--model_type bounded --bound_high 0.1 --lnl_lim -8.0 --lnl_dist -t 6 --verbose


#rebuild the catalog
echo "re-running cstacks"
sh pop2_cstacks.sh

#rematch against the catalog
echo "re-running sstacks"
sh sstacks.sh

#calculate population statistics
echo "running populations on corrected data"
populations -b 1 -P ./cor_stacks/ -t 6 -M ./popgen2_map.txt -s -r 75 \
--fstats -k --structure --genomic

#load data into database to compare results
echo "loading both datasets into the database and web interface"
mysql -e "create database cor_rxstacks_radtags" 
mysql -e "create database unc_rxstacks_radtags" 
mysql cor_rxstacks_radtags < /usr/local/share/stacks/sql/stacks.sql 
mysql unc_rxstacks_radtags < /usr/local/share/stacks/sql/stacks.sql 
load_radtags.pl -D unc_rxstacks_radtags -b 1 -p ./stacks/ -B -e "Uncorrected data" -c
index_radtags.pl -D unc_rxstacks_radtags -c -t
load_radtags.pl -D cor_rxstacks_radtags -b 1 -p ./cor_stacks/ -B -e "Corrected data" -c
index_radtags.pl -D cor_rxstacks_radtags -c -t
