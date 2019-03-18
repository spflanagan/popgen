#!/bin/bash

# Run from within the results directory

############## First create a clust file ##################
CLUST_NAME="$1" #bayenv/fwsw_sub75.clust
PED_NAME="$2" #stacks/populations_subset75/batch_2.plink.ped
MAP_NAME="$3" #stacks/populations_subset75/batch_2.plink.map
OUT_DIR="$4" #bayenv
BAYENV_DIR="$5" # ~/Programs/bayenv/

############## 		Then run plink 		 ##################
~/Programs/plink-1.07-x86_64/plink --ped "$PED_NAME" --map "$MAP_NAME" \
--out "$OUT_DIR"/bayenv --noweb --allow-no-sex --recode --freq --within "$CLUST_NAME"


#run this in the directory where you want the matrix files.
cd ${OUT_DIR}

SNPSFILE="$1"
NUMPOPS="$2"

for i in `seq 1 10`;
do
	echo "Running rep $i"
	${BAYENV_DIR}/bayenv2 -i "$SNPSFILE" -p "$NUMPOPS" -k 100000 -r 628398 > matrix.$i.out
done
