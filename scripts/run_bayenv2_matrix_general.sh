#!/bin/bash

# Run from within the results directory

############## First create a clust file ##################
# Do this manually

############## 			Settings	 	 ##################
FILEMANIP=false
MATRIX=true
BAYENV=false

############## 	 Setup and estimate matrices   ##################
if [ "$FILEMANIP" = true ]; then

	CLUST_NAME="$1" #bayenv/sub75.pruned.clust
	PED_NAME="$2" #stacks/populations_subset75/batch_2.pruned.ped
	MAP_NAME="$3" #stacks/populations_subset75/batch_2.pruned.map
	OUT_DIR="$4" #bayenv
	BAYENV_DIR="$5" # ~/Programs/bayenv/
	NUMPOPS="$7" # 7

	# get plink frq output
	~/Programs/plink-1.07-x86_64/plink --ped "$PED_NAME" --map "$MAP_NAME" \
	--out "$OUT_DIR"/bayenv --noweb --allow-no-sex --recode --freq --within "$CLUST_NAME"

	Rscript ../R/SNPSFILEfromPLINKfrq.R "${OUT_DIR}/bayenv.frq.strat" "$OUT_DIR"

fi

if [ "$MATRIX" = true ]; then

	OUT_DIR="$1" #bayenv
	BAYENV_DIR="$2" # ~/Programs/bayenv/
	NUMPOPS="$3" 
	############## 	   Create matrix files	 ##################
	#run this in the directory where you want the matrix files.
	cd ${OUT_DIR}


	for i in `seq 1 10`;
	do
		echo "Running rep $i"
		${BAYENV_DIR}bayenv2 -i SNPSFILE -p "$NUMPOPS" -k 100000 -r 628398 > matrix.$i.out
	done

	# AFTER THIS IT IS IMPORTANT TO COMPARE MATRICES #
fi

############## 	   		Run bayenv2	 	 ##################
if [ "$BAYENV" = true ]; then
	MATRIXFILE="$1"
	ENVIRONFILE="$2"
	NUMPOPS="$3"
	NUMENVIRON="$4"
	SNPFILEDIR="$5"
	FILES="$SNPFILEDIR/*"

	echo "Matrix: $MATRIXFILE"
	echo "Environment: $ENVIRONFILE"
	echo "Num Pops: $NUMPOPS"
	echo "Num Environmental Variables: $NUMENVIRON"
	echo "Directory with SNPs files: $SNPFILEDIR"

	#run from bayenv folder
	for filename in $FILES; do
		echo "Running bayenv2 on $filename"
		~/Programs/bayenv_2/bayenv2 -i $filename -m "$MATRIXFILE" \
			-e "$ENVIRONFILE" -p "$NUMPOPS" -k 100000 -n "$NUMENVIRON" -t -c -f -X \
			-r 628398 -o $filename.freq

	done
fi