#!/bin/bash

# Run from within the bayenv directory

############## First create a clust file ##################
# Do this manually

############## 			Settings	 	 ##################
ANALYSIS="$1"


############## 	 Setup and estimate matrices   ##################
if [ "$ANALYSIS" == "FILEMANIP" ]; then

	CLUST_NAME="$2" #bayenv/sub75.pruned.clust
	PED_NAME="$3" #stacks/populations_subset75/batch_2.pruned.ped
	MAP_NAME="$4" #stacks/populations_subset75/batch_2.pruned.map
	OUT_DIR="$5" #bayenv
	BAYENV_DIR="$6" # ~/Programs/bayenv/
	NUMPOPS="$7" # 7

	# get plink frq output
	~/Programs/plink-1.07-x86_64/plink --ped "$PED_NAME" --map "$MAP_NAME" \
	--out "$OUT_DIR"/bayenv --noweb --allow-no-sex --recode --freq --within "$CLUST_NAME"

	Rscript ../R/SNPSFILEfromPLINKfrq.R "${OUT_DIR}/bayenv.frq.strat" "$OUT_DIR"

fi

############## 	   Create matrix files	 ##################
if [ "$ANALYSIS" == "MATRIX" ]; then

	OUT_DIR="$2" #bayenv
	BAYENV_DIR="$3" # ~/Programs/bayenv/
	NUMPOPS="$4" 
	
	#run this in the directory where you want the matrix files.
	cd ${OUT_DIR}


	for i in `seq 1 10`;
	do
		echo "Running rep $i"
		${BAYENV_DIR}bayenv2 -i SNPSFILE -p "$NUMPOPS" -k 100000 > matrix.$i.out
		grep "ITER = 100000" matrix.$i.out -A ${NUMPOPS} > matrix.$i.out.last

	done

	# AFTER THIS IT IS IMPORTANT TO COMPARE MATRICES #
	Rscript ../../R/plot_bayenv_matrices.R fwsw75_pruned.png

fi


############## 	   	Create SNPFILES		 ##################
if [ "$ANALYSIS" == "SNPFILES" ]; then
	SNPSFILE="$2"
	SNPFILEDIR="$3"
	# see if the directory exists and if it doesn't make it
	if [ ! -d "$SNPFILEDIR" ];
	then
		mkdir -p "$SNPFILEDIR"
	fi
	
	Rscript ../../R/SNPSfromSNPSFILE.R "$SNPSFILE" ${SNPFILEDIR}	
fi

############## 	   		Run bayenv2	 	 ##################
if [ "$ANALYSIS" == "BAYENV" ]; then
	BAYENV_DIR="$2" # ~/Programs/bayenv/
	MATRIXFILE="$3"
	ENVIRONFILE="$4"
	NUMPOPS="$5"
	NUMENVIRON="$6"
	SNPFILEDIR="$7"
	FILES="$SNPFILEDIR/*"

	echo "Matrix: $MATRIXFILE"
	echo "Environment: $ENVIRONFILE"
	echo "Num Pops: $NUMPOPS"
	echo "Num Environmental Variables: $NUMENVIRON"
	echo "Directory with SNPs files: $SNPFILEDIR"

	#run from bayenv folder
	for filename in $FILES; do
		echo "Running bayenv2 on $filename"
		${BAYENV_DIR}bayenv2 -i $filename -m "$MATRIXFILE" \
			-e "$ENVIRONFILE" -p "$NUMPOPS" -k 100000 -n "$NUMENVIRON" -t -c -f -X \
			-r 628398 -o $filename.freq

	done
fi
