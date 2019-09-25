#!/bin/bash

# Run from within the results directory

############## First create a clust file ##################
# Do this manually

############## 			Settings	 	 ##################
PLINK=false
MATRIX=true


############## 	Then run this script	 ##################
CLUST_NAME="$1" #bayenv/fwsw_sub75.clust
PED_NAME="$2" #stacks/populations_subset75/batch_2.plink.ped
MAP_NAME="$3" #stacks/populations_subset75/batch_2.plink.map
OUT_DIR="$4" #bayenv
BAYENV_DIR="$5" # ~/Programs/bayenv/

SNPSFILE="$6"
NUMPOPS="$7"

############## 		Then run plink 		 ##################
if [ "$PLINK" = true ]; then
	~/Programs/plink-1.07-x86_64/plink --ped "$PED_NAME" --map "$MAP_NAME" \
	--out "$OUT_DIR"/bayenv --noweb --allow-no-sex --recode --freq --within "$CLUST_NAME"
fi

############## 	   Create matrix files	 ##################
#run this in the directory where you want the matrix files.
cd ${OUT_DIR}

if [ "$MATRIX" = true ]; then
	for i in `seq 1 10`;
	do
		echo "Running rep $i"
		${BAYENV_DIR}/bayenv2 -i "$SNPSFILE" -p "$NUMPOPS" -k 100000 -r 628398 > matrix.$i.out
	done
fi

############## 	   		Run bayenv2	 	 ##################

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
