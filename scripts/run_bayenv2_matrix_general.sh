#!/bin/bash

#run this in the directory where you want the matrix files.
SNPSFILE="$1"
NUMPOPS="$2"

for i in {1..10}
do
	echo "Running rep $i"
	~/Programs/bayenv_2/bayenv2 -i "$SNPSFILE" -p "$NUMPOPS" -k 100000 -r 628398 > matrix.$i.out
done
