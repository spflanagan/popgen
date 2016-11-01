#!/bin/bash

MATRIXFILE="$1"
ENVIRONFILE="$2"
NUMPOPS="$3"
NUMENVIRON="$4"
SNPFILEDIR="$5"
FILES="$SNPFILEDIR/*"

#run from bayenv folder
for filename in $FILES; do
	echo "Running bayenv2 on $filename"
	~/Programs/bayenv_2/bayenv2 -i $filename -m "$1" \
		-e "$2" -p "$NUMPOPS" -k 100000 -n "$NUMENVIRON" -t -c -f -X \
		-r 628398 -o $filename.freq

done
