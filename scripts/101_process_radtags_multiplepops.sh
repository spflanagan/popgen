#!/bin/bash

process_radtags -p ./raw/popgen_1 -o ./samples/popgen_1 -b ./popgen1_barcodes.txt -e pstI -r -q -c

process_radtags -p ./raw/popgen_3 -o ./samples/popgen_3 -b ./popgen3_barcodes.txt -e pstI -r -q -c

process_radtags -p ./raw/popgen_4 -o ./samples/popgen_4 -b ./popgen4_barcodes.txt -e pstI -r -q -c
