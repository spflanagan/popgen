#!/bin/bash

process_radtags -p ./raw/popgen_5 -o ./samples/popgen_5 -b ./popgen5_barcodes.txt -e pstI -r -q -c

process_radtags -p ./raw/popgen_6 -o ./samples/popgen_6 -b ./popgen6_barcodes.txt -e pstI -r -q -c

process_radtags -p ./raw/popgen_7 -o ./samples/popgen_7 -b ./popgen7_barcodes.txt -e pstI -r -q -c

process_radtags -p ./raw/popgen_8 -o ./samples/popgen_8 -b ./popgen8_barcodes.txt -e pstI -r -q -c
