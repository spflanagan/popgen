#!/bin/bash

#run from popgen/nerophis directory
process_radtags -p ./nop1 -o ./nop1 -b ../96barcodes.txt -e pstI -r -q -c
sh ./nop1_barcodes.txt

process_radtags -p ./nop2 -o ./nop2 -b ../96barcodes.txt -e pstI -r -q -c
sh ./nop2_barcodes.txt

process_radtags -p ./nop3 -o ./nop3 -b ../96barcodes.txt -e pstI -r -q -c
sh ./nop3_barcodes.txt

wc -l ./samples/*.fq > nop.nseqs.txt
