#!/bin/bash

populations -P ./rxstacks/ -M ./population_map_aligned.txt -b 1 -k -p 2 -r 0.75 -t 6 --structure --genepop --vcf --phase --plink --beagle
