#!/bin/bash
src="/home/joneslab/Documents/Sarahs/PstI_Data"



##do this still:
#load_radtags
echo "load_radtags"
load_radtags.pl -D rx_radtags -p $src/rxstacks/ \
	-b 1 -e "rx test, 10X cov" -c -B
	
#index
index_radtags.pl -D rx_radtags -c
index_radtags.pl -D rx_radtags

