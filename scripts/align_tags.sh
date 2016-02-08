#!/bin/bash

#align rad tags (fasta file converted from batch_x.catalog.matches.tsv) 
#to the genome using bowtie2



bowtie2 --sensitive-local -x scovelli_allpaths -S tags.sensitivelocal.sam -f -U catalog.tags.fasta -p 3 --no-hd

bowtie2 --sensitive -x scovelli_allpaths -S tags.sensitive.sam -f -U catalog.tags.fasta -p 3 --no-hd

bowtie2 -x scovelli_allpaths -S tags.N1ete.sam -f -U catalog.tags.fasta -p 3 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --end-to-end --no-hd

bowtie2 -x scovelli_allpaths -S tags.N1local.sam -f -U catalog.tags.fasta -p 3 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --local --no-hd
