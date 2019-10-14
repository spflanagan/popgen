#!/usr/bin/env Rscript

# need gwscaR tools
require(devtools)
devtools::install_github("https://github.com/spflanagan/gwscaR.git")
library(gwscaR)

# get the input
args = commandArgs(trailingOnly = TRUE)
vcf_name<-args[1]
poplist_file<-args[2]
treemix_name<-args[3]

# read the files and convert
vcf<-parse.vcf(vcf_name)
poplist<-read.delim(args[2])
tm<-treemix.from.vcf(vcf,poplist)
write.table(tm,treemix.name,col.names=TRUE,row.names=FALSE,quote=F,sep=' ')