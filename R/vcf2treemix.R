#!/usr/bin/env Rscript

# need gwscaR tools
# require(devtools)
# devtools::install_github("https://github.com/spflanagan/gwscaR.git")
# library(gwscaR)


# get the input
args = commandArgs(trailingOnly = TRUE)
vcf_name<-args[1]
poplist_file<-args[2]
treemix_name<-args[3]
# need gwscaR - this is the hacky way
gwscaR_dir<-args[4]

source(paste(gwscaR_dir,"R/gwscaR.R",sep="/"))
source(paste(gwscaR_dir,"R/gwscaR_plot.R",sep="/"))
source(paste(gwscaR_dir,"R/gwscaR_utility.R",sep="/"))
source(paste(gwscaR_dir,"R/gwscaR_fsts.R",sep="/"))
source(paste(gwscaR_dir,"R/gwscaR_popgen.R",sep="/"))


# read the files and convert
vcf<-parse.vcf(vcf_name)
poplist<-read.delim(poplist_file)
tm<-treemix.from.vcf(vcf,poplist)
write.table(tm,treemix.name,col.names=TRUE,row.names=FALSE,quote=F,sep=' ')