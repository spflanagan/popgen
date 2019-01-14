#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (snpsfile).n", call.=FALSE)
} else if (length(args)==1) {
  # current directory is 
  directory <- getwd()
}

## convert the args to variables
snpsfile<-read.delim(args[1],row.names = 1)
directory<-args[2]

## pull out each snp.
for(i in 1:(nrow(snpsfile)-1)){
  write.table(snpsfile[i:(i+1),],paste(directory,rownames(snpsfile)[i],sep=""),
              quote=F,col.names=F,row.names=F,sep='\t')
}