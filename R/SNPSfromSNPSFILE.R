#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# Argument 1 (Required): snpsfile
# Argument 2 (Optional): directory for output

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (snpsfile).n", call.=FALSE)
} else if (length(args)==1) {
  # current directory is where we're working 
  directory <- getwd() 
}

## convert the args to variables
snpsfile<-read.delim(args[1],header = FALSE) #first argument is the snpsfile
directory<-args[2] # otherwise the second argument is the directory
if(length(grep("\\/$",directory))>0){
  directory<-gsub("\\/$","\\1",directory) # replace the trailing /
}

## pull out each snp.
for(i in 1:(nrow(snpsfile)-1)){
  write.table(snpsfile[i:(i+1),],paste(directory,rownames(snpsfile)[i],sep="/"),
              quote=F,col.names=F,row.names=F,sep='\t')
}