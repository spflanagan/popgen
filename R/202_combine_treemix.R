#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  stop("You must provide the stem name and outname.n", call.=FALSE)
} else if (length(args)==2) {
  args[3] = 10
}

combine_boots_treemix<-function(stem,out,boots=10){
  if(file.exists(paste0(stem, ".treeout.gz"))){
    tr<- read.tree(gzfile(paste0(stem, ".treeout.gz")))    
    write.tree(tr,paste0(stem, ".treeout.tre"))
  }
  for(i in 1:boots){
    tr<-paste0(stem,"_",i, ".treeout.gz")
    tr<- read.tree(gzfile(tr))    
    write.tree(tr,out,append=TRUE)
  }
  invisible(out)
}

library(ape)
combine_boots_treemix(args[1],args[2],args[3])
