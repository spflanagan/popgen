
# GOAL: subset ped file for power analysis with pcadapt

library(pcadapt)
setwd("Research/popgen/fwsw_results/")

ped<-read.delim("subset.ped",header=FALSE,sep=" ", stringsAsFactors = FALSE)

ns<-c(10,20,30)

pcadapt_subset<-function(ped, num_inds){
  keepinds<-unlist(tapply(ped[,2],ped[,1],sample,size=num_inds,replace=FALSE))
  pedsub<-ped[ped[,2] %in% keepinds,]
  subname<-paste("subset",num_inds,"inds.ped",sep="")
  write.table(pedsub,subname,col.names = FALSE,row.names = FALSE,quote=FALSE,sep=" ")
  filename<-read.pcadapt(subname,type="ped")
  x<-pcadapt(filename, K=6)
  saveRDS(x,paste("subset",num_inds,"inds.RDS",sep=""))
  plot(x,option="screeplot")#K=6
  return(x)
}

par(mfrow=c(2,2))
pcasubs<-lapply(ns,pcadapt_subset,ped = ped)

