
# GOAL: subset ped file for power analysis with pcadapt

library(pcadapt)
setwd("Research/popgen/fwsw_results/")

ped<-read.delim("subset.ped",header=FALSE,sep=" ", stringsAsFactors = FALSE)

ns<-c(10,20,30,0)

pcadapt_subset<-function(ped, num_inds){
  if(num_inds==0){
    keepinds<-ped[,2]
  }else{
    keepinds<-unlist(tapply(ped[,2],ped[,1],sample,size=num_inds,replace=FALSE))
  }
  pedsub<-ped[ped[,2] %in% keepinds,]
  subname<-paste("subset",num_inds,"inds.ped",sep="")
  write.table(pedsub,subname,col.names = FALSE,row.names = FALSE,quote=FALSE,sep=" ")
  filename<-read.pcadapt(subname,type="ped")
  x<-pcadapt(filename,K=20)
  saveRDS(x,paste("subset",num_inds,"inds.RDS",sep=""))
  png(paste("subset",num_inds,"inds.png"),height=7,width=7,units="in",res=300)
  plot(x,option="scores",pop=pedsub[,1])#K=6
  dev.off()
  return(x)
}

pcasubs<-lapply(ns,pcadapt_subset,ped = ped)

