
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



#' Simulate genetic data for populations with different allele frequencies
#' @param npops The number of populations to simulate
#' @param inds The number of individuals per population to simulate
#' @param nloc The number of loci to simulate
#' @param p A vector of average minor allele frequencies per locus to simulate
#' @return A data.frame with ped data
popgen.sim<-function(npops=8,inds=10,nloc=10000,p=rep(0.01,8),outname="simulated.ped",analyze=TRUE){

  #install and load pcadapt if you haven't already
  if("pcadapt" %in% rownames(installed.packages())){
    do.call('library',list("pcadapt"))
  }else{
    install.packages(package,dependencies = TRUE)
    do.call("library",list("pcadapt"))
  }
  #create ped first columns
  simped<-data.frame(cbind(rep(1:npops,each=inds),rep(paste("Ind",rep(1:inds),sep=""),npops),
                           rep(0,inds*npops),rep(0,inds*npops),rep(0,inds*npops),rep(-9,inds*npops)),
                     stringsAsFactors = FALSE)
  
  #simulate allele frequencies for each locus
  frqs<-lapply(p,function(f){
    x<-rnorm(f,n=nloc,sd=f)
    z<-unlist(lapply(x,function(y){ #keeps it between 0 and 1
      while(y < 0 | y > 1){
        y<-rnorm(f,n=1,sd=f)
      }
      return(y)
    }))
    return(z)
  })
  
  #for each locus, assign genotypes for each population
  for(n in seq(1,nloc*2,by=2)){
    #designate columns
    col1<-n+6
    col2<-n+7

    #simulate alleles
    al1<-unlist(lapply(frqs,function(f){ rbinom(f[(n+1)/2],n=inds,size=1) } ))
    al2<-unlist(lapply(frqs,function(f){ rbinom(f[(n+1)/2],n=inds,size=1) } ))
    if(rbinom(1,1,0.5)==0){ #randomly decide if it's going to be G/C or A/T locus
      if(rbinom(1,1,0.5)==0){ #randomly decide if major allele is G
        al1[al1==0]<-"G"
        al1[al1==1]<-"C"
        al2[al2==0]<-"G"
        al2[al2==1]<-"C"
      }else{ #or if it's C
        al1[al1==0]<-"C"
        al1[al1==1]<-"G"
        al2[al2==0]<-"C"
        al2[al2==1]<-"G"    
      }
    }else{
      if(rbinom(1,1,0.5)==0){ #randomly decide if major allele is A
        al1[al1==0]<-"A"
        al1[al1==1]<-"T"
        al2[al2==0]<-"A"
        al2[al2==1]<-"T"
      }else{ #or if it's T
        al1[al1==0]<-"T"
        al1[al1==1]<-"A"
        al2[al2==0]<-"T"
        al2[al2==1]<-"A"
      }
    }
    #put this in the ped object
    simped[,col1]<-al1
    simped[,col2]<-al2
  }
  write.table(simped,outname,col.names = FALSE,row.names = FALSE,quote=FALSE,sep=" ")

  if(isTRUE(analyze)){
    # now run pcadapt
    filename<-read.pcadapt(outname,type="ped")
    x<-pcadapt(filename,K=20,min.maf=0.001)
    png(paste(outname,".png",sep=""),height = 7, width = 7, units="in",res=300)
    plot(x,option="scores",pop=simped[,1])
    dev.off()
  }
  return(simped)
}

# Run the simulations
setwd("~/Research/popgen/simulations/")
ns<-c(10,20,30,40,50)
sims<-lapply(ns,function(n){
  sim<-popgen.sim(npops=8,inds=n,nloc=10000,p=rep(0.05,8),outname=paste("sim",n,"i.ped",sep=""))
  return(sim)
})




