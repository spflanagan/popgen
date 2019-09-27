#!/usr/bin/env Rscript
require(gdata)
args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (plink.frq.strat).n", call.=FALSE)
} else if (length(args)==1) {
  # current directory is 
  directory <- getwd()
} else if (length(args)==2){
  directory<-args[2]
}

# remove trailing slash from directory name
if(length(grep("\\/$",directory))>0){ 
  directory<-gsub("\\/$","",directory)
}

#####CONVERT PLINK TO BAYENV2
freq<-read.table(args[1],header=T, stringsAsFactors=F)
#want to get $MAC for every snp at every pop 
#and NCHROBS-MAC for every stnp at every pop
freq<-cbind(freq,freq$NCHROBS-freq$MAC)
colnames(freq)[ncol(freq)]<-"NAC"
pop.order<-levels(as.factor(freq$CLST))
snp.names<-split(freq$SNP,freq$CLST)[[1]]

mac.by.pop<-as.data.frame(split(freq$MAC,freq$CLST))
rownames(mac.by.pop)<-snp.names
nac.by.pop<-as.data.frame(split(freq$NAC,freq$CLST))
rownames(nac.by.pop)<-snp.names
snpsfile<-interleave(mac.by.pop,nac.by.pop)

write.table(snpsfile, paste(directory,"SNPSFILE",sep="/"), 
            col.names=F,row.names=F,quote=F,sep="\t",eol="\n") #bayenv SNPSFILE
