#!/usr/bin/env Rscript

# need gwscaR tools
# require(devtools)
# devtools::install_github("https://github.com/spflanagan/gwscaR.git")
# library(gwscaR)

# Generate a treemix file from vcf
# from gwscaR
treemix.from.vcf<-function(vcf,pop.list){
  tm.df<-matrix(nrow=nrow(vcf),ncol=length(pop.list))
  for(i in 1:nrow(vcf)){
    vcf.row<-vcf[i,]
    #tm.df<-as.matrix(t(apply(vcf,1,function(vcf.row){ #freaking apply was giving me weird results
    all.alleles<-names(table(vcf.alleles(vcf.row))) 
    if(length(all.alleles)==2){
      this.loc<-do.call("cbind",  
                        lapply(pop.list,function(pop){
                          pop.vcf.row<-cbind(vcf.row[1:9],vcf.row[grep(pop,colnames(vcf.row))])
                          pop.alleles<-vcf.alleles(pop.vcf.row)
                          pop.counts<-table(pop.alleles)
                          tm.pop<-"0,0"
                          if(length(pop.counts)==1){
                            allele<-names(pop.counts)
                            if(grep(allele,all.alleles)==1){ #maintain order
                              tm.pop<-paste(pop.counts[[1]],0,sep=",") 
                            }else{
                              tm.pop<-paste(0,pop.counts[[1]],sep=",") }
                          }
                          if(length(pop.counts)==2){
                            tm.pop<-paste(pop.counts[[all.alleles[1]]],pop.counts[[all.alleles[2]]],sep=",") 
                          }
                          return(tm.pop)
                        }))
    }
    #return(this.loc)
    tm.df[i,]<-as.vector(this.loc)
  }#)))
  colnames(tm.df)<-pop.list
  return(tm.df)
}

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