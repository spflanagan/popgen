
#### CONVERT VCF TO SWEEPFINDER ####
#one file for each chromosome
#first column is position
#second column is count of derived alleles
#third column is the total number of observed alleles
#fourth column is optional and indicates whether it's derived or ancestral. 

parse.vcf<-function(filename){
  #if(substr(filename,nchar(filename)-3,nchar(filename)) != ".vcf") { filename<-paste(filename,"vcf",sep=".") }
  vcf<-read.delim(filename,comment.char="#",sep='\t',header=F,stringsAsFactors = F,strip.white = T)
  header.start<-grep("#CHROM",scan(filename,what="character"))
  header<-scan(filename,what="character")[header.start:(header.start+ncol(vcf)-1)]
  colnames(vcf)<-header
  return(vcf)
}

extract.gt.vcf<-function(vcf){
  if(length(strsplit(as.character(vcf[1,10]),":")[[1]])>1){
    new<-vcf[,1:9]
    for(i in 10:ncol(vcf)){
      new<-cbind(new,
                 sapply(vcf[,i],function(x) {
                   strsplit(as.character(x),":")[[1]][1]})
      )
    }
    colnames(new)<-colnames(vcf[,c(1:9,10:ncol(vcf))])
    vcf<-new
  }
  return(vcf)
}

gts2sfaf<-function(row){
  gt<-row[10:length(row)]
  ref<-(length(gt[gt=="0/0"])*2)+
    length(gt[gt=="0/1"])+length(gt[gt=="1/0"])
  alt<-(length(gt[gt=="1/1"])*2)+
    length(gt[gt=="0/1"])+length(gt[gt=="1/0"])
  return(data.frame(position=row["POS"],x=ref,n=(alt+ref)))
}

vcf2sfaf<-function(vcf,lgs){
  chrs<-lapply(lgs,function(lg){
    gts<-extract.gt.vcf(vcf[vcf$`#CHROM` == lg,])
    sf.af<-do.call(rbind,apply(gts,1,gts2sfaf))
    write.table(sf.af,paste("SF2/",lg,".AF.txt",sep=""),
                col.names = TRUE,row.names=FALSE,sep='\t',
                quote=FALSE,eol='\n')
    print(paste("Writing to file: SF2/",lg,".AF.txt",sep=""))
    return(sf.af)
  })
}

vcf<-parse.vcf("p4.upd.vcf")
lgs<-unique(vcf$`#CHROM`)
#frequencies per LG
af<-vcf2sfaf(vcf,lgs[1:22])
#whole genome frequencies
gts<-extract.gt.vcf(vcf[vcf$`#CHROM` %in% lgs,])
gaf<-do.call(rbind,apply(gts,1,gts2sfaf))
write.table(gaf,"SF2/GenomeWideSpectrum.txt",
            col.names = TRUE,row.names=FALSE,sep='\t',
            quote=FALSE,eol='\n')