source("../../gwscaR/R/gwscaR.R")

pop.list<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLLG")
pop.labs<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLFW")
fw.list<-c("TXFW","LAFW","ALFW","FLLG")
sw.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB",
           "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)
all.colors<-c(rep("black",2),"#2166ac","black","#2166ac","black","#2166ac",
              rep("black",8),"#2166ac")
grp.colors<-c('#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837')

vcf<-parse.vcf("stacks/fw-sw_populations/batch_2.vcf")
vcf$SNP<-paste(vcf$`#CHROM`,as.numeric(as.character(vcf$POS)),sep=".")
scaffs<-levels(as.factor(vcf[,1]))
scaffs[1:22]<-lgs
scaff.starts<-tapply(vcf$POS,vcf$`#CHROM`,max)
scaff.starts<-data.frame(rbind(cbind(names(scaff.starts),scaff.starts)),stringsAsFactors = F)
locus.info<-c(colnames(vcf)[1:9],"SNP")

sexlinked<-NULL
fsts<-list()
uppers<-NULL
for(i in 1:length(sw.list)){
  pop.vcf<-cbind(vcf[,locus.info],vcf[,grep(sw.list[i],colnames(vcf))])
  mal<-cbind(vcf[,1:9],
             pop.vcf[,grep(paste("sample_",sw.list[i],"[PN]\\w+",sep=""),
                           colnames(pop.vcf))])
  fem<-cbind(vcf[,1:9],
             pop.vcf[,grep(paste("sample_",sw.list[i],"[FD]\\w+",sep=""),
                           colnames(pop.vcf))])
  outliers<-gwsca(pop.vcf,locus.info,colnames(mal),colnames(fem))
  outliers$SNP<-paste(outliers$Chrom,as.numeric(as.character(outliers$Pos)),sep=".")
  outliers$locus<-paste(outliers$Chrom,pop.vcf$ID[pop.vcf$SNP %in% outliers$SNP],sep=".")
  fst.thresh<-quantile(outliers$Fst,probs=0.95)
  sexlinked<-rbind(sexlinked,cbind(sw.list[i],outliers[outliers$Chi.p <= 0.05,]))
  #sexlinked<-rbind(sexlinked,cbind(sw.list[i],outliers[outliers$Fst >= fst.thresh,]))
  upp<-outliers[outliers$Fst >= quantile(outliers$Fst,0.95),]
  uppers<-rbind(uppers,cbind(sw.list[i],upp$locus))
  fsts[[i]]<-as.data.frame(outliers)
}
sexlinked$SNP<-paste(sexlinked$Chrom,as.numeric(as.character(sexlinked$Pos)),sep=".")

data.frame(table(sexlinked$locus)) #LG7.8032058    2
dups<-sexlinked[sexlinked$locus %in% sexlinked[duplicated(sexlinked$locus),"locus"],]
dups<-dups[order(dups$SNP),]
dupsnps<-dups$SNP

numerous.dups<-data.frame(table(uppers[,2]))
shared<-numerous.dups[numerous.dups[,2]>=12,]

pops<-list()
for(i in 1:nrow(shared)){
  this.pop<-""
  for(ii in 1:length(sw.list)){
         if(length(fsts[[ii]][fsts[[ii]]$locus %in% shared[i,1],])>0){
           pops[[i]]<-paste(this.pop,sw.list[i],sep="")
         }
  }
}

par(mfrow=c(3,4),oma=c(2,2,2,2),mar=c(2,2,2,2))

for(i in 1:length(sw.list)){
  fst<-fst.plot(fsts[[i]],bp.name="Pos")
  abline(h=quantile(fsts[[i]]$Fst,0.95),col="cornflowerblue")
  fst$SNP<-paste(fst$Chrom,as.numeric(as.character(fst$Pos)),sep=".")
  points(fst[fst$SNP %in% dupsnps,"plot.pos"],fst[fst$SNP %in% dupsnps,"Fst"],col="cornflowerblue")
  points(fst[fst$SNP %in% uppers$SNP,"plot.pos"],fst[fst$SNP %in% uppers$SNP,"Fst"],col="violet",pch=19)
}




#so, now we don't have these things
