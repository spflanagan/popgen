#' ---
#' title: "Looking for sex-linked markers in S. scovelli"
#' subtitle: "Using twelve marine populations"
#' author: "Sarah P. Flanagan"
#' date: "7 June 2017"
#' ---

#' ### Read in the data
setwd("B:/ubuntushare/popgen/fwsw_results") #modify this as needed
source("../../gwscaR/R/gwscaR.R")

#' ```{r,eval=FALSE}
vcf<-parse.vcf("stacks/fw-sw_populations/batch_2.vcf")
vcf$SNP<-paste(vcf$`#CHROM`,as.numeric(as.character(vcf$POS)),sep=".")
#choose one snp per locus
#new.vcf<-choose.one.snp(vcf) #8141 SNPs
#write.table(new.vcf,"stacks/fw-sw_populations/oneSNPperLoc.vcf",
 #           col.names = T,row.names=F,sep='\t',quote=F)
#' ```

#' ```{r,echo=FALSE,eval=FALSE}
new.vcf<-read.table("stacks/fw-sw_populations/oneSNPperLoc.vcf",sep='\t',header=T)
#' ```
#' 
#' #### Define some useful things
sw.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB",
           "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)

scaffs<-levels(as.factor(vcf[,1]))
scaffs[1:22]<-lgs
scaff.starts<-tapply(vcf$POS,vcf$`#CHROM`,max)
scaff.starts<-data.frame(rbind(cbind(names(scaff.starts),scaff.starts)),stringsAsFactors = F)
scaff.widths<-data.frame(levels(as.factor(vcf$`#CHROM`)),tapply(as.numeric(as.character(vcf$POS)),vcf$`#CHROM`,max))
locus.info<-c(colnames(vcf)[1:9],"SNP")

#' ### Run the Fst analyses
sexlinked<-NULL
fsts<-list()
uppers<-NULL
for(i in 1:length(sw.list)){
  pop.vcf<-cbind(new.vcf[,locus.info],new.vcf[,grep(sw.list[i],colnames(new.vcf))])
  #pull out the male and female names
  mal<-grep(paste("sample_",sw.list[i],"[PN]\\w+",sep=""),colnames(pop.vcf))
  fem<-grep(paste("sample_",sw.list[i],"[FD]\\w+",sep=""),colnames(pop.vcf))
  #use gwsca to run a pairwise fst calculation
  outliers<-gwsca(pop.vcf,locus.info,colnames(mal),colnames(fem),maf.cutoff=0.01)
  outliers<-outliers[outliers$Fst>0,] #remove the ones that don't pass the thresholds
  outliers$SNP<-paste(outliers$Chrom,as.numeric(as.character(outliers$Pos)),sep=".")
  outliers$locus<-paste(outliers$Chrom,pop.vcf$ID[pop.vcf$SNP %in% outliers$SNP],sep=".")
  #pull out the ones with sig. chi-squared tests
  sexlinked<-rbind(sexlinked,cbind(sw.list[i],outliers[outliers$Chi.p <= 0.05,]))
  #get the ones above a 95% threshold
  fst.thresh<-quantile(outliers$Fst,probs=0.95)
  upp<-outliers[outliers$Fst >= fst.thresh,]
  uppers<-rbind(uppers,cbind(sw.list[i],upp$locus))
  fsts[[i]]<-as.data.frame(outliers)
}
sexlinked$SNP<-paste(sexlinked$Chrom,as.numeric(as.character(sexlinked$Pos)),sep=".") #add "SNP" to sexlinked

#' ### Looking for shared ones

#' start with ones with significant chi-squared tests
data.frame(table(sexlinked$locus))
dups<-sexlinked[sexlinked$locus %in% sexlinked[duplicated(sexlinked$locus),"locus"],]
dups<-dups[order(dups$locus),]
dupsnps<-dups$locus
max(data.frame(table(sexlinked$locus))[,2])

#' then look at the ones that were above the 95% threshold
numerous.dups<-data.frame(table(uppers[,2]),stringsAsFactors=F)
shared<-numerous.dups[numerous.dups[,2]>=4,]
max(numerous.dups[,2])

#' unfortunately these don't have many overlapping between populations

#' ### Reorganize the data for plotting
pops<-list()
for(i in 1:nrow(numerous.dups)){
  this.pop<-""
  for(ii in 1:length(sw.list)){
         if(length(fsts[[ii]][fsts[[ii]]$locus %in% shared[i,1],])>0){
           this.loc<-fsts[[ii]][fsts[[ii]]$locus %in% shared[i,1],]
           if(nrow(this.loc[this.loc$Fst >= quantile(fsts[[ii]]$Fst,0.95),])>0){
            this.pop<-paste(this.pop,sw.list[ii],sep=",")
           }
         }
  }
  pops[[i]]<-this.pop
}

#' ### Plot the Fsts
#' `#+ fig.width=7, fig.height=5, dpi=50`
png("../popgen_mal-fem.png",height=5,width=7,units="in",res=300)
par(mfrow=c(3,4),oma=c(0.5,0.5,0.5,0.5),mar=c(2,2,2,2))

for(i in 1:length(sw.list)){
  fst<-fst.plot(fsts[[i]],bp.name="Pos",y.lim = c(0,0.3),scaffold.widths=scaff.widths,axis.size=1)
  abline(h=quantile(fsts[[i]]$Fst,0.95),col="cornflowerblue")
  points(fst[fst$locus %in% dupsnps,"plot.pos"],fst[fst$locus %in% dupsnps,"Fst"],col="cornflowerblue") #9 in 2 pops
  points(fst[fst$locus %in% shared[,1],"plot.pos"],fst[fst$locus %in% shared[,1],"Fst"],col="violet",pch=19) #13 in 4 pops
  legend("top",bty='n',legend=sw.list[i])
}

dev.off()