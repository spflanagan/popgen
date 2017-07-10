#Author: Sarah P. Flanagan
#Last Updated: 20 July 2016
#Purpose: Analyze FW-SW dataset

rm(list=ls())

library(ade4)
library(lme4)
library(maps);library(gplots)
library(mapdata)
library(vegan)
library(boot)
library(adegenet)
library(scales)
library(gdata)
library(ape)
library(lattice); library(RColorBrewer); library(grid)
library(devtools)

setwd("B:/ubuntushare/popgen/fwsw_results/")
#source("../scripts/popgen_functions.R")
install_github("spflanagan/gwscaR")
library(gwscaR)
source("../scripts/phenotype_functions.R")

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
#######################**********FILES*********##############################
mar.coor<-read.csv("../sw_results/marine_coordinates_revised.csv", header=T)
fw.coor<-read.csv("fw_coordinates.csv", header=T)
dist<-read.table("fwsw_geographical_distances.txt",header=T,row.names=1,
	sep='\t')
pwise.fst.all<-read.table("stacks/populations/fwsw_fst_summary.txt",header=T,row.names=1,sep='\t')
	pwise.fst.all<-rbind(pwise.fst.all,rep(NA,ncol(pwise.fst.all)))
	rownames(pwise.fst.all)<-colnames(pwise.fst.all)
pwise.fst.sub<-read.table("stacks/fwsw_fst_summary_subset.txt",header=T,row.names=1,sep='\t')
ped.sub<-read.table("stacks/subset.ped",header=F)	
ped.sub$V1<-gsub("sample_(\\w{4})\\w+.*","\\1",ped.sub$V2)
map.sub<-read.table("stacks/subset.map",header = F,stringsAsFactors = F)
map.sub$Locus<-paste(gsub("(\\d+)_\\d+","\\1",map.sub$V2),map.sub$V4,sep=".")
colnames(ped.sub)<-c("Pop","IID","","","","Phenotype","","",map.sub$Locus)
vcf<-parse.vcf("stacks/fw-sw_populations/batch_2.vcf")
vcf$SNP<-paste(vcf$`#CHROM`,vcf$POS,sep=".")
scaffs<-levels(as.factor(vcf[,1]))
scaffs[1:22]<-lgs
scaff.starts<-tapply(vcf$POS,vcf$`#CHROM`,max)
scaff.starts<-data.frame(rbind(cbind(names(scaff.starts),scaff.starts)),stringsAsFactors = F)
locus.info<-c(colnames(vcf[1:9]),"SNP")
#chosen.snps<-choose.one.snp(vcf)$SNP
#write.table(chosen.snps,"chosen.snps.txt",quote=F)
chosen.snps<-unlist(read.table("chosen.snps.txt"))

#######################PLOT THE POINTS ON A MAP##############################
jpeg("all_sites_map.jpg", res=300, height=7,width=14, units="in")
pdf("all_sites_map.pdf",height=7,width=14)
par(oma=c(0,0,0,0),mar=c(0,0,0,0),pin=c(7,7))
map("worldHires", "usa",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray90", mar=c(0,0,0,0),fill=TRUE, res=300,myborder=0)
map("worldHires", "mexico",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray95", fill=TRUE, add=TRUE)
points(mar.coor$lon, mar.coor$lat,  col="black", cex=1.2, pch=19)
points(-1*fw.coor$lon, fw.coor$lat,  col="#2166ac", cex=1.5, pch=18)
abline(h=c(25,30,35),lty=3)
abline(v=c(-80,-85,-90,-95,-100),lty=3)
text(x=c(-99.5,-99.5),y=c(25,30),c("25N","30N"))
text(x=c(-80,-85,-90,-95),y=rep(31.8,4),c("80W","85W","90W","95W"))
text(y=26,x=-90,"Gulf of Mexico")
text(y=25.5,x=-98.5,"Mexico")
text(x=-91,y=31,"USA")
text(x=-78,y=29.5,"Atlantic Ocean")
text(x=-96.5,y=26,"TXSP",font=2)
text(x=-96.9,y=27.2,"TXCC",font=2)
text(x=-96,y=28.3,"TXFW",font=2,col="#2166ac")
text(x=-94.7,y=29,"TXCB",font=2)
text(x=-90.2,y=30.3,"LAFW",font=2,col="#2166ac")
text(x=-88,y=30,"ALST",font=2)
text(x=-87,y=30.75,"ALFW",font=2,col="#2166ac")
text(x=-85,y=29.4,"FLSG",font=2)
text(x=-83.5,y=29.2,"FLKB",font=2)
text(x=-83.2,y=27.6,"FLFD",font=2)
text(x=-82.2,y=26,"FLSI",font=2)
text(x=-80,y=24.8,"FLAB",font=2)
text(x=-79.5,y=26.8,"FLPB",font=2)
text(x=-79.7,y=27.2,"FLHB",font=2)
text(x=-80.2,y=28.2,"FLCC",font=2)
text(x=-80.9,y=29.3,"FLFW",font=2,col="#2166ac")
dev.off()

##################################IBD########################################

#Mantel test using geographical distances and fsts

#read in the subsetted fst summary from running populations with whitelist
ibd.all<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst.all)))
ibd.sub<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst.sub)))

#test
pairwise.fst(ped.sub,9,10,pop.list)

ibd.by.loc<-fst.ibd.byloc(ped.sub,dist,pop.list)  #all NAs
#ignore warnings?  In is.euclid(m1) : Zero distance(s)
rownames(ibd.by.loc)<-sub.map$V2


####PLOT FSTS####
colors<-c("blue","yellow","red")
pal<-colorRampPalette(colors)
ncol=80
cols<-pal(ncol)

png("Fst_fig.png",height=7,width = 8,units="in",res=300)
fst.lv<-levelplot(as.matrix(pwise.fst.all),col.regions=cols,alpha.regions=0.7,
          scales = list(x=list(rot=90),tck = 0),xlab="",ylab="")
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(italic(F)[ST]), 0.2, 0, hjust=0.5, vjust=1.2)
trellis.unfocus()
dev.off()

#compare to treemix
cov<-read.table(gzfile("treemix/fwsw.basic.cov.gz"), as.is = T, head = T, quote = "", comment.char = "")
#reorder
covplot <- data.frame(matrix(nrow = nrow(cov), ncol = ncol(cov)))
for(i in 1:length(pop.list)){
  for( j in 1:length(pop.list)){
    
    covplot[i, j] = cov[which(names(cov)==pop.list[i]), which(names(cov)==pop.list[j])]
    rownames(covplot)[i]<-pop.list[i]
    colnames(covplot)[j]<-pop.list[j]
  }
}
cp<-as.matrix(covplot)
cp[lower.tri(cp)]<-NA
cp[upper.tri(cp)]<-covplot[upper.tri(covplot)]
cp.lv<-levelplot(cp,col.regions=cols,alpha.regions=0.7,
          scales = list(x=list(rot=90),tck = 0),xlab="",ylab="")

#' Plot them together.
png("fst_cov_heatmap.png",height=6,width=11,units="in",res=300)
print(fst.lv,split=c(1,1,2,1),more=TRUE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(italic(F)[ST]), 0.2, 0, hjust=0.5, vjust=1.2)
trellis.unfocus()

print(cp.lv,split=c(2,1,2,1),more=FALSE)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text("covariance", 0.2, 0, hjust=0.5, vjust=1.2)
trellis.unfocus()
dev.off()




####ALLELE FREQUENCY SPECTRA####
#' Calculate AFS from vcf
fw.afs<-lapply(fw.list,function(pop){
  this.vcf<-cbind(vcf[,locus.info],vcf[,grep(pop,colnames(vcf))])
  this.afs<-do.call(rbind,apply(this.vcf,1,calc.afs.vcf))
})
names(fw.afs)<-c("TXFW","LAFW","ALFW","FLFW")
sw.afs<-lapply(sw.list,function(pop){
  this.vcf<-cbind(vcf[,locus.info],vcf[,grep(pop,colnames(vcf))])
  this.afs<-do.call(rbind,apply(this.vcf,1,calc.afs.vcf))
})
names(sw.afs)<-sw.list
all.afs<-c(fw.afs,sw.afs)
png("AFS.png",height=10,width=10,units="in",res=300)
par(mfrow=c(4,4),mar=c(2,2,1,0),oma=c(2,2,0.5,0.5))
for(i in 1:length(pop.labs)){
  if(pop.labs[i] %in% names(fw.afs)){
    color<-"cornflowerblue"
  }else{
    color<-"black"
  }
  hist(all.afs[[pop.labs[i]]]$RefFreq,ylab="",xlab="",main="",
       xlim=c(0,1),ylim=c(0,10000),axes=F,col=color)
  axis(1,pos=0)
  if(i %in% c(1,5,9,13)){
    axis(2,pos=0,las=1)
  }else{
    axis(2,pos=0,las=1,labels = FALSE)
    }
  mtext(pop.labs[i],3,col=color,cex=0.75,line=-1)
}
mtext("Reference Allele Frequency",1,outer=TRUE,cex=0.75)
mtext("Number of SNPs",2,outer = TRUE,cex=0.75,line=1)
dev.off()
#################################OUTLIERS####################################

#### STACKS ####
#Compare neighboring pops.
fwsw.tx<-read.delim("stacks/batch_2.fst_TXCC-TXFW.tsv")
fwsw.la<-read.delim("stacks/batch_2.fst_ALST-LAFW.tsv")
fwsw.al<-read.delim("stacks/batch_2.fst_ALFW-ALST.tsv")
fwsw.fl<-read.delim("stacks/batch_2.fst_FLCC-FLLG.tsv")

tx.sig<-fwsw.tx[fwsw.tx$Fisher.s.P<0.01,"Locus.ID"]
la.sig<-fwsw.la[fwsw.la$Fisher.s.P<0.01,"Locus.ID"]
al.sig<-fwsw.al[fwsw.al$Fisher.s.P<0.01,"Locus.ID"]
fl.sig<-fwsw.fl[fwsw.fl$Fisher.s.P<0.01,"Locus.ID"]
length(tx.sig[(tx.sig %in% c(la.sig,al.sig,fl.sig))])
length(la.sig[(la.sig %in% c(tx.sig,al.sig,fl.sig))])
length(al.sig[(al.sig %in% c(la.sig,tx.sig,fl.sig))])
all.shared<-fl.sig[fl.sig %in% la.sig & fl.sig %in% al.sig & fl.sig %in% tx.sig]
fw.shared.chr<-fwsw.tx[fwsw.tx$Locus.ID %in% all.shared,c("Locus.ID","Chr","BP","Column","Overall.Pi")]
tapply(fw.shared.chr$Locus.ID,factor(fw.shared.chr$Chr),function(x){ length(unique(x)) })
#are they using the same SNPs or different SNPs?
stacks.sig<-data.frame(nrow=nrow(fw.shared.chr),ncol=4)
for(i in 1:nrow(fw.shared.chr)){
  tx.bp<-fwsw.tx[fwsw.tx$Fisher.s.P<0.01 & fwsw.tx$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  la.bp<-fwsw.la[fwsw.la$Fisher.s.P<0.01 & fwsw.la$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  al.bp<-fwsw.al[fwsw.al$Fisher.s.P<0.01 & fwsw.al$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  fl.bp<-fwsw.fl[fwsw.fl$Fisher.s.P<0.01 & fwsw.fl$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  stacks.sig[i,1]<-paste(tx.bp,sep=",",collapse = ",")
  stacks.sig[i,2]<-paste(la.bp,sep=",",collapse = ",")
  stacks.sig[i,3]<-paste(al.bp,sep=",",collapse = ",")
  stacks.sig[i,4]<-paste(fl.bp,sep=",",collapse = ",")
}
colnames(stacks.sig)<-c("TX","LA","AL","FL")
stacks.sig<-data.frame(cbind(fw.shared.chr,stacks.sig))
stacks.sig$SNP<-paste(stacks.sig$Chr,stacks.sig$BP,sep=".")
write.table(stacks.sig,"stacks.sig.snps.txt",sep='\t',row.names=FALSE,col.names=TRUE)

##compare to scovelli genome..
gff<-read.delim(gzfile("../../scovelli_genome/ssc_2016_12_20_chromlevel.gff.gz"),header=F)
colnames(gff)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
genome.blast<-read.csv("../../scovelli_genome/ssc_2016_12_20_cds_nr_blast_results.csv",skip=1,header=T)#I saved it as a csv

#' extract the significant regions from the gff file
fw.sig.reg<-do.call(rbind,apply(fw.shared.chr,1,function(sig){
  this.gff<-gff[as.character(gff$seqname) %in% as.character(unlist(sig["Chr"])),]
  this.reg<-this.gff[this.gff$start <= as.numeric(sig["BP"]) & this.gff$end >= as.numeric(sig["BP"]),]
  if(nrow(this.reg) == 0){
    if(as.numeric(sig["BP"])>max(as.numeric(this.gff$end))){
      new<-data.frame(Locus=sig["Locus.ID"],Chr=sig["Chr"],BP=sig["BP"],SNPCol=sig["Column"],
                      region="beyond.last.contig", description=NA)
    }else{
      new<-data.frame(Locus=sig["Locus.ID"],Chr=sig["Chr"],BP=sig["BP"],SNPCol=sig["Column"],
                      region=NA,description=NA)
    }
  }else{
    if(length(grep("SSCG\\d+",this.reg$attribute))>0){
      geneID<-unique(gsub(".*(SSCG\\d+).*","\\1",this.reg$attribute[grep("SSCG\\d+",this.reg$attribute)]))
      gene<-genome.blast[genome.blast$sscv4_gene_ID %in% geneID,"blastp_hit_description"]
    }else{
      gene<-NA
    }
    new<-data.frame(Locus=sig["Locus.ID"],Chr=sig["Chr"],BP=sig["BP"],SNPCol=sig["Column"],
                    region=paste(this.reg$feature,sep=",",collapse = ","),description=gene)
  }
  return(as.data.frame(new))
}))
write.csv(fw.sig.reg,"StacksFWSWOutliers_annotatedByGenome.csv")

#' graph without highlighted regions
#' first define the scaffold boundaries
all.chr<-data.frame(Chr=c(as.character(fwsw.tx$Chr),as.character(fwsw.la$Chr),as.character(fwsw.al$Chr),as.character(fwsw.fl$Chr)),
                    BP=c(fwsw.tx$BP,fwsw.la$BP,fwsw.al$BP,fwsw.fl$BP),stringsAsFactors = F)
bounds<-data.frame(levels(as.factor(all.chr$Chr)),tapply(as.numeric(as.character(all.chr$BP)),all.chr$Chr,max))
plot.scaffs<-scaffs[scaffs %in% bounds]
colnames(bounds)<-c("Chrom","End")
png("stacks_fsts_fwsw.png",height=8,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(0.85,2,0,0.5),oma=c(1,1,1,0.5))
fwswt.fst<-fst.plot(fwsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, scaffold.widths = bounds,pch=19,y.lim = c(0,1),
                    pt.cols = c(grp.colors[1],grp.colors[2]),pt.cex=1,axis.size = 1)
mtext("TXFW vs. TXCC",3,cex=0.75,line=-1)
fwswl.fst<-fst.plot(fwsw.la,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, scaffold.widths = bounds,pch=19,y.lim = c(0,1),
                    pt.cols=c(grp.colors[2],grp.colors[3]),pt.cex=1,axis.size=1)
mtext("LAFW vs. ALST",3,cex=0.75,line=-1)
fwswa.fst<-fst.plot(fwsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, scaffold.widths = bounds,pch=19,y.lim = c(0,1),
                    pt.cols=c(grp.colors[3],grp.colors[2]),pt.cex=1,axis.size = 1)
mtext("ALFW vs. ALST",3,cex=0.75,line=-1)
fwswf.fst<-fst.plot(fwsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, scaffold.widths =  bounds,pch=19,y.lim = c(0,1),
                    pt.cols=c(grp.colors[6],grp.colors[5]),pt.cex=1,axis.size=1)
mtext("FLFW vs. FLCC",3,cex=0.75,line=-1)
mtext(expression(italic(F)[ST]),2,outer=T,line=-0.8,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(fwswf.fst[fwswf.fst$Chr ==lgs[i],"plot.pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE)
  last<-max(fwswf.fst[fwswf.fst$Chr ==lgs[i],"plot.pos"])
}
dev.off()

#' plot with the outlier regions
png("stacks_fsts_fwsw_withSig.png",height=8,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(0.85,2,0,0.5),oma=c(1,1,1,0.5))
fwswt.fst<-fst.plot(fwsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, y.lim = c(0,1),scaffold.widths = bounds,pch=19,
                    pt.cols = c(grp.colors[1],grp.colors[2]),pt.cex=1,axis.size = 1)
#points(fwswt.fst$BP,fwswt.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[1])
points(fwswt.fst$plot.pos[fwswt.fst$Locus.ID %in% all.shared],fwswt.fst$Corrected.AMOVA.Fst[fwswt.fst$Locus.ID %in% all.shared],
       pch=1,col="black",cex=1.3)
mtext("TXFW vs. TXCC",3,cex=0.75,line=-1)
fwswl.fst<-fst.plot(fwsw.la,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, y.lim = c(0,1),scaffold.widths = bounds,pch=19,
                    pt.cols=c(grp.colors[2],grp.colors[3]),pt.cex=1,axis.size=1)
#points(fwswl.fst$BP,fwswl.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[3])
points(fwswl.fst$plot.pos[fwswl.fst$Locus.ID %in% all.shared],fwswl.fst$Corrected.AMOVA.Fst[fwswl.fst$Locus.ID %in% all.shared],
       pch=1,col="black",cex=1.3)
mtext("LAFW vs. ALST",3,cex=0.75,line=-1)
fwswa.fst<-fst.plot(fwsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, y.lim = c(0,1),scaffold.widths = bounds,pch=19,
                    pt.cols=c(grp.colors[3],grp.colors[2]),pt.cex=1,axis.size = 1)
#points(fwswa.fst$BP,fwswa.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[4])
points(fwswa.fst$plot.pos[fwswa.fst$Locus.ID %in% all.shared],fwswa.fst$Corrected.AMOVA.Fst[fwswa.fst$Locus.ID %in% all.shared],
       pch=1,col="black",cex=1.3)
mtext("ALFW vs. ALST",3,cex=0.75,line=-1)
fwswf.fst<-fst.plot(fwsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, y.lim = c(0,1),scaffold.widths = bounds,pch=19,
                    pt.cols=c(grp.colors[6],grp.colors[5]),pt.cex=1,axis.size=1)
#points(fwswf.fst$BP,fwswf.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[6])
points(fwswf.fst$plot.pos[fwswf.fst$Locus.ID %in% all.shared],fwswf.fst$Corrected.AMOVA.Fst[fwswf.fst$Locus.ID %in% all.shared],
       pch=1,col="black",cex=1.3)
mtext("FLFW vs. FLCC",3,cex=0.75,line=-1)
mtext(expression(italic(F)[ST]),2,outer=T,line=-1,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(fwswf.fst[fwswf.fst$Chr ==lgs[i],"plot.pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE)
  last<-max(fwswf.fst[fwswf.fst$Chr ==lgs[i],"plot.pos"])
}
dev.off()

##For comparison, neighboring sw pops
swsw.tx<-read.delim("stacks/batch_2.fst_TXCB-TXCC.tsv")
swsw.al<-read.delim("stacks/batch_2.fst_ALST-FLSG.tsv")
swsw.fl<-read.delim("stacks/batch_2.fst_FLCC-FLHB.tsv")

tx.sw.sig<-swsw.tx[swsw.tx$Fisher.s.P<0.01,"Locus.ID"]
al.sw.sig<-swsw.al[swsw.al$Fisher.s.P<0.01,"Locus.ID"]
fl.sw.sig<-swsw.fl[swsw.fl$Fisher.s.P<0.01,"Locus.ID"]
length(tx.sw.sig[(tx.sw.sig %in% c(al.sw.sig,fl.sw.sig))])
length(al.sw.sig[(al.sw.sig %in% c(tx.sw.sig,fl.sw.sig))])
length(fl.sw.sig[(fl.sw.sig %in% c(tx.sw.sig,al.sw.sig))])
sw.shared<-fl.sw.sig[fl.sw.sig %in% al.sw.sig & fl.sw.sig %in% tx.sw.sig]
fw.shared.chr<-fwsw.tx[fwsw.tx$Locus.ID %in% all.shared,c("Locus.ID","Chr","BP","Column")]
tapply(fw.shared.chr$Locus.ID,factor(fw.shared.chr$Chr),function(x){ length(unique(x)) })

swswt.fst<-fst.plot(swsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot = scaffs,scaffold.widths = bounds,pt.cex = 1,pch=19)
swswa.fst<-fst.plot(swsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot = scaffs,scaffold.widths = bounds,pt.cex = 1,pch=19)
swswf.fst<-fst.plot(swsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot = scaffs,scaffold.widths = bounds,pt.cex = 1,pch=19)

png("stacks_fsts_swsw.png",height=6,width=7.5,units="in",res=300)
par(mfrow=c(3,1),mar=c(0.85,2,0,0.5),oma=c(1,1,1,0.5))
swswt.fst<-fst.plot(swsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, y.lim = c(0,1),scaffold.widths = bounds,axis.size = 1,
                    pt.cols=c(grp.colors[1],grp.colors[2]),pt.cex = 1,pch=19)
mtext("TXCC vs. TXCB",3,line=-1,cex=0.75)
swswa.fst<-fst.plot(swsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", y.lim=c(0,1),
                    scaffs.to.plot=plot.scaffs, scaffold.widths = bounds,axis.size = 1,
                    pt.cols=c(grp.colors[2],grp.colors[3]),pt.cex = 1,pch=19)
#points(swswa.fst$plot.pos,swswa.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[3])
mtext("ALST vs. FLSG",3,line=-1,cex=0.75)
swswf.fst<-fst.plot(swsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    scaffs.to.plot=plot.scaffs, y.lim = c(0,1),scaffold.widths = bounds,axis.size = 1,
                    pt.cols=c(grp.colors[6],grp.colors[5]),pt.cex = 1,pch=19)
#points(swswf.fst$plot.pos,swswf.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[6])
mtext("FLCC vs. FLHB",3,line=-1,cex=0.75)
mtext(expression(italic(F)[ST]),2,outer=T,line=-0.75,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(swswf.fst[swswf.fst$Chr ==lgs[i],"plot.pos"]),y=-0.07,
       labels=lgn[i], adj=1, xpd=TRUE)
  last<-max(swswf.fst[swswf.fst$Chr ==lgs[i],"plot.pos"])
}
dev.off()

#' What are the allele frequencies of the significant FWSW outliers in FW pops and in SW pops?
#' And what are the SNPs?
#get the RAD loci from the vcf
stacks.sig$SNP<-paste(stacks.sig$Chr,as.character(as.numeric(as.numeric(stacks.sig$BP) + 1)),sep=".")
stacks.sig.vcf<-vcf[vcf$SNP %in% stacks.sig$SNP,] #two don't get added, but I can add them manually
stacks.sig.vcf<-rbind(stacks.sig.vcf,vcf[vcf$ID == 28200,])
fw.sig.vcf<-stacks.sig.vcf[,c(locus.info,unlist(lapply(fw.list,grep,x=colnames(stacks.sig.vcf),value=TRUE)))]
sw.sig.vcf<-stacks.sig.vcf[,c(locus.info,unlist(lapply(sw.list,grep,x=colnames(stacks.sig.vcf),value=TRUE)))]
fw.sig.afs<-do.call(rbind,apply(fw.sig.vcf,1,calc.afs.vcf))
fw.sig.afs$ID<-fw.sig.vcf$ID
fw.sig.afs$SNP<-fw.sig.vcf$SNP
sw.sig.afs<-do.call(rbind,apply(sw.sig.vcf,1,calc.afs.vcf))
sw.sig.afs$ID<-sw.sig.vcf$ID
sw.sig.afs$SNP<-fw.sig.vcf$SNP

#reorder to be in chrom order
fw.sig.afs<-fw.sig.afs[order(factor(fw.sig.afs$Chrom,levels=scaffs)),]
sw.sig.afs<-sw.sig.afs[order(factor(sw.sig.afs$Chrom,levels=scaffs)),]

png("OutlierAlleleFreqs.png",height=5,width=10,units="in",res=300)
par(mar=c(3,2,1,1),oma=c(2,2,0,1))
plot(1:nrow(fw.sig.afs),fw.sig.afs$RefFreq,ylim=c(0,1),type='n',
    xlab="",ylab="",axes=FALSE)
axis(2,ylim=c(0,1),pos=0,las=1,cex.axis=0.75)
mtext("Reference Allele Frequency",2,line=1.5,cex=0.75)
start<-0.1
for(i in 1:length(unique(fw.sig.afs$Chrom))){
  end<-start+nrow(fw.sig.afs[fw.sig.afs$Chrom %in% unique(fw.sig.afs$Chrom)[i],])+0.1
  if(i%%2){
    rect(start,0,end,1,col="lightgrey",border="lightgrey")
  }
  text(x=mean(c(start,end)),y=0,labels=unique(fw.sig.afs$Chrom)[i],
       xpd=TRUE,adj=1,srt=45,cex=0.75)
  start<-end
}
#text(x=1:nrow(fw.sig.afs),y=rep(-0.05,nrow(fw.sig.afs)),
#     labels=fw.sig.afs$Chrom,xpd=TRUE,adj=1,srt=90)
points(1:nrow(fw.sig.afs),fw.sig.afs$RefFreq, 
       pch=19,col="cornflowerblue")
points(1:nrow(fw.sig.afs),sw.sig.afs$RefFreq)
arrows(x0=1:nrow(fw.sig.afs),x1=1:nrow(fw.sig.afs),
       y0=sw.sig.afs$RefFreq,y1=fw.sig.afs$RefFreq,
       angle=45,length=0.1)
legend(x=50,y=0.15,bty='n',pch=c(19,1),col=c("cornflowerblue","black"),
       c("Freshwater Populations","Saltwater Populations"),cex=0.75)
dev.off()
####### Using gwscaR code #####
loci.info<-c(colnames(vcf[1:9]),"SNP")
txcb<-grep("TXCB",colnames(vcf),value = T)
txfw<-grep("TXFW",colnames(vcf),value = T)
alst<-grep("ALST",colnames(vcf),value = T)
alfw<-grep("ALFW",colnames(vcf),value = T)
lafw<-grep("LAFW",colnames(vcf),value = T)
flcc<-grep("FLCC",colnames(vcf),value = T)
fllg<-grep("FLLG",colnames(vcf),value = T)
txcc<-grep("TXCC",colnames(vcf),value = T)
flsg<-grep("FLSG",colnames(vcf),value = T)
flhb<-grep("FLHB",colnames(vcf),value = T)

tfs<-gwsca(vcf=vcf,locus.info=loci.info,group1=txcb,group2=txfw)
tfs$SNP<-paste(tfs$Chrom,as.numeric(as.character(tfs$Pos)),sep=".")
lfs<-gwsca(vcf=vcf,locus.info=loci.info,group1=alst,group2=lafw)
lfs$SNP<-paste(lfs$Chrom,as.numeric(as.character(lfs$Pos)),sep=".")
afs<-gwsca(vcf=vcf,locus.info=loci.info,group1=alst,group2=alfw)
afs$SNP<-paste(afs$Chrom,as.numeric(as.character(afs$Pos)),sep=".")
ffs<-gwsca(vcf=vcf,locus.info=loci.info,group1=flcc,group2=fllg)
ffs$SNP<-paste(ffs$Chrom,as.numeric(as.character(ffs$Pos)),sep=".")

#' Shared fst outliers from gwsca
tfs.sig<-tfs[tfs$Chi.p <= 0.05,"SNP"]
lfs.sig<-lfs[lfs$Chi.p <= 0.05,"SNP"]
afs.sig<-afs[afs$Chi.p <= 0.05,"SNP"]
ffs.sig<-ffs[ffs$Chi.p <= 0.05,"SNP"]
shared.sig<-tfs.sig[tfs.sig %in% lfs.sig & tfs.sig %in% afs.sig & tfs.sig %in% ffs.sig]

tss<-gwsca(vcf=vcf,locus.info=loci.info,group1=txcb,group2=txcc)
ass<-gwsca(vcf=vcf,locus.info=loci.info,group1=alst,group2=flsg)
fss<-gwsca(vcf=vcf,locus.info=loci.info,group1=flcc,group2=flhb)



#### Delta-divergence ####
#' only use chosen SNPs
vcf<-vcf[vcf$SNP %in% chosen.snps,]
#Calculate delta-divergence
#' ``` {r, eval = FALSE}
#' swfw.mu<-calc.mean.fst(vcf = vcf,pop.list1 = sw.list,pop.list2 = fw.list,maf.cutoff=0.01)
#' fwfw.mu<-calc.mean.fst(vcf = vcf,pop.list1 = fw.list,pop.list2 = fw.list, maf.cutoff=0.01)
#' deltad<-merge(swfw.mu,fwfw.mu,by="SNP")
#' deltad<-deltad[,c("SNP","Chrom.x","Pos.x","Mean.Fst.x","Mean.Fst.y")]
#' colnames(deltad)<-c("SNP","Chrom","Pos","MeanSWFW.Fst","MeanFWFW.Fst")
#' deltad$deltad<-deltad$MeanSWFW.Fst - deltad$MeanFWFW.Fst
#' deltad<-deltad[!is.na(deltad$deltad),]#remove NAs
#' ```
#' ```{r, echo=FALSE}
#' deltad<-read.delim("deltadivergence.txt")
#' ```
#' Plot it
png("delta-divergence.png",height=5,width=7,units="in",res=300)
par(mar=c(1,1,1,1),oma=c(1,2,1,1))
dd<-fst.plot(fst.dat = deltad,fst.name = "deltad",bp.name = "Pos",axis=1,pch=19,scaffs.to.plot=scaffs)
mtext(expression(paste(delta,"-divergence")),2,line=1.5)
smooth.par<-data.frame()
smooth.div<-data.frame()
for(i in 1:length(lgs)){#scaffolds are too short
  this.chrom<-dd[dd$Chrom %in% lgs[i],]
  #span<-nrow(this.chrom)/5000
  this.smooth<-loess.smooth(this.chrom$plot.pos,this.chrom$deltad,span=0.1,degree=2) 
  this.plot.pos<-unlist(lapply(this.smooth$x,function(x) { this.chrom[which.min(abs(this.chrom$plot.pos-x)),"plot.pos"] }))
  points(this.smooth$x,this.smooth$y,col="cornflowerblue",type="l",lwd=2)
  # this.div<-cbind(x= this.smooth$x[this.smooth$y>=quantile(this.smooth$y,0.95)],#choosing the upper outliers - smoothed
  #                 smooth.dd=this.smooth$y[this.smooth$y>=quantile(this.smooth$y,0.95)],
  #                 plot.pos=this.plot.pos[this.smooth$y>=quantile(this.smooth$y,0.95)])
  # this.par<-cbind(x=this.smooth$x[this.smooth$y<=quantile(this.smooth$y,0.05)],#choosing the lower outliers - smoothed
  #                 smooth.dd=this.smooth$y[this.smooth$y<=quantile(this.smooth$y,0.05)],
  #                 plot.pos=this.plot.pos[this.smooth$y<=quantile(this.smooth$y,0.05)])
  this.div<-cbind(dd=this.chrom$deltad[this.chrom$deltad>=quantile(this.chrom$deltad,0.95)],#choosing the upper outliers
                  plot.pos=this.chrom$plot.pos[this.chrom$deltad>=quantile(this.chrom$deltad,0.95)])
  this.par<-cbind(dd=this.chrom$deltad[this.chrom$deltad<=quantile(this.chrom$deltad,0.05)],#choosing the lower outliers
                  plot.pos=this.chrom$plot.pos[this.chrom$deltad<=quantile(this.chrom$deltad,0.05)])
  smooth.par<-rbind(smooth.par,this.par)
  smooth.div<-rbind(smooth.div,this.div)
}
dev.off()
smooth.div<-merge(smooth.div,dd,by="plot.pos") 
smooth.par<-merge(smooth.par,dd,by="plot.pos")
smooth.par<-smooth.par[!duplicated(smooth.par$SNP),]#a couple of SNPs had multiple smoothed delat d values
smooth.par$FWpi<-apply(cbind(vcf[vcf$SNP %in% smooth.par$SNP,locus.info],
                             vcf[vcf$SNP %in% smooth.par$SNP,unlist(lapply(fw.list,grep,colnames(vcf)))]),1,calc.pi)
smooth.par$SWpi<-apply(cbind(vcf[vcf$SNP %in% smooth.par$SNP,locus.info],
                             vcf[vcf$SNP %in% smooth.par$SNP,unlist(lapply(sw.list,grep,colnames(vcf)))]),1,calc.pi)
smooth.div$FWpi<-apply(cbind(vcf[vcf$SNP %in% smooth.div$SNP,locus.info],
                             vcf[vcf$SNP %in% smooth.div$SNP,unlist(lapply(fw.list,grep,colnames(vcf)))]),1,calc.pi)
smooth.div$SWpi<-apply(cbind(vcf[vcf$SNP %in% smooth.div$SNP,locus.info],
                             vcf[vcf$SNP %in% smooth.div$SNP,unlist(lapply(sw.list,grep,colnames(vcf)))]),1,calc.pi)
#points(smooth.out$plot.pos,smooth.out$smooth.deltad,col="orchid4")
dim(smooth.par[smooth.par$SNP %in% stacks.sig$SNP,])
write.table(rbind(smooth.par,smooth.div),"smoothed.deltad.out.txt",col.names=T,row.names=F,quote=F,sep='\t')

#' Compare to Stacks Fsts
# Get the significant Fst loci if they're not already here
if(!("stacks.sig" %in% ls())){
  stacks.sig<-read.delim("stacks.sig.snps.txt")
} 
if(is.null(stacks.sig$SNP)){ #make sure the SNP column is there.
  stacks.sig$SNP<-paste(stacks.sig$Chr,as.numeric(stacks.sig$BP)+1,sep=".")
}

#what if I use deltad quantiles, not smoothed quantiles?
dd.div<-dd[dd$deltad >= quantile(dd$deltad,0.95),]
dd.par<-dd[dd$deltad <= quantile(dd$deltad,0.05),]


#### Other Pop gen statistics ####
#pi
avg.pi<-do.call("rbind",sliding.window(vcf,scaffs))
avg.pi.adj<-fst.plot(avg.pi,scaffold.widths=scaff.starts,pch=19,
                    fst.name = "Avg.Pi",chrom.name = "Chr",bp.name = "Avg.Pos")
all.pi<-data.frame(Chrom=vcf$`#CHROM`,Pos=vcf$POS,Pi=unlist(apply(vcf,1,calc.pi)))
all.pi$SNP<-paste(all.pi$Chrom,as.numeric(as.character(all.pi$Pos)),sep=".")
write.table(all.pi,"all.pi.txt",col.names = TRUE,row.names=FALSE, quote=FALSE,sep='\t')
all.pi<-read.table("all.pi.txt",header=T)
#rho
avg.rho<-do.call("rbind",sliding.window(vcf,scaffs,stat = "rho",pop.list=pop.list))
avg.rho.adj<-fst.plot(avg.rho,scaffold.widths=scaff.starts,
                     fst.name = "Avg.Pi",chrom.name = "Chr",bp.name = "Avg.Pos")
all.rho<-data.frame(Chrom=vcf$`#CHROM`,Pos=vcf$POS,Rho=unlist(apply(vcf,1,calc.rho,pop.list=pop.list)))
all.rho$SNP<-paste(all.rho$Chrom,as.numeric(as.character(all.rho$Pos)),sep=".")

#plot
png("FWSWpi.png",height=5,width=7,units="in",res=300)
par(oma=c(2,2,2,2),mar=c(2,2,2,2))
pi.plot<-fst.plot(all.pi,scaffold.widths=scaff.starts,y.lim=c(0,0.5),axis.size = 0.5,pch=19,
                     fst.name = "Pi",chrom.name = "Chrom",bp.name = "Pos",pch=19)
points(x=avg.pi.adj$plot.pos,y=avg.pi.adj$Avg.Pi,col="cornflowerblue",type="l",lwd=2)
dev.off()

png("deltad_pi.png",height=7,width=5,units="in",res=300)
par(mfrow=c(2,1),mar=c(1,1,1,1),oma=c(1,2,1,1),cex=0.75)
dd<-fst.plot(fst.dat = deltad,fst.name = "deltad",bp.name = "Pos",axis=0.75,
             y.lim=c(-0.5,1),scaffold.widths=scaff.starts,pch=19,scaffs.to.plot = scaffs)
mtext(expression(paste(delta,"-divergence")),2,line=2,cex=0.75)
smooth.out<-data.frame()
smoothed.dd<-data.frame()
for(i in 1:length(lgs)){#scaffolds are too short
  this.chrom<-dd[dd$Chrom %in% lgs[i],]
  this.smooth<-loess.smooth(this.chrom$plot.pos,this.chrom$deltad,span=0.1,degree=2) 
  these.sig<-stacks.sig[stacks.sig$Chr %in% lgs[i],]
  this.pos<-unlist(lapply(these.sig$BP,function(x) { this.chrom[which.min(abs(this.chrom$Pos-x)),"Pos"] }))
  points(this.smooth$x,this.smooth$y,col="cornflowerblue",type="l",lwd=2)
  points(dd[dd$Pos %in% this.pos,c("plot.pos","deltad")],col="darkorchid")
  smoothed.dd<-rbind(smoothed.dd,data.frame(Chrom=rep(lgs[i],length(this.smooth$x)),x=this.smooth$x,y=this.smooth$y))
}
points(dd[dd$SNP %in% stacks.sig$SNP,c("plot.pos","deltad")],col="darkorchid")
pi.plot<-fst.plot(all.pi,scaffold.widths=scaff.starts,y.lim=c(0,0.5),axis.size = 0.75,
                  scaffs.to.plot = scaffs,fst.name = "Pi",chrom.name = "Chrom",bp.name = "Pos",pch=19)
points(x=avg.pi.adj$plot.pos[avg.pi.adj$Chr %in% lgs],y=avg.pi.adj$Avg.Pi[avg.pi.adj$Chr %in% lgs],
       col="cornflowerblue",type="l",lwd=2)
points(pi.plot[pi.plot$SNP %in% stacks.sig$SNP,c("plot.pos","Pi")],col="darkorchid")
mtext(expression(pi),2,line=2,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(pi.plot[pi.plot$Chrom ==lgs[i],"plot.pos"]),y=0,
       labels=lgn[i], adj=1, xpd=TRUE,cex=0.75)
  last<-max(pi.plot[pi.plot$Chrom ==lgs[i],"plot.pos"])
}
dev.off()

#### Looking for directional selection
#' High Fst, low pi
pi.sig.fst<-all.pi[all.pi$SNP %in% stacks.sig$SNP,] #70 of 72
points(pi.plot[pi.plot$SNP %in% stacks.sig$SNP,"plot.pos"],pi.plot[pi.plot$SNP %in% stacks.sig$SNP,"Pi"],col="orchid4")

#add snps
fwsw.tx$SNP<-paste(fwsw.tx$Chr,fwsw.tx$BP+1,sep=".")
fwsw.la$SNP<-paste(fwsw.la$Chr,fwsw.la$BP+1,sep=".")
fwsw.al$SNP<-paste(fwsw.al$Chr,fwsw.al$BP+1,sep=".")
fwsw.fl$SNP<-paste(fwsw.fl$Chr,fwsw.fl$BP+1,sep=".")
#plot each chrom with sig. ones
sig.chroms<-unique(pi.sig.fst$Chrom)
png("FstOutliersStats.png",height=10,width=10,units="in",res=300)
par(mfrow=c(4,4),mar=c(1,1,1,1),oma=c(1,2.5,1,1))
for(i in 1:length(sig.chroms)){ #I could turn this into a function
  if(sig.chroms[i] %in% lgs){
    this.tx<-fwsw.tx[fwsw.tx$Chr %in% sig.chroms[i],]
    tx.smooth<-do.call("rbind",lapply(seq(1,nrow(this.tx),(nrow(this.tx)*0.15)/5),sliding.avg,
                                      dat=data.frame(Pos=this.tx$BP,Fst=this.tx$Corrected.AMOVA.Fst),width=nrow(this.tx)*0.15))
      #loess.smooth(this.tx$BP,this.tx$Corrected.AMOVA.Fst,span=0.1,degree=2) 
    this.la<-fwsw.la[fwsw.la$Chr %in% sig.chroms[i],]
    la.smooth<-do.call("rbind",lapply(seq(1,nrow(this.la),(nrow(this.la)*0.15)/5),sliding.avg,
                                      dat=data.frame(Pos=this.la$BP,Fst=this.la$Corrected.AMOVA.Fst),width=nrow(this.la)*0.15))
      #loess.smooth(this.la$BP,this.la$Corrected.AMOVA.Fst,span=0.1,degree=2) 
    this.al<-fwsw.al[fwsw.al$Chr %in% sig.chroms[i],]
    al.smooth<-do.call("rbind",lapply(seq(1,nrow(this.al),(nrow(this.al)*0.15)/5),sliding.avg,
                                      dat=data.frame(Pos=this.al$BP,Fst=this.al$Corrected.AMOVA.Fst),width=nrow(this.al)*0.15))
      #loess.smooth(this.al$BP,this.al$Corrected.AMOVA.Fst,span=0.1,degree=2) 
    this.fl<-fwsw.fl[fwsw.fl$Chr %in% sig.chroms[i],]
    fl.smooth<-do.call("rbind",lapply(seq(1,nrow(this.fl),(nrow(this.fl)*0.15)/5),sliding.avg,
                                      dat=data.frame(Pos=this.fl$BP,Fst=this.fl$Corrected.AMOVA.Fst),width=nrow(this.fl)*0.15))
      #loess.smooth(this.fl$BP,this.fl$Corrected.AMOVA.Fst,span=0.1,degree=2) 
    plot(tx.smooth,col=alpha(grp.colors[1],0.5),type="l",ylim=c(0,1),lwd=2,
         bty="L",xlab="",ylab="",xaxt='n',yaxt='n')
    if(i %in% c(1,5,9,13)){
      axis(2,las=1)  
    }else{
      axis(2,labels=FALSE) }
    mtext(paste("Position on ",sig.chroms[i],sep=""),1,cex=0.75)
    points(la.smooth,col=alpha(grp.colors[2],0.5),type="l",lwd=2)
    points(al.smooth,col=alpha(grp.colors[3],0.5),type="l",lwd=2)
    points(fl.smooth,col=alpha(grp.colors[6],0.5),type="l",lwd=2)
    lapply(pi.sig.fst[pi.sig.fst$Chrom %in% sig.chroms[i],"Pos"],points,y=0.6,pch="-",cex=5)#add bars to sig. fsts
    #add delta divergence
    this.chrom<-dd[dd$Chrom %in% sig.chroms[i],]
    this.smooth<-loess.smooth(this.chrom$Pos,this.chrom$deltad,span=0.1,degree=2)
    points(this.smooth$x,this.smooth$y,col="cornflowerblue",type="l",lwd=2)
    lapply(smooth.out[smooth.out$Chrom %in% sig.chroms[i],"Pos"],points,y=0.5,pch="-",cex=5,col="cornflowerblue")
    #add pi
    fw.inds<-unlist(lapply(fw.list,grep,colnames(vcf)))
    sw.inds<-unlist(lapply(sw.list,grep,colnames(vcf)))
    nloci<-nrow(vcf[vcf$`#CHROM` %in% sig.chroms[i],])
    fw.pi<-do.call("rbind",
                   sliding.window(cbind(vcf[,locus.info],vcf[,fw.inds]),
                                  sig.chroms[i],
                                  nsteps=round((nloci*0.15)/5),
                                  width=nloci*0.15))
    sw.pi<-do.call("rbind",
                   sliding.window(cbind(vcf[,locus.info],vcf[,sw.inds]),
                                  sig.chroms[i],
                                  nsteps=round((nloci*0.15)/5),
                                  width=nloci*0.15))
    points(fw.pi,col="lightgrey",type="l",lwd=3)
    points(sw.pi,col="darkgrey",type="l",lwd=3,lty=2)
  }
  # else{
  #   this.tx<-fwsw.tx[fwsw.tx$Chr %in% sig.chroms[i],]
  #   this.la<-fwsw.la[fwsw.la$Chr %in% sig.chroms[i],]
  #   this.al<-fwsw.al[fwsw.al$Chr %in% sig.chroms[i],]
  #   this.fl<-fwsw.fl[fwsw.fl$Chr %in% sig.chroms[i],]
  #   plot(this.tx$BP,this.tx$Corrected.AMOVA.Fst,col=grp.colors[1],type="l",ylim=c(0,1),lwd=2,xaxt='n',
  #        xlab="",ylab="",xaxt='n',yaxt='n',bty="L")
  #   axis(2, labels=FALSE)
  #   mtext(paste("Position on ",sig.chroms[i],sep=""),1,cex=0.75)
  #   points(this.la$BP,this.la$Corrected.AMOVA.Fst,col=grp.colors[2],type="l",lwd=2)
  #   points(this.al$BP,this.al$Corrected.AMOVA.Fst,col=grp.colors[3],type="l",lwd=2)
  #   points(this.fl$BP,this.fl$Corrected.AMOVA.Fst,col=grp.colors[6],type="l",lwd=2)
  #   #add pi
  #   points(avg.pi$Avg.Pos[avg.pi$Chr %in% sig.chroms[i]],avg.pi$Avg.Pi[avg.pi$Chr %in% sig.chroms[i]],
  #          col="darkgrey",type="l",lwd=2)
  #   lapply(pi.sig.fst[pi.sig.fst$Chrom %in% sig.chroms[i],"Pos"],points,y=0.6,pch="-",cex=5)#add bars to sig. fsts
  #   lapply(smooth.out[smooth.out$Chrom %in% sig.chroms[i],"Pos"],points,y=0.5,pch="-",cex=5,col="cornflowerblue")
  # }
}#this is kind of a mess.
plot(x=c(0,1),y=c(0,1),type="n",axes=FALSE)
legend("left",
  legend=c(expression(TX~FWSW~italic(F)[ST]),
                expression(LA~FWSW~italic(F)[ST]),
         expression(AL~FWSW~italic(F)[ST]),expression(FL~FWSW~italic(F)[ST]),
         expression(delta~-divergence),expression(Average~pi[FW]),
         expression(Average~pi[SW])),bty='n',lwd=2,lty=c(1,1,1,1,1,1,2),
       col=c(alpha(grp.colors[1],0.5),alpha(grp.colors[2],0.5),
             alpha(grp.colors[3],0.5),alpha(grp.colors[6],0.5),
             "cornflowerblue","lightgrey","darkgrey"))
mtext(expression(italic(F)[ST]),2,cex=0.75,outer=T,line=1)
dev.off()


for(i in 1:nrow(stacks.sig)){
  png(paste(stacks.sig$SNP[i],"png",sep="."),height=5,width=7,units="in",res=300)
  this.tx<-fwsw.tx[fwsw.tx$Chr %in% stacks.sig$Chr[i] & 
                     fwsw.tx$BP >= stacks.sig$BP[i]-2500 & fwsw.tx$BP >= stacks.sig$BP[i]-2500,]
  #tx.smooth<-loess.smooth(this.tx$BP,this.tx$Corrected.AMOVA.Fst,span=0.1,degree=2) 
  this.al<-fwsw.al[fwsw.al$Chr %in% stacks.sig$Chr[i] & 
                     fwsw.al$BP >= stacks.sig$BP[i]-2500 & fwsw.al$BP >= stacks.sig$BP[i]-2500,]
#  al.smooth<-loess.smooth(this.al$BP,this.al$Corrected.AMOVA.Fst,span=0.1,degree=2)
  this.la<-fwsw.la[fwsw.la$Chr %in% stacks.sig$Chr[i] & 
                     fwsw.la$BP >= stacks.sig$BP[i]-2500 & fwsw.la$BP >= stacks.sig$BP[i]-2500,]
  #la.smooth<-loess.smooth(this.la$BP,this.la$Corrected.AMOVA.Fst,span=0.1,degree=2) 
  this.fl<-fwsw.fl[fwsw.fl$Chr %in% stacks.sig$Chr[i] & 
                     fwsw.fl$BP >= stacks.sig$BP[i]-2500 & fwsw.fl$BP >= stacks.sig$BP[i]-2500,]
 # fl.smooth<-loess.smooth(this.fl$BP,this.fl$Corrected.AMOVA.Fst,span=0.1,degree=2)
  
  plot(this.tx$BP,this.tx$Corrected.AMOVA.Fst,col=alpha(grp.colors[1],0.5),type="l",ylim=c(0,1),lwd=2,xaxt='n',
      xlab=paste("Position on ",stacks.sig$Chr[i],sep=""),ylab="",bty="L")
  # plot(tx.smooth$x,tx.smooth$y,col=alpha(grp.colors[1],0.5),type="l",
  #      ylim=c(0,1),
  #      lwd=2,xaxt='n',
  #           xlab=paste("Position on ",sig.chroms[i],sep=""),ylab="",bty="L")
  mtext(expression(italic(F)[ST]),2,line=2,cex=0.75)
  points(this.la$BP,this.la$Corrected.AMOVA.Fst,col=alpha(grp.colors[2],0.5),type="l",lwd=2)
  points(this.al$BP,this.al$Corrected.AMOVA.Fst,col=alpha(grp.colors[3],0.5),type="l",lwd=2)
  points(this.fl$BP,this.fl$Corrected.AMOVA.Fst,col=alpha(grp.colors[6],0.5),type="l",lwd=2)
  # points(la.smooth$x,la.smooth$y,col=alpha(grp.colors[2],0.5),type="l",lwd=2)
  # points(al.smooth$x,al.smooth$y,col=alpha(grp.colors[3],0.5),type="l",lwd=2)
  # points(fl.smooth$x,fl.smooth$y,col=alpha(grp.colors[6],0.5),type="l",lwd=2)
  if(stacks.sig$Chr[i] %in% lgs){
    this.chrom<-dd[dd$Chrom %in% stacks.sig$Chr[i],]
    #span<-nrow(this.chrom)/5000
    this.smooth<-loess.smooth(this.chrom$Pos,this.chrom$deltad,span=0.1,degree=2)
   # this.smooth<-smoothed.dd[smoothed.dd$Chrom %in% stacks.sig$Chr[i] &
    #                           smoothed.dd$x > min(this.tx$BP) & smoothed.dd$x < max(this.tx$BP),]
    points(this.smooth$x,this.smooth$y,col="cornflowerblue",type="l",lwd=2)
  }
  #add pi
  #points(avg.pi$Avg.Pos[avg.pi$Chr %in% sig.chroms[i]],avg.pi$Avg.Pi[avg.pi$Chr %in% sig.chroms[i]],
  #       col="darkgrey",type="l",lwd=2)
  
  lapply(stacks.sig[stacks.sig$Chr %in% stacks.sig[i,"Chr"],"BP"],points,y=0.6,pch="-",cex=5,col="darkgrey")#add bars to sig. fsts
  points(stacks.sig[i,"BP"],0.6,cex=5,pch="-")
  lapply(smooth.out[smooth.out$Chrom %in% stacks.sig$Chr[i],"Pos"],points,y=0.5,pch="-",cex=5,col="cornflowerblue")
  dev.off()
}

#which points are in the lower 25%?
put.dir<-pi.sig.fst[pi.sig.fst$Pi <= quantile(pi.plot$Pi,0.25),] #putative directional selection loci.



#### Looking for balancing selection
#' Low Fst, high pi (and het)

#plot fsts


##### Neighbor joining trees #####
#' to get trees and calc gsi (maybe):
#' for each overlapping sliding window (of 33 SNPs, for example), 
#' generate distance matrix (Fsts) using those SNPs
#' ape::nj(as.dist(matrix), "unrooted")
#' genealogicalSorting::gsi(tree, class, assignments, uncertainty)
#' http://molecularevolution.org/software/phylogenetics/gsi/download

library(ape)

fst.trees<-list()
for(vcf.row in 1: nrow(vcf)){
  fst.trees<-c(fst.trees,get.dist(vcf[vcf.row,],pop.list)) #getting an error
}
#Now what? How do I get rid of the numbers in the trees?
ftt<-lapply(fst.trees,gsub,pattern="\\d+",replacement="")
ftt<-lapply(ftt,gsub,pattern="[:.-]",replacement="")
ftt<-lapply(ftt,gsub,pattern="e",replacement="")
ftrees<-do.call(rbind,ftt)
ftrees<-data.frame(Chrom=vcf$`#CHROM`,Pos=vcf$POS,SNP=vcf$SNP,
                   trees=ftrees[,1])
write.table(ftrees,"ftrees.txt",row.names=F,col.names=T,quote=F)
library(combinat)
grps<-permn(fw.list)
fw.grps<-unlist(lapply(grps,function(grp){
  grp.tr<-paste("(((",grp[1],",",grp[2],"),",grp[3],"),",grp[4],")",sep="")
}))

fwtrees<-data.frame() #try to find the ones that match
apply(ftrees,1,function(snp){
  if(length(unlist(lapply(fw.grps,grep,snp["trees"])))>0){
    fwtrees<-rbind(fwtrees,snp)
  }
})

##### TREEMIX #####

tm.fwsw<-treemix.from.vcf(vcf,pop.list)
write.table(tm.fwsw,"fwsw.treemix",col.names=TRUE,row.names=FALSE,quote=F,sep=' ')
#then in unix: gzip -c fwsw.treemix > fwsw.treemix.gz

#In unix: run treemix_analysis.R

##############################POP STRUCTURE##################################

#### ADEGENET ####
dat.plink<-read.PLINK("stacks/subset.raw",parallel=FALSE)
#look at alleles
glPlot(dat.plink, posi="topleft")

#summarize allele freq distribution
mean<-glMean(dat.plink)
myFreq <- c(mean, 1-mean)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
	main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

#summarize missing data
temp<-density(glNA(dat.plink))
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
xlim=c(0,1701))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(glNA(dat.plink), rep(0, nLoc(dat.plink)), pch="|", col="blue")

#PCA!
pca1<-glPca(dat.plink, parallel=FALSE,nf=10)

#discriminant analysis of principal components (DAPC)
dat.clust<-find.clusters(dat.plink, max.n.clust=16) #retain 400 PCs and 6 clusters
                         #parallel=FALSE, n.pca=20, n.clust=NULL,	choose.n.clust=FALSE
dapc1<-dapc(dat.plink, dat.clust$grp, parallel=F) #kept 12 clusters

ind.names<-rownames(dapc1$ind.coord)
pop<-substr(ind.names, 8,11)
colors<-pop
colors[colors %in% sw.list]<-"black"
colors[colors %in% fw.list]<-"#2166ac"
pop.pch<-pop
pop.pch[pop.pch=="TXSP"]<-"0"
pop.pch[pop.pch=="TXCC"]<-"1"
pop.pch[pop.pch=="TXFW"]<-"3"
pop.pch[pop.pch=="TXCB"]<-"5"
pop.pch[pop.pch=="LAFW"]<-"4"
pop.pch[pop.pch=="ALST"]<-"16"
pop.pch[pop.pch=="ALFW"]<-"8"
pop.pch[pop.pch=="FLSG"]<-"17"
pop.pch[pop.pch=="FLKB"]<-"18"
pop.pch[pop.pch=="FLFD"]<-"19"
pop.pch[pop.pch=="FLSI"]<-"21"
pop.pch[pop.pch=="FLAB"]<-"22"
pop.pch[pop.pch=="FLPB"]<-"23"
pop.pch[pop.pch=="FLHB"]<-"24"
pop.pch[pop.pch=="FLCC"]<-"25"
pop.pch[pop.pch=="FLLG"]<-"13"
#pop plot info
ppi<-data.frame(Pop=pop.labs,cols = all.colors,pch=c(0,1,3,5,4,15,8,17,18,19,21,22,23,24,25,13))

#png("adegenet.dapc.png",height=7,width=7,units="in",res=300)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE,cell=0)
#dev.off()
compoplot(dapc1)

da<-data.frame(Individual=rownames(dapc1$ind.coord),Pop=substr(rownames(dapc1$ind.coord),8,11),
               LD1=dapc1$ind.coord[,1],LD2=dapc1$ind.coord[,2],
               LD3=dapc1$ind.coord[,3],LD4=dapc1$ind.coord[,4],LD5=dapc1$ind.coord[,5],
               Group=dapc1$grp, GrpName=names(dapc1$grp),
               stringsAsFactors = F) #include grpname to check
da$Pop[da$Pop == "FLLG"]<-"FLFW"
da$colors<-da$Pop
for(i in 1:nrow(da)){
  da[i,"colors"]<-as.character(ppi[ppi$Pop %in% da[i,"Pop"],"cols"])
}
da$pch<-da$Pop
for(i in 1:nrow(da)){
  da[i,"pch"]<-as.numeric(ppi[ppi$Pop %in% da[i,"Pop"],"pch"])
}

jpeg("subset.adegenet.dapc.jpeg",res=300,height=10.5/3,width=10.5,units="in")
par(mfrow=c(1,3),oma=c(2,2,2,2),mar=c(2,2,2,2))
# plot(pca1$scores[,1], pca1$scores[,2], pch=as.numeric(pop.pch), cex=2,lwd=1.3,
#      col=alpha(colors, 0.5),bg=alpha(colors,0.25), ylab="", xlab="")
# mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
#       1, line = 2,cex=0.75)
# mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
#       2, line = 2,cex=0.75)
# text(x=8,y=-1,"Atlantic",col=grp.colors[6])
# text(x=-5,y=2,"Florida",col="black")
# text(x=-5,y=2,"Florida",col=grp.colors[4])
# text(x=-5,y=-6,"Texas",col=grp.colors[1])
# 
# plot(pca1$scores[,3], pca1$scores[,4], pch=as.numeric(pop.pch), cex=2,lwd=1.3,
#      col=alpha(colors, 0.5),bg=alpha(colors,0.25), ylab="", xlab="")
# mtext(paste("PC3: ", round(pca1$eig[3]/sum(pca1$eig)*100, 2), "%", sep=""), 
#       1, line = 2,cex=0.75)
# mtext(paste("PC4: ", round(pca1$eig[4]/sum(pca1$eig)*100, 2), "%", sep=""), 
#       2, line = 2,cex=0.75)

plot(da$LD1,da$LD2,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="",xlim=c(-20,10),ylim=c(-10,25))
par(lwd=3,bty='n')
s.class(dapc1$ind.coord,fac=factor(da$Group), clabel=0,cstar=0,cellipse=2.5,
        addaxes = F,pch="",grid=F,axesel=F,add.plot = T,col=grp.colors[c(3,5,6,4,2,1)],xlim=c(-20,10),ylim=c(-10,25))
mtext(paste("Discriminant Axis 1 (", round(dapc1$eig[1]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("Discriminant Axis 2 (", round(dapc1$eig[2]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)


plot(da$LD3,da$LD4,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="",xlim=c(-10,10),ylim=c(-5,15))
par(lwd=3,bty='n')
s.class(dapc1$ind.coord[,3:4],fac=factor(da$Group), clabel=0,cstar=0,cellipse=2.5,
        addaxes = F,pch="",grid=F,axesel=F,add.plot = T,col=grp.colors[c(3,5,6,4,2,1)],xlim=c(-10,10),ylim=c(-5,15))
mtext(paste("Discriminant Axis 3 (", round(dapc1$eig[3]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("Discriminant Axis 4 (", round(dapc1$eig[4]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)

plot(da$LD4,da$LD5,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="",xlim=c(-5,15),ylim=c(-10,5))
par(lwd=3,bty='n')
s.class(dapc1$ind.coord[,4:5],fac=factor(da$Group), clabel=0,cstar=0,cellipse=2.5,
        addaxes = F,pch="",grid=F,axesel=F,add.plot = T,col=grp.colors[c(3,5,6,4,2,1)],xlim=c(-5,15),ylim=c(-10,5))
mtext(paste("Discriminant Axis 4 (", round(dapc1$eig[4]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("Discriminant Axis 5 (", round(dapc1$eig[5]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
    cex=1,lwd=1.3)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", legend=ppi$Pop, pch=as.numeric(ppi$pch), pt.cex=1.5,cex=0.85,
       col=alpha(ppi$cols, 0.5),pt.bg=alpha(ppi$cols,0.25), ncol=8,bty='n')
dev.off()

#### PCADAPT ####
library(pcadapt)
filename<-read.pcadapt("stacks/subset.ped",type="ped")
x<-pcadapt("stacks/subset.pcadapt", K=20)
plot(x,option="screeplot")#K=6
pops<-gsub("sample_(\\w{4}).*","\\1",ped.sub$V2)
grp<-pops
grp[grp=="TXFW" | grp=="LAFW" | grp=="ALFW" | grp=="FLLG"]<-"freshwater"
grp[grp!="freshwater"]<-"saltwater"
plot(x,option="scores",pop=pops)
plot(x,option="scores",i=3,j=4,pop=pops)
plot(x,option="scores",i=5,j=6,pop=pops)
plot(x,option="scores",i=7,j=8,pop=pops)#should show no patterns



pa<-pcadapt("stacks/subset.pcadapt",K=6)
pap<-data.frame(Pop=pops,cols=pops,pch=pops,grp=grp,stringsAsFactors = F)
pap$Pop[pap$Pop == "FLLG"]<-"FLFW"
for(i in 1:nrow(pap)){
  pap[i,"cols"]<-as.character(ppi[ppi$Pop %in% pap[i,"Pop"],"cols"])
}
for(i in 1:nrow(pap)){
  pap[i,"pch"]<-as.numeric(ppi[ppi$Pop %in% pap[i,"Pop"],"pch"])
}
pa.props<-round((pa$singular.values/sum(pa$singular.values))*100,2)
png("pcadapt.pc1-6.2017.png",height=8,width=10.5,units="in",res=300)
par(mfrow=c(2,3),oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(pa$scores[,1],pa$scores[,2],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),
     pch=as.numeric(pap$pch),	cex=1.5)
mtext(paste("PC1 (",pa.props[1],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC2 (",pa.props[2],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[,3],pa$scores[,4],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),pch=as.numeric(pap$pch),
	cex=1.5)
mtext(paste("PC3 (",pa.props[3],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC4 (",pa.props[4],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[,5],pa$scores[,6],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),pch=as.numeric(pap$pch),
	cex=1.5)
mtext(paste("PC5 (",pa.props[5],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC6 (",pa.props[6],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[grp=="freshwater",1],pa$scores[grp=="freshwater",2],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
	cex=1.5)
mtext(paste("PC1 (",pa.props[1],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC2 (",pa.props[2],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[grp=="freshwater",3],pa$scores[grp=="freshwater",4],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
	cex=1.5)
mtext(paste("PC3 (",pa.props[3],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC4 (",pa.props[4],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[grp=="freshwater",5],pa$scores[grp=="freshwater",6],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
	cex=1.5)
mtext(paste("PC5 (",pa.props[5],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC6 (",pa.props[2],"%)",sep=""),2,line = 2,cex=0.75)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("top", legend=ppi$Pop, pch=as.numeric(ppi$pch), pt.cex=1.5,cex=0.85,
       col=alpha(ppi$cols, 0.5),pt.bg=alpha(ppi$cols,0.25), ncol=8,bty='n')
dev.off()

###### STRUCTURE #####

#setwd("../../")

structure.k2<-read.table(
  "structure//fwsw//admix//Results//admix_run_2_f_clusters.txt",
  sep='\t', header=F)
structure.k2$V1<-sub('sample_([A-Z]{4})','\\1', structure.k2$V1)
tapply(structure.k2$V2,structure.k2$V1,max) #V2 has TX group

structure.k6<-read.table(
  "structure//fwsw//admix//Results//admix_run_6_f_clusters.txt",
  sep='\t', header=F)
structure.k6$V1<-sub('sample_([A-Z]{4})','\\1', structure.k6$V1)
tapply(structure.k6$V2,structure.k6$V1,max) #V2 has FLAt group
tapply(structure.k6$V3,structure.k6$V1,max) #V3 has TX group
tapply(structure.k6$V4,structure.k6$V1,max) #V4 has AL group
tapply(structure.k6$V5,structure.k6$V1,max) #V5 has FL Gulf group
tapply(structure.k6$V6,structure.k6$V1,max) #V6 has north Gulf group
tapply(structure.k6$V7,structure.k6$V1,max) #V7 has FL atlantic group
str6<-data.frame(structure.k6$V1,structure.k6$V3,structure.k6$V4,structure.k6$V2,
                 structure.k6$V5,structure.k6$V7,structure.k6$V6,stringsAsFactors = F)
str6$structure.k6.V1[str6$structure.k6.V1 == "FLLG"]<-"FLFW"

png("StructureK2K6.png",width=10,height=7.5,units="in",res=300)
par(mfrow=c(2,length(pop.list)),oma=c(1,3.5,1,1),mar=c(1,0,0,0))
plotting.structure(structure.k2,2,pop.list, make.file=FALSE, xlabcol = xcol,plot.new=F,
                   colors=grp.colors[c(1,6)],xlabel=F,
                   ylabel=expression(atop(italic(K)==2,Delta~italic(K)==358.9)))
plotting.structure(str6,2,pop.list, make.file=FALSE, plot.new=F,
                   colors=grp.colors,xlabel=T,xlabcol = xcol,
                   ylabel=expression(atop(italic(K)==6,Delta~italic(K)==326.1)))
dev.off()

##### COMBINED GRAPH ####
npop<-length(pop.list)
pseq<-1:npop
m<-matrix(c(1:32,rep(33,4),rep(34,4),rep(35,4),rep(0,4),
            rep(36,4),rep(37,4),rep(38,4),rep(0,4)),
          nrow=4,ncol=npop,byrow = T)
pdf("pop_structure_comb.pdf",height=8,width=8)
jpeg("pop_structure_comb.jpeg",res=300,height=8,width=8,units="in")
layout(mat=m)
#STRUCTURE
par(oma=c(1,3.5,1,1),mar=c(1,0,0,0))
plotting.structure(structure.k2,2,pop.list, make.file=FALSE, xlabcol = all.colors,plot.new=F,
                   colors=grp.colors[c(1,6)],xlabel=F,
                   ylabel=expression(atop(italic(K)==2,Delta~italic(K)==358.9)))
plotting.structure(str6,2,pop.labs, make.file=FALSE, plot.new=F,
                   colors=grp.colors,xlabel=T,xlabcol = all.colors,
                   ylabel=expression(atop(italic(K)==6,Delta~italic(K)==326.1)))
#PCADAPT
par(mar=c(2,2,2,2))
plot(pa$scores[,1],pa$scores[,2],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),
     pch=as.numeric(pap$pch),	cex=1.5)
points(pa$scores[grp=="freshwater",1],pa$scores[grp=="freshwater",2],
       col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
       bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
       cex=1.5) #to highlight the fw pops
mtext(paste("PC1 (",pa.props[1],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC2 (",pa.props[2],"%)",sep=""),2,line = 2,cex=0.75)

plot(pa$scores[,3],pa$scores[,4],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),pch=as.numeric(pap$pch),
     cex=1.5)
points(pa$scores[grp=="freshwater",3],pa$scores[grp=="freshwater",4],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
     cex=1.5)
mtext(paste("PC3 (",pa.props[3],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC4 (",pa.props[4],"%)",sep=""),2,line = 2,cex=0.75)

plot(pa$scores[,5],pa$scores[,6],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),pch=as.numeric(pap$pch),
     cex=1.5)
points(pa$scores[grp=="freshwater",5],pa$scores[grp=="freshwater",6],
       col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
       bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
       cex=1.5)
mtext(paste("PC5 (",pa.props[5],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC6 (",pa.props[6],"%)",sep=""),2,line = 2,cex=0.75)

# ADEGENET
par(mar=c(2,2,2,2))
plot(da$LD1,da$LD2,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="",xlim=c(-20,10),ylim=c(-10,25))
par(lwd=3,bty='n')
s.class(dapc1$ind.coord,fac=factor(da$Group), clabel=0,cstar=0,cellipse=2.5,
        addaxes = F,pch="",grid=F,axesel=F,add.plot = T,col=grp.colors[c(3,5,6,4,2,1)],xlim=c(-20,10),ylim=c(-10,25))
mtext(paste("DAPC 1 (", round(dapc1$eig[1]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("DAPC 2 (", round(dapc1$eig[2]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)


plot(da$LD3,da$LD4,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="",xlim=c(-10,10),ylim=c(-5,15))
par(lwd=3,bty='n')
s.class(dapc1$ind.coord[,3:4],fac=factor(da$Group), clabel=0,cstar=0,cellipse=2.5,
        addaxes = F,pch="",grid=F,axesel=F,add.plot = T,col=grp.colors[c(3,5,6,4,2,1)],xlim=c(-10,10),ylim=c(-5,15))
mtext(paste("DAPC 3 (", round(dapc1$eig[3]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("DAPC 4 (", round(dapc1$eig[4]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)

plot(da$LD4,da$LD5,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="",xlim=c(-5,15),ylim=c(-10,5))
par(lwd=3,bty='n')
s.class(dapc1$ind.coord[,4:5],fac=factor(da$Group), clabel=0,cstar=0,cellipse=2.5,
        addaxes = F,pch="",grid=F,axesel=F,add.plot = T,col=grp.colors[c(3,5,6,4,2,1)],xlim=c(-5,15),ylim=c(-10,5))
mtext(paste("DAPC 4 (", round(dapc1$eig[4]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("DAPC 5 (", round(dapc1$eig[5]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
    cex=1,lwd=1.3)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(x=0.6,y=-0.1, legend=ppi$Pop, pch=as.numeric(ppi$pch), pt.cex=1.5,cex=0.85,
       col=alpha(ppi$cols, 0.5),pt.bg=alpha(ppi$cols,0.25), ncol=2,bty='n')
dev.off()


##############################BAYENV######################################
ped.pops<-gsub("(sample_)(\\w{4})(\\w+)","\\2",ped.sub[,2])
ped.sex<-sub('(sample_\\w{4})(\\w)(\\w+)','\\2', ped.sub[,2])
ped.sex[ped.sex=="F"]<-2
ped.sex[ped.sex=="D"]<-2
ped.sex[ped.sex=="M"]<-1
ped.sex[ped.sex=="P"]<-1
ped.sex[ped.sex=="N"]<-1
ped.sex[ped.sex=="I"]<-0
ped.sex[ped.sex=="J"]<-0
ped.sub[,1]<-ped.pops
ped.sub[,5]<-ped.sex
write.table(ped.sub,"bayenv/bayenv.plink.ped", 
            row.names=F, col.names=F, quote=F, sep=" ",eol="\n")

clust.plink<-data.frame(FamID=ped.pops, IndID=ped.sub[,2],Pop=ped.pops)
write.table(clust.plink, 
            "bayenv/plink.clust.txt",
            col.names=F, row.names=F, quote=F, sep="\t", eol="\n")
#Then plink --ped bayenv.plink.ped --map subset.map --extract plink.snplist --out bayenv --noweb --allow-no-sex --recode --freq --within plink.clust.txt 
#9820 SNPs in 698 individuals

##### CONVERT PLINK TO BAYENV2
freq<-read.table("stacks/bayenv.frq.strat", 
                 header=T, stringsAsFactors=F)
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

write.table(snpsfile, "bayenv/fwsw.snpsfile", 
            col.names=F,row.names=F,quote=F,sep="\t",eol="\n") #bayenv SNPSFILE

#NOW RUN MATRIX ESTIMATION: run_bayenv2_matrix_general.sh
#../../scripts/run_bayenv2_matrix_general.sh fwsw.snpsfile 16
#last run on 1 May 2017

##### check Bayenv2 matrix
matrix.files<-list.files("bayenv/",pattern="matrix")
matrices<-list()
for(i in 1:length(matrix.files))
{
  matrices[[i]]<-as.matrix(read.table(paste("bayenv/",matrix.files[i],sep=""), skip=3597, header=F))
  rownames(matrices[[i]])<-colnames(matrices[[i]])<-pop.order
  #there are multiple matrices, one every 500 iterations..I'm taking the last one
}
image(matrices[[1]])
image(matrices[[2]])
image(matrices[[3]])
image(matrices[[4]])
image(matrices[[5]])
image(matrices[[6]])
image(matrices[[7]])
image(matrices[[8]])
image(matrices[[9]])
image(matrices[[10]])
#these are all essentially the same-use representative matrix

##### SNPFILEs
#for SNPFILE, need just one file per SNP apparently.
#want to use all of the snps (not just the pruned set)...need to get map with those inds.
all.snps.map<-read.table("stacks/batch_2.plink.map",header=F,stringsAsFactors = F)
all.snps.ped<-read.table("stacks/batch_2.plink.ped", header=F, stringsAsFactors=F)
ped.pop<-sub('sample_(\\w{4}).*','\\1', all.snps.ped[,2])
all.snps.ped[,1]<-ped.pop
all.snps.ped[,5]<-ped.sex
write.table(all.snps.ped,"stacks/fwsw_all.ped",col.names=F,row.names=F,quote=F,sep='\t',eol='\n')
all.snps.clust<-cbind(ped.pop,all.snps.ped[,2],ped.pop)
write.table(all.snps.clust, "stacks/all.clust.txt", sep="\t", eol="\n", quote=F,
            row.names=F, col.names=F)
#then need to run 
#plink --map batch_2.plink.map --ped fwsw_all.ped --freq --within all.clust.txt --allow-no-sex --noweb --out all.bayenv.plink
#57251 markers in 698 individuals

#read in frequency per pop
all.snps.frq<-read.table("stacks/all.bayenv.plink.frq.strat", 
                         header=T, stringsAsFactors=F)
#want to get $MAC for every snp at every pop 
#and NCHROBS-MAC for every stnp at every pop
freq<-cbind(all.snps.frq,all.snps.frq$NCHROBS-all.snps.frq$MAC)
colnames(freq)[ncol(freq)]<-"NAC"
pop.order<-levels(as.factor(freq$CLST))
snp.names<-split(freq$SNP,freq$CLST)[[1]]
mac.by.pop<-as.data.frame(split(freq$MAC,freq$CLST)) #57251 SNPs
rownames(mac.by.pop)<-all.snps.map[,2]
nac.by.pop<-as.data.frame(split(freq$NAC,freq$CLST))
rownames(nac.by.pop)<-all.snps.map[,2]
all.snpsfile<-gdata::interleave(mac.by.pop,nac.by.pop) 

write.table(all.snpsfile, "bayenv/all.fwsw", 
            col.names=F,row.names=T,quote=F,sep="\t",eol="\n")

#Bayenv says to run this:
#./calc bfs.sh SNPSFILE ENVIRONFILE MATRIXFILE NUMPOPS NUMITER NUMENVIRON
#$ ~/Programs/bayenv_2/calc_bf.sh all.fwsw env_data_std.txt representative_matrix.txt 16 100000 3

#INSTEAD:
#sarah@sarah-VirtualBox:~/sf_ubuntushare/popgen/fwsw_results/bayenv$ Rscript --vanilla ../../scripts/SNPSfromSNPSFILE.R all.fwsw ~/sf_ubuntushare/popgen/fwsw_results/bayenv/snpfiles/
#../../scripts/run_bayenv2_general.sh representative_matrix.txt env_data_std.txt 16 3 snpfiles
#####ENVFILE
env.raw<-read.csv("bayenv/env_data_raw.csv",row.names=1) 
#Each environmental variable should be standardized, 
#i.e. subtract the mean and then divided through by the standard deviation 
#of the variable across populations.

env.std<-t(apply(env.raw,1,std.by.mean))
env.std<-env.std[,colnames(snpsfile)] #change the column order to match
write.table(env.std,
            "bayenv/env_data_std.txt",
            sep='\t',quote=F,col.names=F,row.names=F,eol='\n')

##Are they correlated with distance?
colnames(env.raw)[colnames(env.dist) =="FLLG"]<-"FLFW"
env.dist<-as.matrix(vegdist(t(env.raw)))
env.dist<-env.dist[rownames(dist),colnames(dist)]
mantel.rtest(as.dist(t(dist)),as.dist(env.dist),999)
# Monte-Carlo test
# Observation: 0.1329135 
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
# Based on 999 replicates
# Simulated p-value: 0.104


#####GET OUTPUT
bayenv.all<-read.table("bayenv/bf_environ.env_data_std.txt") #59865 rows??
#extras got added in, need to figure out which ones to remove.
dup.snps<-bayenv.all$V1[duplicated(bayenv.all$V1)]
head(bayenv.all[bayenv.all$V1 %in% dup.snps,])
bayenv.all[bayenv.all$V1 %in% dup.snps[1],]
bayenv.all<-bayenv.all[2615:59865,]

colnames(bayenv.all)<-c("SNP",rownames(env.raw))
bayenv.all$SNP<-rownames(nac.by.pop)
bayenv.all$locus<-gsub("(\\d+)_\\d+","\\1",bayenv.all$SNP)
all.snps.map$locus<-gsub("(\\d+)_\\d+","\\1",all.snps.map$V2)

bf<-as.data.frame(do.call("rbind",apply(bayenv.all,1,function(x){
  chrom<-all.snps.map$V1[all.snps.map$locus %in% x["locus"]]
  pos<-all.snps.map$V4[all.snps.map$V2 %in% x["SNP"]]
  if(length(unique(chrom))==1){
    this.chrom<-unique(chrom)
  } else{
    this.chrom<-"BADCHROMNUM"
  }
  new<-data.frame(SNP=x["SNP"],temp=as.numeric(x["temp"]),
                     salinity=as.numeric(x["salinity"]),
                     seagrass=as.numeric(x["seagrass"]),
                     Chrom=as.character(this.chrom),
                     Pos=as.numeric(pos),stringsAsFactors=F)
  return(new)
})),stringsAsFactors=F)
rownames(bf)<-NULL
colnames(bf)<-c("SNP","temp","salinity","seagrass","Chrom","Pos")
bf$logTemp<-log(as.numeric(bf$temp))
bf$logSalt<-log(as.numeric(bf$salinity))
bf$logSeag<-log(as.numeric(bf$seagrass))
write.table(bf,"bayenv/fwsw_environ_corr.txt",
            col.names=T,row.names=F,quote=F,sep='\t')

bf.co<-apply(bf[,2:4],2,quantile,0.95) #get the cutoffs (2863 of them)
temp.sig<-bf[bf[,"temp"]>bf.co["temp"],]
salt.sig<-bf[bf[,"salinity"]>bf.co["salinity"],]
seag.sig<-bf[bf[,"seagrass"]>bf.co["seagrass"],]

tp<-fst.plot(fst.dat = bf,fst.name = "logTemp",chrom.name = "Chrom",bp.name = "Pos")
sp<-fst.plot(fst.dat = bf,fst.name = "logSalt",chrom.name = "Chrom",bp.name = "Pos")
gp<-fst.plot(fst.dat = bf,fst.name = "logSeag",chrom.name = "Chrom",bp.name = "Pos")


#####GET OUTPUT
#' ```{r, eval=FALSE}
setwd("bayenv/snpfiles")
bf.files<-list.files(pattern="bf")
xtx.files<-list.files(pattern="xtx")

bf.dat<-NULL
for(i in 1:length(bf.files)){
  bf<-read.table(bf.files[i])
  bf.dat<-rbind(bf.dat,bf)
}
colnames(bf.dat)<-c("locus", "Temp_BF", "Temp_rho", "Temp_rs", 
                    "Salinity_BF", "Salinity_rho", "Salinity_rs", 
                    "seagrass_BF", "seagrass_rho","seagrass_rs")
bf.dat$locus<-sub("snpfiles/(\\d+.*)","\\1",bf.dat$locus)

xtx<-NULL
for(i in 1:length(xtx.files)){
  x<-read.table(xtx.files[i])
  xtx<-rbind(xtx,x)
}
colnames(xtx)<-c("locus","XtX")
xtx$locus<-sub("snpfiles/(\\d+.*)","\\1",xtx$locus)

setwd("../")
write.table(bf.dat,"BF_summary.txt",sep='\t',quote=F,row.names=F)
write.table(xtx,"XtX_summary.txt",sep='\t',quote=F,row.names=F)

#' ```

xtx<-read.table("bayenv/XtX_summary.txt",header=TRUE)

#' Analyze the environmental associations
bf.dat<-read.table("bayenv/BF_summary.txt",header=TRUE)
bf.scaff<-merge(all.snps.map, bf.dat, by.x="V2", by.y="locus")
colnames(bf.scaff)[1:4]<-c("locus","scaffold","dist","BP")
bf<-cbind(bf.scaff[,1:4],bf.scaff[,grep("BF",colnames(bf.scaff))])
bf$col<-gsub("\\d+_(\\d+)","\\1",bf$locus)
bf$SNP<-paste(bf$scaffold,as.numeric(bf$BP)+1,sep=".")
bf.co<-apply(bf[,5:7],2,quantile,0.99) #focus on Bayes Factors, because of Lotterhos & Whitlock (2015)
temp.bf.sig<-bf[bf$Temp_BF>bf.co["Temp_BF"],c(1,2,4,8,5)]
sal.bf.sig<-bf[bf$Salinity_BF>bf.co["Salinity_BF"],c(1,2,4,8,6)]
grass.bf.sig<-bf[bf$seagrass_BF>bf.co["seagrass_BF"],c(1,2,4,8,7)]

dim(temp.bf.sig[temp.bf.sig$locus %in% sal.bf.sig$locus & 
                  temp.bf.sig$locus %in% grass.bf.sig,])

bf$logSal<-log(bf$Salinity_BF)
bf$logTemp<-log(bf$Temp_BF)
bf$logSeagrass<-log(bf$seagrass_BF)

#' plot bayenv results
png("Bayenv.png",height=10,width=7,units="in",res=300)
par(mfrow=c(3,1),mar=c(2,2,2,1),oma=c(2,2,2,1))
bs.sal<-fst.plot(bf,fst.name="logSal",chrom.name="scaffold",bp.name = "BP",
         scaffs.to.plot = lgs,pch=19,axis.size = 1,pt.cex = 1)
# points(bs.sal[bs.sal$SNP %in% sal.bf.sig$SNP,"plot.pos"],
#        bs.sal[bs.sal$SNP %in% sal.bf.sig$SNP,"logSal"],
#        col="cornflowerblue")
mtext("log(Salinity Bayes Factor)",2,cex=0.75,line=2.1)
points(bs.sal[bs.sal$SNP %in% stacks.sig$SNP,"plot.pos"],
       bs.sal[bs.sal$SNP %in% stacks.sig$SNP,"logSal"],
       col="cornflowerblue",cex=1.3)
clip(min(bs.sal$plot.pos),max(bs.sal$plot.pos),
     min(bs.sal$logSal),max(bs.sal$logSal))
abline(h=log(bf.co["Salinity_BF"]),col="cornflowerblue",lwd=2)

bs.temp<-fst.plot(bf,fst.name="logTemp",chrom.name="scaffold",bp.name = "BP",
                 scaffs.to.plot = lgs,pch=19,axis.size = 1,pt.cex=1)
points(bs.sal[bs.sal$SNP %in% stacks.sig$SNP,"plot.pos"],
       bs.sal[bs.sal$SNP %in% stacks.sig$SNP,"logTemp"],
       col="cornflowerblue",cex=1.3)
clip(min(bs.sal$plot.pos),max(bs.sal$plot.pos),
     min(bs.sal$logTemp),max(bs.sal$logTemp))
abline(h=log(bf.co["Temp_BF"]),col="cornflowerblue",lwd=2)
mtext("log(Temperature Bayes Factor)",2,cex=0.75,line=2.1)

bs.temp<-fst.plot(bf,fst.name="logSeagrass",chrom.name="scaffold",bp.name = "BP",
                  scaffs.to.plot = lgs,pch=19,axis.size = 1,pt.cex=1)
points(bs.sal[bs.sal$SNP %in% stacks.sig$SNP,"plot.pos"],
       bs.sal[bs.sal$SNP %in% stacks.sig$SNP,"logSeagrass"],
       col="cornflowerblue",cex=1.3)
clip(min(bs.sal$plot.pos),max(bs.sal$plot.pos),
     min(bs.sal$logSeagrass),max(bs.sal$logSeagrass))
abline(h=log(bf.co["seagrass_BF"]),col="cornflowerblue",lwd=2)
mtext("log(Seagrass Bayes Factor)",2,cex=0.75,line=2.1)
last<-0
for(i in 1:length(lgs)){
  text(x=median(bs.temp[bs.temp$scaffold ==lgs[i],"plot.pos"]),y=-3,
       labels=lgn[i], adj=1, xpd=TRUE,cex=1)
  last<-max(bs.temp[bs.temp$scaffold ==lgs[i],"plot.pos"])
}
dev.off()

###################BAYESCAN########################
#make the pops file
vcf<-parse.vcf("batch_2.vcf")
inds<-colnames(vcf[10:ncol(vcf)])
pops<-gsub("sample_(\\w{4})\\w+","\\1",inds)
pops.info<-cbind(inds,pops)
write.table(pops.info,"fwsw_pops_info_bayescan.txt",quote=F,row.names=F,col.names=F,sep='\t')
setwd("../bayescan")
bs.fst<-read.table("bayescan_fwsw_fst.txt")
plot_bayescan("bayescan_fwsw_fst.txt")
#plot the fsts
bs.fst<-data.frame(cbind(vcf[,1:2],bs.fst))
bs.fst.plot<-fst.plot(bs.fst,fst.name="fst",chrom.name="X.CHROM",bp.name="POS",axis.size=1)#weird fst peak around 0.2
points(bs.fst.plot$POS[bs.fst.plot$qval<0.05],bs.fst.plot$fst[bs.fst.plot$qval<0.05],pch=19,col="cornflowerblue",cex=0.5)


bs.sel<-read.table("bayescan_fwsw.sel")

#compare to stacks fsts
fst.shared<-paste(fw.shared.chr$Chr,(fw.shared.chr$BP+1),sep=".")
bs.fst.snp<-paste(bs.fst$X.CHROM[bs.fst$qval < 0.05],bs.fst$POS[bs.fst$qval < 0.05],sep=".")
length(bs.fst.snp[bs.fst.snp %in% fst.shared])
