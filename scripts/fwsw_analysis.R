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

setwd("B:/ubuntushare/popgen/fwsw_results/")
#source("../scripts/popgen_functions.R")
source("../../gwscaR/R/gwscaR.R")
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
#############################################################################
#######################**********FILES*********##############################
#############################################################################
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
chosen.snps<-choose.one.snp(vcf)$SNP
#############################################################################
#######################PLOT THE POINTS ON A MAP##############################
#############################################################################
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

#############################################################################
##################################IBD########################################
#############################################################################

#Mantel test using geographical distances and fsts

#read in the subsetted fst summary from running populations with whitelist
ibd.all<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst.all)))
ibd.sub<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst.sub)))

#test
pairwise.fst(ped.sub,9,10,pop.list)

ibd.by.loc<-fst.ibd.byloc(ped.sub,dist,pop.list)  #all NAs
#ignore warnings?  In is.euclid(m1) : Zero distance(s)
rownames(ibd.by.loc)<-sub.map$V2


#############################################################################
#################################OUTLIERS####################################
#############################################################################
###******************************STACKS*********************************####
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
fw.shared.chr<-fwsw.tx[fwsw.tx$Locus.ID %in% all.shared,c("Locus.ID","Chr","BP","Column")]
tapply(fw.shared.chr$Locus.ID,factor(fw.shared.chr$Chr),function(x){ length(unique(x)) })
#are they using the same SNPs or different SNPs?
snps<-data.frame(nrow=nrow(fw.shared.chr),ncol=4)
for(i in 1:nrow(fw.shared.chr)){
  tx.bp<-fwsw.tx[fwsw.tx$Fisher.s.P<0.01 & fwsw.tx$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  la.bp<-fwsw.la[fwsw.la$Fisher.s.P<0.01 & fwsw.la$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  al.bp<-fwsw.al[fwsw.al$Fisher.s.P<0.01 & fwsw.al$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  fl.bp<-fwsw.fl[fwsw.fl$Fisher.s.P<0.01 & fwsw.fl$Locus.ID == fw.shared.chr[i,"Locus.ID"],"BP"]
  snps[i,1]<-paste(tx.bp,sep=",",collapse = ",")
  snps[i,2]<-paste(la.bp,sep=",",collapse = ",")
  snps[i,3]<-paste(al.bp,sep=",",collapse = ",")
  snps[i,4]<-paste(fl.bp,sep=",",collapse = ",")
}
colnames(snps)<-c("TX","LA","AL","FL")
snps<-data.frame(cbind(fw.shared.chr,snps))

##compare to scovelli genome..
gff<-read.delim("../../scovelli_genome/ssc_122016_chromlevel.gff",header=F)
colnames(gff)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")

#agp<-read.table("../../scovelli_genome/ssc_122016_chromlevel.agp")
#colnames(agp)<-c("object","object_beg","object_end","part_number","component_type","component_id","component_beg","orientation")

fw.sig.reg<-do.call(rbind,apply(fw.shared.chr,1,function(sig){
  this.gff<-gff[as.character(gff$seqname) %in% as.character(unlist(sig["Chr"])),]
  this.reg<-this.gff[this.gff$start <= as.numeric(sig["BP"]) & this.gff$end >= as.numeric(sig["BP"]),]
  if(nrow(this.reg) == 0){
    if(as.numeric(sig["BP"])>max(as.numeric(this.gff$end))){
      new<-data.frame(Locus=sig["Locus.ID"],Chr=sig["Chr"],BP=sig["BP"],SNPCol=sig["Column"],
                      region="beyond.last.contig")
    }else{
      new<-data.frame(Locus=sig["Locus.ID"],Chr=sig["Chr"],BP=sig["BP"],SNPCol=sig["Column"],
                      region=NA)
    }
  }else{
    new<-data.frame(Locus=sig["Locus.ID"],Chr=sig["Chr"],BP=sig["BP"],SNPCol=sig["Column"],
                    region=paste(this.reg$feature,sep=",",collapse = ","))
  }
  return(as.data.frame(new))
}))

#graph without highlighted regions
png("stacks_fsts_fwsw.png",height=8,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(0.85,2,0,0.5),oma=c(1,1,1,0.5))
fwswt.fst<-fst.plot(fwsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col = grp.colors[1],pt.cex=0,axis.size = 1)
points(fwswt.fst$BP,fwswt.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[1])
fwswl.fst<-fst.plot(fwsw.la,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[3],pt.cex=0,axis.size=1)
points(fwswl.fst$BP,fwswl.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[3])
fwswa.fst<-fst.plot(fwsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[4],pt.cex=0,axis.size = 1)
points(fwswa.fst$BP,fwswa.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[4])
fwswf.fst<-fst.plot(fwsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[6],pt.cex=0,axis.size=1)
points(fwswf.fst$BP,fwswf.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[6])
mtext(expression(italic(F)[ST]),2,outer=T,line=-1,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(fwswf.fst[fwswf.fst$Chr ==lgs[i],"BP"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE)
  last<-max(fwswf.fst[fwswf.fst$Chr ==lgs[i],"BP"])
}
dev.off()

#plot with the outlier regions
png("stacks_fsts_fwsw_withSig.png",height=8,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(0.85,2,0,0.5),oma=c(1,1,1,0.5))
fwswt.fst<-fst.plot(fwsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col = grp.colors[1],pt.cex=0,axis.size = 1)
points(fwswt.fst$BP,fwswt.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[1])
points(fwswt.fst$BP[fwswt.fst$Locus.ID %in% all.shared],fwswt.fst$Corrected.AMOVA.Fst[fwswt.fst$Locus.ID %in% all.shared],
       pch=0,col="red",cex=1.3)
fwswl.fst<-fst.plot(fwsw.la,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[3],pt.cex=0,axis.size=1)
points(fwswl.fst$BP,fwswl.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[3])
points(fwswl.fst$BP[fwswl.fst$Locus.ID %in% all.shared],fwswl.fst$Corrected.AMOVA.Fst[fwswl.fst$Locus.ID %in% all.shared],
       pch=0,col="red",cex=1.3)
fwswa.fst<-fst.plot(fwsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[4],pt.cex=0,axis.size = 1)
points(fwswa.fst$BP,fwswa.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[4])
points(fwswa.fst$BP[fwswa.fst$Locus.ID %in% all.shared],fwswa.fst$Corrected.AMOVA.Fst[fwswa.fst$Locus.ID %in% all.shared],
       pch=0,col="red",cex=1.3)
fwswf.fst<-fst.plot(fwsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[6],pt.cex=0,axis.size=1)
points(fwswf.fst$BP,fwswf.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[6])
points(fwswf.fst$BP[fwswf.fst$Locus.ID %in% all.shared],fwswf.fst$Corrected.AMOVA.Fst[fwswf.fst$Locus.ID %in% all.shared],
       pch=0,col="red",cex=1.3)
mtext(expression(italic(F)[ST]),2,outer=T,line=-1,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(fwswf.fst[fwswf.fst$Chr ==lgs[i],"BP"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE)
  last<-max(fwswf.fst[fwswf.fst$Chr ==lgs[i],"BP"])
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
                    groups = scaffs,group.boundaries = scaff.starts)
swswa.fst<-fst.plot(swsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts)
swswf.fst<-fst.plot(swsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts)

png("stacks_fsts_swsw.png",height=6,width=7.5,units="in",res=300)
par(mfrow=c(3,1),mar=c(0.85,2,0,0.5),oma=c(1,1,1,0.5))
swswt.fst<-fst.plot(swsw.tx,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col = grp.colors[1],pt.cex=0,axis.size = 1)
points(swswt.fst$BP,swswt.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[1])
swswa.fst<-fst.plot(swsw.al,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", y.lim=c(0,1),
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[3],pt.cex=0,axis.size = 1)
points(swswa.fst$BP,swswa.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[3])
swswf.fst<-fst.plot(swsw.fl,fst.name = "Corrected.AMOVA.Fst", bp.name = "BP",chrom.name = "Chr", 
                    groups = scaffs,group.boundaries = scaff.starts,pt.col=grp.colors[6],pt.cex=0,axis.size=1)
points(swswf.fst$BP,swswf.fst$Corrected.AMOVA.Fst,pch=21,bg=grp.colors[6])
mtext(expression(italic(F)[ST]),2,outer=T,line=-1,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(swswf.fst[swswf.fst$Chr ==lgs[i],"BP"]),y=-0.015,
       labels=lgn[i], adj=1, xpd=TRUE)
  last<-max(swswf.fst[swswf.fst$Chr ==lgs[i],"BP"])
}
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
lfs<-gwsca(vcf=vcf,locus.info=loci.info,group1=alst,group2=lafw)
afs<-gwsca(vcf=vcf,locus.info=loci.info,group1=alst,group2=alfw)
ffs<-gwsca(vcf=vcf,locus.info=loci.info,group1=flcc,group2=fllg)

tss<-gwsca(vcf=vcf,locus.info=loci.info,group1=txcb,group2=txcc)
ass<-gwsca(vcf=vcf,locus.info=loci.info,group1=alst,group2=flsg)
fss<-gwsca(vcf=vcf,locus.info=loci.info,group1=flcc,group2=flhb)

#With ped file
txcb<-grep("TXCB",ped.sub$V2,value = T)
txfw<-grep("TXFW",ped.sub$V2,value = T)
alst<-grep("ALST",ped.sub$V2,value = T)
alfw<-grep("ALFW",ped.sub$V2,value = T)
lafw<-grep("LAFW",ped.sub$V2,value = T)
flcc<-grep("FLCC",ped.sub$V2,value = T)
fllg<-grep("FLLG",ped.sub$V2,value = T)
txcc<-grep("TXCC",ped.sub$V2,value = T)
flsg<-grep("FLSG",ped.sub$V2,value = T)
flhb<-grep("FLHB",ped.sub$V2,value = T)

tfs.fst<-fst.one.plink(ped.sub,txcc,txfw)
  tfs.fst<-fst.sig(tfs.fst)
lfs.fst<-fst.one.plink(ped.sub,alst,lafw)
  lfs.fst<-fst.sig(lfs.fst)
afs.fst<-fst.one.plink(ped.sub,alst,alfw)
afs.fst<-fst.sig(afs.fst)
ffs.fst<-fst.one.plink(ped.sub,flcc,fllg)
  ffs.fst<-fst.sig(ffs.fst)
tss.fst<-fst.one.plink(ped.sub,txcc,txcb)
tss.fst<-fst.sig(tss.fst)
ass.fst<-fst.one.plink(ped.sub,alst,flsg)
ass.fst<-fst.sig(ass.fst)
fss.fst<-fst.one.plink(ped.sub,flcc,flhb)
fss.fst<-fst.sig(fss.fst)
# fwfw.tf<-read.delim("stacks/batch_2.fst_FLLG-TXFW.tsv")
# fwfw.ta<-read.delim("stacks/batch_2.fst_ALFW-TXFW.tsv")
# fwfw.tl<-read.delim("stacks/batch_2.fst_LAFW-TXFW.tsv")
# fwfw.la<-read.delim("stacks/batch_2.fst_ALFW-LAFW.tsv")
# fwfw.lf<-read.delim("stacks/batch_2.fst_FLLG-LAFW.tsv")
# fwfw.af<-read.delim("stacks/batch_2.fst_ALFW-FLLG.tsv")
# 
# swsw.tf<-read.delim("stacks/batch_2.fst_FLCC-TXCB.tsv")
# swsw.ta<-read.delim("stacks/batch_2.fst_ALST-TXCB.tsv")
# swsw.af<-read.delim("stacks/batch_2.fst_ALST-FLCC.tsv")

scaffs<-levels(as.factor(c(as.character(swsw.tx[,"Chr"]),
                           as.character(fwsw.tx[,"Chr"]),
                           as.character(swsw.al[,"Chr"]),
                           as.character(fwsw.al[,"Chr"]),
                           as.character(fwsw.la[,"Chr"]),
                           as.character(swsw.fl[,"Chr"]),
                           as.character(fwsw.fl[,"Chr"]))))
scaffs[1:22]<-lgs
scaff.starts<-tapply(vcf$POS,vcf$`#CHROM`,max)
scaff.starts<-data.frame(rbind(cbind(names(scaff.starts),scaff.starts)),stringsAsFactors = F)
png("FW-SW_Fsts.png")
par(mfrow=c(4,2),oma=c(3,3,1,0),mar=c(1,0.5,1,0))
ss.t<-fst.plot(swsw.tx,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)
mtext("Texas",2,line=1.4)
mtext("Saltwater",3)
fs.t<-fst.plot(fwsw.tx,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)
mtext("Freshwater",3)

#LA
ss.l<-fst.plot(swsw.al,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)
mtext("Louisiana",2,line=1.4)
fs.l<-fst.plot(fwsw.la,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)
#AL
ss.a<-fst.plot(swsw.al,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)
mtext("Alabama",2,line=1.4)
fs.a<-fst.plot(fwsw.al,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),pt.cex =1.25,
               groups=scaffs,group.boundaries = scaff.starts)

#FL
ss.f<-fst.plot(swsw.fl,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(ss.f[ss.f$Chr ==lgs[i],"BP"]),y=-0.03,
       labels=lgs[i], adj=1, xpd=TRUE,srt=90,cex=1.25)
  last<-max(ss.f[ss.f$Chr ==lgs[i],"BP"])
}
mtext("Florida",2,line=1.4)
fs.f<-fst.plot(fwsw.fl,fst.name="Corrected.AMOVA.Fst",
               axis.size=1.25,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),pt.cex=1.25,
               groups=scaffs,group.boundaries = scaff.starts)

last<-0
for(i in 1:length(lgs)){
  text(x=mean(fs.f[fs.f$Chr ==lgs[i],"BP"]),y=-0.03,
       labels=lgs[i], adj=1, xpd=TRUE,srt=90,cex=1.25)
  last<-max(fs.f[fs.f$Chr ==lgs[i],"BP"])
}
mtext(expression(italic(F)[ST]),2,outer=T,line=1.5)
dev.off()
# 
# par(mfrow=c(6,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
# ff.tf<-plotting.fsts.scaffs(fwfw.tf,"Fst",pt.lty=1)
# ff.ta<-plotting.fsts.scaffs(fwfw.ta,"Fst",pt.lty=1)
# ff.tl<-plotting.fsts.scaffs(fwfw.tl,"Fst",pt.lty=1)
# ff.la<-plotting.fsts.scaffs(fwfw.la,"Fst",pt.lty=1)
# ff.lf<-plotting.fsts.scaffs(fwfw.lf,"Fst",pt.lty=1)
# ff.af<-plotting.fsts.scaffs(fwfw.af,"Fst",pt.lty=1)
# 
# par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
# ss.tf<-plotting.fsts.scaffs(swsw.tf,"Fst",pt.lty=1)
# ss.ta<-plotting.fsts.scaffs(swsw.ta,"Fst",pt.lty=1)
# ss.af<-plotting.fsts.scaffs(swsw.af,"Fst",pt.lty=1)

#### Delta-divergence ####
#' only use chosen SNPs
vcf<-vcf[vcf$SNP %in% chosen.snps,]
#sw-fw
swfw.mu<-calc.mean.fst(vcf = vcf,pop.list1 = sw.list,pop.list2 = fw.list,maf.cutoff=0.01)
fwfw.mu<-calc.mean.fst(vcf = vcf,pop.list1 = fw.list,pop.list2 = fw.list, maf.cutoff=0.01)
deltad<-merge(swfw.mu,fwfw.mu,by="SNP")
deltad<-deltad[,c("SNP","Chrom.x","Pos.x","Mean.Fst.x","Mean.Fst.y")]
colnames(deltad)<-c("SNP","Chrom","Pos","MeanSWFW.Fst","MeanFWFW.Fst")
deltad$deltad<-deltad$MeanSWFW.Fst - deltad$MeanFWFW.Fst
deltad<-deltad[!is.na(deltad$deltad),]#remove NAs
dd<-fst.plot(fst.dat = deltad,fst.name = "deltad",bp.name = "Pos",axis=1)
mtext(expression(paste(delta,"-divergence")),2,line=1.5)
smooth.out<-data.frame()
for(i in 1:length(lgs)){#scaffolds are too short
  this.chrom<-dd[dd$Chrom %in% lgs[i],]
  #span<-nrow(this.chrom)/5000
  this.smooth<-loess.smooth(this.chrom$plot.pos,this.chrom$deltad,span=0.1,degree=2) 
  points(this.smooth$x,this.smooth$y,col="cornflowerblue",type="l",lwd=2)
  this.out<-cbind(this.smooth$x[this.smooth$y>=0.2],this.smooth$y[this.smooth$y>=0.2])
  smooth.out<-rbind(smooth.out,this.out)
}


#' sliding window pi and rho - across all snps
#' pi = 1-sum((ni choose 2)/(n choose i)); ni is number of alleles i in sample, n = sum(ni)
#' rho=1 if allele in pop j is only found in that pop and at least one ind was genotyped at that site in each pop; rho = 0 otherwise
#' Jones et al. (2012) used 2500bp sliding windows with a step size 500bp<-more than just SNPs, but I'll just focus on SNPs
#' Hohenlohe did a similar thing and weighted pi by all nt sites (not just SNPs) but rho by SNPs only
#' Not sure how to make this work for me.

#' to get trees and calc gsi (maybe):
#' for each overlapping sliding window (of 33 SNPs, for example), 
#' generate distance matrix (Fsts) using those SNPs
#' ape::nj(as.dist(matrix), "unrooted")
#' genealogicalSorting::gsi(tree, class, assignments, uncertainty)


#############################################################################
##############################POP STRUCTURE##################################
#############################################################################
#******************************ADEGENET*********************************#
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

#*******************************PCADAPT***********************************#
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


###########################################################################
##############################BAYENV######################################
##########################################################################
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

#####CONVERT PLINK TO BAYENV2
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

#####check Bayenv2 matrix--NOT YET DONE
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

#####SNPFILEs
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
mac.by.pop<-as.data.frame(split(freq$MAC,freq$CLST))
rownames(mac.by.pop)<-all.snps.map[,2]
nac.by.pop<-as.data.frame(split(freq$NAC,freq$CLST))
rownames(nac.by.pop)<-all.snps.map[,2]
all.snpsfile<-interleave(mac.by.pop,nac.by.pop)

write.table(all.snpsfile, "bayenv/all.fwsw", 
            col.names=F,row.names=F,quote=F,sep="\t",eol="\n")

#then run this:
#./calc bfs.sh SNPSFILE ENVIRONFILE MATRIXFILE NUMPOPS NUMITER NUMENVIRON
#$ ~/Programs/bayenv_2/calc_bf.sh all.fwsw env_data_std.txt representative_matrix.txt 16 100000 3

## or try running it on your own
for(i in 1:nrow(all.snpsfile)){
  write.table(all.snpsfile[i:(i+1),],paste("bayenv/snpfiles/",rownames(all.snpsfile)[i],sep=""),
              quote=F,col.names=F,row.names=F,sep='\t')
}
#THINGS ARENT WORKING RIGHT - yields segmentation fault. maybe because of eol?

#~/Programs/bayenv_2/bayenv2 -i $f -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t
#~/Programs/bayenv_2/bayenv2 -i ./20645_26 -e ../env_data_std.txt -m ../representative_matrix.txt -k 100000 -r 416 -p 16 -n 3 -t

#####ENVFILE
env.raw<-read.csv("bayenv/env_data_raw.csv",row.names=1)
#Each environmental variable should be standardized, 
#i.e. subtract the mean and then divided through by the standard deviation 
#of the variable across populations.
std.by.mean<-function(x){
  m<-mean(x)
  s<-sd(x)
  newx<-(x-m)/s
  return(newx)
}
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
bayenv.all<-read.table("bayenv/bf_environ.env_data_std.txt")
#why are there 2 rows per snp???
colnames(bayenv.all)<-c("SNP",rownames(env.raw))
bayenv.all$SNP<-rownames(snpsfile)
bayenv.all$locus<-gsub("(\\d+)_\\d+","\\1",bayenv.all$SNP)
map.all$locus<-gsub("(\\d+)_\\d+","\\1",map.all$V2)

bf<-t(as.data.frame(do.call("cbind",apply(bayenv.all,1,function(x){
  chrom<-map.all$V1[map.all$locus %in% x["locus"]]
  pos<-map.all$V4[map.all$V2 %in% x["SNP"]]
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
})),stringsAsFactors=F))
rownames(bf)<-NULL
colnames(bf)<-c("SNP","temp","salinity","seagrass","locus","Chrom","Pos")
bf<-as.data.frame(bf,stringsAsFactors = F)
bf$logTemp<-log(as.numeric(bf$temp))
bf$logSalt<-log(as.numeric(bf$salinity))
bf$logSeag<-log(as.numeric(bf$seagrass))
write.table(bf,"bayenv/fwsw_environ_corr.txt",
            col.names=T,row.names=F,quote=F,sep='\t')

bf.co<-apply(bayenv.all[,2:4],2,quantile,0.95) #get the cutoffs
temp.sig<-bf[bf[,"temp"]>bf.co["temp"],]
salt.sig<-bf[bf[,"salinity"]>bf.co["salinity"],]
seag.sig<-bf[bf[,"seagrass"]>bf.co["seagrass"],]

fst.plot(fst.dat = bf,fst.name = "logTemp",chrom.name = "Chrom")

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
