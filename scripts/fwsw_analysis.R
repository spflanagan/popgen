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


setwd("E:/ubuntushare/popgen/fwsw_results/")
#source("../scripts/popgen_functions.R")
source("../../gwscaR/R/gwscaR.R")
source("../phenotype_functions.R")

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
map.sub<-read.table("stacks/subset.map",header = F)
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
#******************************STACKS*********************************#
#Compare neighboring pops.
fwsw.tx<-read.delim("stacks/batch_2.fst_TXCB-TXFW.tsv")
fwsw.la<-read.delim("stacks/batch_2.fst_ALST-LAFW.tsv")
fwsw.al<-read.delim("stacks/batch_2.fst_ALFW-ALST.tsv")
fwsw.fl<-read.delim("stacks/batch_2.fst_FLCC-FLLG.tsv")

swsw.tx<-read.delim("stacks/batch_2.fst_TXCB-TXCC.tsv")
swsw.al<-read.delim("stacks/batch_2.fst_ALST-FLSG.tsv")
swsw.fl<-read.delim("stacks/batch_2.fst_FLCC-FLHB.tsv")
  
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
                           as.character(fwsw.tx[,"Chr"]))))
png("FW-SW_Fsts.png")
par(mfrow=c(4,2),oma=c(0,0,0,0),mar=c(1,1,1,1))
ss.t<-fst.plot(swsw.tx,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                  levels(factor(swsw.tx$Chr))]))
mtext("Texas",2)
fs.t<-fst.plot(fwsw.tx,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(fwsw.tx$Chr))]))

#LA
ss.l<-fst.plot(swsw.al,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(swsw.al$Chr))]))
mtext("Louisiana",2)
fs.l<-fst.plot(fwsw.la,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(fwsw.la$Chr))]))
#AL
ss.a<-fst.plot(swsw.al,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(swsw.al$Chr))]))
mtext("Alabama",2)
fs.a<-fst.plot(fwsw.al,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(fwsw.al$Chr))]))

#FL
ss.f<-fst.plot(swsw.fl,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="black",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(swsw.fl$Chr))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(ss.f[ss.f$Chr ==lgs[i],"BP"]),y=-0.13,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
  last<-max(ss.f[ss.f$Chr ==lgs[i],"BP"])
}
mtext("Florida",2)
fs.f<-fst.plot(fwsw.fl,fst.name="Corrected.AMOVA.Fst",
               axis.size=1,chrom.name="Chr",pt.col="#2166ac",
               bp.name="BP",y.lim=c(0,1),
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(fwsw.fl$Chr))]))

last<-0
for(i in 1:length(lgs)){
  text(x=mean(fs.f[fs.f$Chr ==lgs[i],"BP"]),y=-0.13,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
  last<-max(fs.f[fs.f$Chr ==lgs[i],"BP"])
}

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

ind.names<-dimnames(pca1$scores)[[1]]

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

#discriminant analysis of principal components (DAPC)
dat.clust<-find.clusters(dat.plink, parallel=FALSE, n.pca=20, n.clust=NULL,
	choose.n.clust=FALSE, max.n.clust=16)
dapc1<-dapc(dat.plink, dat.clust$grp, parallel=F) #kept 12 clusters
#png("adegenet.dapc.png",height=7,width=7,units="in",res=300)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE,cell=0)
#dev.off()
compoplot(dapc1)

da<-data.frame(Individual=rownames(dapc1$ind.coord),Pop=substr(rownames(dapc1$ind.coord),8,11),
               LD1=dapc1$ind.coord[,1],LD2=dapc1$ind.coord[,2],
               LD3=dapc1$ind.coord[,3],LD4=dapc1$ind.coord[,4],
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

jpeg("subset.adegenet.dapc.jpeg",res=300,height=7,width=7,units="in")
par(mfrow=c(2,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(pca1$scores[,1], pca1$scores[,2], pch=as.numeric(pop.pch), cex=2,lwd=1.3,
     col=alpha(colors, 0.5),bg=alpha(colors,0.25), ylab="", xlab="")
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
      1, line = 2,cex=0.75)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
      2, line = 2,cex=0.75)
text(x=8,y=-1,"Atlantic",col=grp.colors[6])
text(x=-5,y=2,"Florida",col="black")
text(x=-5,y=2,"Florida",col=grp.colors[4])
text(x=-5,y=-6,"Texas",col=grp.colors[1])

plot(pca1$scores[,3], pca1$scores[,4], pch=as.numeric(pop.pch), cex=2,lwd=1.3,
     col=alpha(colors, 0.5),bg=alpha(colors,0.25), ylab="", xlab="")
mtext(paste("PC3: ", round(pca1$eig[3]/sum(pca1$eig)*100, 2), "%", sep=""), 
      1, line = 2,cex=0.75)
mtext(paste("PC4: ", round(pca1$eig[4]/sum(pca1$eig)*100, 2), "%", sep=""), 
      2, line = 2,cex=0.75)

plot(da$LD1,da$LD2,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="")
mtext(paste("Discriminant Axis 1 (", round(dapc1$eig[1]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("Discriminant Axis 2 (", round(dapc1$eig[2]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)
points(dapc1$grp.coord[,1],dapc1$grp.coord[,2],cex=5,col="purple")

plot(da$LD3,da$LD4,col=alpha(da$colors,0.5),pch=as.numeric(da$pch),cex=2,lwd=1.3,
     bg=alpha(colors,0.25),xlab="",ylab="")
mtext(paste("Discriminant Axis 3 (", round(dapc1$eig[3]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      1, line = 2,cex=0.75)
mtext(paste("Discriminant Axis 4 (", round(dapc1$eig[4]/sum(dapc1$eig)*100, 2), "%)", sep=""),
      2, line = 2,cex=0.75)
points(dapc1$grp.coord[,3],dapc1$grp.coord[,4],cex=5,col="purple")

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
    cex=1)
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

setwd("../../")

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
                 structure.k6$V5,structure.k6$V7,structure.k6$V6)

png("StructureK2K6.png",width=10,height=7.5,units="in",res=300)
par(mfrow=c(2,length(pop.list)),oma=c(1,3.5,1,1),mar=c(1,0,0,0))
plotting.structure(structure.k2,2,pop.list, make.file=FALSE, xlabcol = xcol,plot.new=F,
                   colors=grp.colors[c(1,6)],xlabel=F,
                   ylabel=expression(atop(italic(K)==2,Delta~italic(K)==358.9)))
plotting.structure(str6,2,pop.list, make.file=FALSE, plot.new=F,
                   colors=grp.colors,xlabel=T,xlabcol = xcol,
                   ylabel=expression(atop(italic(K)==6,Delta~italic(K)==326.1)))
dev.off()
