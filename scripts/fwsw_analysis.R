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
source("../scripts/popgen_functions.R")
source("../phenotype_functions.R")

pop.list<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
	"FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLLG")
fw.list<-c("TXFW","LAFW","ALFW","FLLG")
sw.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB",
	"FLFD","FLSI","FLAB","FLPB","FLHB","FLCC")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
lgn<-seq(1,22)

#############################################################################
#######################**********FILES*********##############################
#############################################################################
mar.coor<-read.csv("F://Docs//PopGen//marine_coordinates_revised.csv", header=T)
fw.coor<-read.csv("F://Docs//PopGen//fw_coordinates.csv", header=T)
dist<-read.table("fwsw_geographical_distances.txt",header=T,row.names=1,
	sep='\t')
pwise.fst.all<-read.table("stacks/populations/fwsw_fst_summary.txt",header=T,row.names=1,sep='\t')
	pwise.fst.all<-rbind(pwise.fst.all,rep(NA,ncol(pwise.fst.all)))
	rownames(pwise.fst.all)<-colnames(pwise.fst.all)
pwise.fst.sub<-read.table("stacks/fwsw_fst_summary_subset.txt",header=T,row.names=1,sep='\t')
ped.sub<-read.table("stacks/subset.ped",header=F)	
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
points(-1*fw.coor$lon, fw.coor$lat,  col="forest green", cex=1.5, pch=18)
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
text(x=-96,y=28.3,"TXFW",font=2,col="forest green")
text(x=-94.7,y=29,"TXCB",font=2)
text(x=-90.2,y=30.3,"LAFW",font=2,col="forest green")
text(x=-88,y=30,"ALST",font=2)
text(x=-87,y=30.75,"ALFW",font=2,col="forest green")
text(x=-85,y=29.4,"FLSG",font=2)
text(x=-83.5,y=29.2,"FLKB",font=2)
text(x=-83.2,y=27.6,"FLFD",font=2)
text(x=-82.2,y=26,"FLSI",font=2)
text(x=-80,y=24.8,"FLAB",font=2)
text(x=-79.5,y=26.8,"FLPB",font=2)
text(x=-79.7,y=27.2,"FLHB",font=2)
text(x=-80.2,y=28.2,"FLCC",font=2)
text(x=-80.9,y=29.3,"FLLG",font=2,col="forest green")
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
fwsw.tx<-read.delim("stacks/batch_2.fst_TXCB-TXFW.tsv")
fwsw.la<-read.delim("stacks/batch_2.fst_ALST-LAFW.tsv")
fwsw.al<-read.delim("stacks/batch_2.fst_ALFW-ALST.tsv")
fwsw.fl<-read.delim("stacks/batch_2.fst_FLCC-FLLG.tsv")

fwfw.tf<-read.delim("stacks/batch_2.fst_FLLG-TXFW.tsv")
fwfw.ta<-read.delim("stacks/batch_2.fst_ALFW-TXFW.tsv")
fwfw.tl<-read.delim("stacks/batch_2.fst_LAFW-TXFW.tsv")
fwfw.la<-read.delim("stacks/batch_2.fst_ALFW-LAFW.tsv")
fwfw.lf<-read.delim("stacks/batch_2.fst_FLLG-LAFW.tsv")
fwfw.af<-read.delim("stacks/batch_2.fst_ALFW-FLLG.tsv")

swsw.tf<-read.delim("stacks/batch_2.fst_FLCC-TXCB.tsv")
swsw.ta<-read.delim("stacks/batch_2.fst_ALST-TXCB.tsv")
swsw.af<-read.delim("stacks/batch_2.fst_ALST-FLCC.tsv")

png("FW-SW_Fsts.png")
par(mfrow=c(4,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
fs.t<-plotting.fsts.scaffs(fwsw.tx,"Fst",pt.lty=1,axis.size=0.75)
legend("top",c("TX FW-SW"),bty='n')
last<-0
for(i in 1:length(lgs)){
	text(x=mean(fs.t[fs.t$Chr ==lgs[i],"BP"]),y=-0.13,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(fs.t[fs.t$Chr ==lgs[i],"BP"])
}
fs.l<-plotting.fsts.scaffs(fwsw.la,"Fst",pt.lty=1,axis.size=0.75)
legend("top",c("LA FW-SW"),bty='n')
last<-0
for(i in 1:length(lgs)){
	text(x=mean(fs.l[fs.l$Chr ==lgs[i],"BP"]),y=-0.03,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(fs.l[fs.l$Chr ==lgs[i],"BP"])
}
fs.a<-plotting.fsts.scaffs(fwsw.al,"Fst",pt.lty=1,axis.size=0.75)
legend("top",c("AL FW-SW"),bty='n')
last<-0
for(i in 1:length(lgs)){
	text(x=mean(fs.a[fs.a$Chr ==lgs[i],"BP"]),y=-0.03,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(fs.t[fs.a$Chr ==lgs[i],"BP"])
}
fs.f<-plotting.fsts.scaffs(fwsw.fl,"Fst",pt.lty=1,axis.size=0.75)
legend("top",c("FL FW-SW"),bty='n')
last<-0
for(i in 1:length(lgs)){
	text(x=mean(fs.f[fs.f$Chr ==lgs[i],"BP"]),y=-0.13,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(fs.f[fs.f$Chr ==lgs[i],"BP"])
}
dev.off()

par(mfrow=c(6,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
ff.tf<-plotting.fsts.scaffs(fwfw.tf,"Fst",pt.lty=1)
ff.ta<-plotting.fsts.scaffs(fwfw.ta,"Fst",pt.lty=1)
ff.tl<-plotting.fsts.scaffs(fwfw.tl,"Fst",pt.lty=1)
ff.la<-plotting.fsts.scaffs(fwfw.la,"Fst",pt.lty=1)
ff.lf<-plotting.fsts.scaffs(fwfw.lf,"Fst",pt.lty=1)
ff.af<-plotting.fsts.scaffs(fwfw.af,"Fst",pt.lty=1)

par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
ss.tf<-plotting.fsts.scaffs(swsw.tf,"Fst",pt.lty=1)
ss.ta<-plotting.fsts.scaffs(swsw.ta,"Fst",pt.lty=1)
ss.af<-plotting.fsts.scaffs(swsw.af,"Fst",pt.lty=1)


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
colors[colors=="TXSP"]<-rainbow(16)[1]
colors[colors=="TXCC"]<-rainbow(16)[2]
colors[colors=="TXFW"]<-rainbow(16)[3]
colors[colors=="TXCB"]<-rainbow(16)[4]
colors[colors=="LAFW"]<-rainbow(16)[5]
colors[colors=="ALST"]<-rainbow(16)[6]
colors[colors=="ALFW"]<-rainbow(16)[7]
colors[colors=="FLSG"]<-rainbow(16)[8]
colors[colors=="FLKB"]<-rainbow(16)[9]
colors[colors=="FLFD"]<-rainbow(16)[10]
colors[colors=="FLSI"]<-rainbow(16)[11]
colors[colors=="FLAB"]<-rainbow(16)[12]
colors[colors=="FLPB"]<-rainbow(16)[13]
colors[colors=="FLHB"]<-rainbow(16)[14]
colors[colors=="FLCC"]<-rainbow(16)[15]
colors[colors=="FLLG"]<-rainbow(16)[16]

jpeg("subset.adegenet.pca1.2.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,lwd=1.3,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("bottomright", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(16), 0.5), ncol=4)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

jpeg("subset.adegenet.pca1.3.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,3], pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("bottomleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(16), 0.5), ncol=4)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC3: ", round(pca1$eig[3]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

#discriminant analysis of principal components (DAPC)
dat.clust<-find.clusters(dat.plink, parallel=FALSE, n.pca=20, n.clust=NULL,
	choose.n.clust=FALSE, max.n.clust=16)
dapc1<-dapc(dat.plink, dat.clust$grp, n.pca=20,n.da=15, parallel=F)
png("adegenet.dapc.png",height=7,width=7,units="in",res=300)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
dev.off()
compoplot(dapc1)

dapc6<-dapc(dat.plink, dat.clust$grp, n.pca=6,n.da=6, parallel=F)


#output k=3 clusters
adegenet.groups<-as.data.frame(cbind(names(dat.clust$grp), dat.clust$grp))
adegenet.groups[,1]<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', adegenet.groups[,1])
adegenet.groups[,1]<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', 
	adegenet.groups[,1])
#get the discriminant analysis loadings
adegenet.da<-merge(adegenet.groups,dapc1$ind.coord,by=0)
adegenet.da$pop<-substr(adegenet.da$V1, 1,4)
adegenet.da$pop[adegenet.da$pop=="TXSP"]<-rainbow(16)[1]
adegenet.da$pop[adegenet.da$pop=="TXCC"]<-rainbow(16)[2]
adegenet.da$pop[adegenet.da$pop=="TXFW"]<-rainbow(16)[3]
adegenet.da$pop[adegenet.da$pop=="TXCB"]<-rainbow(16)[4]
adegenet.da$pop[adegenet.da$pop=="LAFW"]<-rainbow(16)[5]
adegenet.da$pop[adegenet.da$pop=="ALST"]<-rainbow(16)[6]
adegenet.da$pop[adegenet.da$pop=="ALFW"]<-rainbow(16)[7]
adegenet.da$pop[adegenet.da$pop=="FLSG"]<-rainbow(16)[8]
adegenet.da$pop[adegenet.da$pop=="FLKB"]<-rainbow(16)[9]
adegenet.da$pop[adegenet.da$pop=="FLFD"]<-rainbow(16)[10]
adegenet.da$pop[adegenet.da$pop=="FLSI"]<-rainbow(16)[11]
adegenet.da$pop[adegenet.da$pop=="FLAB"]<-rainbow(16)[12]
adegenet.da$pop[adegenet.da$pop=="FLPB"]<-rainbow(16)[13]
adegenet.da$pop[adegenet.da$pop=="FLHB"]<-rainbow(16)[14]
adegenet.da$pop[adegenet.da$pop=="FLCC"]<-rainbow(16)[15]
adegenet.da$pop[adegenet.da$pop=="FLLG"]<-rainbow(16)[16]
adegenet.da$V2<-as.numeric(adegenet.da$V2)


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
shapes<-grp
shapes[shapes=="freshwater"]<-as.numeric(as.character(17))
shapes[shapes=="saltwater"]<-as.numeric(as.character(15))
pa.props<-round((pa$singular.values/sum(pa$singular.values))*100,2)
png("pcadapt.pc1-6.png",height=8,width=10.5,units="in",res=300)
par(mfrow=c(2,3),oma=c(2,2,2,2))
plot(pa$scores[,1],pa$scores[,2],col=alpha(colors,0.5),pch=as.numeric(shapes),
	cex=1.5,xlab=paste("PC1 (",pa.props[1],"%)",sep=""),
	ylab=paste("PC2 (",pa.props[2],"%)",sep=""))
plot(pa$scores[,3],pa$scores[,4],col=alpha(colors,0.5),pch=as.numeric(shapes),
	cex=1.5,xlab=paste("PC3 (",pa.props[3],"%)",sep=""),
	ylab=paste("PC4 (",pa.props[4],"%)",sep=""))
plot(pa$scores[,5],pa$scores[,6],col=alpha(colors,0.5),pch=as.numeric(shapes),
	cex=1.5,xlab=paste("PC5 (",pa.props[5],"%)",sep=""),
	ylab=paste("PC6 (",pa.props[6],"%)",sep=""))
plot(pa$scores[grp=="freshwater",1],pa$scores[grp=="freshwater",2],
	col=alpha(colors[grp=="freshwater"],0.5),pch=17,
	cex=1.5,xlab=paste("PC1 (",pa.props[1],"%)",sep=""),
	ylab=paste("PC2 (",pa.props[2],"%)",sep=""))
plot(pa$scores[grp=="freshwater",3],pa$scores[grp=="freshwater",4],
	col=alpha(colors[grp=="freshwater"],0.5),pch=17,
	cex=1.5,xlab=paste("PC3 (",pa.props[3],"%)",sep=""),
	ylab=paste("PC4 (",pa.props[4],"%)",sep=""))
plot(pa$scores[grp=="freshwater",5],pa$scores[grp=="freshwater",6],
	col=alpha(colors[grp=="freshwater"],0.5),pch=17,
	cex=1.5,xlab=paste("PC5 (",pa.props[5],"%)",sep=""),
	ylab=paste("PC6 (",pa.props[2],"%)",sep=""))

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(x=-0.5,y=0, pop.list, pch=c(15,15,17,15,17,15,17,rep(15,8),17), 
	pt.cex=1.5,	col=alpha(rainbow(16), 0.5), ncol=8,bty='n')
dev.off()

###### STRUCTURE #####

setwd("../../")
xcol<-c(rep("black",2),"#2166ac","black","#2166ac","black","#2166ac",
        rep("black",8),"#2166ac")
all.colors<-c('#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837')

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
                   colors=all.colors[c(1,6)],xlabel=F,
                   ylabel=expression(atop(italic(K)==2,Delta~italic(K)==358.9)))
plotting.structure(str6,2,pop.list, make.file=FALSE, plot.new=F,
                   colors=all.colors,xlabel=T,xlabcol = xcol,
                   ylabel=expression(atop(italic(K)==6,Delta~italic(K)==326.1)))
dev.off()
