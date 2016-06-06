#Author: Sarah P. Flanagan
#Last updated: 8 May 2016
#Date: 8 May 2016
#Purpose: Analyze Nerophis ophidion Population genetics data 

rm(list=ls())

library(ade4)
library(adegenet)

setwd("B:/ubuntushare/popgen/nerophis/")
source("../scripts/plotting_functions.R")

pop.list<-c("SEW","LEM","GEL","STR","GTL","FIN")
#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################

pwise.fst.sub<-read.table("stacks/fst_summary_subset.txt",
	 header=T, row.names=1, sep='\t')


#########################################################################
#***********************************************************************#
########################POPULATION STRUCTURE#############################
#***********************************************************************#
#########################################################################

#******************************ADEGENET*********************************#
dat.plink<-read.PLINK("stacks/subset.raw",parallel=FALSE)
#look at alleles
png("Missingness.png",height=7, width=7,units="in",res=300)
glPlot(dat.plink, posi="topleft")
dev.off()
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
pca1<-glPca(dat.plink, parallel=FALSE,nf=5)
scatter(pca1)

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")

#or my custom pca plotting skillzz
ind.names<-dimnames(pca1$scores)[[1]]

pop<-substr(ind.names, 1,3)
colors<-pop
colors[colors=="FIN"]<-rainbow(6)[1]
colors[colors=="GEL"]<-rainbow(6)[2]
colors[colors=="GTL"]<-rainbow(6)[3]
colors[colors=="LEM"]<-rainbow(6)[4]
colors[colors=="SEW"]<-rainbow(6)[5]
colors[colors=="STR"]<-rainbow(6)[6]
pop.list<-levels(as.factor(pop))

jpeg("subset.adegenet.pca1.2.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,lwd=1.3,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("topleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(6), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

jpeg("subset.adegenet.pca1.3.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,3], pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("topleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(6), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC3: ", round(pca1$eig[3]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

#discriminant analysis of principal components (DAPC)
dat.clust<-find.clusters(dat.plink, parallel=FALSE, n.pca=20, n.clust=NULL,
	choose.n.clust=FALSE, max.n.clust=6)#
dapc1<-dapc(dat.plink, dat.clust$grp, n.pca=20,n.da=3, parallel=F)
png("E:/Docs/PopGen/adegenet.dapc.png",height=7,width=7,units="in",res=300)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
dev.off()
compoplot(dapc1)

#output k=5 clusters
adegenet.groups<-as.data.frame(cbind(names(dat.clust$grp), dat.clust$grp))
#get the discriminant analysis loadings
adegenet.da<-merge(adegenet.groups,dapc1$ind.coord,by=0)
adegenet.da$pop<-substr(adegenet.da$V1, 1,3)
adegenet.da$colors[adegenet.da$pop=="FIN"]<-rainbow(6)[1]
adegenet.da$colors[adegenet.da$pop=="GEL"]<-rainbow(6)[2]
adegenet.da$colors[adegenet.da$pop=="GTL"]<-rainbow(6)[3]
adegenet.da$colors[adegenet.da$pop=="LEM"]<-rainbow(6)[4]
adegenet.da$colors[adegenet.da$pop=="SEW"]<-rainbow(6)[5]
adegenet.da$colors[adegenet.da$pop=="STR"]<-rainbow(6)[6]

adegenet.da$shape<-as.numeric(adegenet.da$V2)
adegenet.da$shape[adegenet.da$V2=="1"]<-21
adegenet.da$shape[adegenet.da$V2=="2"]<-22
adegenet.da$shape[adegenet.da$V2=="3"]<-23
adegenet.da$shape[adegenet.da$V2=="4"]<-24
adegenet.da$shape[adegenet.da$V2=="5"]<-25

plot(adegenet.da$LD1,adegenet.da$LD2,pch=as.numeric(adegenet.da$shape),
	bg=alpha(adegenet.da$colors,0.5),col=alpha("black",0.5),ylab="",xlab="",cex=2)
legend("topright",pch=c(21,22,23,24,25),pt.cex=2,
	c("Group 1","Group 2","Group 3", "Group 4", "Group 5"),
	col=alpha("black",0.5),ncol=3)
mtext("Discriminant Axis 1",1,line=2,cex=1.3)
mtext("Discriminant Axis 2",2,line=1.5,las=0,cex=1.3)
text(x=-8,y=4.5,"Adegenet,\nBest K = 5")

par(fig=c(0.5,1,0,0.45),new=T)
pca1$pops<-substr(rownames(pca1$scores),1,3)
pca1$colors[pca1$pops=="FIN"]<-rainbow(6)[1]
pca1$colors[pca1$pops=="GEL"]<-rainbow(6)[2]
pca1$colors[pca1$pops=="GTL"]<-rainbow(6)[3]
pca1$colors[pca1$pops=="LEM"]<-rainbow(6)[4]
pca1$colors[pca1$pops=="SEW"]<-rainbow(6)[5]
pca1$colors[pca1$pops=="STR"]<-rainbow(6)[6]
plot(pca1$scores[,1], pca1$scores[,2], pch=21, cex=2,lwd=1.3,
	bg=alpha(pca1$colors, 0.5),col=alpha("black",0.5), ylab="", xlab="")
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2,cex=1.3)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 1.5,las=0,cex=1.3)
text(x=1,y=2.8,"Adegenet,\nBest K = 3")

#*****************************STRUCTURE***********************************#
str.path<-"structure/nop_str/admix/Results/"
structure.k2<-read.table(paste(str.path,"admix_run_2_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k2$V1<-sub('([A-Z]{3}).*','\\1', structure.k2$V1)
structure.k3<-read.table(paste(str.path,"admix_run_3_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k3$V1<-sub('([A-Z]{3}).*','\\1', structure.k3$V1)
structure.k4<-read.table(paste(str.path,"admix_run_4_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k4$V1<-sub('([A-Z]{3}).*','\\1', structure.k4$V1)
structure.k5<-read.table(paste(str.path,"admix_run_5_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k5$V1<-sub('([A-Z]{3}).*','\\1', structure.k5$V1)
structure.k6<-read.table(paste(str.path,"admix_run_6_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k6$V1<-sub('([A-Z]{3}).*','\\1', structure.k6$V1)


all.colors<-c("palegreen","goldenrod1","orchid3","tomato","darkblue", "forest green")

tapply(structure.k2$V2,structure.k2$V1,max) #V2 has FLCC group
str2<-data.frame(structure.k2$V1,structure.k2$V3, structure.k2$V2)

tapply(structure.k3$V2,structure.k3$V1,max) #V2 if FLFD, V3 is TX, V4 is FLCC
str3<-data.frame(structure.k3$V1,structure.k3$V3, structure.k3$V2, structure.k3$V4)

tapply(structure.k4$V2,structure.k4$V1,max)#V2=TXCB,V3=FLSG,V4=FLCC,V5=TXCC
str4<-data.frame(structure.k4$V1,structure.k4$V5,structure.k4$V2,structure.k4$V3, 
	structure.k4$V4)

tapply(structure.k5$V2,structure.k5$V1,max)#V2=TX,V3=FLKB,V4=TXCB,V5=FLAB,V6=FLCC
str5<-data.frame(structure.k5$V1,structure.k5$V2, structure.k5$V4,structure.k5$V3,
	structure.k5$V5,structure.k5$V6)

tapply(structure.k6$V2,structure.k6$V1,max)#V2=TX,V3=FLKB,V4=TXCB,V5=FLAB,V6=FLCC
str6<-data.frame(structure.k6$V1,structure.k6$V2, structure.k6$V4,structure.k6$V3,
	structure.k6$V5,structure.k6$V6)

png("structure_k2-6.png",height=10,width=7,units="in",res=300)
par(mfrow=c(5,length(pop.list)),mar=c(0.5,0,1,0),oma=c(1,3,1,0))
plotting.structure(str2,2,pop.list, make.file=FALSE, 
	colors=all.colors[c(1,5)],xlabel=F,ylabel="STRUCTURE\nK=2")
plotting.structure(str3,3,pop.list, make.file=FALSE, 
	colors=all.colors[c(1,3,5)],xlabel=F,ylabel="STRUCTURE\nK=3")
plotting.structure(str4,4,pop.list, make.file=FALSE,
	colors=all.colors[c(1,2,3,5)],xlabel=F,ylabel="STRUCTURE\nK=4")
plotting.structure(str5,5,pop.list, make.file=FALSE,
	colors=all.colors[c(1,2,3,5,6)],xlabel=F,ylabel="STRUCTURE\nK=5")
plotting.structure(str6,6,pop.list, make.file=FALSE, colors=all.colors,
	xlabel=F,ylabel="STRUCTURE\nK=6")
dev.off()
#It looks like k=2 is best...they are kind of panmictic.

