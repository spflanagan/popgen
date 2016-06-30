#Author: Sarah P. Flanagan
#Last updated: 4 May 2016
#Date: 26 March 2015
#Purpose: Analyze Population genetics data
#Re-analyses in response to reviewer comments

rm(list=ls())

library(ade4)
library(lme4)
library(maps);library(gplots)
library(mapdata)
library(vegan)
library(boot)
library(adegenet)
library(scales)
library(gdata);library(matrixcalc)
library(psych)

setwd("E:/ubuntushare/popgen/sw_results/")
source("../scripts/plotting_functions.R")
source("../phenotype_functions.R")

pop.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
	"FLAB","FLPB","FLHB","FLCC")
#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################
summ.dat<-read.table("sw_results/stacks/populations/batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
ld.hwe<-read.table("stacks/populations/ld.hwe.whitelist.txt")

catalog<-read.delim("stacks/batch_1.catalog.tags.tsv", 
	header=F)

mar.coor<-read.csv("F://Docs//PopGen//marine_coordinates_revised.csv", header=T)
m.f.summ.dat<-read.table("stacks//populations_sex//batch_1.sumstats.tsv",
	sep='\t', skip=2, header=T, comment.char="")
dist<-read.table("F://Docs//PopGen//geographical_distances.txt", 
	header=T, row.names=1, sep='\t')
pwise.fst.all<-read.table("stacks/populations/fst_summary_all.txt",
	 header=T, row.names=1, sep='\t')
pwise.fst.sub<-read.table("stacks/populations_subset/fst_summary_subset.txt",
	 header=T, row.names=1, sep='\t')
#####Re-name plink files so that Family ID contains Population ID.
sub.ped<-read.table("stacks/populations_subset/batch_1.plink.ped")
sub.ped$V1<-gsub("sample_(\\w{4})\\w+.*align","\\1",sub.ped$V2)
write.table(sub.ped,"migrate/subset.ped",col.names=F,row.name=F,quote=F)
all.map<-read.table("stacks/populations/batch_1.plink.map")
sub.map<-read.table("stacks/populations_subset/batch_1.plink.map")
sub.scaffs<-all.map[all.map$V2 %in% sub.map$V2,]#not in the correct order!

raw.pheno<-read.table("popgen.pheno.txt", sep="\t", header=T)
fem.pheno<-read.table("fem.pheno.txt", sep="\t", header=T)
	fem.pheno$TailLength<-fem.pheno$std.length-fem.pheno$SVL
	fem.pheno$HeadLength<-fem.pheno$HeadLength-fem.pheno$SnoutLength
	fem.pheno<-fem.pheno[,c(1,2,3,11,5,6,7,8,9,10)]
mal.pheno<-read.table("mal.pheno.txt", sep="\t", header=T)
	mal.pheno$TailLength<-mal.pheno$std.length-mal.pheno$SVL
	mal.pheno$HeadLength<-mal.pheno$HeadLength-mal.pheno$SnoutLength
	mal.pheno<-mal.pheno[,c(1,2,3,9,5,6,7,8)]

#############################################################################
#######################PLOT THE POINTS ON A MAP##############################
#############################################################################
jpeg("mar_sites_map_again.jpg", res=300, height=7,width=14, units="in")
pdf("marine_sites_map.pdf",height=7,width=14)
par(oma=c(0,0,0,0),mar=c(0,0,0,0),pin=c(7,7))
map("worldHires", "usa",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray90", mar=c(0,0,0,0),fill=TRUE, res=300,myborder=0)
map("worldHires", "mexico",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray95", fill=TRUE, add=TRUE)
points(mar.coor$lon, mar.coor$lat,  col="black", cex=1.2, pch=19)
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
text(x=-94.7,y=29,"TXCB",font=2)
text(x=-88,y=30,"ALST",font=2)
text(x=-85,y=29.4,"FLSG",font=2)
text(x=-83.5,y=29.2,"FLKB",font=2)
text(x=-83.2,y=27.6,"FLFD",font=2)
text(x=-82.2,y=26,"FLSI",font=2)
text(x=-80,y=24.8,"FLAB",font=2)
text(x=-79.5,y=26.8,"FLPB",font=2)
text(x=-79.7,y=27.2,"FLHB",font=2)
text(x=-80.2,y=28.2,"FLCC",font=2)
dev.off()

##############################################################################
#****************************************************************************#
#########################TEST FOR ISOLATION BY DISTANCE#######################
#****************************************************************************#
##############################################################################
#Mantel test using geographical distances and fsts

#read in the subsetted fst summary from running populations with whitelist
ibd.all<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst.all)))
ibd.sub<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst.sub)))

######BY LOCUS#############
pairwise.fst<-function(ped,allele1,allele2,pop.order){
	#V1 of ped should be pop index
	ped.split<-split(ped[,c(allele1,allele2)], factor(ped[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(pop.order),numeric(0), simplify = F), pop.order))
	for(i in 1:(length(pop.order)-1)){
	  for(j in (i+1):length(pop.order)){
		pop1<-factor(ped.split[[pop.order[i]]][ped.split[[pop.order[i]]]!="0"])
		pop2<-factor(ped.split[[pop.order[j]]][ped.split[[pop.order[j]]]!="0"])
		freq1<-summary(pop1)/sum(summary(pop1))	
		freq2<-summary(pop2)/sum(summary(pop2))	
		freqall<-summary(as.factor(c(pop1,pop2)))/
			sum(summary(as.factor(c(pop1,pop2))))
		if(length(freq1)>1){ hs1<-2*freq1[1]*freq1[2] 
		} else {
			hs1<-0
		}
		if(length(freq2)>1){ hs2<-2*freq2[1]*freq2[2] 
		} else {
			hs2<-0
		}
		if(length(freqall)>1){
			hs<-mean(c(hs1,hs2))
			ht<-2*freqall[1]*freqall[2]
			fst<-(ht-hs)/ht
		}
		if(length(freqall)<=1){ fst<-1 }
		dat.var[pop.order[i],pop.order[j]]<-fst
	  }
	}
	dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
	rownames(dat.var)<-colnames(dat.var)
	return(as.matrix(dat.var))
}


pairwise.fst(sub.ped,9,10,pop.list)

fst.ibd.byloc<-function(ped.file,dist.mat,pop.order){
	results.mantel<-data.frame()
	for(i in seq(7,ncol(ped.file),2)){
		res<-mantel.rtest(
			as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
			as.dist(t(dist.mat)), nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	colnames(results.mantel)<-c("Obs","P")
	return(results.mantel)
}

ibd.by.loc<-fst.ibd.byloc(sub.ped,dist,pop.list) 
#ignore warnings?  In is.euclid(m1) : Zero distance(s)
rownames(ibd.by.loc)<-sub.map$V2

#fst.test<-as.matrix(pairwise.fst(ped,7,8,pop.order))
#fst.test[which.min(abs(fst.test))]<-0.0001
#mantel.rtest(as.dist(t(fst.test)),
#			as.dist(t(dist)), nrepet=9999)
#is.euclid(as.dist(t(fst.test)))

#########################################################################
#***********************************************************************#
########################POPULATION STRUCTURE#############################
#***********************************************************************#
#########################################################################

#******************************ADEGENET*********************************#
dat.plink<-read.PLINK("stacks/populations/subsetA.raw",
	parallel=FALSE)
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
pca1<-glPca(dat.plink, parallel=FALSE,nf=5)
scatter(pca1)

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")

#or my custom pca plotting skillzz
ind.names<-dimnames(pca1$scores)[[1]]

pop<-substr(ind.names, 8,11)
colors<-pop
colors[colors=="TXSP"]<-rainbow(12)[1]
colors[colors=="TXCC"]<-rainbow(12)[2]
colors[colors=="TXCB"]<-rainbow(12)[3]
colors[colors=="ALST"]<-rainbow(12)[4]
colors[colors=="FLSG"]<-rainbow(12)[5]
colors[colors=="FLKB"]<-rainbow(12)[6]
colors[colors=="FLFD"]<-rainbow(12)[7]
colors[colors=="FLSI"]<-rainbow(12)[8]
colors[colors=="FLAB"]<-rainbow(12)[9]
colors[colors=="FLPB"]<-rainbow(12)[10]
colors[colors=="FLHB"]<-rainbow(12)[11]
colors[colors=="FLCC"]<-rainbow(12)[12]

jpeg("E://Docs//PopGen//subset.pca1.2.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,lwd=1.3,
	col=alpha(col, 0.5), ylab="", xlab="")
legend("topleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

jpeg("E://Docs//PopGen//subset.pca1.3.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,3], pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="")
legend("topleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC3: ", round(pca1$eig[3]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

#discriminant analysis of principal components (DAPC)
dat.clust<-find.clusters(dat.plink, parallel=FALSE, n.pca=20, n.clust=NULL,
	choose.n.clust=FALSE, max.n.clust=12)#created 3 clusters
dapc1<-dapc(dat.plink, dat.clust$grp, n.pca=20,n.da=3, parallel=F)
png("E:/Docs/PopGen/adegenet.dapc.png",height=7,width=7,units="in",res=300)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
dev.off()
compoplot(dapc1)

#output k=3 clusters
adegenet.groups<-as.data.frame(cbind(names(dat.clust$grp), dat.clust$grp))
adegenet.groups[,1]<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', adegenet.groups[,1])
adegenet.groups[,1]<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', 
	adegenet.groups[,1])
#get the discriminant analysis loadings
adegenet.da<-merge(adegenet.groups,dapc1$ind.coord,by=0)
adegenet.da$pop<-substr(adegenet.da$V1, 1,4)
adegenet.da$pop[adegenet.da$pop=="TXSP"]<-rainbow(12)[1]
adegenet.da$pop[adegenet.da$pop=="TXCC"]<-rainbow(12)[2]
adegenet.da$pop[adegenet.da$pop=="TXCB"]<-rainbow(12)[3]
adegenet.da$pop[adegenet.da$pop=="ALST"]<-rainbow(12)[4]
adegenet.da$pop[adegenet.da$pop=="FLSG"]<-rainbow(12)[5]
adegenet.da$pop[adegenet.da$pop=="FLKB"]<-rainbow(12)[6]
adegenet.da$pop[adegenet.da$pop=="FLFD"]<-rainbow(12)[7]
adegenet.da$pop[adegenet.da$pop=="FLSI"]<-rainbow(12)[8]
adegenet.da$pop[adegenet.da$pop=="FLAB"]<-rainbow(12)[9]
adegenet.da$pop[adegenet.da$pop=="FLPB"]<-rainbow(12)[10]
adegenet.da$pop[adegenet.da$pop=="FLHB"]<-rainbow(12)[11]
adegenet.da$pop[adegenet.da$pop=="FLCC"]<-rainbow(12)[12]
adegenet.da$V2<-as.numeric(adegenet.da$V2)
adegenet.da$V2[adegenet.da$V2=="1"]<-as.numeric(15)
adegenet.da$V2[adegenet.da$V2=="2"]<-16
adegenet.da$V2[adegenet.da$V2=="3"]<-17

#*******************************PCADAPT***********************************#
#K=4 WAS BEST
setwd("E:/ubuntushare/popgen/sw_results/pcadapt/pruned")
scores.files<-list.files(pattern="4_.*.scores")
loadings.files<-list.files(pattern="4_.*.loadings")
snps.files<-list.files(pattern="4_.*.topBF")

snp.list<-list()
for(i in 1: length(snps.files)){
	#read in files
	snp.list[[i]]<-read.table(snps.files[i],header=T)
}
#compare lists of snps in all of the runs

all.snps<-as.vector(sapply(snp.list, "[[","snp"))
all.snps.dup<-all.snps[duplicated(all.snps)]
rep.snps<-all.snps.dup[!duplicated(all.snps.dup)]

scores<-read.table("pcadapt.4.scores")
loading<-read.table("pcadapt.4", header=T, sep="\t")
#bf.log<-log10(snp.list[[1]]$BF)
#loading[round(loading$logBF, 4) %in% round(bf.log),] #doesn't work..

#the snp is the row number in the map file for the snps
sub.map<-read.table("../../stacks/populations/subset.map",header=F)
pcadapt.outliers<-sub.map[rep.snps,]
pa.out.radloc<-sub('(\\d+)_\\d+','\\1',pcadapt.outliers$V2)

summ.dat<-read.table("../../stacks/populations/batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
pa.out.dat<-summ.dat[summ.dat$Locus.ID %in% pa.out.radloc,]
pa.out.dat$Chr<-factor(pa.out.dat$Chr)
#are any on the linkage map?
#length(levels(as.factor(use.contigs[use.contigs$Scaffold %in% pa.out.dat$Chr,1])))

#plot individual scores
jpeg("E:/Docs/PopGen/pcadapt.scores1.2.jpeg", height=12, width=12, units="in", res=300)
plot(as.numeric(scores[1,]),as.numeric(scores[2,]),pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="",lwd=1.3)
legend("bottomleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0117)", 1, line = 2)
mtext("PC2 (rho2: 0.0091)",2, line = 2)
dev.off()

jpeg("E:/Docs/PopGen/pcadapt.scores1.3.jpeg", height=12, width=12, units="in", res=300)
plot(as.numeric(scores[1,]),as.numeric(scores[3,]),pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="")
legend("bottomleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0117)", 1, line = 2)
mtext("PC3 (rho2: 0.0018)",2, line = 2)
dev.off()

######################MAKE FIGURE 3################################
jpeg("F:/Docs/PopGen/Figure3_revised.jpeg", height=8, width=8, units="in",
	 res=300)
par(lwd=1.3,cex=1.3,mfrow=c(2,2),oma=c(2,2,2,1),mar=c(2,2,2,1),las=1)
#PCAdapt
plot(as.numeric(scores[1,]),as.numeric(scores[2,]),pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="",lwd=1.3)
#legend("bottomleft", pop.list, pch=19, pt.cex=2,
#	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0175)", 1, line = 2)
mtext("PC2 (rho2: 0.0131)",2, line = 1.5,las=0)
text(x=-3,y=2.5,"PCAdapt,\nBest K = 4")

plot(as.numeric(scores[1,]),as.numeric(scores[3,]),pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="")
#legend("bottomleft", pop.list, pch=19, pt.cex=2,
#	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0175)", 1, line = 2)
mtext("PC3 (rho2: 0.0019)",2, line = 1.5,las=0)
text(x=-3,y=3,"PCAdapt,\nBest K = 4")

#Adegenet
plot(adegenet.da$LD1,adegenet.da$LD2,pch=as.numeric(adegenet.da$V2),
	col=alpha(adegenet.da$pop,0.5),ylab="",xlab="",cex=2)
legend("bottomleft",pch=c(15,16,17),pt.cex=2,c("Group 1","Group 2","Group 3"),
	col=alpha("black",0.5))
mtext("Discriminant Axis 1",1,line=2)
mtext("Discriminant Axis 2",2,line=1.5,las=0)
text(x=-8,y=4.5,"Adegenet,\nBest K = 3")

plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,lwd=1.3,
	col=alpha(col, 0.5), ylab="", xlab="")
legend("bottomleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 1.5,las=0)
text(x=-3.5,y=2.8,"Adegenet,\nBest K = 3")
dev.off()


#*****************************STRUCTURE***********************************#
setwd("../../")
structure.k2<-read.table(
	"structure//popgen//admixture//Results//admixture_run_11_f_clusters.txt",
	sep='\t', header=F)
structure.k2$V1<-sub('sample_([A-Z]{4})','\\1', structure.k2$V1)
structure.k3<-read.table(
	"structure//popgen//admixture//Results//admixture_run_21_f_clusters.txt",
	sep='\t', header=F)
structure.k3$V1<-sub('sample_([A-Z]{4})','\\1', structure.k3$V1)
structure.k4<-read.table(
	"structure//popgen//admixture//Results//admixture_run_31_f_clusters.txt",
	sep='\t', header=F)
structure.k4$V1<-sub('sample_([A-Z]{4})','\\1', structure.k4$V1)
structure.k5<-read.table(
	"structure//popgen//admixture//Results//admixture_run_41_f_clusters.txt",
	sep='\t', header=F)
structure.k5$V1<-sub('sample_([A-Z]{4})','\\1', structure.k5$V1)



all.colors<-c("palegreen","goldenrod1","orchid3","tomato","darkblue")

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


##################################MAKE FIGURE 2#############################
jpeg("Figure2_revisions.jpeg",height=7,width=7.5,units="in",res=300)
par(mfrow=c(4,length(pop.list)),mar=c(0.5,0,1,0),oma=c(1,3,1,0))

plotting.structure(str2,2,pop.list, make.file=FALSE, 
	colors=all.colors[c(1,5)],xlabel=F,ylabel="STRUCTURE\nK=2")
#plotting.structure(fstr2,2,pop.list, make.file=FALSE, 
#	colors=all.colors[c(1,5)],xlabel=F,ylabel="FAST\nK=2")
plotting.structure(str3,3,pop.list, make.file=FALSE, 
	colors=all.colors[c(1,3,5)],xlabel=F,ylabel="STRUCTURE\nK=3")
#plotting.structure(fstr3,3,pop.list, make.file=FALSE, 
#	colors=all.colors[c(1,3,5)],xlabel=F,ylabel="FAST\nK=3")
plotting.structure(str4,4,pop.list, make.file=FALSE,
	colors=all.colors[c(1,2,3,5)],xlabel=F,ylabel="STRUCTURE\nK=4")
#plotting.structure(fstr4,4,pop.list, make.file=FALSE,
#	colors=all.colors[c(1,2,3,5)],xlabel=F,ylabel="FAST\nK=4")
plotting.structure(str5,5,pop.list, make.file=FALSE, colors=all.colors,
	xlabel=F,ylabel="STRUCTURE\nK=5")
#plotting.structure(fstr5,5,pop.list, make.file=FALSE, colors=all.colors
#	,xlabel=T,ylabel="FAST\nK=5")
dev.off()
########################MAKE POP STRUCTURE FIGURE#############################
jpeg("PopStructure_revised.jpeg", height=12, width=12, units="in",
	 res=300)
pdf("PopStructure_revised.pdf", height=12, width=12)
par(lwd=1.3,cex=1.3,fig=c(0,0.5,0.45,0.9),oma=c(2,2,2,1),mar=c(2,2,2,1),
	las=1,new=F)
#PCAdapt
plot(as.numeric(scores[1,]),as.numeric(scores[2,]),pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="",lwd=1.3)
mtext("PC1 (rho2: 0.0175)", 1, line = 2,cex=1.3)
mtext("PC2 (rho2: 0.0131)",2, line = 2,las=0,cex=1.3)
text(x=-3,y=2.5,"PCAdapt,\nBest K = 4")
legend("bottomleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12),0.5), ncol=2)

par(fig=c(0.5,1,0.45,0.9),new=T)
plot(as.numeric(scores[1,]),as.numeric(scores[3,]),pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="")
mtext("PC1 (rho2: 0.0175)", 1, line = 2,cex=1.3)
mtext("PC3 (rho2: 0.0019)",2, line = 2,las=0,cex=1.3)
text(x=-3,y=3,"PCAdapt,\nBest K = 4")

#Adegenet
par(fig=c(0,0.5,0,0.45),new=T)
plot(adegenet.da$LD1,adegenet.da$LD2,pch=as.numeric(adegenet.da$V2),
	col=alpha(adegenet.da$pop,0.5),ylab="",xlab="",cex=2)
legend("bottomleft",pch=c(15,16,17),pt.cex=2,c("Group 1","Group 2","Group 3"),
	col=alpha("black",0.5),ncol=3)
mtext("Discriminant Axis 1",1,line=2,cex=1.3)
mtext("Discriminant Axis 2",2,line=1.5,las=0,cex=1.3)
text(x=-8,y=4.5,"Adegenet,\nBest K = 3")

par(fig=c(0.5,1,0,0.45),new=T)
plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,lwd=1.3,
	col=alpha(pop.colors, 0.5), ylab="", xlab="")
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2,cex=1.3)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 1.5,las=0,cex=1.3)
text(x=1,y=2.8,"Adegenet,\nBest K = 3")

par(mar=c(1,0,0,0), 
		oma=c(1,2,2,2),cex=0.5)
fig.start<-seq(0,11,1)/12
fig.end<-seq(1,12,1)/12
str.split<-split(str5,str5[,1])
for(i in 1:length(str.split)){
	pop.index<-pop.list[i]
	par(fig=c(fig.start[i],fig.end[i],0.9,1),mar=c(1,0,0,0), 
		oma=c(1,2,2,2),cex=0.5,new=T)
	barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
		beside=FALSE, space=0,	border=NA, col=all.colors,
		xlab="", ylab="", xaxt='n', yaxt='n')
	mtext(pop.index, 1, line=1, cex=1.3, outer=F)	
}
dev.off()


##############################################################################
#****************************************************************************#
##############################OUTLIER ANALYSES################################
#****************************************************************************#
##############################################################################

#**********************************BAYENV2***********************************#
#####STARTING WITH PLINK FILES
ped<-read.table("/results/stacks/populations/subset.ped", 
	stringsAsFactors=F, colClasses="character")
ped.pops<-substr(ped[,2],8,11)
ped.sex<-sub('sample_\\w{4}(\\w+).*[_.].*','\\1', ped[,2])
ped.sex[substr(ped.sex,1,2)=="DP"]<-"M"
ped.sex[substr(ped.sex,1,2)=="DB"]<-"F"
ped.sex[substr(ped.sex,1,2)=="NP"]<-"M"
ped.sex[substr(ped.sex,1,1)=="P"]<-"M"
ped.sex<-substr(ped.sex,1,1)
ped.sex[ped.sex=="F"]<-2
ped.sex[ped.sex=="M"]<-1

ped[,1]<-ped.pops
ped[,5]<-ped.sex

write.table(ped,"results/bayenv/bayenv.plink.ped", 
	row.names=F, col.names=F, quote=F, sep="\t",eol="\n")

clust.plink<-cbind(ped.pops, ped[,2],ped.pops)
write.table(clust.plink, 
	"results/stacks/populations/plink.clust.txt",
	col.names=F, row.names=F, quote=F, sep="\t", eol="\n")


#plink.map -> numbers instead of scaffold_#
sub.map<-read.table("stacks/populations/subset.map", skip = 1)
chr.nums<-sub('scaffold_','',sub.map[,1])
map[,1]<-chr.nums

write.table(map, "bayenv.plink.map", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n")

##get any snps *not* in plink.ld
plink.ld<-read.table("results/stacks/populations/plink.ld",
	header=T)
all.snps<-c(as.character(plink.ld$SNP_A), as.character(plink.ld$SNP_B))
above.2<-map[map$V2 %in% all.snps,]
below.2<-map[!(map$V2 %in% all.snps),]

write.table(below.2$V2, "results/stacks/populations/ld.subset.list",
	col.names=F, row.names=F, quote=F, sep='\t', eol='\n')

#####CONVERT PLINK TO BAYENV2
freq<-read.table("results/bayenv/ld.hwe.bayenv.plink.frq.strat", 
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

write.table(snpsfile, "results/bayenv/null.1833", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n") #bayenv SNPSFILE

#####check Bayenv2 matrix
matrix.files<-list.files("results/environmental_assoc/new_bayenv/",pattern="matrix")
matrices<-list()
for(i in 1:length(matrix.files))
{
	matrices[[i]]<-as.matrix(read.table(
		paste("results/environmental_assoc/new_bayenv/",matrix.files[i],sep=""), 
		skip=2801, header=F))
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

#I took a representative matrix rather than averaging.

#####SNPFILEs
#for SNPFILE, need just one file per SNP apparently.
#want to use all of the snps...need to get map with those inds.
all.snps.ped<-read.table("results/stacks/populations/batch_1.plink.ped", header=F, stringsAsFactors=F)
ped.pop<-sub('sample_(\\w{4})\\w+.*[_.].*','\\1', all.snps.ped[,2])
all.snps.clust<-cbind(ped.pop,all.snps.ped[,2],ped.pop)
write.table(all.snps.clust, "all.6348.clust.txt", sep="\t", eol="\n", quote=F,
	row.names=F, col.names=F)
#then need to run plink --file batch_1.plink --freq --within all.6348.clust.txt \
	#--allow-no-sex --noweb --out all.bayenv.plink

#read in frequency per pop
all.snps.frq<-read.table("all.bayenv.plink.frq.strat", 
	header=T, stringsAsFactors=F)
freq<-cbind(freq,freq$NCHROBS-freq$MAC)
colnames(freq)[ncol(freq)]<-"NAC"
pop.order<-levels(as.factor(freq$CLST))
snp.names<-split(freq$SNP,freq$CLST)[[1]]

mac.by.pop<-as.data.frame(split(freq$MAC,freq$CLST))
rownames(mac.by.pop)<-snp.names

write.table(mac.by.pop, "all.6348", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n")

#####ENVFILE
env.raw<-read.table("bayenv2/new_bayenv/env_data_bayenv_raw.txt")
#Each environmental variable should be standardized, 
#i.e. subtract the mean and then divided through by the standard deviation 
#of the variable across populations.
std.by.mean<-function(x){
	m<-mean(x)
	s<-sd(x)
	newx<-(x-m)/s	
}
env.std<-t(apply(env.raw,1,std.by.mean))
write.table(env.std,
	"bayenv2/new_bayenv/env_data_bayenv_std.txt",
	sep='\t',quote=F,col.names=F,row.names=F,eol='\n')

##Are they correlated with distance?
colnames(env.std)<-colnames(env.raw)
rownames(env.std)<-rownames(env.raw)
env.dist<-as.matrix(vegdist(t(env.raw)))
env.dist<-env.dist[rownames(dist),colnames(dist)]
mantel.rtest(as.dist(t(dist)),as.dist(env.dist),999)
#Monte-Carlo test
#Observation: 0.2614484 
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
#Based on 999 replicates
#Simulated p-value: 0.061 

#####GET OUTPUT
setwd("bayenv2/new_bayenv/pruned_snps")
bf.files<-list.files(pattern="bf")
xtx.files<-list.files(pattern="xtx")

bf.dat<-NULL
for(i in 1:length(bf.files)){
	bf<-read.table(bf.files[i])
	bf.dat<-rbind(bf.dat,bf)
}
colnames(bf.dat)<-c("locus", "Temp_BF", "Temp_rho", "Temp_rs", 
	"Salinity_BF", "Salinity_rho", "Salinity_rs", "coll.temp_BF", 
	"coll.temp_rho", "coll.temp_rs", "coll.sal_BF", "coll.sal_rho", 
	"coll.sal_rs", "seagrass_BF", "seagrass_rho","seagrass_rs")
bf.dat$locus<-sub("./pruned_snps/(\\d+.*)","\\1",bf.dat$locus)
xtx<-NULL
for(i in 1:length(xtx.files)){
	x<-read.table(xtx.files[i])
	xtx<-rbind(xtx,x)
}
colnames(xtx)<-c("locus","XtX")
xtx$locus<-sub("./pruned_snps/(\\d+.*)","\\1",xtx$locus)

setwd("../")
write.table(bf.dat,"BF_summary.pruned_snps.txt",sep='\t',quote=F,row.names=F)
write.table(xtx,"XtX_summary.pruned_snps.txt",sep='\t',quote=F,row.names=F)

setwd("/new_bayenv/temp_var")
bf.files<-list.files(pattern="bf")
bf.tempvar<-NULL
for(i in 1:length(bf.files)){
	bf<-read.table(bf.files[i])
	bf.tempvar<-rbind(bf.tempvar,bf)
}
colnames(bf.tempvar)<-c("SNP","BFtempvar","rtempvar","rStempvar")
bf.tempvar$SNP<-sub("./temp_var/(\\d+.*)","\\1",bf.tempvar$SNP)
write.table(bf.tempvar,"../tempvar_bf.txt",quote=F)
setwd("../../")

bf.dat<-merge(bf.dat,bf.tempvar,by.x="locus",by.y="SNP")
#####BAYENV: ENV
bf.scaff<-merge(sub.map, bf.dat, by.x="V2", by.y="locus")
colnames(bf.scaff)[1:4]<-c("locus","scaffold","dist","BP")
#focus on Bayes Factors, because of Lotterhos & Whitlock (2015)
bf<-bf.scaff[,c(1,2,4,5,8,11,14,17,20)]
bf.co<-apply(bf[,4:9],2,quantile,0.95)
temp.bf.sig<-bf[bf$Temp_BF>bf.co["Temp_BF"],c(1,2,3,4)]
sal.bf.sig<-bf[bf$Salinity_BF>bf.co["Salinity_BF"],c(1,2,3,5)]
ctemp.bf.sig<-bf[bf$coll.temp_BF>bf.co["coll.temp_BF"],c(1,2,3,6)]
csal.bf.sig<-bf[bf$coll.sal_BF>bf.co["coll.sal_BF"],c(1,2,3,7)]
grass.bf.sig<-bf[bf$seagrass_BF>bf.co["seagrass_BF"],c(1,2,3,8)]
tvar.bf.sig<-bf[bf$BFtempvar>bf.co["BFtempvar"],c(1,2,3,9)]

dim(tvar.bf.sig[tvar.bf.sig$locus %in% temp.bf.sig$locus & 
	tvar.bf.sig$locus %in% ctemp.bf.sig,])
dim(sal.bf.sig[sal.bf.sig$locus %in% csal.bf.sig$locus,])
dim(grass.bf.sig[grass.bf.sig$locus %in% temp.bf.sig$locus &
	grass.bf.sig$locus %in% sal.bf.sig$locus,])
dim(grass.bf.sig[grass.bf.sig$locus %in% ctemp.bf.sig$locus &
	grass.bf.sig$locus %in% csal.bf.sig$locus,])
out.venn<-venn( list("MeanSalinity"=sal.bf.sig$locus,
	"MeanTemp"=temp.bf.sig$locus,"Seagrass"=grass.bf.sig$locus))

temp.bf.sig$start<-temp.bf.sig$BP-2500
temp.bf.sig$end<-temp.bf.sig$BP+2500
sal.bf.sig$start<-sal.bf.sig$BP-2500
sal.bf.sig$end<-sal.bf.sig$BP+2500
ctemp.bf.sig$start<-ctemp.bf.sig$BP-2500
ctemp.bf.sig$end<-ctemp.bf.sig$BP+2500
csal.bf.sig$start<-csal.bf.sig$BP-2500
csal.bf.sig$end<-csal.bf.sig$BP+2500
grass.bf.sig$start<-grass.bf.sig$BP-2500
grass.bf.sig$end<-grass.bf.sig$BP+2500
tvar.bf.sig$start<-tvar.bf.sig$BP-2500
tvar.bf.sig$end<-tvar.bf.sig$BP+2500
chroms<-c(as.character(temp.bf.sig$scaffold),as.character(sal.bf.sig$scaffold),
	as.character(ctemp.bf.sig$scaffold),as.character(csal.bf.sig$scaffold),
	as.character(grass.bf.sig$scaffold),as.character(tvar.bf.sig$scaffold))
chroms<-chroms[!duplicated(chroms)]
write.table(chroms,
	"bayenv2/new_bayenv/outliers/all_chroms.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')
write.table(csal.bf.sig[,c("scaffold","start","end")],
	"bayenv2/new_bayenv/outliers/csal.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')
write.table(sal.bf.sig[,c("scaffold","start","end")],
	"bayenv2/new_bayenv/outliers/sal.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')
write.table(grass.bf.sig[,c("scaffold","start","end")],
	"bayenv2/new_bayenv/outliers/grass.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')
write.table(temp.bf.sig[,c("scaffold","start","end")],
	"bayenv2/new_bayenv/outliers/temp.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')
write.table(ctemp.bf.sig[,c("scaffold","start","end")],
	"bayenv2/new_bayenv/outliers/ctemp.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')
write.table(tvar.bf.sig[,c("scaffold","start","end")],
	"bayenv2/new_bayenv/outliers/tvar.txt",col.names=F,row.names=F,
	quote=F,eol='\n',sep='\t')

tmp.sig<-tvar.bf.sig[tvar.bf.sig$locus %in% temp.bf.sig$locus,]#4
write.table(tmp.sig[,1:3],"shared_temp_outliers.txt",sep='\t',quote=F,
	row.names=F)
tmp.region<-data.frame(tmp.sig$scaffold,tmp.sig$BP-2500,tmp.sig$BP+2500)
tmp.region[,2][tmp.region[,2]<0]<-0
write.table(tmp.region,
	"bayenv2/new_bayenv/outliers/extract_shared_tmp_region.sh",
	row.names=F,col.names=F,quote=F,sep='\t')
#non-outliers: col = black, pch=19
#Temp Variance BF outliers: col=blue, pch=19 or 1
#Temp BF outliers: col=purple, pch=2

neutral.col<-"black"
bfv.out.col<-"mediumvioletred"
bfm.out.col<-"orangered"
shared.out.col<-"red"
neutral.pch<-19
bf.scaff$Temp_BF<-log(bf.scaff$Temp_BF)
bf.scaff$BFtempvar<-log(bf.scaff$BFtempvar)

#############PLOTTING
jpeg("temp_revisions.jpeg", height=9,width=7.5,units="in",res=300)
pdf("Temp_Outliers.pdf",height=9,width=7.5)
par(mfrow=c(2,1),oma=c(1,1,0,0),mar=c(0,1,1,0),cex=2,mgp=c(3,0.5,0))

###GLOBAL FSTS

bf.scaff$scaffold<-factor(bf.scaff$scaffold)
all.scaff<-split(bf.scaff, bf.scaff$scaffold)
last.max<-0
rect.xs<-NULL
addition.values<-0
for(i in 1:length(all.scaff)){
	new.max<-last.max+round(max(all.scaff[[i]]$BP), -2)
	rect.xs<-rbind(rect.xs,c(last.max, new.max))
	addition.values<-c(addition.values, new.max)
	last.max<-new.max
}
#change BP to plot
for(i in 1:length(all.scaff)){
	all.scaff[[i]]$BP<-all.scaff[[i]]$BP+addition.values[i]
}
x.max<-max(addition.values)
x.min<-min(all.scaff[[1]]$BP)
y.max<-max(bf.scaff$Temp_BF)+0.1*max(bf.scaff$Temp_BF)
y.min<-min(bf.scaff$Temp_BF)-0.1*min(bf.scaff$Temp_BF)
if(min(bf.scaff$Temp_BF) < 0) {
	y.min<-min(bf.scaff$Temp_BF) - 0.1*min(bf.scaff$Temp_BF)
} else {
	y.min<-0
}

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray75"
	}
	rect(rect.xs[i,1],-2,rect.xs[i,2],4, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plotting.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$Temp_BF,plot.rect=FALSE,
		y.max=3.5,x.max, rect.xs[i,],y.min=-1.5,x.min=x.min, 
		pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in% temp.bf.sig$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col=bfm.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in%
		tmp.sig$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col="red", pch=8,
		cex=0.75)
}
axis(2, at = seq(-2,4,1),ylim = c(-1.5, 3.5), pos=0,
	labels=seq(-2,4,1),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)

mtext(side=2, "log(Bayes Factor)", 
	outer=FALSE, line=1.2,cex=1)


#######Temp variance
y.max<-max(bf.scaff$BFtempvar)+0.1*max(bf.scaff$BFtempvar)
y.min<-min(bf.scaff$BFtempvar)-0.1*min(bf.scaff$BFtempvar)
if(min(bf.scaff$BFtempvar) < 0) {
	y.min<-min(bf.scaff$BFtempvar) - 0.1*min(bf.scaff$BFtempvar)
} else {
	y.min<-0
}

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray75"
	}
	rect(rect.xs[i,1],-3.5,rect.xs[i,2],4, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plotting.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$BFtempvar,plot.rect=FALSE,
		y.max=3.5,x.max, rect.xs[i,],y.min=-3.5,x.min=x.min, 
		pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in% tvar.bf.sig$locus,]
	points(temp.sig$BP, temp.sig$BFtempvar, col=bfv.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in%
		tmp.sig$locus,]
	points(temp.sig$BP, temp.sig$BFtempvar, col="red", pch=8,
		cex=0.75)
}
axis(2, at = seq(-3,4,1),
	ylim = c(-4,3.5), pos=0,
	labels=seq(-3,4,1),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=2, "log(Bayes Factor)", 
	outer=FALSE, line=1.2,cex=1)
mtext("Position in Genome",1,line=1.2,cex=1,outer=F)

#PLOT THE LEGEND
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("topright", ncol=3,
	col=c(bfm.out.col,bfv.out.col, "red"),
	c("Mean Temp Outliers","Temp Variance Outliers",
		"Shared Outliers"),
	pch=c(19,19,8), box.lty=0)

dev.off()




#####BAYENV: XTX
xtx$rad.loc<-sub('(\\d+)(_\\d+)','\\1',xtx$locus)
xtx.scaff<-merge(sub.map, xtx, by.x="V2", by.y="locus")
colnames(xtx.scaff)<-c("locus","scaffold","dist","BP","XtX", "rad.loc")
xtx.1<-xtx.scaff[xtx.scaff$XtX >= quantile(xtx.scaff$XtX,0.99),] #18
xtx.5<-xtx.scaff[xtx.scaff$XtX >= quantile(xtx.scaff$XtX,0.95),] #88

#****************************GLOBAL FST OUTLIERS*****************************#
#Import global values generated by C++ program calculate_global_fsts
#pruned loci found in all 12 populations, 1753 of them.
global.fsts<-read.table("stacks//populations//subset.globalstats.txt",
	header=T, sep='\t')
global.fsts<-merge(global.fsts,sub.map,by.x="Locus",by.y="V2")
global.fsts$Chrom<-global.fsts$V1
global.fsts$radloc<-sub('(\\d+)(_\\d+)','\\1',global.fsts$Locus)
#outliers:
fst.ci.99<-mean(global.fsts$Fst) + 2.57583*sd(global.fsts$Fst)
glob.fst.out<-global.fsts[global.fsts$Fst >= fst.ci.99,]#44
glob.fst.99<-global.fsts[global.fsts$Fst >=quantile(global.fsts$Fst,0.99),]#18
glob.fst.95<-global.fsts[global.fsts$Fst >=quantile(global.fsts$Fst,0.95),]#88

#**********************************PCADAPT**********************************#
#K=4 WAS BEST
all.map<-read.table("stacks//populations//all.6348.plink.map",header=F)
pca.files<-list.files("pcadapt//6348//",pattern="4_")
scores.files<-list.files("pcadapt//6348//", pattern="scores")
loadings.files<-list.files("pcadapt//6348//", pattern="loadings")
snps.files<-list.files("pcadapt//6348//", pattern="topBF")
stats.files<-list.files("pcadapt//6348//", pattern="stats")
bf.files<-pca.files[!(pca.files %in% scores.files) & 
	!(pca.files %in% loadings.files) & !(pca.files %in% snps.files) &
	!(pca.files %in% stats.files)]


bf.list<-list()
snp.list<-list()
scores.list<-list()
loadings.list<-list()
for(i in 1: length(snps.files)){
	#read in files
	snp.list[[i]]<-read.table(
		paste("pcadapt//6348//", snps.files[i], sep=""), 
		header=T)
	scores.list[[i]]<-read.table(
		paste("pcadapt//6348//", scores.files[i], sep=""), 
		header=T)
	loadings.list[[i]]<-read.table(
		paste("pcadapt//6348//", loadings.files[i], sep=""), 
		header=F)
	bf.dat<-read.table(
		paste("pcadapt//6348//", bf.files[i], sep=""), 
		skip=1)
	bf.dat<-cbind(bf.dat, all.map[rownames(bf.dat),])
	
	bf.list[[i]]<-bf.dat
}
#compare lists of snps in all of the runs

pca.scaff<-bf.list[[1]]
colnames(pca.scaff)<-c("logBF","logPO","PZ1","PZ2","PZ3","PZ4","chrom","loc",
	"dist","BP")
pca.scaff$radloc<-sub('(\\d+)_\\d+','\\1',pca.scaff$loc)

all.snps<-as.vector(sapply(snp.list, "[[","snp"))
all.snps.dup<-all.snps[duplicated(all.snps)]
pca.snps<-all.snps.dup[!duplicated(all.snps.dup)]
#the snp is the row number in the map file for the snps

pca.out<-all.map[pca.snps,]
colnames(pca.out)<-c("chr","loc","dist","bp")
pca.out$radloc<-sub('(\\d+)_\\d+','\\1',pca.out$loc)
pca.out$index<-rownames(pca.out)



#*********************COMPARE ALL OF THE OUTLIER ANALYSES********************#
glob.fst.99$radloc<-factor(glob.fst.99$radloc)
pca.out$radloc<-factor(pca.out$radloc)
#temp.bf.sig$radloc<-factor(sub('(\\d+)_\\d+','\\1',temp.bf.sig$locus))
xtx.5$rad.loc<-factor(xtx.5$rad.loc)
locad.out<-c(as.character(pca.out$loc),as.character(glob.fst.99$Locus),
	as.character(xtx.5$locus))


pca.loc<-pca.out$radloc[!duplicated(pca.out$radloc)]
fst.loc<-glob.fst.99$radloc[!duplicated(glob.fst.99$radloc)]
xtx.loc<-xtx.5$rad.loc[!duplicated(xtx.5$rad.loc)]
###base=fst
fst.pca<-glob.fst.out[glob.fst.out$Locus %in% pca.out$radloc,]
fst.bft<-glob.fst.out[glob.fst.out$Locus %in% temp.bf.sig$rad.loc,]
fst.xtx<-glob.fst.out[glob.fst.out$Locus %in% xtx.5$rad.loc,]

###base=xtx
xtx.pca<-xtx.5[xtx.5$rad.loc %in% pca.out$radloc,]
xtx.bft<-xtx.5[xtx.5$rad.loc %in% temp.bf.sig$rad.loc,]
xtx.fst<-xtx.5[xtx.5$rad.loc %in% glob.fst.out$radloc,]

#compare to environment
dim(temp.bf.sig[temp.bf.sig$locus %in% locad.out,])
dim(ctemp.bf.sig[ctemp.bf.sig$locus %in% locad.out,])
dim(tvar.bf.sig[tvar.bf.sig$locus %in% locad.out,])
dim(sal.bf.sig[sal.bf.sig$locus %in% locad.out,])
dim(csal.bf.sig[csal.bf.sig$locus %in% locad.out,])
dim(grass.bf.sig[grass.bf.sig$locus %in% locad.out,])

temp.bf.sig$analysis<-"MeanTemp"
ctemp.bf.sig$analysis<-"CollectionTemp"
tvar.bf.sig$analysis<-"TempVariance"
sal.bf.sig$analysis<-"Salinity"
csal.bf.sig$analysis<-"CollectionSalinity"
grass.bf.sig$analysis<-"Seagrass"
xtx.5$analysis<-"XtX"
xtx<-xtx.5[,c("locus","scaffold","BP","analysis")]
glob.fst.99$analysis<-"Fst"
fst.99<-glob.fst.99[,c("Locus","Chrom","BP","analysis")]
colnames(fst.99)<-c("locus","scaffold","BP","analysis")
pca.out$analysis<-"PCAdapt"
pca.o<-pca.out[,c("loc","chr","bp","analysis")]
colnames(pca.o)<-c("locus","scaffold","BP","analysis")
outliers<-rbind(temp.bf.sig[,c("locus","scaffold","BP","analysis")],
	ctemp.bf.sig[,c("locus","scaffold","BP","analysis")],
	tvar.bf.sig[,c("locus","scaffold","BP","analysis")],	
	sal.bf.sig[,c("locus","scaffold","BP","analysis")],
	csal.bf.sig[,c("locus","scaffold","BP","analysis")],
	grass.bf.sig[,c("locus","scaffold","BP","analysis")],
	fst.99,xtx,pca.o)
outliers<-outliers[sort(outliers$locus),]
write.table(outliers,"AllOutliers.txt",sep='\t',row.names=F,col.names=T,
	quote=F)

########################PLOT APPENDIX 2: VENN DIAGRAM#########################
out.venn<-venn( list("."=pca.loc, "."=fst.loc,"."=xtx.loc))

jpeg("Fig4_Venn_revisions.jpeg", height=7,width=7,units="in", res=300)
pdf("local_adaptation_venn.pdf",height=7,width=7)
plot(out.venn, small=0.9)
text(x=40,y=120,"PCAdapt")
text(x=200,y=365,expression("X"^T*"X"))
text(x=345,y=120, expression(italic(F)[italic(ST)]))
dev.off()

#which one is shared among all 4?
all.sig<-pca.out[pca.out$radloc %in% fst.loc & pca.out$radloc %in% xtx.loc,]
write.table(all.sig[,c(1,2,4)],"shared_outlier.txt",quote=F,row.names=F)

sig.region<-data.frame(all.sig$chr,all.sig$bp-2500,all.sig$bp+2500)
write.table(sig.region,"shared_outlier_region.sh",quote=F,row.names=F,
	col.names=F,sep='\t',eol='\n')
##########################PLOT FIGURE 3: OUTLIERS#############################
#non-outliers: col = grey53, pch=19
#Fst outliers: col=black
#XtX outliers: col=blue, pch=19 or 1
#pcadapt=dark green, pch=0
neutral.col<-"black"
fst.out.col<-"purple"
xtx.out.col<-"blue"
pca.out.col<-"forestgreen"
neutral.pch<-19
fst.out.pch<-1
xtx.out.pch<-2
pca.out.pch<-0

#############FIG 3
jpeg("LocalAdaptationOutliers_revisions.jpeg", height=9,width=7.5,units="in",res=300)
pdf("LocalAdaptationOutliers_revisions.pdf", height=9,width=7.5)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),cex=2,mgp=c(3,0.5,0))

###GLOBAL FSTS
global.fsts$Chrom<-factor(global.fsts$Chrom)
all.scaff<-split(global.fsts, global.fsts$Chrom)
last.max<-0
rect.xs<-NULL
addition.values<-0
for(i in 1:length(all.scaff)){
	new.max<-last.max+round(max(all.scaff[[i]]$BP), -2)
	rect.xs<-rbind(rect.xs,c(last.max, new.max))
	addition.values<-c(addition.values, new.max)
	last.max<-new.max
}
#change BP to plot
for(i in 1:length(all.scaff)){
	all.scaff[[i]]$BP<-all.scaff[[i]]$BP+addition.values[i]
}
x.max<-max(addition.values)
x.min<-min(all.scaff[[1]]$BP)
y.max<-0.35
y.min<-0
plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray75"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plotting.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$Fst,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=y.min,x.min=x.min, pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$radloc %in% glob.fst.99$radloc,]
	points(temp.sig$BP, temp.sig$Fst, col=fst.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$radloc %in%
		all.sig$radloc,]
	points(temp.sig$BP, temp.sig$Fst, col="red", pch=8,
		cex=0.75)
}
axis(2, at = seq(0,0.4,0.1),ylim = c(y.min, y.max), pos=0,
	labels=seq(0,0.4,0.1),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=2, expression(italic(F)[italic(ST)]), 
	outer=FALSE, line=1.2,cex=1)

#######PCADAPT
#plot a representative run
pca.scaff$chrom<-factor(pca.scaff$chrom)
pca.scaff$logBF<-as.numeric(pca.scaff$logBF)
all.scaff<-split(pca.scaff, factor(pca.scaff$chrom))
last.max<-0
rect.xs<-NULL
addition.values<-0
for(i in 1:length(all.scaff)){
	new.max<-last.max+round(max(all.scaff[[i]]$BP), -2)
	rect.xs<-rbind(rect.xs,c(last.max, new.max))
	addition.values<-c(addition.values, new.max)
	last.max<-new.max
}
#change BP to plot
for(i in 1:length(all.scaff)){
	all.scaff[[i]]$BP<-all.scaff[[i]]$BP+addition.values[i]
}
x.max<-max(addition.values)
x.min<-min(all.scaff[[1]]$BP)
y.max<-30
y.min<--1

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray75"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plotting.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$logBF,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=0,x.min=0, pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$loc %in% pca.out$loc,]
	points(temp.sig$BP, temp.sig$logBF, col=pca.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$radloc %in%
		all.sig$radloc,]
	points(temp.sig$BP, temp.sig$logBF, col="red", pch=8,
		cex=0.75)
}
axis(2, at = c(0,15,30),
	ylim = c(y.min, y.max), pos=0,
	labels = c(0,15,30),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=2, expression(PCAdapt~log[10](BF)), outer=FALSE,line=1.2,cex=1)


#######BAYENV2 POPULATION DIFFERENTIATION
xtx.scaff$scaffold<-factor(xtx.scaff$scaffold)
all.scaff<-split(xtx.scaff, xtx.scaff$scaffold)
last.max<-0
rect.xs<-NULL
addition.values<-0
for(i in 1:length(all.scaff)){
	new.max<-last.max+round(max(all.scaff[[i]]$BP), -2)
	rect.xs<-rbind(rect.xs,c(last.max, new.max))
	addition.values<-c(addition.values, new.max)
	last.max<-new.max
}
#change BP to plot
for(i in 1:length(all.scaff)){
	all.scaff[[i]]$BP<-all.scaff[[i]]$BP+addition.values[i]
}
x.max<-max(addition.values)
x.min<-min(all.scaff[[1]]$BP)
y.max<-60
y.min<-0

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray75"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plotting.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$XtX,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=0,x.min=0, pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in% xtx.5$locus,]
	points(temp.sig$BP, temp.sig$XtX, col=xtx.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$rad.loc %in%
		all.sig$radloc,]
	points(temp.sig$BP, temp.sig$XtX, col="red", pch=8,
		cex=0.75)
}
axis(2, at = seq(0,60,10),
	ylim = c(y.min, y.max), pos=0,
	labels = seq(0,60,10),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=1, "Genomic Location", outer = FALSE, line=-0.5,cex=1)
mtext(side=2,  expression("X"^T*"X"), outer=FALSE,line=1.2,cex=1)

#PLOT THE LEGEND
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", ncol=4,
	col=c(fst.out.col,xtx.out.col,pca.out.col, "red"),
	c("Global Fst Outliers","XtX Outliers","PCAdapt Outliers",
		"Shared"),
	pch=c(19,19,19,8), box.lty=0)

dev.off()



##############################################################################
#****************************************************************************#
####################################PST-FST###################################
#****************************************************************************#
##############################################################################
pairwise.pst<-function(dat, pop.order){
	#first column must be pop id/grouping factor
	library(nlme)
	dat.split<-split(dat, factor(dat[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(pop.order),numeric(0), simplify = F), pop.order))
	for(i in 1:(length(pop.order)-1)){
	  for(j in (i+1):length(pop.order)){
		temp.data<-rbind(as.data.frame(dat.split[[pop.order[i]]]),
			as.data.frame(dat.split[[pop.order[j]]]))
		colnames(temp.data)<-c("PopID","Var")
		temp.data$PopID<-factor(temp.data$PopID)
		anv <- lme(fixed=Var ~ 1, random=~1|PopID,data=temp.data)
		varcomp <- VarCorr(anv)
		v.btwn<- as.numeric(varcomp[1])
		v.wthn <- as.numeric(varcomp[2])
		pst <- v.btwn/(v.btwn+2*v.wthn)
		dat.var[pop.order[i],pop.order[j]]<-pst
		#aov.var<-summary.aov(
		#	aov(temp.data[,2]~temp.data[,1]))[[1]]$`Sum Sq`
		#aov.df<-summary.aov(
		#	aov(temp.data[,2]~temp.data[,1]))[[1]]$`Df`
		#dat.var[pop.order[i],pop.order[j]]<-aov.var[2]/(aov.var[2]+
		#	(2*(aov.var[1]/(aov.df[2]-1))))	
	  }
	}
	dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
	rownames(dat.var)<-colnames(dat.var)
	return(dat.var)
}
all.traits.pst.mantel<-function(trait.df,comp.df,id.index){
	results.mantel<-data.frame()
	for(i in 3:ncol(trait.df)){
		res<-mantel.rtest(
			as.dist(t(pairwise.pst(trait.df[,c(id.index,i)],pop.order))),
			as.dist(t(comp.df)), nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	rownames(results.mantel)<-colnames(trait.df)[3:ncol(trait.df)]
	colnames(results.mantel)<-c("Obs","P")
	return(results.mantel)
}

pop.order<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
	"FLAB","FLPB","FLHB","FLCC")

#**********PST COMPARISONS**********#
fem.pheno.sep<-split(fem.pheno, fem.pheno$PopID)
fem.unstd.new<-rbind(fem.pheno.sep$TXSP,fem.pheno.sep$TXCC,fem.pheno.sep$TXCB,
	fem.pheno.sep$ALST,fem.pheno.sep$FLSG,fem.pheno.sep$FLKB,
	fem.pheno.sep$FLFD,fem.pheno.sep$FLSI,fem.pheno.sep$FLAB,
	fem.pheno.sep$FLPB,fem.pheno.sep$FLHB,fem.pheno.sep$FLCC)
mal.pheno.sep<-split(mal.pheno, mal.pheno$PopID)
mal.unstd.new<-rbind(mal.pheno.sep$TXSP,mal.pheno.sep$TXCC,mal.pheno.sep$TXCB,
	mal.pheno.sep$ALST,mal.pheno.sep$FLSG,mal.pheno.sep$FLKB,
	mal.pheno.sep$FLFD,mal.pheno.sep$FLSI,mal.pheno.sep$FLAB,
	mal.pheno.sep$FLPB,mal.pheno.sep$FLHB,mal.pheno.sep$FLCC)

fem.psts<-apply(fem.unstd.new[,3:10],2,function(x){
	pst<-pairwise.pst(data.frame(fem.unstd.new[,1],x),pop.list)
	return(pst)
})
for(i in 1:length(fem.psts)){
	write.table(fem.psts[[i]],paste(names(fem.psts)[i],".fem.pst.txt",sep=""),
		sep='\t',quote=F)
}
mal.psts<-apply(mal.unstd.new[,3:8],2,function(x){
	pst<-pairwise.pst(data.frame(mal.unstd.new[,1],x),pop.list)
	return(pst)
})
for(i in 1:length(mal.psts)){
	write.table(mal.psts[[i]],paste(names(mal.psts)[i],".mal.pst.txt",sep=""),
		sep='\t',quote=F)
}

fem.fst.upst<-all.traits.pst.mantel(fem.unstd.new,pwise.fst.sub,1)
mal.fst.upst<-all.traits.pst.mantel(mal.unstd.new,pwise.fst.sub,1)

fem.upst.dist<-all.traits.pst.mantel(fem.unstd.new,dist,1)
mal.upst.dist<-all.traits.pst.mantel(mal.unstd.new,dist,1)

fst.pst.byloc<-function(ped.file,trait.df,pop.order,trait.ind){
	results.list<-list()
	for(j in 3:ncol(trait.df)){
	results.mantel<-data.frame()
	for(i in seq(7,ncol(ped.file),2)){
		res<-mantel.rtest(
			as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
			as.dist(t(pairwise.pst(trait.df[,c(trait.ind,j)],pop.order))),
			 nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	colnames(results.mantel)<-c("Obs","P")
	results.list<-append(results.list,data.frame(results.mantel))
	}
	#names(results.list)<-colnames(trait.df)[3:ncol(trait.df)]
	return(results.list)
}

fem.pst.fst.loc<-fst.pst.byloc(sub.ped,fem.unstd.new,pop.list,1)
fpf<-data.frame(SVL.Obs=fem.pst.fst.loc[[1]],SVL.P=fem.pst.fst.loc[[2]],
	TailLength.Obs=fem.pst.fst.loc[[3]],TailLength.P=fem.pst.fst.loc[[4]],
	BodyDepth.Obs=fem.pst.fst.loc[[5]],BodyDepth.P=fem.pst.fst.loc[[6]],
	SnoutLength.Obs=fem.pst.fst.loc[[7]],SnoutLength.P=fem.pst.fst.loc[[8]],
	SnoutDepth.Obs=fem.pst.fst.loc[[9]],SnoutDepth.P=fem.pst.fst.loc[[10]],
	HeadLength.Obs=fem.pst.fst.loc[[11]],HeadLength.P=fem.pst.fst.loc[[12]],
	BandArea.Obs=fem.pst.fst.loc[[13]],BandArea.P=fem.pst.fst.loc[[14]],
	BandNum.Obs=fem.pst.fst.loc[[15]],BandNum.P=fem.pst.fst.loc[[16]])
row.names(fpf)<-sub.map$V2
mal.pst.fst.loc<-fst.pst.byloc(sub.ped,mal.unstd.new,pop.list,1)
mpf<-data.frame(SVL.Obs=mal.pst.fst.loc[[1]],SVL.P=mal.pst.fst.loc[[2]],
	TailLength.Obs=mal.pst.fst.loc[[3]],TailLength.P=mal.pst.fst.loc[[4]],
	BodyDepth.Obs=mal.pst.fst.loc[[5]],BodyDepth.P=mal.pst.fst.loc[[6]],
	SnoutLength.Obs=mal.pst.fst.loc[[7]],SnoutLength.P=mal.pst.fst.loc[[8]],
	SnoutDepth.Obs=mal.pst.fst.loc[[9]],SnoutDepth.P=mal.pst.fst.loc[[10]],
	HeadLength.Obs=mal.pst.fst.loc[[11]],HeadLength.P=mal.pst.fst.loc[[12]])
row.names(mpf)<-sub.map$V2

fpf.sig<-fpf[fpf$SVL.P <= 0.05 | fpf$TailLength.P <= 0.05 | 
	fpf$BodyDepth.P <= 0.05 | fpf$SnoutLength.P <= 0.05 | 
	fpf$SnoutDepth.P <= 0.05 | fpf$HeadLength.P <= 0.05 | 
	fpf$BandArea.P <= 0.05 | fpf$BandNum.P <= 0.05,]
mpf.sig<-mpf[mpf$SVL.P <= 0.05 | mpf$TailLength.P <= 0.05 | 
	mpf$BodyDepth.P <= 0.05 | mpf$SnoutLength.P <= 0.05 | 
	mpf$SnoutDepth.P <= 0.05 | mpf$HeadLength.P <= 0.05,]

create.extract.sh<-function(df){
	#the df needs chrom, start, end as columns
	commands<-paste(
		"../SCA/programs/extract_sequence_part/extract_sequence_part",
		" -f ./sw_results/pstfst/sig_regions/",df[,1],".fasta -s ",
		df[,2], " -e ", df[,3],sep="")
	return(commands)
}
svl.sig<-rownames(fpf)[rownames(fpf[fpf$SVL.P <= 0.05,]) %in%
	rownames(mpf[mpf$SVL.P <= 0.05,])]
write.table(svl.sig,"pstfst/SVL_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",svl.sig),
	"pstfst/SVL_radloc.txt",col.names=F,row.names=F,quote=F)
svl.5kb<-all.map[all.map$V2 %in% svl.sig,c(1,4)]
svl.5kb$start<-svl.5kb$V4-2500
svl.5kb$stop<-svl.5kb$V4+2500
write.table(create.extract.sh(svl.5kb[,-2]),
	"pstfst/SVL_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

tail.sig<-rownames(fpf)[rownames(fpf[fpf$TailLength.P <= 0.05,]) %in%
	rownames(mpf[mpf$TailLength.P <= 0.05,])]
write.table(tail.sig,"pstfst/TailLength_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",tail.sig),
	"pstfst/TailLength_radloc.txt",col.names=F,row.names=F,quote=F)
tail.5kb<-all.map[all.map$V2 %in% tail.sig,c(1,4)]
tail.5kb$start<-tail.5kb$V4-2500
tail.5kb$stop<-tail.5kb$V4+2500
write.table(create.extract.sh(tail.5kb[,-2]),
	"pstfst/TailLength_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

body.sig<-rownames(fpf)[rownames(fpf[fpf$BodyDepth.P <= 0.05,]) %in%
	rownames(mpf[mpf$BodyDepth.P <= 0.05,])]
write.table(body.sig,"pstfst/BodyDepth_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",body.sig),
	"pstfst/TailLength_radloc.txt",col.names=F,row.names=F,quote=F)
body.5kb<-all.map[all.map$V2 %in% body.sig,c(1,4)]
body.5kb$start<-body.5kb$V4-2500
body.5kb$stop<-body.5kb$V4+2500
write.table(create.extract.sh(body.5kb[,-2]),
	"pstfst/BodyDepth_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

sntl.sig<-rownames(fpf)[rownames(fpf[fpf$SnoutLength.P <= 0.05,]) %in%
	rownames(mpf[mpf$SnoutLength.P <= 0.05,])]
write.table(sntl.sig,"pstfst/SnoutLength_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",sntl.sig),
	"pstfst/SnoutLength_radloc.txt",col.names=F,row.names=F,quote=F)
sntl.5kb<-all.map[all.map$V2 %in% sntl.sig,c(1,4)]
sntl.5kb$start<-sntl.5kb$V4-2500
sntl.5kb$stop<-sntl.5kb$V4+2500
write.table(create.extract.sh(sntl.5kb[,-2]),
	"pstfst/SnoutLength_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

head.sig<-rownames(fpf)[rownames(fpf[fpf$HeadLength.P <= 0.05,]) %in%
	rownames(mpf[mpf$HeadLength.P <= 0.05,])]
write.table(head.sig,"pstfst/HeadLength_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+)","\\1",head.sig),
	"pstfst/HeadLength_radloc.txt",col.names=F,row.names=F,quote=F)
head.5kb<-all.map[all.map$V2 %in% head.sig,c(1,4)]
head.5kb$start<-head.5kb$V4-2500
head.5kb$stop<-head.5kb$V4+2500
write.table(create.extract.sh(head.5kb[,-2]),
	"pstfst/HeadLength_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

sntd.sig<-rownames(fpf)[rownames(fpf[fpf$SnoutDepth.P <= 0.05,]) %in%
	rownames(mpf[mpf$SnoutDepth.P <= 0.05,])]
write.table(sntd.sig,"pstfst/SnoutDepth_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",sntd.sig),
	"pstfst/SnoutDepth_radloc.txt",col.names=F,row.names=F,quote=F)
sntd.5kb<-all.map[all.map$V2 %in% sntd.sig,c(1,4)]
sntd.5kb$start<-sntd.5kb$V4-2500
sntd.5kb$stop<-sntd.5kb$V4+2500
write.table(create.extract.sh(sntd.5kb[,-2]),
	"pstfst/SnoutDepth_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

band.sig<-rownames(fpf[fpf$BandArea.P <= 0.05 | 
	fpf$BandNum.P <=0.05,])
write.table(band.sig,"pstfst/Bands_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",band.sig),
	"pstfst/Bands_radloc.txt",col.names=F,row.names=F,quote=F)
band.5kb<-all.map[all.map$V2 %in% band.sig,c(1,4)]
band.5kb$start<-band.5kb$V4-2500
band.5kb$stop<-band.5kb$V4+2500
write.table(create.extract.sh(band.5kb[,-2]),
	"pstfst/Bands_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

sig.all<-rownames(fpf)[rownames(fpf) %in% svl.sig & 
	rownames(fpf) %in% tail.sig & rownames(fpf) %in% body.sig & 
	rownames(fpf) %in% sntd.sig & rownames(fpf) %in% sntl.sig & 
	rownames(fpf) %in% head.sig & rownames(fpf) %in% band.sig] # "7875_9"   "24777_60" "5317_10" 

#WRITE TO FILE, after other files have been written
svl.sig<-read.table("pstfst/SVL_pstfst.txt")
tail.sig<-read.table("pstfst/TailLength_pstfst.txt")
body.sig<-read.table("pstfst/BodyDepth_pstfst.txt")
sntl.sig<-read.table("pstfst/SnoutLength_pstfst.txt")
head.sig<-read.table("pstfst/HeadLength_pstfst.txt")
sntd.sig<-read.table("pstfst/SnoutDepth_pstfst.txt")
band.sig<-read.table("pstfst/Bands_pstfst.txt")

svl.sig<-all.map[all.map$V2 %in% svl.sig$V1,]
colnames(svl.sig)<-c("scaffold","SNP","Dist","BP")
svl.sig$locus<-gsub("(\\d+)_\\d+","\\1",svl.sig$SNP)
tail.sig<-all.map[all.map$V2 %in% tail.sig$V1,]
colnames(tail.sig)<-c("scaffold","SNP","Dist","BP")
tail.sig$locus<-gsub("(\\d+)_\\d+","\\1",tail.sig$SNP)
body.sig<-all.map[all.map$V2 %in% body.sig$V1,]
colnames(body.sig)<-c("scaffold","SNP","Dist","BP")
body.sig$locus<-gsub("(\\d+)_\\d+","\\1",body.sig$SNP)
sntl.sig<-all.map[all.map$V2 %in% sntl.sig$V1,]
colnames(sntl.sig)<-c("scaffold","SNP","Dist","BP")
sntl.sig$locus<-gsub("(\\d+)_\\d+","\\1",sntl.sig$SNP)
head.sig<-all.map[all.map$V2 %in% head.sig$V1,]
colnames(head.sig)<-c("scaffold","SNP","Dist","BP")
head.sig$locus<-gsub("(\\d+)_\\d+","\\1",head.sig$SNP)
sntd.sig<-all.map[all.map$V2 %in% sntd.sig$V1,]
colnames(sntd.sig)<-c("scaffold","SNP","Dist","BP")
sntd.sig$locus<-gsub("(\\d+)_\\d+","\\1",sntd.sig$SNP)
band.sig<-all.map[all.map$V2 %in% band.sig$V1,]
colnames(band.sig)<-c("scaffold","SNP","Dist","BP")
band.sig$locus<-gsub("(\\d+)_\\d+","\\1",band.sig$SNP)

tags<-read.table("stacks/batch_1.catalog.tags.tsv",sep='\t',header=F)
seqs<-tags[,c("V3","V10")]
svl.sig<-merge(svl.sig,seqs, by.x="locus",by.y="V3")
svl.sig$Trait<-"SVL"
tail.sig<-merge(tail.sig,seqs, by.x="locus",by.y="V3")
tail.sig$Trait<-"Tail Length"
body.sig<-merge(body.sig,seqs, by.x="locus",by.y="V3")
body.sig$Trait<-"Body Depth"
sntl.sig<-merge(sntl.sig,seqs, by.x="locus",by.y="V3")
sntl.sig$Trait<-"Snout Length"
head.sig<-merge(head.sig,seqs, by.x="locus",by.y="V3")
head.sig$Trait<-"Head Length"
sntd.sig<-merge(sntd.sig,seqs, by.x="locus",by.y="V3")
sntd.sig$Trait<-"Snout Depth"
band.sig<-merge(band.sig,seqs, by.x="locus",by.y="V3")
band.sig$Trait<-"Bands"
pstfst.sig<-data.frame(rbind(svl.sig,tail.sig,body.sig,sntl.sig,head.sig,
	sntd.sig,band.sig))
write.csv(pstfst.sig,"pstfst/pstfst_loci_summary.csv")


chroms<-levels(factor(c(as.character(svl.5kb$V1),as.character(body.5kb$V1),
	as.character(sntl.5kb$V1),as.character(sntd.5kb$V1),
	as.character(tail.5kb$V1),as.character(head.5kb$V1),
	as.character(band.5kb$V1))))
write.table(chroms,"pstfst/pstfst.sig_scaffolds.txt",col.names=F,row.names=F,
	quote=F)

sig.fst.ibd<-ibd.by.loc[ibd.by.loc$P <= 0.05,]
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(svl.sig)])
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(tail.sig)])
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(body.sig)])
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(sntl.sig)])
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(sntd.sig)])
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(head.sig)])
length(rownames(sig.fst.ibd)[rownames(sig.fst.ibd) %in% rownames(band.sig)])

length(outliers$SNP[outliers$SNP %in% svl.sig$SNP])
length(outliers$SNP[outliers$SNP %in% tail.sig$SNP])
length(outliers$SNP[outliers$SNP %in% body.sig$SNP])
length(outliers$SNP[outliers$SNP %in% sntl.sig$SNP])
length(outliers$SNP[outliers$SNP %in% sntd.sig$SNP])
length(outliers$SNP[outliers$SNP %in% head.sig$SNP])
length(outliers$SNP[outliers$SNP %in% band.sig$SNP])

#merge with outlier and BF analyses
pstfst.sig$analysis<-pstfst.sig$Trait
pst.out.bf<-data.frame(rbind(pstfst.sig[,c("locus","scaffold","BP","analysis")],
	outliers))
#********************************PCA*****************************************#
fem.pheno$PopID<-factor(fem.pheno$PopID)
mal.pheno$PopID<-factor(mal.pheno$PopID)
bands.pcdat<-fem.pheno[,c(1,2,9,10)]
#pca per pop
band.pca<-rda(bands.pcdat[,3:4])
fem.pheno.pca<-rda(fem.pheno[,3:8])
mal.pheno.pca<-rda(mal.pheno[,3:8])


fem.pop<-bands.pcdat$PopID
fem.colors<-as.character(fem.pop)
fem.colors[fem.colors=="TXSP"]<-rainbow(12)[1]
fem.colors[fem.colors=="TXCC"]<-rainbow(12)[2]
fem.colors[fem.colors=="TXCB"]<-rainbow(12)[3]
fem.colors[fem.colors=="ALST"]<-rainbow(12)[4]
fem.colors[fem.colors=="FLSG"]<-rainbow(12)[5]
fem.colors[fem.colors=="FLKB"]<-rainbow(12)[6]
fem.colors[fem.colors=="FLFD"]<-rainbow(12)[7]
fem.colors[fem.colors=="FLSI"]<-rainbow(12)[8]
fem.colors[fem.colors=="FLAB"]<-rainbow(12)[9]
fem.colors[fem.colors=="FLPB"]<-rainbow(12)[10]
fem.colors[fem.colors=="FLHB"]<-rainbow(12)[11]
fem.colors[fem.colors=="FLCC"]<-rainbow(12)[12]

mal.pop<-mal.pheno$PopID
mal.colors<-as.character(mal.pop)
mal.colors[mal.colors=="TXSP"]<-rainbow(12)[1]
mal.colors[mal.colors=="TXCC"]<-rainbow(12)[2]
mal.colors[mal.colors=="TXCB"]<-rainbow(12)[3]
mal.colors[mal.colors=="ALST"]<-rainbow(12)[4]
mal.colors[mal.colors=="FLSG"]<-rainbow(12)[5]
mal.colors[mal.colors=="FLKB"]<-rainbow(12)[6]
mal.colors[mal.colors=="FLFD"]<-rainbow(12)[7]
mal.colors[mal.colors=="FLSI"]<-rainbow(12)[8]
mal.colors[mal.colors=="FLAB"]<-rainbow(12)[9]
mal.colors[mal.colors=="FLPB"]<-rainbow(12)[10]
mal.colors[mal.colors=="FLHB"]<-rainbow(12)[11]
mal.colors[mal.colors=="FLCC"]<-rainbow(12)[12]

#*********************************D-MATRIX***********************************#
y<-pops.pheno[,3:8]
pops.pheno$sex<-gsub("\\w{4}(\\w{1}).*","\\1",pops.pheno$ID)
pops.pheno$sex[pops.pheno$sex=="P" | pops.pheno$sex=="N"]<-"Male"
pops.pheno$sex[pops.pheno$sex=="F"]<-"Female"
pops.pheno$sex[pops.pheno$sex=="D"]<-"Juvenile"

pops.pca<-rda(y)
y.pca<-pops.pca$CA$u

man<-manova(y.pca~pops.pheno$sex*pops.pheno$PopID)
summary.manova(man)
#                                 Df  Pillai approx F num Df den Df    Pr(>F)
#pops.pheno$sex                    2 0.73400   46.479     12    962 < 2.2e-16
#pops.pheno$PopID                 11 1.60513   16.103     66   2910 < 2.2e-16
#pops.pheno$sex:pops.pheno$PopID  15 0.36243    2.079     90   2910 1.988e-08
#Residuals                       485                                         
#                                   
#pops.pheno$sex                  ***
#pops.pheno$PopID                ***
#pops.pheno$sex:pops.pheno$PopID ***
#Residuals      

C<-pops.pca$CA$v
          
D.sex<-solve(summary.manova(man)$SS$Residuals) %*% 
	summary.manova(man)$SS$`pops.pheno$sex`
c.sex<-rda(D.sex)$CA$u
d.sex<-C%*%c.sex
D.pop<-solve(summary.manova(man)$SS$Residuals) %*% 
	summary.manova(man)$SS$`pops.pheno$PopID`
c.pop<-rda(D.pop)$CA$u
d.pop<-C%*%c.pop
D.pox<-solve(summary.manova(man)$SS$Residuals) %*% 
	summary.manova(man)$SS$`pops.pheno$sex:pops.pheno$PopID`
c.pox<-rda(D.pox)$CA$u
d.pox<-C%*%c.pox

###BANDS
y.bands<-pops.pheno[!is.na(pops.pheno$MBandArea),9:10]
bands.popid<-pops.pheno[!is.na(pops.pheno$MBandArea),"PopID"]
bands.sex<-pops.pheno[!is.na(pops.pheno$MBandArea),"sex"]
bands.pca<-rda(y.bands)
y.bands.pca<-bands.pca$CA$u

bands.man<-manova(y.bands.pca~bands.popid*bands.sex) #includes DB & females
summary.manova(bands.man)
#                       Df  Pillai approx F num Df den Df    Pr(>F)    
#bands.popid            11 0.91355  17.8108     22    466 < 2.2e-16 ***
#bands.sex               1 0.18476  26.2900      2    232 5.115e-11 ***
#bands.popid:bands.sex   3 0.03041   1.1992      6    466    0.3054    
#Residuals             233     

C.bands<-bands.pca$CA$v
          
D.bands<-solve(summary.manova(bands.man)$SS$Residuals) %*% 
	summary.manova(bands.man)$SS$bands.popid
c.bands<-rda(as.matrix(D.bands))$CA$u
d.bands<-C.bands%*%c.bands

d.out<-data.frame(dp = c(d.pop[,1],NA,NA),ds=c(d.sex[,1],NA,NA),
	dps=c(d.pox[,1],NA,NA),db=c(rep(NA,6),d.bands[,1]),
	row.names=c(rownames(d.pop),rownames(d.bands)))
write.csv(d.out,"dx.csv")

sem<-function(x){
	sem<-sd(x)/sqrt(length(x))
	return(sem)
}

fem.means<-rbind(tapply(fem.pheno$SVL,fem.pheno$PopID,mean),
	tapply(fem.pheno$TailLength,fem.pheno$PopID,mean),
	tapply(fem.pheno$depth,fem.pheno$PopID,mean),
	tapply(fem.pheno$SnoutLength,fem.pheno$PopID,mean),
	tapply(fem.pheno$SnoutDepth,fem.pheno$PopID,mean),
	tapply(fem.pheno$HeadLength,fem.pheno$PopID,mean))

fem.means<-fem.means[,pop.list]
rownames(fem.means)<-colnames(fem.pheno[,3:8])
fem.means<-data.frame(t(fem.means))
fem.sem<-data.frame(cbind(tapply(fem.pheno$SVL,fem.pheno$PopID,sem),
	tapply(fem.pheno$TailLength,fem.pheno$PopID,sem),
	tapply(fem.pheno$depth,fem.pheno$PopID,sem),
	tapply(fem.pheno$SnoutLength,fem.pheno$PopID,sem),
	tapply(fem.pheno$SnoutDepth,fem.pheno$PopID,sem),
	tapply(fem.pheno$HeadLength,fem.pheno$PopID,sem)))
colnames(fem.sem)<-colnames(fem.pheno[,3:8])


mal.means<-rbind(tapply(mal.pheno$SVL,mal.pheno$PopID,mean),
	tapply(mal.pheno$TailLength,mal.pheno$PopID,mean),
	tapply(mal.pheno$depth,mal.pheno$PopID,mean),
	tapply(mal.pheno$SnoutLength,mal.pheno$PopID,mean),
	tapply(mal.pheno$SnoutDepth,mal.pheno$PopID,mean),
	tapply(mal.pheno$HeadLength,mal.pheno$PopID,mean))
mal.sem<-data.frame(cbind(tapply(mal.pheno$SVL,mal.pheno$PopID,sem),
	tapply(mal.pheno$TailLength,mal.pheno$PopID,sem),
	tapply(mal.pheno$depth,mal.pheno$PopID,sem),
	tapply(mal.pheno$SnoutLength,mal.pheno$PopID,sem),
	tapply(mal.pheno$SnoutDepth,mal.pheno$PopID,sem),
	tapply(mal.pheno$HeadLength,mal.pheno$PopID,sem)))
colnames(mal.sem)<-colnames(mal.pheno[,3:8])
mal.means<-mal.means[,pop.list]
rownames(mal.means)<-colnames(mal.pheno[,3:8])
mal.means<-data.frame(t(mal.means))

band.means<-rbind(tapply(fem.pheno$MBandArea,fem.pheno$PopID,mean),
	tapply(fem.pheno$BandNum,fem.pheno$PopID,mean))
band.sem<-data.frame("BandArea"=tapply(fem.pheno$MBandArea,fem.pheno$PopID,sem),
	"BandNum"=tapply(fem.pheno$BandNum,fem.pheno$PopID,sem))
band.means<-band.means[,pop.list]
rownames(band.means)<-c("BandArea","BandNum")
band.means<-data.frame(t(band.means))


fem.upp<-fem.means+fem.sem
fem.low<-fem.means-fem.sem
mal.upp<-mal.means+mal.sem
mal.low<-mal.means-mal.sem
band.upp<-band.means+band.sem
band.low<-band.means-band.sem
fem.upp<-log(fem.upp)
fem.low<-log(fem.low)
mal.upp<-log(mal.upp)
mal.low<-log(mal.low)
band.upp<-log(band.upp)
band.low<-log(band.low)
mal.means<-log(mal.means)
fem.means<-log(fem.means)
band.means<-log(band.means)

#***********************************PLOT*************************************#
png("PhenotypePCA_Dmatrix.png",height=8,width=10,units="in",res=300)
pdf("PhenotypePCA_Dmatrix.pdf",height=8,width=10)
par(mfrow=c(2,3),oma=c(3,2,2,2),mar=c(3,2,2,2),lwd=1.3)
plot(mal.pheno.pca,type="n",xlim=c(-3,3),ylim=c(-8.2,4)
	,xlab="",ylab="",las=1,cex.axis=1.5)
points(mal.pheno.pca,col=alpha(mal.colors,0.5),cex=1.5,pch=19)
mtext("PC1 (95.75%)",1,line=2)
mtext("PC2 (3.77%)",2,line=2.5)
legend("top",bty='n',c("Male Body Traits"),cex=1.5)

plot(fem.pheno.pca,type="n",xlab="",ylab="",las=1,cex.axis=1.5,ylim=c(-4,12),
	xlim=c(-3,3))
points(fem.pheno.pca,col=alpha(fem.colors,0.5),cex=1.5,pch=19)
mtext("PC1 (91.39%)",1,line=2)
mtext("PC2 (7.54%)",2,line=2.5)
legend("top",bty='n',c("Female Body Traits"),cex=1.5)

plot(band.pca,type="n",xlab="",ylab="",las=1,cex.axis=1.5,xlim=c(-2,2))
points(band.pca,pch=19,col=alpha(fem.colors,0.5),cex=1.5)
mtext("PC1 (98.23%)",1,line=2)
mtext("PC2 (1.77%)",2,line=2.5)
legend("top",bty='n',c("Female Band Traits"),cex=1.5)

par(mar=c(3,2,3,2))
plot(fem.means$TailLength,ylim=c(0,5),pch=15,axes=F,xlab="",ylab="",cex=1.5)
arrows(x0=seq(1,12,1),y0=fem.upp$TailLength,
	x1=seq(1,12,1),y1=fem.low$TailLength,
	angle=90,code=3,length=0.1,cex=1.5)
points(fem.means$SVL,pch=16,cex=1.5)
arrows(x0=seq(1,12,1),y0=fem.upp$SVL,
	x1=seq(1,12,1),y1=fem.low$SVL,
	angle=90,code=3,length=0.1,cex=1.5)
points(fem.means$depth,pch=17,cex=1.5)
arrows(x0=seq(1,12,1),y0=fem.upp$depth,
	x1=seq(1,12,1),y1=fem.low$depth,
	angle=90,code=3,length=0.1,cex=1.5)
points(fem.means$SnoutLength,pch=0,col="dodgerblue",cex=1.5)
arrows(x0=seq(1,12,1),y0=fem.upp$SnoutLength,
	x1=seq(1,12,1),y1=fem.low$SnoutLength,
	angle=90,code=3,length=0.1,col="dodgerblue",cex=1.5)
points(fem.means$SnoutDepth,pch=1,col="dodgerblue",cex=1.5)
arrows(x0=seq(1,12,1),y0=fem.upp$SnoutDepth,
	x1=seq(1,12,1),y1=fem.low$SnoutDepth,
	angle=90,code=3,length=0.1,col="dodgerblue",cex=1.5)
points(fem.means$HeadLength,pch=2,col="dodgerblue",cex=1.5)
arrows(x0=seq(1,12,1),y0=fem.upp$HeadLength,
	x1=seq(1,12,1),y1=fem.low$HeadLength,
	angle=90,code=3,length=0.1,col="dodgerblue",cex=1.5)
axis(2,pos=0.5,las=1,cex.axis=1.5)
mtext("log(mm)",2,outer=F,line=2.2)
axis(1, at=seq(0,13,1), labels=F,pos=0,cex=1.5)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-1.1,
	srt=-90, xpd=TRUE,cex=1.5)

plot(mal.means$TailLength,ylim=c(0,5),pch=15,axes=F,xlab="",ylab="",cex=1.5)
arrows(x0=seq(1,12,1),y0=mal.upp$TailLength,
	x1=seq(1,12,1),y1=mal.low$TailLength,
	angle=90,code=3,length=0.1,cex=1.5)
points(mal.means$SVL,pch=16,cex=1.5)
arrows(x0=seq(1,12,1),y0=mal.upp$SVL,
	x1=seq(1,12,1),y1=mal.low$SVL,
	angle=90,code=3,length=0.1,cex=1.5)
points(mal.means$depth,pch=17,cex=1.5)
arrows(x0=seq(1,12,1),y0=mal.upp$depth,
	x1=seq(1,12,1),y1=mal.low$depth,
	angle=90,code=3,length=0.1,cex=1.5)
points(mal.means$SnoutLength,pch=0,col="dodgerblue",cex=1.5)
arrows(x0=seq(1,12,1),y0=mal.upp$SnoutLength,
	x1=seq(1,12,1),y1=mal.low$SnoutLength,
	angle=90,code=3,length=0.1,col="dodgerblue",cex=1.5)
points(mal.means$SnoutDepth,pch=1,col="dodgerblue",cex=1.5)
arrows(x0=seq(1,12,1),y0=mal.upp$SnoutDepth,
	x1=seq(1,12,1),y1=mal.low$SnoutDepth,
	angle=90,code=3,length=0.1,col="dodgerblue",cex=1.5)
points(mal.means$HeadLength,pch=2,col="dodgerblue",cex=1.5)
arrows(x0=seq(1,12,1),y0=mal.upp$HeadLength,
	x1=seq(1,12,1),y1=mal.low$HeadLength,
	angle=90,code=3,length=0.1,col="dodgerblue",cex=1.5)
axis(2,pos=0.5,las=1,cex.axis=1.5,at=seq(0,5,1))
axis(1, at=seq(0,13,1), labels=F,pos=0,cex=1.5)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-1.1,
	srt=-90, xpd=TRUE,cex=1.5)

plot(band.means$BandNum,ylim=c(0,20),pch=18,axes=F,xlab="",ylab="",cex=1.5)
arrows(x0=seq(1,12,1),y0=band.means$BandNum-band.sem$BandNum,
	x1=seq(1,12,1),y1=band.means$BandNum+band.sem$BandNum,
	angle=90,code=3,length=0.1,cex=1.5)
axis(2,pos=0.5,las=1,cex.axis=1.5)
mtext("Count",2,outer=F,line=2.3)
par(new=TRUE)
plot(band.means$BandArea,pch=5,cex=1.5,ylim=c(0,1.1),axes=F,xlab="",ylab="")
arrows(x0=seq(1,12,1),y0=band.means$BandArea-band.sem$BandArea,
	x1=seq(1,12,1),y1=band.means$BandArea+band.sem$BandArea,
	angle=90,code=3,length=0.1,cex=1.5)
axis(4,las=1,cex.axis=1.5,at=seq(0,1.5,0.5))
mtext(expression(mm^2),4,outer=F,line=2.7)
axis(1, at=seq(0,13,1), labels=F,pos=0,cex=1.5)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.68,
	srt=-90, xpd=TRUE,cex=1.5)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", pop.list, pch=19, pt.cex=1,bty='n',
	col=alpha(rainbow(12), 0.5), ncol=12)
legend(x=-.5,y=0,c("Tail Length","SVL","Body Depth","Snout Length",
	"Snout Depth","Head Length","Band Area", "Band Number"), 
	col=c(rep("black",3),rep("dodgerblue",3),"black","black"),
	bty='n',ncol=4,pch=c(15,16,17,0,1,2,18,5))

dev.off()

#*********************************PCA ON PST*********************************#
#extract eigenvalue
band.eig<-band.pca$CA$eig

#extract PC scores
band.u<-data.frame(bands.pcdat[,1:2],"BandPC1"=band.pca$CA$u[,1],stringsAsFactors=F)
band.u.sep<-split(band.u, band.u[,1])
band.u.new<-rbind(band.u.sep$TXSP,band.u.sep$TXCC,band.u.sep$TXCB,
	band.u.sep$ALST,band.u.sep$FLSG,band.u.sep$FLKB,
	band.u.sep$FLFD,band.u.sep$FLSI,band.u.sep$FLAB,
	band.u.sep$FLPB,band.u.sep$FLHB,band.u.sep$FLCC)

fem.pheno.eig<-fem.pheno.pca$CA$eig

#extract PC scores
fem.pheno.u<-data.frame(fem.pheno[,1:2],
	"FemBodyPC1"=fem.pheno.pca$CA$u[,1],stringsAsFactors=F)
fem.u.sep<-split(fem.pheno.u, fem.pheno.u[,1])
fem.u.new<-rbind(fem.u.sep$TXSP,fem.u.sep$TXCC,fem.u.sep$TXCB,
	fem.u.sep$ALST,fem.u.sep$FLSG,fem.u.sep$FLKB,
	fem.u.sep$FLFD,fem.u.sep$FLSI,fem.u.sep$FLAB,
	fem.u.sep$FLPB,fem.u.sep$FLHB,fem.u.sep$FLCC)

mal.pheno.eig<-mal.pheno.pca$CA$eig

#extract PC scores
mal.u<-data.frame(mal.pheno[,1:2],"MalBodyPC1"=mal.pheno.pca$CA$u[,1],
	stringsAsFactors=F)
mal.u.sep<-split(mal.u, mal.u[,1])
mal.u.new<-rbind(mal.u.sep$TXSP,mal.u.sep$TXCC,mal.u.sep$TXCB,
	mal.u.sep$ALST,mal.u.sep$FLSG,mal.u.sep$FLKB,
	mal.u.sep$FLFD,mal.u.sep$FLSI,mal.u.sep$FLAB,
	mal.u.sep$FLPB,mal.u.sep$FLHB,mal.u.sep$FLCC)

band.pst<-pairwise.pst(band.u.new[,c(1,3)],pop.order)
fem.pst<-pairwise.pst(fem.u.new[,c(1,3)],pop.order)
mal.pst<-pairwise.pst(mal.u.new[,c(1,3)],pop.order)

mantel.rtest(as.dist(t(band.pst)),as.dist(t(dist)), nrepet=9999)
mantel.rtest(as.dist(t(fem.pst)),as.dist(t(dist)), nrepet=9999)
mantel.rtest(as.dist(t(mal.pst)),as.dist(t(dist)), nrepet=9999)


mantel.rtest(as.dist(t(band.pst)),as.dist(t(pwise.fst.sub)), nrepet=9999)
mantel.rtest(as.dist(t(fem.pst)),as.dist(t(pwise.fst.sub)), nrepet=9999)
mantel.rtest(as.dist(t(mal.pst)),as.dist(t(pwise.fst.sub)), nrepet=9999)

env.dat<-read.table("bayenv2//env_data_bayenv_raw.txt")
env.u<-rda(t(env.dat))$CA$u
env.u.new<-env.u[match(pop.order,rownames(env.u)),1]
env.dist<-dist(env.u.new)

###************************************PLOT********************************###
jpeg("Fig5.pst.fst.dist.jpeg",height=7,width=7, units="in", res=300)
pdf("pst.fst.dist.pdf",height=7,width=7)
par(las=1, oma=c(1,1,2.5,1), mar=c(3,3,1,3))
plot(dist[upper.tri(dist)], pwise.fst.sub[upper.tri(pwise.fst.sub)], pch=19,
	ylim=c(0,1),xlab="",ylab="")
points(dist[upper.tri(dist)],band.pst[upper.tri(band.pst)], pch=6,col="darkgreen")
points(dist[upper.tri(dist)],fem.pst[upper.tri(fem.pst)],pch=4,col="red")
points(dist[upper.tri(dist)],mal.pst[upper.tri(mal.pst)],pch=5,col="blue")
#points(dist[upper.tri(dist)],env.dist[lower.tri(env.dist)],pch=15,col="purple")
axis(4)
mtext("Distance (miles)",1, outer=T, line=-0.5)
mtext(expression(Smoothed~Pairwise~italic(F)[ST]),2, las=0, outer=T, line=-0.5)
mtext(expression(Pairwise~italic(P)[ST]),4, outer=T, las=0,line=-0.5)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", ncol=2, col=c("black","blue","red","darkgreen"),pch=c(19,5,4,6),
	c(expression(italic(F)[ST]),expression(Male~PCA~italic(P)[ST]), 
		expression(Female~PCA~italic(P)[ST]), 
		expression(Female~Bands~PCA~italic(P)[ST])),bty="n")

dev.off()


#*****************************STANDARDIZED**********************************#
#************create male and female files for pst analysis******************#
#standardize pop.pheno by population.
std.by.mean<-function(trait){
	trait.trans <- trait/mean(trait)
}
fun.over.dfcol<-function(df,vec,index,fun){
	tout<-tapply(df[,vec],INDEX=index,fun)
	out<-data.frame()
	for(i in 1:length(tout)){
		out<-rbind(out,
			cbind(rep(names(tout)[i],length(tout[[i]])),tout[[i]]))
	}
	return(out)
}
pops.pheno<-mal.pheno
pops.pheno$MBandArea<-NA
pops.pheno$BandNum<-NA
pops.pheno<-rbind(pops.pheno,fem.pheno)
pheno.split<-split(pops.pheno,pops.pheno$PopID)
pheno.split.std<-lapply(pheno.split,function(x){
	std<-data.frame(x[,1:2])
	for(i in 3:ncol(x)){	
		newx<-x[!is.na(x[,i]),] #this removes any empty ones in this trait
		tout<-std.by.mean(newx[,i])
		tout<-data.frame(newx$ID,tout)
		colnames(tout)<-c("ID",colnames(newx)[i])
		std<-merge(std,tout,by="ID",all.x=T) #gives ones that were empty NA
	}
	return(std)
})
pops.pheno.std<-data.frame(do.call("rbind",pheno.split.std))
colnames(pops.pheno.std)<-c("ID","PopID",colnames(pops.pheno)[3:10])

fem.pheno.std<-pops.pheno.std[pops.pheno.std$ID %in% fem.pheno$ID,]
mal.pheno.std<-pops.pheno.std[pops.pheno.std$ID %in% mal.pheno$ID,]
mal.pheno.std<-mal.pheno.std[,-9]
mal.pheno.std<-mal.pheno.std[,-9]

#reorder them to match dist file
fem.pheno.sep<-split(fem.pheno.std, fem.pheno.std$PopID)
fem.pheno.new<-rbind(fem.pheno.sep$TXSP,fem.pheno.sep$TXCC,fem.pheno.sep$TXCB,
	fem.pheno.sep$ALST,fem.pheno.sep$FLSG,fem.pheno.sep$FLKB,
	fem.pheno.sep$FLFD,fem.pheno.sep$FLSI,fem.pheno.sep$FLAB,
	fem.pheno.sep$FLPB,fem.pheno.sep$FLHB,fem.pheno.sep$FLCC)
mal.pheno.sep<-split(mal.pheno.std, mal.pheno.std$PopID)
mal.pheno.new<-rbind(mal.pheno.sep$TXSP,mal.pheno.sep$TXCC,mal.pheno.sep$TXCB,
	mal.pheno.sep$ALST,mal.pheno.sep$FLSG,mal.pheno.sep$FLKB,
	mal.pheno.sep$FLFD,mal.pheno.sep$FLSI,mal.pheno.sep$FLAB,
	mal.pheno.sep$FLPB,mal.pheno.sep$FLHB,mal.pheno.sep$FLCC)

fem.fst.pst<-all.traits.pst.mantel(fem.pheno.new,pwise.fst,2)
mal.fst.pst<-all.traits.pst.mantel(mal.pheno.new,pwise.fst,2)

fem.pst.dist<-all.traits.pst.mantel(fem.pheno.new,dist,2)
mal.pst.dist<-all.traits.pst.mantel(mal.pheno.new,dist,2)


#***************PCA*****************#
fem.pheno.std$PopID<-factor(fem.pheno.std$PopID)
mal.pheno.std$PopID<-factor(mal.pheno.std$PopID)
bands.std.pcdat<-fem.pheno.std[,c(1,2,9,10)]
#pca per pop
band.std.pca<-rda(bands.std.pcdat[,3:4])
fem.std.pca<-rda(fem.pheno.std[,3:8])
mal.std.pca<-rda(mal.pheno.std[,3:8])

#extract eigenvalue
band.eig<-band.std.pca$CA$eig

#extract PC scores
sband.u<-data.frame(bands.std.pcdat[,1:2],"BandPC1"=band.std.pca$CA$u[,1],
	stringsAsFactors=F)
sband.u.sep<-split(sband.u, sband.u[,2])
sband.u.new<-rbind(sband.u.sep$TXSP,sband.u.sep$TXCC,sband.u.sep$TXCB,
	sband.u.sep$ALST,sband.u.sep$FLSG,sband.u.sep$FLKB,
	sband.u.sep$FLFD,sband.u.sep$FLSI,sband.u.sep$FLAB,
	sband.u.sep$FLPB,sband.u.sep$FLHB,sband.u.sep$FLCC)

fem.std.eig<-fem.std.pca$CA$eig

#extract PC scores
fem.std.u<-data.frame(fem.pheno.std[,1:2],
	"FemBodyPC1"=fem.std.pca$CA$u[,1],stringsAsFactors=F)
sfem.u.sep<-split(fem.std.u, fem.std.u[,2])
sfem.u.new<-rbind(sfem.u.sep$TXSP,sfem.u.sep$TXCC,sfem.u.sep$TXCB,
	sfem.u.sep$ALST,sfem.u.sep$FLSG,sfem.u.sep$FLKB,
	sfem.u.sep$FLFD,sfem.u.sep$FLSI,sfem.u.sep$FLAB,
	sfem.u.sep$FLPB,sfem.u.sep$FLHB,sfem.u.sep$FLCC)

mal.std.eig<-mal.std.pca$CA$eig

#extract PC scores
smal.u<-data.frame(mal.pheno.std[,1:2],"MalBodyPC1"=mal.std.pca$CA$u[,1],
	stringsAsFactors=F)
smal.u.sep<-split(smal.u, smal.u[,2])
smal.u.new<-rbind(smal.u.sep$TXSP,smal.u.sep$TXCC,smal.u.sep$TXCB,
	smal.u.sep$ALST,smal.u.sep$FLSG,smal.u.sep$FLKB,
	smal.u.sep$FLFD,smal.u.sep$FLSI,smal.u.sep$FLAB,
	smal.u.sep$FLPB,smal.u.sep$FLHB,smal.u.sep$FLCC)

sband.pst<-pairwise.pst(sband.u.new[,c(2,3)],pop.order)
sfem.pst<-pairwise.pst(sfem.u.new[,c(2,3)],pop.order)
smal.pst<-pairwise.pst(smal.u.new[,c(2,3)],pop.order)

mantel.rtest(as.dist(t(sband.pst)),as.dist(t(dist)), nrepet=9999)
mantel.rtest(as.dist(t(sfem.pst)),as.dist(t(dist)), nrepet=9999)
mantel.rtest(as.dist(t(smal.pst)),as.dist(t(dist)), nrepet=9999)


mantel.rtest(as.dist(t(sband.pst)),as.dist(t(pwise.fst)), nrepet=9999)
mantel.rtest(as.dist(t(sfem.pst)),as.dist(t(pwise.fst)), nrepet=9999)
mantel.rtest(as.dist(t(smal.pst)),as.dist(t(pwise.fst)), nrepet=9999)

env.dat<-read.table("bayenv2//env_data_bayenv_raw.txt")
env.u<-rda(t(env.dat))$CA$u
env.u.new<-env.u[match(pop.order,rownames(env.u)),1]
env.dist<-dist(env.u.new)

###PLOT###
#jpeg("Fig6_standardized.jpeg",height=7,width=7, units="in", res=300)
#par(las=1, oma=c(1,1,2.5,1), mar=c(3,3,1,3))
#plot(dist[upper.tri(dist)], pwise.fst.sub[upper.tri(pwise.fst.sub)], pch=19,#
#	ylim=c(0,1),xlab="",ylab="")
#points(dist[upper.tri(dist)],sband.pst[upper.tri(sband.pst)], pch=6,col="darkgreen")
#points(dist[upper.tri(dist)],sfem.pst[upper.tri(sfem.pst)],pch=4,col="red")
#points(dist[upper.tri(dist)],smal.pst[upper.tri(smal.pst)],pch=5,col="blue")
#points(dist[upper.tri(dist)],env.dist[lower.tri(env.dist)],pch=15,col="purple")
#axis(4)
#mtext("Distance (miles)",1, outer=T, line=-0.5)
#mtext("Smoothed Pairwise Fst",2, las=0, outer=T, line=-0.5)
#mtext("Pairwise Pst",4, outer=T, las=0,line=-0.5)
#par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,#
#	cex=1)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("top", ncol=2, col=c("black","blue","red","darkgreen"),pch=c(19,5,4,6),
#	c("Fst","Male PCA Pst", "Female PCA Pst", "Female Bands PCA Pst"),bty="n")

#dev.off()

##############################################################################
#****************************************************************************#
################################P-MATRIX######################################
#****************************************************************************#
##############################################################################
#BETWEEN POPS
#females
#STANDARDIZED: fem.pheno.new
#UNSTANDARDIZED: fem.unstd.new
pops.fem.std.dat<-split(fem.pheno.new, fem.pheno.new$PopID)
pops.fem.uns.dat<-split(fem.unstd.new, fem.unstd.new$PopID)

pmat.fem.std.pops<-calc.pmat(pops.fem.std.dat,3,10)
pmat.fem.uns.std.pops<-calc.pmat(pops.fem.std.dat,3,8)
pmat.fem.uns.std.pops<-calc.pmat(pops.fem.std.dat,9,10)
pmat.fem.uns.pops<-calc.pmat(pops.fem.uns.dat,3,10)
pmat.fem.uns.body.pops<-calc.pmat(pops.fem.uns.dat,3,8)
pmat.fem.uns.band.pops<-calc.pmat(pops.fem.uns.dat,9,10)

#males
pops.mal.std.dat<-split(mal.pheno.new, mal.pheno.new$PopID)
pops.mal.uns.dat<-split(mal.unstd.new, mal.unstd.new$PopID)

pmat.mal.std.pops<-calc.pmat(pops.mal.std.dat,3,8)
pmat.mal.uns.pops<-calc.pmat(pops.mal.uns.dat,3,8)

#write pmatrices to file
#following format of presentation in Bertram et al 2011
#variance on diagonal, covariance lower, correlations upper, 
#with p1,p2,and p3 to right
setwd("pmatrix")
for(i in 1:length(pmat.fem.uns.pops)){
	dat<-pmat.fem.uns.pops[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.fem.uns.pops[[i]])[
		upper.tri(cor(pmat.fem.uns.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.fem.uns.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(fem.unstd.new)[3:10],
		"p1","p2","p3")
	rownames(dat)<-colnames(fem.unstd.new)[3:10]
	write.table(dat, 
		paste("pop.",names(pmat.fem.uns.pops)[i], ".pmat_unstd.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}
for(i in 1:length(pmat.fem.std.pops)){
	dat<-pmat.fem.std.pops[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.fem.std.pops[[i]])[
			upper.tri(cor(pmat.fem.std.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.fem.std.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(fem.pheno.new)[3:10], "p1","p2","p3")
	rownames(dat)<-colnames(fem.pheno.new)[3:10]
	write.table(dat, 
		paste("pop.",names(pmat.fem.std.pops)[i], ".pmat.std.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}


for(i in 1:length(pmat.mal.uns.pops)){
	dat<-pmat.mal.uns.pops[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.mal.uns.pops[[i]])[
		upper.tri(cor(pmat.mal.uns.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.mal.uns.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(mal.unstd.new)[3:8],
		"p1","p2","p3")
	rownames(dat)<-colnames(mal.unstd.new)[3:8]
	write.table(dat, 
		paste("pop.",names(pmat.mal.uns.pops)[i], 
			".male.pmat_unstd.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}
for(i in 1:length(pmat.mal.std.pops)){
	dat<-pmat.mal.std.pops[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.mal.std.pops[[i]])[
			upper.tri(cor(pmat.mal.std.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.mal.std.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(mal.pheno.new)[3:8], "p1","p2","p3")
	rownames(dat)<-colnames(mal.pheno.new)[3:8]
	write.table(dat, 
		paste("pop.",names(pmat.mal.std.pops)[i],
			".male.pmat.std.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}


################################STANDARDIZED P-MATRIX ANALYSES################
#####################MULTIPLE SUBSPACES###################
h.fem<-calc.h(pmat.fem.std.pops)
write.csv(h.fem,"HFemales_Std.csv")
h.mal<-calc.h(pmat.mal.std.pops)
write.csv(h.mal,"HMales_Std.csv")

h.fem.ang<-data.frame(FemAngle1=do.call("rbind",
	lapply(pmat.fem.std.pops,pop.h.angle,H=h.fem,h.eig=1)),
	FemAngle2=do.call("rbind",
	lapply(pmat.fem.std.pops,pop.h.angle,H=h.fem,h.eig=2)),
	FemAngle3=do.call("rbind",
	lapply(pmat.fem.std.pops,pop.h.angle,H=h.fem,h.eig=3)))
rownames(h.fem.ang)<-names(pmat.fem.std.pops)
h.fem.ang<-h.fem.ang[match(pop.order,rownames(h.fem.ang)),]
write.csv(h.fem.ang,"CommonSubspaceAngles_Fem_Std.csv")

h.mal.ang<-data.frame(MalAngle1=do.call("rbind",
	lapply(pmat.mal.std.pops,pop.h.angle,H=h.mal,h.eig=1)),
	MalAngle2=do.call("rbind",
	lapply(pmat.mal.std.pops,pop.h.angle,H=h.mal,h.eig=2)),
	MalAngle3=do.call("rbind",
	lapply(pmat.mal.std.pops,pop.h.angle,H=h.mal,h.eig=3)))
rownames(h.mal.ang)<-names(pmat.mal.std.pops)
h.mal.ang<-h.mal.ang[match(pop.order,rownames(h.mal.ang)),]
write.csv(h.mal.ang,"CommonSubspaceAngles_Mal_Std.csv")

write.csv(h.fem,"CommonSubspaceFemalesStd.csv")
write.csv(h.mal,"CommonSubspaceMalesStd.csv")

write.csv(rbind(eigen(h.fem)$values,eigen(h.fem)$vectors),
	"HfemalesEigenvectorsStd.csv")
write.csv(rbind(eigen(h.mal)$values,eigen(h.mal)$vectors),
	"HmalesEigenvectorsStd.csv")


#####################TENSORS############################
#Adapted from Aguirre et al. 2013 supplemental material
f.n.traits<-8
m.n.traits<-6
pop.names<-names(pmat.fem.std.pops)
fem.tensor<-covtensor(pmat.fem.std.pops)
mal.tensor<-covtensor(pmat.mal.std.pops)
#eigenvalues for the nonzero eigentensors
fem.tensor$s.alpha[,1:fem.tensor$nonzero]
fem.tensor$tensor.summary[
	1:((ncol(fem.tensor$tensor.summary)-2)*fem.tensor$nonzero),]
#plot coordinates of female p matrix in the space of e1 and e2 
plot(fem.tensor$p.coord[,1], ylim=c(-0.5,0.5), xaxt="n", las=1,
	xlab="Population", ylab="alpha")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.names, par("usr")[1]-230,
	srt=-45, xpd=TRUE)
lines(fem.tensor$p.coord[,1], lty=2)
points(fem.tensor$p.coord[,2], pch=19)
lines(fem.tensor$p.coord[,2], lty=1)
points(fem.tensor$p.coord[,3], pch=15)
lines(fem.tensor$p.coord[,3],lty=4)
legend("bottomright", pch=c(1,19,15),lty=c(2,1,4), c("e1","e2","e3")) 
#trait combinations for e1, e2, and e3
round(fem.tensor$tensor.summary[1:(f.n.traits*3),
	2:dim(fem.tensor$tensor.summary)[2]], 3)

#determine which eigenvector explains the most variation in eigentensors
e1.max<-max.eig.val(fem.tensor$tensor.summary$eT.val[1:f.n.traits])
e2.max<-max.eig.val(fem.tensor$tensor.summary$eT.val[(f.n.traits+1):(2*f.n.traits)])
e3.max<-max.eig.val(fem.tensor$tensor.summary$eT.val[((f.n.traits*2)+1):(3*f.n.traits)])

e1.max<-max.eig.val(mal.tensor$tensor.summary$eT.val[1:m.n.traits])
e2.max<-max.eig.val(mal.tensor$tensor.summary$eT.val[(m.n.traits+1):(2*m.n.traits)])
e3.max<-max.eig.val(mal.tensor$tensor.summary$eT.val[((m.n.traits*2)+1):(3*m.n.traits)])


#project eigentensors on the observed array
f.e1.L <- c(as.numeric(fem.tensor$tensor.summary[1,
	3:dim(fem.tensor$tensor.summary)[2]]))
f.e2.L <- c(as.numeric(fem.tensor$tensor.summary[(f.n.traits+1),
	3:dim(fem.tensor$tensor.summary)[2]]))
f.e3.L <- c(as.numeric(fem.tensor$tensor.summary[((f.n.traits*2)+1),
	3:dim(fem.tensor$tensor.summary)[2]]))

m.e1.L <- c(as.numeric(mal.tensor$tensor.summary[1,
	3:dim(mal.tensor$tensor.summary)[2]]))
m.e2.L <- c(as.numeric(mal.tensor$tensor.summary[(m.n.traits+1),
	3:dim(mal.tensor$tensor.summary)[2]]))
m.e3.L <- c(as.numeric(mal.tensor$tensor.summary[((m.n.traits*2)+1),
	3:dim(mal.tensor$tensor.summary)[2]]))

#variance along e1 and e2
f.e11.proj <- lapply(pmat.fem.std.pops, proj, b = f.e1.L)
f.e21.proj <- lapply(pmat.fem.std.pops, proj, b = f.e2.L)
f.e31.proj <- lapply(pmat.fem.std.pops, proj, b = f.e3.L)

m.e11.proj <- lapply(pmat.mal.std.pops, proj, b = m.e1.L)
m.e21.proj <- lapply(pmat.mal.std.pops, proj, b = m.e2.L)
m.e31.proj <- lapply(pmat.mal.std.pops, proj, b = m.e3.L)

###Write info to files
write.csv(fem.tensor$s.mat,"FemaleSmatrixStd.csv")
write.csv(mal.tensor$s.mat,"MaleSmatrixStd.csv")

write.csv(fem.tensor$tensor.summary,"FemaleTensorSummaryStd.csv")
write.csv(mal.tensor$tensor.summary,"MaleTensorSummaryStd.csv")
#******************************FIGURE 6***********************************#
jpeg("Pmatrix_analyses_std.jpeg", width=14, height=14, units="in", res=300)
pdf("Pmatrix_analyses_std.pdf",width=14,height=14)
#layout(matrix(c(0,1,2,3),2,2,byrow=F))
#layout.show(3)
par(mfrow=c(2,2),lwd=1.3,cex=1.3,oma=c(2,2,2,2))

#common subspace females
par(mar=c(3,5,2,1))
plot(x=seq(1,12,1),h.fem.ang[,1],col="red",
	 xaxt="n", las=1,ylim=c(0,1.3),xlab="Population", ylab="angle")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.71,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),h.fem.ang[,1], lty=2,col="red")
points(x=seq(1,12,1),h.fem.ang[,2], pch=19,col="red")
lines(x=seq(1,12,1),h.fem.ang[,2], lty=1,col="red")
points(x=seq(1,12,1),h.fem.ang[,3], pch=15,col="red")
lines(x=seq(1,12,1),h.fem.ang[,3],lty=4,col="red")
legend("topleft", pch=c(1,19,15), lty=c(2,1,4), c(expression(bolditalic(h)[1]),
	expression(bolditalic(h)[2]),expression(bolditalic(h)[3])),col="red",
	ncol=3)
mtext("Females",3,line=1.5,cex=1.3)
text(x=12,y=1.2, "A", font=2)

#common subspace males
par(mar=c(3,5,2,1))
plot(x=seq(1,12,1),h.mal.ang[,1],col="blue",
	 xaxt="n", las=1,ylim=c(0,1.2),xlab="Population", ylab="angle")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.7,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),h.mal.ang[,1], lty=2,col="blue")
points(x=seq(1,12,1),h.mal.ang[,2], pch=19,col="blue")
lines(x=seq(1,12,1),h.mal.ang[,2], lty=1,col="blue")
points(x=seq(1,12,1),h.mal.ang[,3], pch=15,col="blue")
lines(x=seq(1,12,1),h.mal.ang[,3],lty=4,col="blue")
legend("topleft", pch=c(1,19,15), lty=c(2,1,4), c(expression(bolditalic(h)[1]),
	expression(bolditalic(h)[2]),expression(bolditalic(h)[3])),col="blue",
	ncol=3)
mtext("Males",3,line=1.5,cex=1.3)
text(x=12,y=1.1, "B", font=2)

#plot the variance in each population in the direction of e11,e21, and e31
par(mar=c(3,5,2,1))
f.e11.proj<-f.e11.proj[match(pop.order,names(f.e11.proj))]
f.e21.proj<-f.e21.proj[match(pop.order,names(f.e21.proj))]
f.e31.proj<-f.e31.proj[match(pop.order,names(f.e31.proj))]
plot(x=seq(1,12,1),f.e11.proj, col="red",
	xaxt="n", las=1,ylim=c(-.5,.5),xlab="", ylab="lambda")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-1.2,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),f.e11.proj, lty=2,col="red")
points(x=seq(1,12,1),f.e21.proj, pch=19,col="red")
lines(x=seq(1,12,1),f.e21.proj, lty=1,col="red")
#legend("bottomright", pch=c(1,19), lty=c(2,1),col="red",
#	 c(expression(bolditalic(e)[11]),expression(bolditalic(e)[21])))
points(x=seq(1,12,1),f.e31.proj, pch=15,col="red")
lines(x=seq(1,12,1),f.e31.proj,lty=4,col="red")
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4),
	 c(expression(bolditalic(e)[11]),expression(bolditalic(e)[21]),
	expression(bolditalic(e)[31])),col="red")

text(x=12,y=0.42, "C", font=2)

par(mar=c(3,5,2,1))
m.e11.proj<-m.e11.proj[match(pop.order,names(m.e11.proj))]
m.e21.proj<-m.e21.proj[match(pop.order,names(m.e21.proj))]
m.e31.proj<-m.e31.proj[match(pop.order,names(m.e31.proj))]
plot(x=seq(1,12,1),m.e11.proj,col="blue",
	 xaxt="n", las=1,ylim=c(-0.25,0.25),xlab="Population", ylab="lambda")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-.865,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),m.e11.proj, lty=2,col="blue")
points(x=seq(1,12,1),m.e21.proj, pch=19,col="blue")
lines(x=seq(1,12,1),m.e21.proj, lty=1,col="blue")
#legend("bottomright", pch=c(1,19,15), lty=c(2,1,4),col="blue", 
#	 c(expression(bolditalic(e)[11]),expression(bolditalic(e)[21])))
points(x=seq(1,12,1),m.e31.proj, pch=15,col="blue")
lines(x=seq(1,12,1),m.e31.proj,lty=4,col="blue")
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4), 
	 c(expression(bolditalic(e)[11]),expression(bolditalic(e)[21]),
	expression(bolditalic(e)[31])),col="blue")
#text(x=2,y=100, "Males")
text(x=12,y=0.22, "D", font=2)
mtext("Population",1,outer=T,cex=1.3)
dev.off()


#**************************SUPPLEMENTAL FIGURES*****************************#
#summary plots
jpeg("Subspace_eigenvalues_std.jpeg", width=7, height=7, units="in", res=300)
pdf("Subspace_eigenvalues_std.pdf", width=7, height=7)
par(lwd=1.3,cex=1.3,oma=c(2,2,2,2))
#common subspace
par(mar=c(3,5,2,1))
plot(seq(1,8,1),eigen(h.fem)$values,pch=25,col="red",bg="red",xaxt='n',
	yaxt='n',xlab="",ylab="")
axis(1,at=seq(1.1,8.1,1),labels=c(expression(bolditalic(h)[1]),
	expression(bolditalic(h)[2]),expression(bolditalic(h)[3]),
	expression(bolditalic(h)[4]),expression(bolditalic(h)[5]),
	expression(bolditalic(h)[6]),expression(bolditalic(h)[7]),
	expression(bolditalic(h)[8])))
mtext(expression(Eigenvectors~of~bold(H)),1,line=2,cex=1.3)
axis(2,las=1)
mtext(expression(Eigenvalues~of~bold(H)),2,line=2,cex=1.3)
points(seq(1.2,6.2,1),eigen(h.mal)$values,pch=15,col="blue")
legend("top",c("Female","Male"),pch=c(25,15),
	col=c("red","blue"),pt.bg=c("red","blue"))
dev.off()

jpeg("Eigentensor_alphasStd.jpeg", width=7, height=7, units="in", res=300)
pdf("Eigentensor_alphasStd.pdf", width=7, height=7)
par(lwd=1.3,cex=1.3,oma=c(2,2,2,2))
#plot eigenvalues of non-zero eigentensors for S
par(mar=c(3,5,2,1))
plot(fem.tensor$s.alpha[,1:fem.tensor$nonzero], ylab="alpha",
	xlab="", xaxt="n", las=1,pch=25,col="red",bg="red")#3 until it levels off.
axis(1, at=seq(1.1,fem.tensor$nonzero+0.1,1), 
	labels=c(expression(bolditalic(E)[1]),
	expression(bolditalic(E)[2]),expression(bolditalic(E)[3]),
	expression(bolditalic(E)[4]),expression(bolditalic(E)[5]),
	expression(bolditalic(E)[6]),expression(bolditalic(E)[7]),
	expression(bolditalic(E)[8]),expression(bolditalic(E)[9]),
	expression(bolditalic(E)[10]),expression(bolditalic(E)[11])))
points(seq(1.2,mal.tensor$nonzero+0.2,1),
	mal.tensor$s.alpha[,1:mal.tensor$nonzero],pch=15,col="blue")
legend("top", ,c("Female","Male"),pch=c(25,15),
	col=c("red","blue"),pt.bg=c("red","blue"))
mtext("Eigentensors of fourth-order covariance tensor",1,line=2,cex=1.3)
text(x=11,y=395, "B", font=2)
dev.off()

fc1<-fem.tensor$p.coord[,1][match(pop.order,names(fem.tensor$p.coord[,1]))]
fc2<-fem.tensor$p.coord[,2][match(pop.order,names(fem.tensor$p.coord[,2]))]
fc3<-fem.tensor$p.coord[,3][match(pop.order,names(fem.tensor$p.coord[,3]))]
jpeg("FemaleCoordinatesInEspace.jpeg",height=7,width=21,units="in",res=300)
pdf("FemaleCoordinatesInEspace.pdf",height=7,width=21)
par(mfrow=c(1,3),oma=c(3,2,1,1),mar=c(4,2,2,2),cex=1.3,lwd=1.3)
plot(fc1,pch=6,xaxt='n',ylab="Coordinates",xlab="")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.5405,
	srt=-90, xpd=TRUE)
text(x=12,y=0.27,"E1")
plot(fc2,pch=6,xaxt='n',ylab="Coordinates",xlab="")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-.605,
	srt=-90, xpd=TRUE)
text(x=12,y=0.08,"E2")
plot(fc3,pch=6,xaxt='n',ylab="Coordinates",xlab="")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.6205,
	srt=-90, xpd=TRUE)
text(x=12,y=0.03,"E3")
dev.off()



mc1<-mal.tensor$p.coord[,1][match(pop.order,names(mal.tensor$p.coord[,1]))]
mc2<-mal.tensor$p.coord[,2][match(pop.order,names(mal.tensor$p.coord[,2]))]
mc3<-mal.tensor$p.coord[,3][match(pop.order,names(mal.tensor$p.coord[,3]))]
jpeg("MaleCoordinatesInEspaceStd.jpeg",height=7,width=21,units="in",res=300)
pdf("MaleCoordinatesInEspaceStd.pdf",height=7,width=21)
par(mfrow=c(1,3),oma=c(3,2,1,1),mar=c(4,2,2,2),cex=1.3,lwd=1.3)
plot(mc1,xaxt='n',ylab="Coordinates",xlab="")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.73,
	srt=-90, xpd=TRUE)
text(x=12,y=-0.03,"E1")
plot(mc2,xaxt='n',ylab="Coordinates",xlab="")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.583,
	srt=-90, xpd=TRUE)
text(x=12,y=0.019,"E2")
plot(mc3,xaxt='n',ylab="Coordinates",xlab="")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.order, par("usr")[1]-0.588,
	srt=-90, xpd=TRUE)
text(x=12,y=0,"E3")
dev.off()



