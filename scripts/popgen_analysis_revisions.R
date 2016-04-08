#Author: Sarah P. Flanagan
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

setwd("E:/ubuntushare/popgen/sw_results/")
source("../scripts/plotting_functions.R")

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
sub.map<-read.table("stacks/populations/subset.map")
sub.scaffs<-all.map[all.map$V2 %in% sub.map$V2,]#not in the correct order!

raw.pheno<-read.table("popgen.pheno.txt", sep="\t", header=T)
fem.pheno<-read.table("fem.pheno.txt", sep="\t", header=T)
mal.pheno<-read.table("mal.pheno.txt", sep="\t", header=T)


> fem.pheno$TailLength<-fem.pheno$std.length-fem.pheno$SVL
> fem.pheno$HeadLength<-fem.pheno$HeadLength-fem.pheno$SnoutLength

mal.pheno$TailLength<-mal.pheno$std.length-mal.pheno$SVL
mal.pheno$HeadLength<-mal.pheno$HeadLength-mal.pheno$SnoutLength
#############################################################################
#######################PLOT THE POINTS ON A MAP##############################
#############################################################################
jpeg("mar_sites_map_again.jpg", res=300, height=7,width=14, units="in")
par(oma=c(0,0,0,0),mar=c(0,0,0,0),pin=c(7,7))
map("worldHires", "usa",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray90", mar=c(0,0,0,0),fill=TRUE, res=300,myborder=0)
map("worldHires", "mexico",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray95", fill=TRUE, add=TRUE)
points(mar.coor$lon, mar.coor$lat,  col="black", cex=1, pch=19)
abline(h=c(25,30,35),lty=3)
abline(v=c(-80,-85,-90,-95,-100),lty=3)
text(x=c(-99.5,-99.5,-99.5),y=c(25,30,35),c("25N","30N","35N"))
text(x=c(-80,-85,-90,-95),y=rep(35.3,4),c("80W","85W","90W","95W"))
text(y=26,x=-90,"Gulf of Mexico")
text(y=26,x=-90,"Mexico")
text(x=-88,y=32,"USA")
text(x=-78,y=29.5,"Atlantic Ocean")
abline(h=c(25,30),lty=3)
abline(v=c(-80,-85,-90,-95,-100),lty=3)
text(x=c(-99.5,-99.5,-99.5),y=c(25,30,35),c("25N","30N"))
text(x=c(-80,-85,-90,-95),y=rep(31.7,4),c("80W","85W","90W","95W"))
text(x=-96,y=26,"TXSP",font=2)
text(x=-96.4,y=27,"TXCC",font=2)
text(x=-94,y=29,"TXCB",font=2)
text(x=-88,y=29.7,"ALST",font=2)
text(x=-85.7,y=29.2,"FLSG",font=2)
text(x=-84,y=28.8,"FLKB",font=2)
text(x=-83.6,y=27.6,"FLFD",font=2)
text(x=-82.8,y=26,"FLSI",font=2)
text(x=-79.5,y=25,"FLAB",font=2)
text(x=-79,y=26.2,"FLPB",font=2)
text(x=-79,y=27.2,"FLHB",font=2)
text(x=-79.4,y=28.2,"FLCC",font=2)
dev.off()



###############SUMMARY STATS PER POP######################################


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
dat.plink<-read.PLINK("sw_results/stacks/populations/subsetA.raw",
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

#***************************FASTSTRUCTURE*********************************#
str.in<-read.table("faststructure/subset.structure.recode.str")
inds<-str.in[,1]
inds<-sub('.*_([ATF]\\w+)[_.].*','\\1', inds)
inds<-inds[c(TRUE,FALSE)]
pop.id<-substr(inds,1,4)
pop.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
	"FLAB","FLPB","FLHB","FLCC")

#process a faststructure file
faststr.barplot<-function(meanQ.file, k, plot.order, to.file=TRUE){
	rownames(meanQ.file)<-inds
	meanQ.file<-cbind(meanQ.file,pop.id)
	colnames(meanQ.file)<-c(seq(1,k,1), "pop.id")
	pop.means<-rowsum(meanQ.file[,1:k],
		meanQ.file$pop.id)/summary(meanQ.file$pop.id)
	pop.means.plot<-pop.means[match(plot.order,rownames(pop.means)),]
	filename<-paste("structure.barplot.",k,".jpeg",sep="")
	if(to.file==TRUE){
		jpeg(filename,res=300,height=7,width=7, units="in")
		par(mar=c(5,4,4,1),oma=c(2,2,2,1),xpd=TRUE)
	}
	
	bp<-barplot(as.matrix(t(pop.means.plot)), col=rainbow(k, 0.5), 
		legend=FALSE, axes=FALSE, axisnames=FALSE)
	legend("top", inset=c(0,-0.15), box.lty=0,title="Group",
		legend=seq(1,k,1),
		pch=15, col=rainbow(k,0.5), horiz=TRUE)
	axis(2,las=1, pos=0)
	axis(1,labels=FALSE, at=bp,pos=0)
	text(bp, labels=plot.order,srt=30, xpd=TRUE, pos=1,
		par("usr")[3])
	if(to.file==TRUE){
		dev.off()}
}
#simple
fstr2<-read.table("faststructure/pruned_out_simple.2.meanQ",header=F)
fstr2<-cbind(pop.id,fstr2)
fstr2<-data.frame(fstr2$pop.id,fstr2$V2,fstr2$V1)
fstr3<-read.table("faststructure/pruned_out_simple.3.meanQ",header=F)
fstr3<-cbind(pop.id,fstr3)
fstr3<-data.frame(fstr3$pop.id,fstr3$V3,fstr3$V2,fstr3$V1)
fstr4<-read.table("faststructure/pruned_out_simple.4.meanQ",header=F)
fstr4<-cbind(pop.id,fstr4)
fstr4<-data.frame(fstr4$pop.id,fstr4$V1,fstr4$V2,fstr4$V3,fstr4$V4)
fstr5<-read.table("faststructure/pruned_out_simple.5.meanQ",header=F)
fstr5<-cbind(pop.id,fstr5)
fstr5<-data.frame(fstr5$pop.id,fstr5$V1,fstr5$V3,fstr5$V2,fstr5$V4,fstr5$V5)

#res
fstr2res<-read.table("faststructure/res.k2.output_log.2.meanQ",header=F)
fstr2res<-cbind(pop.id,fstr2)
fstr2<-data.frame(fstr2$pop.id,fstr2$V2,fstr2$V1)
fstr3<-read.table("faststructure/pruned_out_simple.3.meanQ",header=F)
fstr3<-cbind(pop.id,fstr3)
fstr3<-data.frame(fstr3$pop.id,fstr3$V3,fstr3$V2,fstr3$V1)
fstr4<-read.table("faststructure/pruned_out_simple.4.meanQ",header=F)
fstr4<-cbind(pop.id,fstr4)
fstr4<-data.frame(fstr4$pop.id,fstr4$V3,fstr4$V2,fstr4$V1)

fstr5<-read.table("faststructure/pruned_out_simple.5.meanQ",header=F)
fstr5<-cbind(pop.id,fstr5)


#assign each ind to a group
fast.groups.3<-as.data.frame(cbind(inds, apply(stru3, 1, which.max)))
fast.groups.3$inds<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', 
	group.inds.3$inds)

fast.groups.5<-as.data.frame(cbind(inds, apply(stru5, 1, which.max)))
fast.groups.5$inds<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', 
	fast.groups.5$inds)

#make a populations map for stacks
#first need to get names from marine_map
marine.map<-read.table("E://ubuntushare//stacks//marine_map.txt") 
marine.map$V1<-sub('([A-Z]{5,7})([[:digit:]]{1}[[:punct:]])', '\\10\\2', 
	marine.map$V1)

fast.groups.5$inds<-marine.map$V1[match(
	sub(
	'(sample_)([A-Z]{5,7}[0-9]?[0-9]?[0-9]?)(.fq)?(_align)',
	'\\2',marine.map$V1), 
	fast.groups.5$inds)]

write.table(fast.groups.5, "E://ubuntushare//stacks//fstru.groups.popmap.txt", 
	col.names=F, sep="\t", eol="\n", quote=F, row.names=F)

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
setwd("../../")

#####BAYENV: ENV
bf.scaff<-merge(sub.map, bf.dat, by.x="V2", by.y="locus")
colnames(bf.scaff)[1:4]<-c("locus","scaffold","dist","BP")
#focus on Bayes Factors, because of Lotterhos & Whitlock (2015)
bf<-bf.dat[,c(1,2,5,8,11,14)]
bf.co<-apply(bf.dat[,-1],2,quantile,0.95)
temp.bf.sig<-bf[bf$Temp_BF>bf.co["Temp_BF"],c(1,2)]#88
sal.bf.sig<-bf[bf$Sal_BF>bf.co["Salinity_BF"],c(1,3)]#0
ctemp.bf.sig<-bf[bf$ctemp_bf>bf.co["coll.temp_bf"],c(1,4)]#0
csal.bf.sig<-bf[bf$csal_bf>bf.co["coll.sal_bf"],c(1,5)]#0
grass.bf.sig<-bf[bf$seagrass_bf>bf.co["seagrass_bf"],c(1,6)]#0


#####BAYENV: XTX
xtx$rad.loc<-sub('(\\d+)(_\\d+)','\\1',xtx$locus)
xtx.scaff<-merge(all.map, xtx, by.x="V2", by.y="locus")
colnames(xtx.scaff)<-c("locus","scaffold","dist","BP","radloc","XtX", "rad.loc")
xtx.1<-xtx.scaff[xtx.scaff$XtX >= quantile(xtx.scaff$XtX,0.99),] #18
xtx.5<-xtx.scaff[xtx.scaff$XtX >= quantile(xtx.scaff$XtX,0.95),] #88

#****************************GLOBAL FST OUTLIERS*****************************#
#Import global values generated by C++ program calculate_global_fsts
#pruned loci found in all 12 populations, 1753 of them.
global.fsts<-read.table("stacks//populations//subset.globalstats.txt",
	header=T, sep='\t')
#missing useful chrom info. let's get some.
all.map<-read.table("stacks//populations//all.6348.plink.map",header=F)
all.map$radloc<-sub("(\\d+)(_\\d+)","\\1",all.map$V2)
global.fsts<-merge(global.fsts,all.map,by.x="Locus",by.y="V2")
global.fsts$Chrom<-global.fsts$V1
#outliers:
fst.ci.99<-mean(global.fsts$Fst) + 2.57583*sd(global.fsts$Fst)
glob.fst.out<-global.fsts[global.fsts$Fst >= fst.ci.99,]


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
glob.fst.out$Locus<-factor(glob.fst.out$radloc)
pca.out$radloc<-factor(pca.out$radloc)
temp.bf.sig$radloc<-factor(sub('(\\d+)_\\d+','\\1',temp.bf.sig$locus))
xtx.1$rad.loc<-factor(xtx.1$rad.loc)


pca.loc<-pca.out$radloc[!duplicated(pca.out$radloc)]
fst.loc<-glob.fst.out$radloc[!duplicated(glob.fst.out$radloc)]
xtx.loc<-xtx.5$rad.loc[!duplicated(xtx.5$rad.loc)]
bft.loc<-temp.bf.sig$radloc[!duplicated(temp.bf.sig$radloc)]
###base=fst
fst.pca<-np.glob.fst.out[np.glob.fst.out$Locus %in% pca.out$radloc,]
fst.bft<-np.glob.fst.out[np.glob.fst.out$Locus %in% temp.bf.sig$rad.loc,]
fst.xtx<-np.glob.fst.out[np.glob.fst.out$Locus %in% xtx.1$rad.loc,]

###base=xtx
xtx.pca<-xtx.1[xtx.1$rad.loc %in% pca.out$radloc,]
xtx.bft<-xtx.1[xtx.1$rad.loc %in% temp.bf.sig$rad.loc,]
xtx.fst<-xtx.1[xtx.1$rad.loc %in% np.glob.fst.out$Locus,]

###base=BF (focus on temp)
bft.pca<-temp.bf.sig[temp.bf.sig$rad.loc %in% pca.out$radloc,]
bft.fst<-temp.bf.sig[temp.bf.sig$rad.loc %in% np.glob.fst.out$Locus,]
bft.xtx<-temp.bf.sig[temp.bf.sig$rad.loc %in% xtx.1$rad.loc,]


##########################PLOT FIGURE 4: VENN DIAGRAM#########################
out.venn<-venn( list("."=pca.loc, "."=fst.loc,"."=xtx.loc))

jpeg("Fig4_revisions.jpeg", height=7,width=7,units="in", res=300)
plot(out.venn, small=0.9)
text(x=40,y=120,"PCAdapt")
text(x=200,y=370,expression("X"^T*"X"))
text(x=345,y=120, expression(italic(F)[italic(ST)]))
dev.off()

#which two is shared among all 4?
all.sig<-pca.out[pca.out$radloc %in% fst.loc & pca.out$radloc %in% xtx.loc,]


##########################PLOT FIGURE 5: OUTLIERS#############################
#non-outliers: col = grey53, pch=19
#Fst outliers: col=black
#XtX outliers: col=blue, pch=19 or 1
#Temp BF outliers: col=purple, pch=2

#pcadapt=dark green, pch=0
neutral.col<-"grey53"
fst.out.col<-"purple"
xtx.out.col<-"blue"
#bft.out.col<-"black"
pca.out.col<-"forestgreen"
neutral.pch<-19
fst.out.pch<-1
xtx.out.pch<-2
#bft.out.pch<-2
pca.out.pch<-0

#############FIG 5
jpeg("Fig5_revisions.jpeg", height=9,width=7.5,units="in",res=300)
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
y.max<-max(global.fsts$Fst)+0.1*max(global.fsts$Fst)
y.min<-min(global.fsts$Fst)-0.1*min(global.fsts$Fst)
if(min(global.fsts$Fst) < 0) {
	y.min<-min(global.fsts$Fst) - 0.1*min(global.fsts$Fst)
} else {
	y.min<-0
}

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray96"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plotting.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$Fst,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=y.min,x.min=x.min, pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$Fst >= fst.ci.99,]
	points(temp.sig$BP, temp.sig$Fst, col=fst.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$radloc %in%
		all.sig$radloc,]
	points(temp.sig$BP, temp.sig$Fst, col="red", pch=8,
		cex=0.75)
}
axis(2, at = c(y.min,y.max/2,y.max),ylim = c(y.min, y.max), pos=0,
	labels=c(round(y.min,2),round(y.max/2,2),round(y.max,2)),
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
y.max<-max(pca.scaff$logBF)+0.2*max(pca.scaff$logBF)
y.min<-min(pca.scaff$logBF)-0.1*min(pca.scaff$logBF)
if(min(pca.scaff$logBF) < 0) {
	y.min<-min(pca.scaff$logBF) - 0.1*min(pca.scaff$logBF)
} else {
	y.min<-0
}

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray96"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plot.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$logBF,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=0,x.min=0, pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$radloc %in% pca.out$radloc,]
	points(temp.sig$BP, temp.sig$logBF, col=pca.out.col, pch=19, cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$radloc %in%
		all.sig$radloc,]
	points(temp.sig$BP, temp.sig$logBF, col="red", pch=8,
		cex=0.75)
}
axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
	ylim = c(y.min, y.max), pos=0,
	labels = seq(round(y.min,2),round(y.max,2),
		round((y.max-y.min)/2, digits=2)),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=2, "PCAdapt log10(BF)", outer=FALSE,line=1.2,cex=1)


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
y.max<-max(xtx.scaff$XtX)+0.1*max(xtx.scaff$XtX)
y.min<-min(xtx.scaff$XtX)-0.1*min(xtx.scaff$XtX)
if(min(xtx.scaff$XtX) < 0) {
	y.min<-min(xtx.scaff$XtX) - 0.1*min(xtx.scaff$XtX)
} else {
	y.min<-0
}

plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray96"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plot.genome.wide(all.scaff[[i]]$BP, 
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
axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
	ylim = c(y.min, y.max), pos=0,
	labels = seq(round(y.min,2),round(y.max,2),
		round((y.max-y.min)/2, digits=2)),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=1, "Genomic Location", outer = FALSE, line=-0.5,cex=1)
mtext(side=2,  expression("X"^T*"X"), outer=FALSE,line=1.2,cex=1)

#PLOT THE LEGEND
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", ncol=3,
	col=c(neutral.col,fst.out.col,xtx.out.col,pca.out.col, "red"),
	c("Neutral","Global Fst Outliers","XtX Outliers","PCAdapt Outliers",
		"Shared"),
	pch=c(19,19,19,19,8), box.lty=0)

dev.off()



#########BAYENV2 TEMPERATURE
bf.plot<-merge(bf.scaff,all.map,by.x="locus",by.y="V2")
bf.plot$scaffold<-factor(bf.plot$V1)
bf.plot$Temp_BF<-log10(bf.plot$Temp_BF)
all.scaff<-split(bf.plot, bf.plot$scaffold)
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
y.max<-max(bf.plot$Temp_BF)+0.1*max(bf.plot$Temp_BF)
y.min<-min(bf.plot$Temp_BF)-0.1*min(bf.plot$Temp_BF)
if(min(bf.plot$Temp_BF) < 0) {
	y.min<-min(bf.plot$Temp_BF) - 0.1*min(bf.plot$Temp_BF)
} else {
	y.min<-0
}

jpeg("bf.temp.assoc.jpeg", width=8, height=4.5, units="in", res=300)
par(mar=c(0,1,1,0), oma=c(1,1,0,0),cex=2)
plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), ylim=c(y.min, y.max), 
	bty="n",type="n",	axes=F, xlab="", ylab="")
for(i in 1:nrow(rect.xs)){
	if(i%%2 == 0) {
		rect.color<-"white"
	} else {
		rect.color<-"gray96"
	}
	rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
		col=rect.color, border=NA)
}
for(i in 1:length(all.scaff)){
	plot.genome.wide(all.scaff[[i]]$BP, 
		all.scaff[[i]]$Temp_BF,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=y.min,x.min=x.min, pt.col=neutral.col,
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in% temp.bf.sig$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col=bft.out.col, pch=19, cex=0.5)
	#also allsig
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$rad.loc %in%
		all.sig$radloc,]
	points(temp.sig$BP, temp.sig$Temp_BF, col="red", pch=8,
		cex=0.75)
	
	
}
axis(2, at = seq(round(y.min,2),round(y.max+0.02,2),
		round((y.max-y.min)/2, digits=8)),
	ylim = c(y.min, y.max), pos=0,
	labels = seq(round(y.min,2),round(y.max+0.02,2),
		round((y.max-y.min)/2, digits=2)),
	las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
mtext(side=1, "Genomic Location", outer = FALSE, cex=1, line=-0.5)
mtext(side=2, "Temp log10(BF)", line=1.2,outer=FALSE,cex=1, las=0)
dev.off()


##############################################################################
#****************************************************************************#
####################################PST-FST###################################
#****************************************************************************#
##############################################################################
pairwise.pst<-function(dat, pop.order){
	#first column must be pop id/grouping factor
	dat.split<-split(dat, factor(dat[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(pop.order),numeric(0), simplify = F), pop.order))
	for(i in 1:(length(pop.order)-1)){
	  for(j in (i+1):length(pop.order)){
		temp.data<-rbind(as.data.frame(dat.split[[pop.order[i]]]),
			as.data.frame(dat.split[[pop.order[j]]]))
		aov.var<-summary.aov(
			aov(temp.data[,2]~temp.data[,1]))[[1]]$`Sum Sq`
		aov.df<-summary.aov(
			aov(temp.data[,2]~temp.data[,1]))[[1]]$`Df`
		dat.var[pop.order[i],pop.order[j]]<-aov.var[2]/(aov.var[2]+
			(2*(aov.var[1]/(aov.df[2]-1))))	
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

svl.sig<-rownames(fpf)[rownames(fpf[fpf$SVL.P <= 0.05,]) %in%
	rownames(mpf[mpf$SVL.P <= 0.05,])]
write.table(svl.sig,"pstfst/SVL_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",svl.sig),
	"pstfst/SVL_radloc.txt",col.names=F,row.names=F,quote=F)
svl.5kb<-sub.scaffs[sub.scaffs$V2 %in% svl.sig,c(1,4)]
svl.5kb$start<-svl.5kb$V4-2500
svl.5kb$stop<-svl.5kb$V4+2500
write.table(svl.5kb[,-2],"pstfst/SVL_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

tail.sig<-rownames(fpf)[rownames(fpf[fpf$TailLength.P <= 0.05,]) %in%
	rownames(mpf[mpf$TailLength.P <= 0.05,])]
write.table(tail.sig,"pstfst/TailLength_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",tail.sig),
	"pstfst/TailLength_radloc.txt",col.names=F,row.names=F,quote=F)
tail.5kb<-sub.scaffs[sub.scaffs$V2 %in% tail.sig,c(1,4)]
tail.5kb$start<-tail.5kb$V4-2500
tail.5kb$stop<-tail.5kb$V4+2500
write.table(tail.5kb[,-2],"pstfst/TailLength_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

body.sig<-rownames(fpf)[rownames(fpf[fpf$BodyDepth.P <= 0.05,]) %in%
	rownames(mpf[mpf$BodyDepth.P <= 0.05,])]
write.table(body.sig,"pstfst/BodyDepth_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",body.sig),
	"pstfst/TailLength_radloc.txt",col.names=F,row.names=F,quote=F)
body.5kb<-sub.scaffs[sub.scaffs$V2 %in% body.sig,c(1,4)]
body.5kb$start<-body.5kb$V4-2500
body.5kb$stop<-body.5kb$V4+2500
write.table(body.5kb[,-2],"pstfst/BodyDepth_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

sntl.sig<-rownames(fpf)[rownames(fpf[fpf$SnoutLength.P <= 0.05,]) %in%
	rownames(mpf[mpf$SnoutLength.P <= 0.05,])]
write.table(sntl.sig,"pstfst/SnoutLength_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",sntl.sig),
	"pstfst/SnoutLength_radloc.txt",col.names=F,row.names=F,quote=F)
sntl.5kb<-sub.scaffs[sub.scaffs$V2 %in% sntl.sig,c(1,4)]
sntl.5kb$start<-sntl.5kb$V4-2500
sntl.5kb$stop<-sntl.5kb$V4+2500
write.table(sntl.5kb[,-2],"pstfst/SnoutLength_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

head.sig<-rownames(fpf)[rownames(fpf[fpf$HeadLength.P <= 0.05,]) %in%
	rownames(mpf[mpf$HeadLength.P <= 0.05,])]
write.table(head.sig,"pstfst/HeadLength_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+)","\\1",head.sig),
	"pstfst/HeadLength_radloc.txt",col.names=F,row.names=F,quote=F)
head.5kb<-sub.scaffs[sub.scaffs$V2 %in% head.sig,c(1,4)]
head.5kb$start<-head.5kb$V4-2500
head.5kb$stop<-head.5kb$V4+2500
write.table(head.5kb[,-2],"pstfst/HeadLength_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

sntd.sig<-rownames(fpf)[rownames(fpf[fpf$SnoutDepth.P <= 0.05,]) %in%
	rownames(mpf[mpf$SnoutDepth.P <= 0.05,])]
write.table(sntd.sig,"pstfst/SnoutDepth_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",sntd.sig),
	"pstfst/SnoutDepth_radloc.txt",col.names=F,row.names=F,quote=F)
sntd.5kb<-sub.scaffs[sub.scaffs$V2 %in% sntd.sig,c(1,4)]
sntd.5kb$start<-sntd.5kb$V4-2500
sntd.5kb$stop<-sntd.5kb$V4+2500
write.table(sntd.5kb[,-2],"pstfst/SnoutDepth_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

band.sig<-rownames(fpf[fpf$BandArea.P <= 0.05 | 
	fpf$BandNum.P <=0.05,])
write.table(band.sig,"pstfst/Bands_pstfst.txt",col.names=F,row.names=F,quote=F)
write.table(gsub("(\\d+)_\\d+","\\1",band.sig),
	"pstfst/Bands_radloc.txt",col.names=F,row.names=F,quote=F)
band.5kb<-sub.scaffs[sub.scaffs$V2 %in% band.sig,c(1,4)]
band.5kb$start<-band.5kb$V4-2500
band.5kb$stop<-band.5kb$V4+2500
write.table(band.5kb[,-2],"pstfst/Bands_extract.sh",col.names=F, row.names=F,
	quote=F,eol='\n')

sig.all<-rownames(fpf)[rownames(fpf) %in% svl.sig & 
	rownames(fpf) %in% tail.sig & rownames(fpf) %in% body.sig & 
	rownames(fpf) %in% sntd.sig & rownames(fpf) %in% sntd.sig & 
	rownames(fpf) %in% head.sig & rownames(fpf) %in% band.sig] #17514_61

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

#***************PST PCA*****************#
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

png("PhenotypePCA.png",height=4,width=10,units="in",res=300)
par(mfrow=c(1,3),oma=c(2,2,2,2),mar=c(2,2,2,2),lwd=1.3)
plot(mal.pheno.pca,type="n",xlim=c(-3,3),ylim=c(-8.2,4)
	,xlab="",ylab="",las=1,cex.axis=1.5)
points(mal.pheno.pca,col=alpha(mal.colors,0.5),cex=1.5,pch=19)
mtext("PC1 (95.75%)",1,line=2)
mtext("PC2 (3.77%)",2,line=2)
legend("top",bty='n',c("Male Body Traits"),cex=1.5)

plot(fem.pheno.pca,type="n",xlab="",ylab="",las=1,cex.axis=1.5,ylim=c(-4,12),
	xlim=c(-3,3))
points(fem.pheno.pca,col=alpha(fem.colors,0.5),cex=1.5,pch=19)
mtext("PC1 (91.39%)",1,line=2)
mtext("PC2 (7.54%)",2,line=2)
legend("top",bty='n',c("Female Body Traits"),cex=1.5)

plot(band.pca,type="n",xlab="",ylab="",las=1,cex.axis=1.5,xlim=c(-2,2))
points(band.pca,pch=19,col=alpha(fem.colors,0.5),cex=1.5)
mtext("PC1 (98.23%)",1,line=2)
mtext("PC2 (1.77%)",2,line=2)
legend("top",bty='n',c("Female Band Traits"),cex=1.5)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", pop.list, pch=19, pt.cex=1,bty='n',
	col=alpha(rainbow(12), 0.5), ncol=12)
dev.off()

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


mantel.rtest(as.dist(t(band.pst)),as.dist(t(pwise.fst)), nrepet=9999)
mantel.rtest(as.dist(t(fem.pst)),as.dist(t(pwise.fst)), nrepet=9999)
mantel.rtest(as.dist(t(mal.pst)),as.dist(t(pwise.fst)), nrepet=9999)

env.dat<-read.table("bayenv2//env_data_bayenv_raw.txt")
env.u<-rda(t(env.dat))$CA$u
env.u.new<-env.u[match(pop.order,rownames(env.u)),1]
env.dist<-dist(env.u.new)

###PLOT###
jpeg("Fig6.pst.fst.dist.jpeg",height=7,width=7, units="in", res=300)
par(las=1, oma=c(1,1,2.5,1), mar=c(3,3,1,3))
plot(dist[upper.tri(dist)], pwise.fst[upper.tri(pwise.fst)], pch=19,
	ylim=c(0,1),xlab="",ylab="")
points(dist[upper.tri(dist)],band.pst[upper.tri(band.pst)], pch=6,col="darkgreen")
points(dist[upper.tri(dist)],fem.pst[upper.tri(fem.pst)],pch=4,col="red")
points(dist[upper.tri(dist)],mal.pst[upper.tri(mal.pst)],pch=5,col="blue")
#points(dist[upper.tri(dist)],env.dist[lower.tri(env.dist)],pch=15,col="purple")
axis(4)
mtext("Distance (miles)",1, outer=T, line=-0.5)
mtext("Smoothed Pairwise Fst",2, las=0, outer=T, line=-0.5)
mtext("Pairwise Pst",4, outer=T, las=0,line=-0.5)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", ncol=2, col=c("black","blue","red","darkgreen"),pch=c(19,5,4,6),
	c("Fst","Male PCA Pst", "Female PCA Pst", "Female Bands PCA Pst"),bty="n")

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
jpeg("Fig6_standardized.jpeg",height=7,width=7, units="in", res=300)
par(las=1, oma=c(1,1,2.5,1), mar=c(3,3,1,3))
plot(dist[upper.tri(dist)], pwise.fst[upper.tri(pwise.fst)], pch=19,
	ylim=c(0,1),xlab="",ylab="")
points(dist[upper.tri(dist)],sband.pst[upper.tri(sband.pst)], pch=6,col="darkgreen")
points(dist[upper.tri(dist)],sfem.pst[upper.tri(sfem.pst)],pch=4,col="red")
points(dist[upper.tri(dist)],smal.pst[upper.tri(smal.pst)],pch=5,col="blue")
#points(dist[upper.tri(dist)],env.dist[lower.tri(env.dist)],pch=15,col="purple")
axis(4)
mtext("Distance (miles)",1, outer=T, line=-0.5)
mtext("Smoothed Pairwise Fst",2, las=0, outer=T, line=-0.5)
mtext("Pairwise Pst",4, outer=T, las=0,line=-0.5)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", ncol=2, col=c("black","blue","red","darkgreen"),pch=c(19,5,4,6),
	c("Fst","Male PCA Pst", "Female PCA Pst", "Female Bands PCA Pst"),bty="n")

dev.off()

##############################################################################
#****************************************************************************#
################################P-MATRIX######################################
#****************************************************************************#
##############################################################################

#GLOBAL
pmat.fem.global<-cov(fem.unstd.new[,3:10])

#pca of traits
fem.phen.uns.rda<-rda(fem.unstd.new[,3:10])
fem.phen.std.rda<-rda(fem.pheno.new[,3:10])
#pca of p-matrices
fem.pmat.uns.rda<-rda(data.frame(cov(fem.unstd.new[,3:10])))
fem.pmat.std.rda<-rda(data.frame(pmat.fem.global))

fem.phen.uns.rda$CA$v #shows loadings per trait

pmat.mal.global<-cov(mal.unstd.new[,3:8])
pmat.mal.uns.rda<-rda(data.frame(pmat.mal.global))

#BETWEEN POPS
calc.pmat<-function(phen.dat.list, dim1, dim2){
	pmat<-lapply(phen.dat.list, function(x){
		cov(x[,dim1:dim2])
	})
	return(pmat)
}

#females
#STANDARDIZED: fem.pheno.new
#UNSTANDARDIZED: fem.unstd.new
pops.fem.std.dat<-split(fem.pheno.new, fem.pheno.new$PopID)
pops.fem.uns.dat<-split(fem.unstd.new, fem.unstd.new$PopID)

pmat.fem.std.pops<-calc.pmat(pops.fem.std.dat,3,10)
pmat.fem.uns.pops<-calc.pmat(pops.fem.uns.dat,3,10)
pmat.fem.uns.pops<-calc.pmat(pops.fem.uns.dat,3,8)

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
	colnames(dat)<-c(colnames(fem.pheno.new)[3:10],
		"p1","p2","p3")
	rownames(dat)<-colnames(fem.pheno.new)[3:10]
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
	colnames(dat)<-c(colnames(mal.pheno.new)[3:8],
		"p1","p2","p3")
	rownames(dat)<-colnames(mal.pheno.new)[3:8]
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



###############################BLOWS METHOD#################################
#A=first three eigenvectors of P for first species
#B=first three eigenvectors of P for second species
#S=t(A)Bt(B)A
#sim=sum(eigen(S))/3

#####FUNCTIONS#####
calc.sim<-function(P.one, P.two){
	A<-eigen(P.one)$vectors[,1:3]
	B<-eigen(P.two)$vectors[,1:3]
	S<-t(A)%*%B%*%t(B)%*%A
	sim<-sum(eigen(S)$values)/3
	return(sim)
}

generate.sim.mat<-function(list.pmatrices){
	sim.mat<-matrix(nrow=length(list.pmatrices), ncol=length(list.pmatrices))
	for(i in 1:length(list.pmatrices)){
		for(j in 1:length(list.pmatrices)){
			sim.mat[i,j]<-
				calc.sim(list.pmatrices[[i]], list.pmatrices[[j]])
		}
	}
	colnames(sim.mat)<-names(list.pmatrices)
	rownames(sim.mat)<-names(list.pmatrices)
	return(sim.mat)
}

#vector correlation of pmax
find.pmax<-function(pmatrix){
	pmax<-eigen(pmatrix)$vectors[,which.max(eigen(pmatrix)$values)]
	pmax<-pmax/sqrt(sum(pmax*pmax))
	return(as.vector(pmax))
}
vector.correlations<-function(list.pmatrices){
	pmax.mat<-matrix(nrow=length(list.pmatrices), ncol=length(list.pmatrices))
	for(i in 1:(length(list.pmatrices)-1)){
		for(j in (i+1):length(list.pmatrices)){
			p1<-find.pmax(list.pmatrices[[i]])
			p2<-find.pmax(list.pmatrices[[j]])
			pmax.mat[i,j]<-p1%*%p2
		}
	}
	colnames(pmax.mat)<-names(list.pmatrices)
	rownames(pmax.mat)<-names(list.pmatrices)
	return(pmax.mat)

}

calc.h<-function(P.list){
	H<-0
	for(i in 1:length(P.list)){
		A<-eigen(P.list[[i]])$vectors[,1:3]
		A.out<-A%*%t(A)
		H<-H+A.out
	}
	colnames(H)<-colnames(P.list[[1]])
	rownames(H)<-colnames(P.list[[1]])
	return(H)
}

pop.h.angle<-function(Pmatrix,H,h.eig){
	b<-eigen(H)$vectors[,h.eig]
	A<-eigen(Pmatrix)$vectors[,1:3]
	trans<-sqrt(t(b)%*%A%*%t(A)%*%b)
	delta<-acos(trans^0.5)
}

#####ANALYSIS#####
#unstandardized
#females
pmat.f.u.p<-pmat.fem.uns.pops[match(pop.list, names(pmat.fem.uns.pops))]
fem.sim.unstd.mat<-generate.sim.mat(pmat.f.u.p)
fem.sim.unstd.fst<-mantel.rtest(as.dist(fem.sim.unstd.mat), 
	as.dist(t(pwise.fst.sub)),	nrepet=9999)
pmax.unstd.mat<-vector.correlations(pmat.fem.uns.pops)
blows.sim.unstd.mat<-fem.sim.unstd.mat
blows.sim.unstd.mat[upper.tri(blows.sim.unstd.mat)]<-
	pmax.unstd.mat[upper.tri(pmax.unstd.mat)]
write.table(blows.sim.unstd.mat, "unstd_p_sim_matrix.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)
#males
pmat.m.u.p<-pmat.mal.uns.pops[match(pop.list, names(pmat.mal.uns.pops))]
mal.sim.unstd.mat<-generate.sim.mat(pmat.m.u.p)
mal.sim.unstd.fst<-mantel.rtest(as.dist(mal.sim.unstd.mat), 
	as.dist(t(pwise.fst.sub)),	nrepet=9999)
pmax.mal.unstd.mat<-vector.correlations(pmat.mal.uns.pops)
blows.mal.sim.unstd.mat<-mal.sim.unstd.mat
blows.mal.sim.unstd.mat[upper.tri(blows.mal.sim.unstd.mat)]<-
	pmax.mal.unstd.mat[upper.tri(pmax.mal.unstd.mat)]
write.table(blows.mal.sim.unstd.mat, "unstd_p_sim_matrix_male.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)
#standardized
#females
pmat.f.s.p<-pmat.fem.std.pops[match(pop.list, names(pmat.fem.std.pops))]
fem.sim.std.mat<-generate.sim.mat(pmat.f.s.p)
fem.sim.std.fst<-mantel.rtest(as.dist(fem.sim.std.mat), as.dist(t(pwise.fst)),
	nrepet=9999)
pmax.fem.std.mat<-vector.correlations(pmat.fem.std.pops)
blows.sim.std.mat<-fem.sim.std.mat
blows.sim.std.mat[upper.tri(blows.sim.std.mat)]<-
	pmax.fem.std.mat[upper.tri(pmax.fem.std.mat)]
write.table(blows.sim.std.mat, "standardized_p_sim_matrix_female.txt", 
	col.names=T,row.names=T, eol='\n', sep='\t', quote=F)
#males
pmat.m.s.p<-pmat.mal.std.pops[match(pop.list, names(pmat.mal.std.pops))]
mal.sim.std.mat<-generate.sim.mat(pmat.m.s.p)
mal.sim.std.fst<-mantel.rtest(as.dist(mal.sim.std.mat), as.dist(t(pwise.fst)),
	nrepet=9999)
pmax.mal.std.mat<-vector.correlations(pmat.mal.std.pops)
blows.mal.sim.std.mat<-mal.sim.std.mat
blows.mal.sim.std.mat[upper.tri(blows.mal.sim.std.mat)]<-
	pmax.mal.std.mat[upper.tri(pmax.std.mat)]
write.table(blows.mal.sim.std.mat, "standardized_p_sim_matrix_female.txt", 
	col.names=T,row.names=T, eol='\n', sep='\t', quote=F)


#######PLOT##########
#make a heatmap, males above the diagonal, females below
#unstandardized
plot.pmat<-mal.sim.unstd.mat
plot.pmat[lower.tri(plot.pmat)]<-
	fem.sim.unstd.mat[lower.tri(fem.sim.unstd.mat)]
jpeg("blows.pmatrix.heatmap.jpeg",height=7,width=7,units="in",res=300)
par(oma=c(3,1,2,3),mar=c(3,1,2,3),las=1)
heatmap.2(plot.pmat, Rowv=NA, Colv=NA, dendrogram="none",col=bluered(100),
	trace="none", margins,key.title=NA, key.ylab=NA,key.xlab=NA,srtCol=25, 
	lmat=rbind(c(2,1),c(3,4)), lwid=c(1,4),lhei=c(4,0.5),
	keysize=0.5, margins=c(2,2),key.par=list(yaxt="n"))
mtext("Females", 2,outer=T, line=-4)
mtext("Males", 3, outer=T,line=1)
dev.off()

#standardized
plot.pmat<-mal.sim.std.mat
plot.pmat[lower.tri(plot.pmat)]<-
	fem.sim.std.mat[lower.tri(fem.sim.std.mat)]
jpeg("blows.pmatrix.std.heatmap.jpeg",height=7,width=7,units="in",res=300)
par(oma=c(3,1,2,3),mar=c(3,1,2,3),las=1)
heatmap.2(plot.pmat, Rowv=NA, Colv=NA, dendrogram="none",col=bluered(100),
	trace="none", margins,key.title=NA, key.ylab=NA,key.xlab=NA,srtCol=25, 
	lmat=rbind(c(2,1),c(3,4)), lwid=c(1,4),lhei=c(4,0.5),
	keysize=0.5, margins=c(2,2),key.par=list(yaxt="n"))
mtext("Females", 2,outer=T, line=-4)
mtext("Males", 3, outer=T,line=1)
dev.off()

#####################MULTIPLE SUBSPACES###################
h.fem<-calc.h(pmat.fem.uns.pops)
h.mal<-calc.h(pmat.mal.uns.pops)

h.fem.ang<-data.frame(FemAngle1=do.call("rbind",
	lapply(pmat.fem.uns.pops,pop.h.angle,H=h.fem,h.eig=1)),
	FemAngle2=do.call("rbind",
	lapply(pmat.fem.uns.pops,pop.h.angle,H=h.fem,h.eig=2)),
	FemAngle3=do.call("rbind",
	lapply(pmat.fem.uns.pops,pop.h.angle,H=h.fem,h.eig=3)))

h.mal.ang<-data.frame(MalAngle1=do.call("rbind",
	lapply(pmat.mal.uns.pops,pop.h.angle,H=h.mal,h.eig=1)),
	MalAngle2=do.call("rbind",
	lapply(pmat.mal.uns.pops,pop.h.angle,H=h.mal,h.eig=2)),
	MalAngle3=do.call("rbind",
	lapply(pmat.mal.uns.pops,pop.h.angle,H=h.mal,h.eig=3)))

#####################TENSORS############################
#Adapted from Aguirre et al. 2013 supplemental material
pmat.mal.unst.pops<-calc.pmat(pops.mal.unstd.dat,3,8)
pmat.mal.unst.pops<-pmat.mal.unst.pops[match(pop.list, names(pmat.mal.unst.pops))]
pmat.fem.unst.pops<-calc.pmat(pops.fem.unstd.dat,3,10)
pmat.fem.unst.pops<-pmat.fem.unst.pops[match(pop.list, names(pmat.fem.unst.pops))]

covtensor<-function(matrix.list){
#Adapted from Aguirre et al. 2013 supplemental material
	n.traits<-dim(matrix.list[[1]])[1]
	n.pops<-length(matrix.list)
	trait.names<-colnames(matrix.list[[1]])
	pop.names<-names(matrix.list)
	n.eigten<-n.traits*(n.traits+1)/2
	mat.var<-matrix(nrow=n.pops,ncol=n.traits)
	mat.cov<-matrix(nrow=n.pops,ncol=length(lowerTriangle(matrix.list[[1]])))
	for(i in 1:n.pops){
		mat.var[i,]<-diag(matrix.list[[i]])
		mat.cov[i,]<-lowerTriangle(matrix.list[[i]])
	}
	#create the S matrix
	s.mat<-matrix(nrow=n.eigten,ncol=n.eigten)
	colnames(s.mat)<-paste("E", 1:n.eigten, sep="")
	rownames(s.mat)<-paste("E", 1:n.eigten, sep="")
	#fill in upper left quadrant
	s.mat[1:n.traits,1:n.traits]<-cov(mat.var,mat.var)
	#fill in lower right quadrant
	s.mat[(n.traits+1):n.eigten,(n.traits+1):n.eigten]<-2*cov(mat.cov,mat.cov)
	#fill in upper right quadrant
	s.mat[1:n.traits,(n.traits+1):n.eigten]<-sqrt(2)*cov(mat.var, mat.cov)
	#fill in lower left quadrant
	s.mat[(n.traits+1):n.eigten,1:n.traits]<-sqrt(2)*cov(mat.cov, mat.var)
	
	#eigenvalues and vectors of S
	s.val<-eigen(s.mat)$values
	s.vec<-eigen(s.mat)$vectors

	#make the eigentensors matrix
	et.mat<-array(, c(n.traits, n.traits, n.eigten))
	dimnames(et.mat) <- list(trait.names, trait.names, 
		paste("E", 1:n.eigten, sep=""))  
	for(i in 1:n.eigten){
		e.mat <- matrix(0, n.traits, n.traits) 
		lowerTriangle(e.mat) <- 1/sqrt(2)*s.vec[(n.traits+1):n.eigten,i]
		e.mat <- e.mat + t(e.mat)
		diag(e.mat) <- s.vec[1:n.traits,i]
		et.mat[,,i] <- e.mat 
	}
	
	#eigenvectors of eigentensors
	et.eigen<-array(,c((n.traits+1),n.traits, n.eigten))
	for(i in 1:n.eigten){
		#eigenvalues of ith eigentensor
		et.eigen[1,,i]<-t(eigen(et.mat[,,i])$values)
		#eigenvectors of the ith eigentensor
		et.eigen[2:(n.traits+1),,i]<-eigen(et.mat[,,i])$vectors
		et.eigen [,,i]<-et.eigen[,order(abs(et.eigen[1,,i]),
			decreasing = T), i]
	}

	#project eigenvectors of s onto s to determine alpha
	s.alpha<-matrix(ncol=n.eigten)
	for(i in 1:n.eigten){
		s.alpha[,i]<-t(s.vec[,i]) %*% s.mat %*% s.vec[,i]
	}
	#distribution of phenotypic variance for eigenvectors of S
	p.coord<-matrix(nrow=n.pops, ncol=n.eigten)
	for(i in 1:n.eigten){
		p.coord[,i]<-unlist(lapply(matrix.list, frobenius.prod, y = et.mat[,,i]))
	}
	rownames(p.coord)<-pop.names
	colnames(p.coord)<-paste("E", 1:n.eigten, sep="")
	tensor.summary <- data.frame(rep(s.val,each=n.traits), 
		t(data.frame(et.eigen)))
	colnames(tensor.summary) <- c("S.eigval", "eT.val", trait.names)
	rownames(tensor.summary)<- paste(
		paste("e", rep(1:n.eigten, each=n.traits), sep=""), 
		rep(1:n.traits,n.eigten), sep=".")
	#maximum number of nonzero eigentensors
	nonzero<-min(n.traits*(n.traits+1)/2,n.pops-1)
	

	return(list(tensor.summary = tensor.summary, s.mat = s.mat, 
		et.mat = et.mat, p.coord = p.coord, s.val = s.val, 
		s.alpha=s.alpha, nonzero=nonzero))
	#tensor.summary: sumamry of the covariance tensor for S.
		#contains eigenvalues of the tensor and 
		#the eigenvalues of the eigentensors. And the eigenvectors of 
		#eigentensors with trait loadings identified by the column names.
	#s.mat: the S matrix
	#et.mat: the eigentensors of S. Rows and columns identify elements 
		#of an eigentensor. 3rd dimension identifies the eigentensor
	#p.coord: Coordinates of the P matrix in the space of the eigentensors 
		#of S
	#s.val: eigenvalues of S
	#s.alpha: projection of the eigenvectors of S on S.
	#nonzero: The maximum number of nonzero eigentensors.
}
f.n.traits<-8
m.n.traits<-6
pop.names<-names(pmat.fem.unst.pops)
fem.tensor<-covtensor(pmat.fem.unst.pops)
mal.tensor<-covtensor(pmat.mal.unst.pops)
#eigenvalues for the nonzero eigentensors
fem.tensor$s.alpha[,1:fem.tensor$nonzero]
fem.tensor$tensor.summary[
	1:((ncol(fem.tensor$tensor.summary)-2)*fem.tensor$nonzero),]
#plot coordinates of female p matrix in the space of e1 and e2 
plot(fem.tensor$p.coord[,1], ylim=c(-200,50), xaxt="n", las=1,
	xlab="Population", ylab="alpha")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.names, par("usr")[1]-230,
	srt=-45, xpd=TRUE)
lines(fem.tensor$p.coord[,1], lty=2)
points(fem.tensor$p.coord[,2], pch=19)
lines(fem.tensor$p.coord[,2], lty=1)
points(fem.tensor$p.coord[,3], pch=15)
lines(fem.tensor$p.coord[,3],lty=4)
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4), c("e1","e2","e3"))
#trait combinations for e1, e2, and e3
round(fem.tensor$tensor.summary[1:(f.n.traits*3),
	2:dim(fem.tensor$tensor.summary)[2]], 3)

#determine which eigenvector explains the most variation in eigentensors
max.eig.val<-function(eig.vec){
	max.perc<-max(abs(eig.vec))/sum(abs(eig.vec))
	max.loc<-which.max(abs(eig.vec))
	return(list(max.loc=max.loc, max.perc=max.perc))
}
e1.max<-max.eig.val(fem.tensor$tensor.summary$eT.val[1:f.n.traits])
e2.max<-max.eig.val(fem.tensor$tensor.summary$eT.val[(f.n.traits+1):(2*f.n.traits)])
e3.max<-max.eig.val(fem.tensor$tensor.summary$eT.val[((f.n.traits*2)+1):(3*f.n.traits)])


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

#Function to do projection
proj<- function(G, b) t(b) %*% G %*% (b)

#variance along e1 and e2
f.e11.proj <- lapply(pmat.fem.unst.pops, proj, b = f.e1.L)
f.e21.proj <- lapply(pmat.fem.unst.pops, proj, b = f.e2.L)
f.e31.proj <- lapply(pmat.fem.unst.pops, proj, b = f.e3.L)

m.e11.proj <- lapply(pmat.mal.unst.pops, proj, b = m.e1.L)
m.e21.proj <- lapply(pmat.mal.unst.pops, proj, b = m.e2.L)
m.e31.proj <- lapply(pmat.mal.unst.pops, proj, b = m.e3.L)

#summary plots
jpeg("blows.method3.summary.jpeg", width=7, height=7, units="in", res=300)
layout(matrix(c(0,1,2,3),2,2,byrow=F))
layout.show(3)

#plot eigenvalues of non-zero eigentensors for S
par(mar=c(5,4,2,0))
plot(fem.tensor$s.alpha[,1:fem.tensor$nonzero], ylab="alpha",
	xlab="", xaxt="n", las=1)#3 until it levels off.
axis(1, at=seq(1,fem.tensor$nonzero,1), 
	labels=paste("E",1:fem.tensor$nonzero, sep=""))
points(mal.tensor$s.alpha[,1:mal.tensor$nonzero],pch=6)
legend("top", pch=c(1,6), c("Males", "Females"))
text(x=11,y=2000, "B", font=2)

#plot the variance in each population in the direction of e11,e21, and e31
par(mar=c(3,5,2,1))
plot(x=seq(1,12,1),f.e11.proj, xaxt="n", las=1,ylim=c(-50,200),
	xlab="", ylab="lambda")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.names, par("usr")[1]-80,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),f.e11.proj, lty=2)
points(x=seq(1,12,1),f.e21.proj, pch=19)
lines(x=seq(1,12,1),f.e21.proj, lty=1)
points(x=seq(1,12,1),f.e31.proj, pch=15)
lines(x=seq(1,12,1),f.e31.proj,lty=4)
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4), c("e1","e2","e3"))
text(x=2,y=200, "Females")
text(x=12,y=200, "C", font=2)

par(mar=c(5,5,0,1))
plot(x=seq(1,12,1),m.e11.proj, xaxt="n", las=1,ylim=c(-50,200),
	xlab="Population", ylab="lambda")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.names, par("usr")[1]-80,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),m.e11.proj, lty=2)
points(x=seq(1,12,1),m.e21.proj, pch=19)
lines(x=seq(1,12,1),m.e21.proj, lty=1)
points(x=seq(1,12,1),m.e31.proj, pch=15)
lines(x=seq(1,12,1),m.e31.proj,lty=4)
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4), c("e1","e2","e3"))
text(x=2,y=200, "Males")
text(x=12,y=200, "D", font=2)
dev.off()


###Try to plot pc1 and pc2 for males and females
plot.eigenvectors<-function(dat.cov, col.num, eig.lty=1, eig.col="black"){
	#Plots the eigenvectors (scaled by the eigenvalues) 
	j = 1.96*sqrt(eigen(dat.cov)$values[col.num])
	delx1 = 0 + j*eigen(dat.cov)$vectors[1,col.num]
	delx2 = 0 + j*eigen(dat.cov)$vectors[2,col.num]
	delx11 = 0 - j*eigen(dat.cov)$vectors[1,col.num]
	delx22 = 0 - j*eigen(dat.cov)$vectors[2,col.num]
	x.values=c(0 , delx1, delx11)
	y.values=c(0 , delx2, delx22)
	lines(x.values, y.values, lwd=2, col=eig.col, lty=eig.lty)
}
plot.pmatrix<-function(dat, plot.new=TRUE, eig.col="black", eig.lty=1){
	#get leading eigenvector for bands
	#and leading eigenvector for length
	#Plots the P-matrix 95% confidence ellipse and eigenvectors.
	if(plot.new==TRUE){
	plot(c(-1,1), c(-1,1) , asp=1, type="n", xlim=c(-1,1), ylim=c(-1,1), 
		cex.lab=1.5, xlab="",ylab="", las=1)
	}
	dat.cov<-cov(dat)
	ellipse(center=c(0, 0), shape=dat.cov, lty=eig.lty,
		radius=1.96, center.cex=0.1, 
		lwd=2, col=eig.col, add=TRUE )
	plot.eigenvectors(dat.cov, 1, eig.col=eig.col, eig.lty=eig.lty)
	plot.eigenvectors(dat.cov, 2, eig.col=eig.col, eig.lty=eig.lty)
}

band.u<-list()
fem.pheno.u<-list()
mal.pheno.u<-list()
for(i in 1:length(band.pca)){
	band.u[[i]]<-data.frame(cbind(as.numeric(band.pca[[i]]$CA$u[,1]),
		as.numeric(band.pca[[i]]$CA$u[,2])))
	fem.pheno.u[[i]]<-data.frame(cbind(as.numeric(fem.pheno.pca[[i]]$CA$u[,1]),
		as.numeric(fem.pheno.pca[[i]]$CA$u[,2])))
	mal.pheno.u[[i]]<-data.frame(cbind(as.numeric(mal.pheno.pca[[i]]$CA$u[,1]),
		as.numeric(mal.pheno.pca[[i]]$CA$u[,2])))
}
names(band.u)<-names(band.pca)
names(fem.pheno.u)<-names(band.pca)
names(mal.pheno.u)<-names(band.pca)

jpeg("P-matrices.female.jpeg", res=300, height=14, width=14, units="in")
par(mfrow=c(3,4), mar=c(2,2.5,2,2.5), oma=c(2,2.5,2,2.5))
for(i in 1:length(fem.pheno.u)){
	plot.pmatrix(fem.pheno.u[[pop.list[i]]])
	plot.pmatrix(cbind(fem.pheno.u[[pop.list[i]]][,1],
		band.u[[pop.list[i]]][,1]),F, "blue", 2)
	axis(4, las=1)
	legend("top", pop.list[i], bty="n", cex=1.5)
}
par(fig = c(0, 1, 0, 1), oma=c(2,2.5,1,2.5), mar = c(0, 0, 0, 0), 
	new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", col=c("black","blue"),cex=1.5,
	c("Body","Bands"),lty=c(1,2), box.lty=0,ncol=4)
mtext("Body PC1",1, outer=T, line=0, cex=1)
mtext("Body PC2",2, outer=T, line=0, cex=1)
mtext("Band PC1",4, outer=T, line=1, cex=1)
dev.off()

jpeg("P-matrices.male.jpeg", res=300, height=14, width=14, units="in")
par(mfrow=c(3,4), mar=c(2,2.5,2,2), oma=c(2,2.5,2,2))
for(i in 1:length(mal.pheno.u)){
	plot.pmatrix(mal.pheno.u[[pop.list[i]]])
	legend("top", pop.list[i], bty="n", cex=1.5)
}
par(fig = c(0, 1, 0, 1), oma=c(2,2.5,1,2), mar = c(0, 0, 0, 0), 
	new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
mtext("Body PC1",1, outer=T, line=0.5, cex=1)
mtext("Body PC2",2, outer=T, line=0, cex=1)
dev.off()







