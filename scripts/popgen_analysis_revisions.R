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

setwd("B:/ubuntushare/popgen/")


#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################
summ.dat<-read.table("sw_results/stacks/populations/batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
ld.hwe<-read.table("sw_results/stacks/populations/ld.hwe.whitelist.txt")

catalog<-read.delim("results/stacks/batch_1.catalog.tags.tsv", 
	header=F)

mar.coor<-read.csv("E://Docs//PopGen//marine_coordinates.csv", header=T)
fw.coor<-read.csv("E://Docs//PopGen//fw_coordinates.csv", header=T)
m.f.summ.dat<-read.table("results//stacks//populations_sex//batch_1.sumstats.tsv",
	sep='\t', skip=2, header=T, comment.char="")
dist<-read.table("E://Docs//PopGen//geographical_distances.txt", 
	header=T, row.names=1, sep='\t')
pwise.fst<-read.table("results\\stacks\\populations\\batch_1.fst_summary.tsv",
	 header=T, row.names=1, sep='\t')
	pwise.fst<-rbind(pwise.fst,rep(NA, ncol(pwise.fst)))
	rownames(pwise.fst)<-colnames(pwise.fst)

#####Re-name plink files so that Family ID contains Population ID.
ped<-read.table("sw_results/stacks/populations/subset.ped")
ped$V1<-gsub("sample_(\\w{4})\\w+.*align","\\1",ped$V2)
write.table(ped,"sw_results/migrate/subset.ped",col.names=F,row.name=F,quote=F)


#############################################################################
#######################PLOT THE POINTS ON A MAP##############################
#############################################################################
mar.coor$lon<-(-1*mar.coor$lon)
fw.coor$lon<-(-1*fw.coor$lon)

jpeg("mar_sites_map.jpg", res=300, width=9, height=7, units="in")
map("worldHires", "usa",xlim=c(-100,-76), ylim=c(23,35.5), 
	col="gray90", fill=TRUE)
map("worldHires", "mexico",xlim=c(-100,-76), ylim=c(23,35.5), 
	col="gray95", fill=TRUE, add=TRUE)
map("worldHires", "cuba",xlim=c(-100,-76), ylim=c(23,35), 
	col="gray95", fill=TRUE, add=TRUE)
points(mar.coor$lon, mar.coor$lat,  col="black", cex=1, pch=19)
abline(h=c(25,30,35),lty=3)
abline(v=c(-80,-85,-90,-95,-100),lty=3)
text(x=c(-99.5,-99.5,-99.5),y=c(25,30,35),c("25N","30N","35N"))
text(x=c(-80,-85,-90,-95),y=rep(35.3,4),c("80W","85W","90W","95W"))
text(y=26,x=-90,"Gulf of Mexico")
text(x=-88,y=32,"USA")
text(x=-78,y=29.5,"Atlantic Ocean")
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
pwise.fst<-read.table("sw_results\\stacks\\populations\\batch_1.fst_summary.tsv",
	 header=T, row.names=1, sep='\t')
	pwise.fst<-rbind(pwise.fst,rep(NA, ncol(pwise.fst)))
	rownames(pwise.fst)<-colnames(pwise.fst)

ibd<-mantel.rtest(as.dist(dist),as.dist(pwise.fst))



#########################################################################
#***********************************************************************#
########################POPULATION STRUCTURE#############################
#***********************************************************************#
#########################################################################

#******************************ADEGENET*********************************#
dat.plink<-read.PLINK(
	"E:/ubuntushare/popgen/sw_results/stacks/populations/subsetA.raw",
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
pop.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
	"FLAB","FLPB","FLHB","FLCC")
pop<-substr(ind.names, 8,11)
col<-pop
col[col=="TXSP"]<-rainbow(12)[1]
col[col=="TXCC"]<-rainbow(12)[2]
col[col=="TXCB"]<-rainbow(12)[3]
col[col=="ALST"]<-rainbow(12)[4]
col[col=="FLSG"]<-rainbow(12)[5]
col[col=="FLKB"]<-rainbow(12)[6]
col[col=="FLFD"]<-rainbow(12)[7]
col[col=="FLSI"]<-rainbow(12)[8]
col[col=="FLAB"]<-rainbow(12)[9]
col[col=="FLPB"]<-rainbow(12)[10]
col[col=="FLHB"]<-rainbow(12)[11]
col[col=="FLCC"]<-rainbow(12)[12]

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
#try another number
dat.clust.5<-find.clusters(dat.plink, parallel=FALSE, n.pca=20, n.clust=5)
dapc5<-dapc(dat.plink, dat.clust.5$grp, n.pca=20,n.da=3, parallel=F)
scatter(dapc5, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
compoplot(dapc5)

#output k=3 clusters
adegenet.groups<-as.data.frame(cbind(names(dat.clust$grp), dat.clust$grp))
adegenet.groups[,1]<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', adegenet.groups[,1])
adegenet.groups[,1]<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', 
	adegenet.groups[,1])


#***************************FASTSTRUCTURE*********************************#
str.in<-read.table("E://ubuntushare//stacks//populations//ld.hwe.pgd.str")
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
#K=2
stru2<-read.table("E://ubuntushare//ld.hwe_out_simple.2.meanQ",header=F)
structure.barplot(stru2,2, pop.list)
stru3<-read.table("E://ubuntushare//pop_structure//faststructure//ld.hwe_out_simple.3.meanQ",header=F)
stru3<-cbind(pop.id,stru3)
stru4<-read.table("E://ubuntushare//pop_structure//faststructure//ld.hwe_out_simple.4.meanQ",header=F)
stru4<-cbind(pop.id,stru4)
stru5<-read.table("E://ubuntushare//pop_structure//faststructure//ld.hwe_out_simple.5.meanQ",header=F)
stru5<-cbind(pop.id,stru5)
stru6<-read.table("E://ubuntushare//pop_structure//faststructure//ld.hwe_out_simple.6.meanQ",header=F)
stru6<-cbind(pop.id,stru6)

#plot pop averages
jpeg("faststructure.k3-6.summ.jpeg", height=12, width=12, units="in", res=300)
par(mfrow=c(2,2),mar=c(5,4,4,1),oma=c(2,2,2,1))
#K=3
structure.barplot(stru3,3, pop.list, FALSE)
#K=4
structure.barplot(stru4,4, pop.list, FALSE)
#K=5
structure.barplot(stru5,5, pop.list, FALSE)
#K=6
structure.barplot(stru6,6, pop.list, FALSE)
dev.off()

#plot structure-like plot
plot.structure(stru3, 3, pop.list ,"faststructure.k3.jpeg")
plot.structure(stru4, 4, pop.list ,"faststructure.k4.jpeg")
plot.structure(stru5, 5, pop.list ,"faststructure.k5.jpeg")
plot.structure(stru6, 6, pop.list ,"faststructure.k6.jpeg")

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




#*****************************STRUCTURE***********************************#
structure.k3<-read.table(
	"E://ubuntushare//stacks//populations//structure//popstr//admixture_set//Results//admixture_k3//admixture_set_run_13_f_clusters.txt",
	sep='\t', header=F)
structure.k3$V1<-sub('sample_([A-Z]{4})','\\1', structure.k3$V1)
structure.k4<-read.table(
	"E://ubuntushare//stacks//populations//structure//popstr//admixture_set//Results//admixture_k4//admixture_set_run_23_f_clusters.txt",
	sep='\t', header=F)
structure.k4$V1<-sub('sample_([A-Z]{4})','\\1', structure.k4$V1)
structure.k5<-read.table(
	"E://ubuntushare//stacks//populations//structure//popstr//admixture_set//Results//admixture_k5//admixture_set_run_33_f_clusters.txt",
	sep='\t', header=F)
structure.k5$V1<-sub('sample_([A-Z]{4})','\\1', structure.k5$V1)
structure.k6<-read.table(
	"E://ubuntushare//stacks//populations//structure//popstr//admixture_set//Results//admixture_k6//admixture_set_run_43_f_clusters.txt",
	sep='\t', header=F)
structure.k6$V1<-sub('sample_([A-Z]{4})','\\1', structure.k6$V1)

pop.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
		"FLAB","FLPB","FLHB","FLCC")

plot.structure<-function(structure.out, k, pop.order, 
	filename=paste("str.k",k,".jpeg",sep=""),make.file=TRUE,
	plot.new=TRUE){
	str.split<-split(structure.out,structure.out[,1])
	bar.colors<-rainbow(k,s=0.5)
	if(make.file==TRUE){
		jpeg(filename,width=7, height=1.25, units="in", res=300)
		par(mfrow=c(1,length(str.split)))
	} 
	par(mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
	for(i in 1:length(str.split)){
		pop.index<-pop.order[i]
		barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
			beside=FALSE, space=0,	border=NA, col=bar.colors,
			xlab="", ylab="", xaxt='n', yaxt='n', new=plot.new)
		mtext(pop.index, 1, line=0.5, cex=0.5, outer=F)
	}
	if(make.file==TRUE) {dev.off()}
}
par(mfrow=c(4,length(pop.list)),mar=c(1,0,1,0),oma=c(1,0,1,0))
plot.structure(structure.k3,3,pop.list)#, make.file=FALSE, plot.new=FALSE)
plot.structure(structure.k4,4,pop.list)#, make.file=FALSE, plot.new=TRUE)
plot.structure(structure.k5,5,pop.list)#, make.file=FALSE, plot.new=TRUE)
plot.structure(structure.k6,6,pop.list)#, make.file=FALSE, plot.new=TRUE)


#*******************************PCADAPT***********************************#
#K=4 WAS BEST
setwd("E:/ubuntushare/popgen/sw_results/pcadapt")
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
sub.map<-read.table("../stacks/populations/subset.map",header=F)
pcadapt.outliers<-sub.map[rep.snps,]
pa.out.radloc<-sub('(\\d+)_\\d+','\\1',pcadapt.outliers$V2)

summ.dat<-read.table("../stacks/populations/batch_1.sumstats.tsv",
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


#########################################################################
#***********************************************************************#
################################BAYENV2##################################
#***********************************************************************#
#########################################################################


#**************************STARTING WITH PLINK FILES********************#
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
map<-read.table("results/stacks/populations/ld.hwe.sub.map", skip = 1)
chr.nums<-sub('scaffold_','',map[,1])
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

#***************************CONVERT PLINK TO BAYENV2************************#
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

#*******************************check Bayenv2 matrix*************************#
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

#************************************SNPFILEs******************************#
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

#************************************ENVFILE******************************#
env.raw<-read.table("sw_results/environmental_assoc/new_bayenv/env_data_bayenv_raw.txt")
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
	"sw_results/environmental_assoc/new_bayenv/env_data_bayenv_std.txt",
	sep='\t',quote=F,col.names=F,row.names=F,eol='\n')


##############################################################################
#****************************************************************************#
####################################PST-FST###################################
#****************************************************************************#
##############################################################################
raw.pheno<-read.table("sw_results/pstfst/popgen.pheno.txt", sep="\t", header=T)

#**********************************RAW**************************************#
#************create male and female files for pst analysis******************#
fem.keep<-raw.pheno[!is.na(raw.pheno$BandNum) & substr(raw.pheno$ID,5,5)=="F",]
db.keep<-raw.pheno[substr(raw.pheno$ID,5,6)=="DB",]
db.keep<-db.keep[(db.keep$Side=="LEFT" || 
	db.keep$Side == "Left"),]
not.fem.keep<-raw.pheno[substr(raw.pheno$ID,5,5)!="F" & 
	substr(raw.pheno$ID,5,6)!="DB",]
not.fem.keep<-not.fem.keep[(not.fem.keep$Side=="LEFT" || 
	not.fem.keep$Side == "Left"),]
pheno.dat<-rbind(fem.keep, not.fem.keep, db.keep)

#read in ped file from whitelisted populations run
ped<-read.table("sw_results/stacks/populations/batch_1.plink.ped", 
	skip = 1, stringsAsFactors=F, colClasses="character") 
ped.names<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', ped[,2])
ped.names<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', ped.names)
pheno.dat$ID<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', pheno.dat$ID)
pops.pheno<-pheno.dat[pheno.dat$ID %in% ped.names,] #this does not have the juveniles
pops.pheno<-pops.pheno[match(ped.names,pops.pheno$ID),]
#pops.pheno<-replace(pops.pheno, is.na(pops.pheno),as.numeric(-9))
pops.pheno<-pops.pheno[,-8]
pops.pheno<-cbind(substr(pops.pheno$ID,1,4), pops.pheno)
colnames(pops.pheno)[1]<-"PopID"
fem.pheno<-pops.pheno[pops.pheno$BandNum!=-9,]
fem.pheno<-replace(fem.pheno, fem.pheno==-9,NA)
write.table(fem.pheno, 
	"sw_results/pstfst/fem.pheno.txt",
	sep="\t", quote=F, col.names=T, row.names=F)
mal.pheno<-pops.pheno[pops.pheno$BandNum==-9,]
mal.pheno<-replace(mal.pheno, mal.pheno==-9, NA)
mal.pheno<-mal.pheno[,-9]
mal.pheno<-mal.pheno[,-9]
mal.pheno<-mal.pheno[complete.cases(mal.pheno),]
write.table(mal.pheno, 
	"sw_results/pstfst/mal.pheno.txt",
	sep="\t", quote=F, col.names=T, row.names=F)
trait.names<-c("SVL", "Standard Length", "Body Depth", 
	"Snout Length", "Snout Depth", "Head Length")

#*****************************STANDARDIZED**********************************#
#************create male and female files for pst analysis******************#
#standardize pop.pheno by population.
std.by.mean<-function(trait){
	trait.trans <- trait/mean(trait,na.rm=F) -1
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
#pops.pheno<-pops.pheno[complete.cases(pops.pheno),]
pheno.split<-split(pops.pheno,pops.pheno$PopID)
pheno.split.std<-lapply(pheno.split,function(x){
	std<-data.frame(x[,1:2])
	for(i in 3:8){	
		tout<-std.by.mean(x[,i])
		std<-cbind(std,tout)
	}
	return(std)
})
pops.pheno.std<-data.frame(do.call("rbind",pheno.split.std))
colnames(pops.pheno.std)<-colnames(pops.pheno)[1:8]

fem.pheno.std<-pops.pheno.std[pops.pheno.std$ID %in% fem.pheno$ID,]
mal.pheno.std<-pops.pheno.std[pops.pheno.std$ID %in% mal.pheno$ID,]


fem.split<-split(fem.pheno,fem.pheno$PopID)
fem.split.std<-lapply(fem.split,function(x){
	std<-data.frame(x[,1:2])
	for(i in 9:10){	
		tout<-std.by.mean(x[,i])
		std<-cbind(std,tout)
	}
	return(std)
})
fem.pheno.std<-merge(fem.pheno.std,do.call("rbind",fem.split.std)[,2:4],by="ID")
colnames(fem.pheno.std)<-c("ID","PopID",colnames(fem.pheno)[3:10])

#reorder them to match dist file
fem.pheno.sep<-split(fem.pheno.std, fem.pheno$PopID)
fem.pheno.new<-rbind(fem.pheno.sep$TXSP,fem.pheno.sep$TXCC,fem.pheno.sep$TXCB,
	fem.pheno.sep$ALST,fem.pheno.sep$FLSG,fem.pheno.sep$FLKB,
	fem.pheno.sep$FLFD,fem.pheno.sep$FLSI,fem.pheno.sep$FLAB,
	fem.pheno.sep$FLPB,fem.pheno.sep$FLHB,fem.pheno.sep$FLCC)
mal.pheno.sep<-split(mal.pheno.std, mal.pheno$PopID)
mal.pheno.new<-rbind(mal.pheno.sep$TXSP,mal.pheno.sep$TXCC,mal.pheno.sep$TXCB,
	mal.pheno.sep$ALST,mal.pheno.sep$FLSG,mal.pheno.sep$FLKB,
	mal.pheno.sep$FLFD,mal.pheno.sep$FLSI,mal.pheno.sep$FLAB,
	mal.pheno.sep$FLPB,mal.pheno.sep$FLHB,mal.pheno.sep$FLCC)


pairwise.pst<-function(dat){
	#first column must be pop id/grouping factor
	names<-levels(factor(dat[,1]))
	dat.split<-split(dat, factor(dat[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(names),numeric(0), simplify = F), names))
	for(i in 1:(length(dat.split)-1)){
	  for(j in (i+1):length(dat.split)){
		temp.data<-rbind(as.data.frame(dat.split[[i]]),
			as.data.frame(dat.split[[j]]))
		aov.var<-summary.aov(
			aov(temp.data[,2]~temp.data[,1]))[[1]]$`Sum Sq`
		aov.df<-summary.aov(
			aov(temp.data[,2]~temp.data[,1]))[[1]]$`Df`
		dat.var[i,j]<-aov.var[2]/(aov.var[2]+
			(2*(aov.var[1]/(aov.df[2]-1))))	
	  }
	}
	dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
	rownames(dat.var)<-colnames(dat.var)
	return(dat.var)
}

all.traits.pst.mantel<-function(trait.df,comp.df,id.index){
	results.mantel<-NULL
	for(i in 3:ncol(trait.df)){
		res<-mantel.rtest(as.dist(t(pairwise.pst(
			trait.df[,c(id.index,i)]))),
			as.dist(t(comp.df)), nrepet=9999)
		results.mantel<-c(results.mantel,res)
	}
	return(results.mantel$pvalue)
}
fem.fst.pst<-all.traits.pst.mantel(fem.pheno.new,pwise.fst,2)

