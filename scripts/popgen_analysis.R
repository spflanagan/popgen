#Author: Sarah P. Flanagan
#Date: 27 July 2015
#Purpose: Analyze Population genetics data
#made from used analyses from fst_cline_analysis and population_structure

rm(list=ls())

library(ade4)
library(lme4)
library(maps)
library(mapdata)
library(vegan)
library(boot)
library(car)
library(adegenet)
library(scales)
library(gdata)

#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################
summ.dat<-read.table("E://ubuntushare//stacks//populations//batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")

contigs<-read.table("E://ubuntushare//linkage_map//ordered_contigs_denovo50.txt",
	sep='\t', header=T)
	cont.by.scaf<-split(contigs, contigs$Scaffold)
	n.lg.per.scaff<-lapply(cont.by.scaf, 
		function(x){length(levels(as.factor(x$LG)))})
	mult.lg<-unlist(n.lg.per.scaff[n.lg.per.scaff!=1])
	cont.mult<-contigs[contigs$Scaffold %in% as.factor(names(mult.lg)),]
	cont.mult$Scaffold<-factor(cont.mult$Scaffold)
	cont.mult.s<-split(cont.mult, cont.mult$Scaffold)
	dup.scaff<-contigs$Scaffold %in% as.factor(names(mult.lg))
use.contigs<-contigs[!dup.scaff,]

catalog<-read.delim("E://ubuntushare//stacks//batch_1.catalog.tags.tsv", 
	header=F)
sig.diff<-read.delim("E://ubuntushare//fst_geographical_clines//fst_geographical_clines//fst_sig_diff_4cline.txt", 
	header=T, sep='\t')
mar.coor<-read.csv("E://Docs//PopGen//marine_coordinates.csv", header=T)
fw.coor<-read.csv("E://Docs//PopGen//fw_coordinates.csv", header=T)
m.f.summ.dat<-read.table("E://ubuntushare//stacks//populations_sex//batch_1.sumstats.tsv",
	sep='\t', skip=2, header=T, comment.char="")
sig.diff<-read.delim("E://ubuntushare//fst_geographical_clines//fst_geographical_clines//fst_sig_diff.txt", header=T, sep='\t')

dist<-read.table("C:\\Users\\sflanagan\\Dropbox\\RAD\\goegraphical_distances.txt", 
	header=T, row.names=1, sep='\t')
pwise.fst<-read.table("E:\\ubuntushare\\stacks\\populations\\batch_1.fst_summary.tsv",
	 header=T, row.names=1, sep='\t')
	pwise.fst<-rbind(pwise.fst,rep(NA, ncol(pwise.fst)))
	rownames(pwise.fst)<-colnames(pwise.fst)

setwd("E://Docs//PopGen")
#############################################################################
#***************************************************************************#
#################################FUNCTIONS###################################
#***************************************************************************#
#############################################################################

#***************************************************************************#
#PLOT ANY GENOME-WIDE STATISTIC
#***************************************************************************#
plot.genome.wide<-function(bp,var,y.max,x.max, rect.xs,y.min=0,x.min=0, 
	plot.new=FALSE, plot.axis=TRUE, rect.color="white", pt.cex=1){
	par(new=new)
	plot(bp, var,xlab="",ylab="", 
		type="n", bg="transparent", axes=F, bty="n", 
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	num.rect<-nrow(rect.xs)
	if(is.null(num.rect)) {
		rect(rect.xs[1],y.min,rect.xs[2],y.max, 
				col=rect.color, border=NA)
	} else {
		for(i in 1:nrow(rect.xs)){
			rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
				col=rect.color, border=NA)
		}
	}
	if(plot.axis){
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
	points(bp, var, pch=19, cex=pt.cex,
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

#***************************************************************************#
##AUTOMATED FST PLOTTING
#***************************************************************************#
plot.fsts.scaffs<-function(dat, dat.name, ci.dat=NULL, 
	col.pts=NULL, col.pt.col="dark green", col.pt.pch=8, col.pt.cex=2){
	byscaff<-split(dat, factor(dat$Chr))#dat is a fst file from stacks_genomewideCIs
	if(!is.null(col.pts)){
		col.pt.byscaff<-split(col.pts, factor(col.pts$Chr))}
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	for(i in 1:length(byscaff)){
		new.max<-last.max+round(max(byscaff[[i]]$BP), -2)
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	
	new.x<-list()
	col.pt.x<-NULL
	col.pt.y<-NULL
	for(i in 1:length(byscaff)){
		new.x[[i]]<-byscaff[[i]]$BP+addition.values[i]
		if(!is.null(col.pts)){
		j<-match(names(byscaff[i]),names(col.pt.byscaff))
		if(!is.na(j)){
			col.pt.x<-c(col.pt.x, byscaff[[i]][(
				byscaff[[i]]$Locus %in% 
				col.pt.byscaff[[j]]$Locus) &
				(byscaff[[i]]$BP %in% col.pt.byscaff[[j]]$BP)
				,"BP"]+addition.values[i])
			col.pt.y<-c(col.pt.y, byscaff[[i]][(
				byscaff[[i]]$Locus
				 %in% col.pt.byscaff[[j]]$Locus) &
				(byscaff[[i]]$BP %in% col.pt.byscaff[[j]]$BP),
				"SmoothFst"])
		}}
	}

	x.min<-min(addition.values)
	x.max<-max(addition.values)
	y.max<-max(dat$SmoothFst)+0.5*max(dat$SmoothFst)
	if(min(dat$SmoothFst) < 0) {
		y.min<-min(dat$SmoothFst) + 0.5*min(dat$SmoothFst)
	} else {
		y.min<-0
	}

	plot(new.x[[1]], byscaff[[1]]$SmoothFst, 
		xlim=c(x.min,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(byscaff)){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(new.x[[j]], 
			byscaff[[j]]$SmoothFst,
			y.max,x.max, rect.xs[j,],y.min=y.min,x.min=x.min, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
		
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=10)),
		xlim=c(x.min,x.max),ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
	mtext(dat.name,2,las=1, line=4,cex=.75)
	if(!is.null(ci.dat)){
		clip(x.min,x.max,y.min,y.max)
		abline(h=ci.dat$CI99smooth,col="red")
		if(!is.null(col.pts)){
			clip(x.min,x.max,ci.dat$CI99smooth, y.max)
			points(col.pt.x, col.pt.y,
				col=col.pt.col, pch=col.pt.pch, cex=col.pt.cex)
		}
	}
}


#############################################################################
#######################PLOT THE POINTS ON A MAP##############################
#############################################################################
mar.coor$lon<-(-1*mar.coor$lon)
fw.coor$lon<-(-1*fw.coor$lon)

jpeg("mar_sites_map.jpg", res=300, width=9, height=7, units="in")
map("worldHires", "usa",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray90", fill=TRUE)
map("worldHires", "mexico",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray95", fill=TRUE, add=TRUE)
map("worldHires", "cuba",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray95", fill=TRUE, add=TRUE)
#points(mar.coor$lon, mar.coor$lat,  col="red", cex=1, pch=19)
text(mar.coor$lon, mar.coor$lat, labels=mar.coor$site, col="red")
dev.off()

#fst cline analysis graphic
jpeg("mar_sites_clines.jpg", res=300, width=9, height=7, units="in")
map("worldHires", "usa",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray95", fill=TRUE)
map("worldHires", "mexico",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray100", fill=TRUE, add=TRUE)
map("worldHires", "cuba",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray100", fill=TRUE, add=TRUE)
mar.north<-rbind(mar.coor[3,],mar.coor[6,],mar.coor[12,])
mar.south<-rbind(mar.coor[1,],mar.coor[8,],mar.coor[10,])

arrows(mar.north$lon[1]+1, mar.north$lat[1],
	mar.north$lon[2]-1, mar.north$lat[2],
	col="darkgreen", lwd=2,length=0.15, code=3)
arrows(mar.north$lon[2]+1, mar.north$lat[2],
	mar.north$lon[3], mar.north$lat[3]+.25,
	col="darkgreen", lwd=2,length=0.15, code=3)

arrows(mar.south$lon[1]+1, mar.south$lat[1],
	mar.south$lon[2]-1, mar.south$lat[2],
	col="darkblue", lwd=2,length=0.15, code=3)
arrows(mar.south$lon[2]+.75, mar.south$lat[2]-.25,
	mar.south$lon[3], mar.south$lat[3]-.35,
	col="darkblue", lwd=2,length=0.15, code=3)
arrows(mar.north$lon+.25, mar.north$lat-.25,
	mar.south$lon+.25,mar.south$lat+.25,
	col="darkcyan", lwd=2,length=0.15, code=3)
text(mar.north$lon, mar.north$lat, labels=mar.north$site, 
	col="dark green", cex=1.25)
text(mar.south$lon, mar.south$lat, labels=mar.south$site, 
	col="dark blue", cex=1.25)
dev.off()


###############SUMMARY STATS PER POP######################################

ss<-split(summ.dat, summ.dat$Pop.ID)#create a list of data frames for each pop
y<-split(ss, ss$Locus.ID)
sum.stats<-as.data.frame(cbind(
	unlist(lapply(ss, function(x){length(levels(as.factor(x$Locus.ID)))})),#num RAD loci
	unlist(lapply(ss, nrow)), #num SNPS
	unlist(lapply(ss, function(x){nrow(x[x$Private==1,])})),#num private
	unlist(lapply(ss, function(x){(nrow(x[x$P!=1,])/nrow(x))*100})),#%SNPS polymorphic
	unlist(lapply(ss, function(x){ 
			den<-length(levels(as.factor(x$Locus.ID)))
			num<-length(levels(as.factor(x[x$P!=1,]$Locus.ID)))
			return(num/den*100)
	})),#%RAD loci polymorphic
	unlist(lapply(ss, function(x){mean(x$P)})),#average frequencuy of major allele
	unlist(lapply(ss, function(x){mean(x$Obs.Het)})),#average observed heterozygosity
	unlist(lapply(ss, function(x){mean(x$Pi)})),#average nucleotide diversity
	unlist(lapply(ss, function(x){mean(x$Fis)}))#average Fis
))
colnames(sum.stats)<-c("N.RAD","N.SNP","Private","%Poly.SNP","%Poly.RAD","P","Hobs","Pi","Fis")	
write.table(sum.stats, "E://Docs//PopGen//fst_summary_stats.txt", sep='\t', )

setwd("E://Docs//PopGen")
##Allele frequency spectrum
jpeg("marine_AFS.jpeg", width=480*8, height=480*6, res=300)
par(mfrow=c(3,4), mar=c(2,2,1,1), oma=c(4,4,0,0))
hist(ss$TXSP$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Texas South Padre",, font=2)
hist(ss$TXCC$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Texas Corpus Christi", font=2)
hist(ss$TXCB$P, xlim=c(0.5,1),xlab="",ylab="", main=
	"Texas Christmas Bay", font=2)
hist(ss$ALST$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Alabama Saltwater", font=2)
hist(ss$FLSI$P, xlim=c(0.5,1),xlab="",ylab="", 
	main= "Florida Sanibel Island", font=2)
hist(ss$FLFD$P, xlim=c(0.5,1),xlab="",ylab="", 
	main= "Florida Fort Desoto", font=2)
hist(ss$FLKB$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Florida Keaton Beach", font=2)
hist(ss$FLSG$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Florida St. George", font=2)
hist(ss$FLAB$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Florida Anne's Beach", font=2)
hist(ss$FLPB$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Florida Palm Beach", font=2)
hist(ss$FLHB$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Florida Harbor Branch", font=2)
hist(ss$FLCC$P, xlim=c(0.5,1),xlab="",ylab="", 
	main="Florida Cape Canaveral", font=2)
mtext(side=1, "Allele Frequency", outer = TRUE, line=1)
mtext(side=2, "Frequency", outer=TRUE, line=1)
dev.off()

##FIS
jpeg(file="E:\\Docs\\PopGen\\marine_Fis_hist.jpeg", width=480*8, height=480*6, res=300)
par(mfrow=c(3,4), mar=c(2,2,2,1), oma=c(4,4,0,0))
hist(ss$TXSP$Fis, xlim=c(-1,1),xlab="",ylab="", main="Texas South Padre")
hist(ss$TXCC$Fis, xlim=c(-1,1),xlab="",ylab="", main="Texas Corpus Christi")
hist(ss$TXCB$Fis, xlim=c(-1,1),xlab="",ylab="", main="Texas Christmas Bay")
hist(ss$ALST$Fis, xlim=c(-1,1),xlab="",ylab="", main="Alabama Saltwater")
hist(ss$FLSI$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida Sanibel Island")
hist(ss$FLFD$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida Fort Desoto")
hist(ss$FLKB$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida Keaton Beach")
hist(ss$FLSG$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida St. George")
hist(ss$FLAB$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida Anne's Beach")
hist(ss$FLPB$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida Palm Beach")
hist(ss$FLHB$Fis, xlim=c(-1,1),xlab="",ylab="", main="Florida Harbor Branch")
hist(ss$FLCC$Fis, xlim=c(-1,1),xlab="",ylab="",
	main="Florida Cape Canaveral")
mtext(side=1, "Fis value", outer = TRUE, line=1)
mtext(side=2, "Frequency", outer=TRUE, line=1)
dev.off()


##############################################################################
#****************************************************************************#
#########################TEST FOR ISOLATION BY DISTANCE#######################
#****************************************************************************#
##############################################################################
#Date: 3 April 2015
#Purpose: Analyze fst data
#Mantel test using geographical distances and fsts

fsts<-t(fsts)
TXCC<-rep(NA,nrow(fsts))
fsts<-cbind(fsts,TXCC)

ibd<-mantel.rtest(as.dist(dist),as.dist(fsts))
#Monte-Carlo test
#Observation: 0.675948 
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
#Based on 99 replicates
#Simulated p-value: 0.02 


#########################################################################
#***********************************************************************#
################################ADEGENET#################################
#***********************************************************************#
#########################################################################
#dat.plink<-read.PLINK("E://ubuntushare//stacks//populations//plinkA.raw",
#	parallel=FALSE)
dat.plink<-read.PLINK("E://ubuntushare//stacks//populations//ld.hwe.A.raw",
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
pca1<-glPca(dat.plink, parallel=FALSE)
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
plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,
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
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
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


#########################################################################
#***********************************************************************#
##############################FASTSTRUCTURE##############################
#***********************************************************************#
#########################################################################
#get individual names in the correct order from original structure input file
#str.in<-read.table("E://ubuntushare//stacks//populations//batch_1.structure.str")
str.in<-read.table("E://ubuntushare//stacks//populations//ld.hwe.pgd.str")
inds<-str.in[,1]
inds<-sub('.*_([ATF]\\w+)[_.].*','\\1', inds)
inds<-inds[c(TRUE,FALSE)]
pop.id<-substr(inds,1,4)

#process a structure file
structure.barplot<-function(meanQ.file, k){
	rownames(meanQ.file)<-inds
	meanQ.file<-cbind(meanQ.file,pop.id)
	colnames(meanQ.file)<-c(seq(1,k,1), "pop.id")
	pop.means<-rowsum(meanQ.file[,1:k],
		meanQ.file$pop.id)/summary(meanQ.file$pop.id)
	pop.means.plot<-pop.means[match(pop.list,rownames(pop.means)),]
	filename<-paste("structure.barplot.",k,".jpeg",sep="")
	jpeg(filename,res=300,height=7,width=7, units="in")
	par(mar=c(5,4,4,1),oma=c(2,2,2,1),xpd=TRUE)
	bp<-barplot(as.matrix(t(pop.means.plot)), col=rainbow(k, 0.5), 
		legend=FALSE, axes=FALSE, axisnames=FALSE)
	legend("top", inset=c(0,-0.15), box.lty=0,title="Group",
		legend=seq(1,k,1),
		pch=15, col=rainbow(k,0.5), horiz=TRUE)
	axis(2,las=1, pos=0)
	axis(1,labels=FALSE, at=bp,pos=0)
	text(bp, labels=pop.list,srt=30, xpd=TRUE, pos=1,
		par("usr")[3])
	dev.off()
}
#K=2
stru2<-read.table("E://ubuntushare//ld.hwe_out_simple.2.meanQ",header=F)
structure.barplot(stru2,2)

#K=3
stru3<-read.table("E://ubuntushare//ld.hwe_out_simple.3.meanQ",header=F)
structure.barplot(stru3,3)

#K=4
stru4<-read.table("E://ubuntushare//ld.hwe_out_simple.4.meanQ",header=F)
structure.barplot(stru4,4)

#K=5
stru5<-read.table("E://ubuntushare//ld.hwe_out_simple.5.meanQ",header=F)
structure.barplot(stru5,5)

#K=6
stru6<-read.table("E://ubuntushare//ld.hwe_out_simple.6.meanQ",header=F)
structure.barplot(stru6,6)

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

####RUN STRUCTURE ON ONLY ONE POP..TXSP is largest
subset.txsp<-inds[pop.id=="TXSP"]
write.table(subset.txsp,"E://ubuntushare//stacks//populations//txsp.indlist", 
	quote=FALSE, col.names=F, row.names=F, eol="\n", sep="\t")


#########################################################################
#***********************************************************************#
################################BAYENV2##################################
#***********************************************************************#
#########################################################################

#**************************STARTING WITH PLINK FILES********************#
ped<-read.table("B:\\ubuntushare\\current_analysis\\ld.hwe.sub.ped", 
	skip = 1, stringsAsFactors=F, colClasses="character")
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

write.table(ped,"bayenv.plink.ped", 
	row.names=F, col.names=F, quote=F, sep="\t",eol="\n")

clust.plink<-cbind(ped.pops, ped[,2],ped.pops)
write.table(clust.plink, 
	"E://ubuntushare//stacks//populations//plink.clust.txt",
	col.names=F, row.names=F, quote=F, sep="\t", eol="\n")


#plink.map -> numbers instead of scaffold_#
map<-read.table("ld.hwe.sub.map", skip = 1)
chr.nums<-sub('scaffold_','',map[,1])
map[,1]<-chr.nums

write.table(map, "bayenv.plink.map", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n")

##get any snps *not* in plink.ld
plink.ld<-read.table("E://ubuntushare//stacks//populations//plink.ld",
	header=T)
all.snps<-c(as.character(plink.ld$SNP_A), as.character(plink.ld$SNP_B))
above.2<-map[map$V2 %in% all.snps,]
below.2<-map[!(map$V2 %in% all.snps),]

write.table(below.2$V2, "E://ubuntushare//stacks//populations//ld.subset.list",
	col.names=F, row.names=F, quote=F, sep='\t', eol='\n')

#***************************CONVERT PLINK TO BAYENV2************************#
freq<-read.table("ld.hwe.bayenv.plink.frq.strat", 
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

write.table(snpsfile, "null.1833", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n") #bayenv SNPSFILE

#*******************************check Bayenv2 matrix*************************#
matrix.files<-list.files("B://ubuntushare//environmental_assoc//",pattern="matrix")
matrices<-list()
for(i in 1:length(matrix.files))
{
	matrices[[i]]<-as.matrix(read.table(
		paste("B://ubuntushare//environmental_assoc//",matrix.files[i],sep=""), 
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

#************************************SNPFILEs******************************#
#for SNPFILE, need just one file per SNP apparently.
#want to use all of the snps...need to get map with those inds.
all.snps.ped<-read.table("batch_1.plink.ped", header=F, stringsAsFactors=F)
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

#########################################################################
#***********************************************************************#
################################PCADAPT##################################
#***********************************************************************#
#########################################################################
#K=4 WAS BEST
scores.files<-list.files("E://ubuntushare//pcadapt//", pattern="scores")
loadings.files<-list.files("E://ubuntushare//pcadapt//", pattern="loadings")
snps.files<-list.files("E://ubuntushare//pcadapt//", pattern="topBF")

snp.list<-list()
for(i in 1: length(snps.files)){
	#read in files
	snp.list[[i]]<-read.table(
		paste("E://ubuntushare//pcadapt//", snps.files[i], sep=""), 
		header=T)
}
#compare lists of snps in all of the runs

all.snps<-as.vector(sapply(snp.list, "[[","snp"))
all.snps.dup<-all.snps[duplicated(all.snps)]
rep.snps<-all.snps.dup[!duplicated(all.snps.dup)]

scores<-read.table("E://ubuntushare//pcadapt//pcadapt.4.scores")
loading<-read.table("E://ubuntushare//pcadapt//pcadapt.4", header=T, sep="\t")
#bf.log<-log10(snp.list[[1]]$BF)
#loading[round(loading$logBF, 4) %in% round(bf.log),] #doesn't work..

#the snp is the row number in the map file for the snps
ld.map<-read.table("E://ubuntushare//stacks//populations//ld.subset.map",
	header=F)
pcadapt.outliers<-ld.map[rep.snps,]
pa.out.radloc<-sub('(\\d+)_\\d+','\\1',pcadapt.outliers$V2)

summ.dat<-read.table("E://ubuntushare//stacks//populations//batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
pa.out.dat<-summ.dat[summ.dat$Locus.ID %in% pa.out.radloc,]
pa.out.dat$Chr<-factor(pa.out.dat$Chr)
#are any on the linkage map?
length(levels(as.factor(use.contigs[use.contigs$Scaffold %in% pa.out.dat$Chr,1])))

#plot individual scores
jpeg("pcadapt.scores1.2.jpeg", height=12, width=12, units="in", res=300)
plot(as.numeric(scores[1,]),as.numeric(scores[2,]),pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="")
legend("bottomright", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0117)", 1, line = 2)
mtext("PC2 (rho2: 0.0091)",2, line = 2)
dev.off()

jpeg("pcadapt.scores1.3.jpeg", height=12, width=12, units="in", res=300)
plot(as.numeric(scores[1,]),as.numeric(scores[3,]),pch=16, cex=2,
	col=alpha(col, 0.5), ylab="", xlab="")
legend("topright", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0117)", 1, line = 2)
mtext("PC3 (rho2: 0.0018)",2, line = 2)
dev.off()

scores5.files<-list.files("E://ubuntushare//pcadapt//k5//", pattern="scores")
loadings5.files<-list.files("E://ubuntushare//pcadapt//k5//", pattern="loadings")
snps5.files<-list.files("E://ubuntushare//pcadapt//k5//", pattern="topBF")

snp5.list<-list()
for(i in 1: length(snps5.files)){
	#read in files
	snp5.list[[i]]<-read.table(
		paste("E://ubuntushare//pcadapt//k5//", snps5.files[i], sep=""), 
		header=T)
}
#compare lists of snps in all of the runs

all5.snps<-as.vector(sapply(snp5.list, "[[","snp"))
all5.snps.dup<-all5.snps[duplicated(all5.snps)]
rep5.snps<-all5.snps.dup[!duplicated(all5.snps.dup)]
#recover locus IDs from plink map
ld.hwe.map<-read.table("E://ubuntushare//stacks//populations//ld.hwe.sub.map")
pcadapt5.outliers<-ld.hwe.map[rep5.snps,]
pa5.out.radloc<-sub('(\\d+)_\\d+','\\1',pcadapt5.outliers$V2)
pa5.out.dat<-summ.dat[summ.dat$Locus.ID %in% pa5.out.radloc,]
pa5.out.dat$Chr<-factor(pa5.out.dat$Chr)

#retrieve bayes factors 
all5.bf<-as.data.frame(cbind(as.vector(sapply(snp5.list, "[[","snp")),
	as.vector(sapply(snp5.list, "[[","BF"))))
colnames(all5.bf)<-c("snp", "BF")
rep5.bf<-all5.bf[all5.bf$snp %in% rep5.snps,]
rep5.bf<-cbind(rep5.bf, sub('(\\d+)(_\\d+)','\\1',ld.hwe.map[rep5.bf$snp,2]))
colnames(rep5.bf)[3]<-"Locus"


##############################LOSITAN######################################
lositan.loci<-read.delim("E://ubuntushare//ld.hwe.lositan.loci", header=T)
lositan.ci<-read.delim("E://ubuntushare//ld.hwe.lositan.ci", header=T, sep='\t')

plot_ci <- function(df.ci, bcolor, mcolor, tcolor) {
  lines(df.ci[,1],df.ci[,2], type='l', col=bcolor)#low
  lines(df.ci[,1],df.ci[,3], type='l', col=mcolor, lty=3)#med
  lines(df.ci[,1],df.ci[,4], type='l', col=tcolor)#upp
}

plot_loci <- function(df.loc, color) {
  points(df.loc[,2],df.loc[,3], col=color, pch=19)
}

lositan.loci<-cbind(lositan.loci, sub('(\\d+)(_\\d+)','\\1',lositan.loci[,1])
colnames(lositan.loci)[5]<-"RAD.loc"

lositan.sig<-lositan.loci[(lositan.loci[,4] <= 0.05) | 
	(lositan.loci[,4]>=0.95),]
#compare to outlier fst analysis
lositan.allfsts<-lositan.sig[lositan.sig$RAD.loc 
	%in% no.dups.outliers$Locus,]

#compare to pcadapt with k=5
lositan.pcadapt<-lositan.sig[lositan.sig$RAD.loc %in% pa5.out.dat$Locus.ID,]

jpeg("Fst_het_lositan.jpeg", res=300, height=7, width=7, units="in")
plot(-10,ylim=c(0,0.25),xlim=c(0,1), xlab='He', ylab='Fst', las=1)
plot_ci(lositan.ci, 'black', 'gray', 'black')
plot_loci(lositan.loci, 'black')
points(lositan.sig$Het, lositan.sig$Fst, cex=1, col="red", pch=8)
points(lositan.allfsts$Het, lositan.allfsts$Fst, cex=3, col="green")
points(lositan.pcadapt$Het, lositan.pcadapt$Fst, cex=3, col="purple", pch=6)
legend(x=0.6, y=0.25,bty="n",
	c("99% Quantile", "Median"),lty=c(1,3), col=c("black", "gray"))
legend(x=0.65, y=0.23, bty="n",
	c("Lositan", "Lositan & Fst Out", "Lositan & PCAdapt"),
	pch=c(8,1,6), col=c("red", "green", "purple"))

dev.off()

sig.in.three<-lositan.pcadapt[lositan.pcadapt$RAD.loc 
	%in% lositan.allfsts$RAD.loc,]
sig.in3.summ<-pa5.out.dat[pa5.out.dat$Locus.ID %in% sig.in.three$RAD.loc,]
sig.in3.global<-global.fsts[global.fsts$Locus %in% sig.in.three$RAD.loc,]

all.fst.files<-list.files("E://ubuntushare//stacks//populations",
	pattern="fst_[A-Z].*[A-Z].txt")
sig.in3.stacks<-data.frame()
for(i in 1:length(all.fst.files))
{
	fst.file<-read.delim(paste("E://ubuntushare//stacks//populations//",
		all.fst.files[i], sep=""), header=T)
	sig.in3.stacks<-rbind(sig.in3.stacks,
		fst.file[fst.file$Locus %in% sig.in.three$RAD.loc,])
}
sig.in3.stacks.split<-split(sig.in3.stacks, sig.in3.stacks$Locus)
sig.in3.stacks.avg<-lapply(sig.in3.stacks.split, function(x) 
	{ mean(x$SmoothFst) })

sig.in3.stacks.het<-summ.dat[summ.dat$Locus.ID %in% sig.in.three$RAD.loc,]
sig.in3.stacks.h.split<-split(sig.in3.stacks.het, sig.in3.stacks.het$Locus)
sig.in3.stacks.h.avg<-lapply(sig.in3.stacks.h.split, function(x) 
	{ avgs<-cbind(mean(as.numeric(x$Obs.Het)), 
		mean(as.numeric(x$Exp.Het)))
	return(avgs) })

sig.in3.bf<-rep5.bf[rep5.bf$Locus %in% sig.in.three$RAD.loc,]
sig.in3.bf.split<-split(sig.in3.bf, factor(sig.in3.bf$Locus))
sig.in3.bf.avg<-lapply(sig.in3.bf.split, function(x){ mean(x$BF) })

##############################################################################
##############################PLOT ALL SCAFFOLDS##############################
##############################################################################
#####PREP DATA FOR PLOTTING#####
all.scaff<-split(summ.dat, summ.dat$Chr)
plot.scaffs<-lapply(split(summ.dat, summ.dat$Pop.ID), function(x){split(x, x$Chr)})
last.max<-0
rect.xs<-NULL
addition.values<-0
for(i in 1:length(all.scaff)){
	new.max<-last.max+round(max(all.scaff[[i]]$BP), -2)
	rect.xs<-rbind(rect.xs,c(last.max, new.max))
	addition.values<-c(addition.values, new.max)
	last.max<-new.max
}

for(i in 1:length(plot.scaffs)){
	for(j in 1:length(plot.scaffs[[i]])){
		plot.scaffs[[i]][[j]]$BP<-plot.scaffs[[i]][[j]]$BP+addition.values[j]
	}
}

x.max<-max(addition.values)
y.max<-max(summ.dat$Smoothed.Pi)+0.1*max(summ.dat$Smoothed.Pi)
if(min(summ.dat$Smoothed.Pi) < 0) {
	y.min<-min(summ.dat$Smoothed.Pi) - 0.1*min(summ.dat$Smoothed.Pi)
} else {
	y.min<-0
}
pop.names<-c("Alabama Saltwater", "Florida Anne's Beach", 
	"Florida Cape Canaveral", "Florida Fort Desoto","Florida Harbor Branch",
	"Florida Keaton Beach","Florida Palatka Bay", "Florida St. George",
	"Florida Sanibel Island", "Texas Christmas Bay", "Texas Corpus Christi",
	"Texas South Padre")
############################PLOT PI###########################################
##ALL LOCI
jpeg("marine_all_pi.jpeg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,10,1,0), oma=c(4,10,0,0))

for(i in 1:length(plot.scaffs)){
	plot(plot.scaffs[[i]][[1]]$BP, plot.scaffs[[i]][[1]]$Smoothed.Pi, 
		xlim=c(0,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(plot.scaffs[[i]])){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(plot.scaffs[[i]][[j]]$BP, 
			plot.scaffs[[i]][[j]]$Smoothed.Pi,
			y.max,x.max, rect.xs[j,],y.min=0,x.min=0, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
	mtext(pop.names[i],2,las=1, line=2)
}
mtext(side=1, "Genomic Location", outer = TRUE, line=1)
mtext(side=2, "Smoothed pi value for each population", outer=TRUE, line=6)
dev.off()



y.max<-max(sig.scaff$Smoothed.Pi)+0.1*max(sig.scaff$Smoothed.Pi)
x.max<-1481000
rect.pts<-rbind(c(50000,134000), c(1370000,x.max))
jpeg("marine_scaff_pi.jpeg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,8,1,0), oma=c(4,8,0,0))

plot.genome.wide(scaff.ss$TXSP$BP,scaff.ss$TXSP$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Texas South Padre",2,las=1, line=1)

plot.genome.wide(scaff.ss$TXCC$BP,scaff.ss$TXCC$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Texas Corpus Christi",2,las=1, line=1)

plot.genome.wide(scaff.ss$TXCB$BP,scaff.ss$TXCB$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Texas Christmas Bay",2,las=1, line=1)

plot.genome.wide(scaff.ss$ALST$BP,scaff.ss$ALST$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Alabama Salt Water",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLSI$BP,scaff.ss$FLSI$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Sanibel Island",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLFD$BP,scaff.ss$FLFD$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Fort Desoto",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLKB$BP,scaff.ss$FLKB$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Keaton Beach",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLSG$BP,scaff.ss$FLSG$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida St. George",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLAB$BP,scaff.ss$FLAB$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Anne's Beach",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLPB$BP,scaff.ss$FLPB$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Palm Beach",2,las=1, line=1)

plot.genome.wide(scaff.ss$FLHB$BP,scaff.ss$FLHB$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Harbor Branch",2,las=1, line=1)

par(mar=c(4,8,1,0))
plot.genome.wide(scaff.ss$FLCC$BP,scaff.ss$FLCC$Smoothed.Pi,
	y.max,x.max,rect.pts)
mtext("Florida Cape Canaveral",2,las=1, line=1)

axis(1, at=c(0,25000, 92000, 752000, 1425500,1481000), xlim=c(0,1481000), 
	ylim=c(0,y.max), las=1, tck=-0.01, labels=FALSE, pos=0, xpd=TRUE)
text(c(25000, 92000, 752000, 1425500), 
	labels=c("scaffold_1189","scaffold_468","scaffold_7","scaffold_854"),
	par("usr")[3]-0.00051,srt = 38, adj = 1, xpd = TRUE,cex=0.75)

mtext(side=1, "Genomic Location", outer = TRUE, line=1)
mtext(side=2, "Smoothed pi value for each population", outer=TRUE, line=6)
dev.off()

############################PLOT FIS###########################################
##ALL LOCI
y.max<-max(summ.dat$Smoothed.Fis)+0.1*max(summ.dat$Smoothed.Fis)
if(min(summ.dat$Smoothed.Fis) < 0) {
	y.min<-min(summ.dat$Smoothed.Fis) + 0.1*min(summ.dat$Smoothed.Fis)
} else {
	y.min<-0
}

jpeg("marine_all_fis.jpeg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,10,1,0), oma=c(4,10,0,0))

for(i in 1:length(plot.scaffs)){
	plot(plot.scaffs[[i]][[1]]$BP, plot.scaffs[[i]][[1]]$Smoothed.Fis, 
		xlim=c(0,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(plot.scaffs[[i]])){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(plot.scaffs[[i]][[j]]$BP, 
			plot.scaffs[[i]][[j]]$Smoothed.Fis,
			y.max,x.max, rect.xs[j,],y.min,x.min=0, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
	mtext(pop.names[i],2,las=1, line=2)
}
mtext(side=1, "Genomic Location", outer = TRUE, line=1)
mtext(side=2, "Smoothed Fis value for each population", outer=TRUE, line=6)
dev.off()


y.max<-max(sig.scaff$Smoothed.Fis)+0.1*max(sig.scaff$Smoothed.Fis)
x.max<-1481000
y.min<-min(sig.scaff$Smoothed.Fis)-0.1*min(sig.scaff$Smoothed.Fis)

jpeg("marine_scaff_fis.jpg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,10,1,0), oma=c(4,10,0,0))

plot.genome.wide(scaff.ss$TXSP$BP,scaff.ss$TXSP$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Texas South Padre",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$TXCC$BP,scaff.ss$TXCC$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Texas Corpus Christi",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$TXCB$BP,scaff.ss$TXCB$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Texas Christmas Bay",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$ALST$BP,scaff.ss$ALST$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Alabama Salt Water",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLSI$BP,scaff.ss$FLSI$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Sanibel Island",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLFD$BP,scaff.ss$FLFD$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Fort Desoto",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLKB$BP,scaff.ss$FLKB$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Keaton Beach",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLSG$BP,scaff.ss$FLSG$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida St. George",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLAB$BP,scaff.ss$FLAB$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Anne's Beach",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLPB$BP,scaff.ss$FLPB$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Palm Beach",2,las=1, line=3,cex=1)

plot.genome.wide(scaff.ss$FLHB$BP,scaff.ss$FLHB$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Harbor Branch",2,las=1, line=3,cex=1)

par(mar=c(4,10,1,0))
plot.genome.wide(scaff.ss$FLCC$BP,scaff.ss$FLCC$Smoothed.Fis,
	y.max,x.max,rect.pts, y.min=y.min)
mtext("Florida Cape Canaveral",2,las=1, line=3,cex=1)

axis(1, at=c(0,25000, 92000, 752000, 1425500,1481000), xlim=c(0,1481000), 
	ylim=c(y.min,y.max), las=1, tck=-0.01, labels=FALSE, pos=0, xpd=TRUE)
text(c(25000, 92000, 752000, 1425500), 
	labels=c("scaffold_1189","scaffold_468","scaffold_7","scaffold_854"),
	par("usr")[3]-0.00051,srt = 38, adj = 1, xpd = TRUE,cex=0.75)

mtext(side=1, "Genomic Location", outer = TRUE, line=1)
mtext(side=2, "Smoothed Fis value for each population", outer=TRUE, line=7)
dev.off()


##############################################################################
*****************************************************************************#
#################################ALL OUTLIERS#################################
#****************************************************************************#
##############################################################################
all.fst.files<-list.files("E://ubuntushare//stacks//populations",
	pattern="fst_[A-Z].*[A-Z].txt")
all.summ.files<-list.files("E://ubuntushare//stacks//populations",
	pattern="fst_[A-Z].*[A-Z]_summary.txt")
all.outliers<-data.frame()
for(i in 1:length(all.fst.files))
{
	fst.file<-read.delim(paste("E://ubuntushare//stacks//populations//",
		all.fst.files[i], sep=""), header=T)
	summ.file<-read.delim(paste("E://ubuntushare//stacks//populations//",
		all.summ.files[i], sep=""), header=T)
	all.outliers<-rbind(all.outliers,
		fst.file[fst.file$SmoothFst >= summ.file$CI99smooth,])
}

no.dups.outliers<-all.outliers[!duplicated(all.outliers$Locus),]
length(levels(factor(no.dups.outliers$Chrom))) #number scaffolds
length(levels(as.factor(
	contigs[contigs$Scaffold %in% no.dups.outliers$Chrom,1]))) #num lgs


#select those also found in non-outlier analysis
length(levels(as.factor(no.dups.outliers[no.dups.outliers$Locus 
	%in% pa.out.dat$Locus.ID,1])))
#summarize fsts for those (updated 26 May 2015)
summary(all.outliers$SmoothFst[all.outliers$Locus==3149])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08358 0.10720 0.11510 0.10720 0.11510 0.11510 
summary(all.outliers$SmoothFst[all.outliers$Locus==32348])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3460  0.3534  0.3608  0.3608  0.3683  0.3757  
summary(all.outliers$SmoothFst[all.outliers$Locus==37147])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.1231  0.1244  0.1258  0.1258  0.1272  0.1285
summary(all.outliers$SmoothFst[all.outliers$Locus==29715])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06973 0.08017 0.09062 0.09062 0.10110 0.11150 
summary(all.outliers$SmoothFst[all.outliers$Locus==15712])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.2391  0.2498  0.2605  0.2605  0.2712  0.2819 




##############################################################################
#****************************************************************************#
##########################FST-HETEROZYGOSITY PLOT#############################
#****************************************************************************#
##############################################################################
#Import global values generated by C++ program calculate_global_fsts
global.fsts<-read.table(
	"E://ubuntushare//stacks//populations//ld.hwe.sub.globalstats.txt",
	header=T, sep='\t')

fst.ci.99<-mean(global.fsts$Fst) + 2.57583*sd(global.fsts$Fst)
#fst.het.lm<-lm(global.fsts$Fst~global.fsts$Obs.Het)
#fst.het.pr<-predict(fst.het.lm, interval="confidence", level=0.99)

jpeg("fst_het.jpeg", res=300, width=7,height=7, units="in")
plot(global.fsts$Obs.Het, global.fsts$Fst, xlab="", ylab="", 
	xaxt='n', yaxt='n', pch=19)
abline(h=fst.ci.99, lwd=2)
axis(1)
axis(2, las=2)
mtext("Observed Heterozygosity", 1, outer=FALSE, line =3)
mtext("Fst", 2, outer=FALSE, line=3)
dev.off()

#what are the outliers??
glob.fst.out<-global.fsts[global.fsts$Fst >= fst.ci.99,]
 summ.gfo<-summ.dat[summ.dat$Locus.ID %in% glob.fst.out$Locus,]
 summ.gfo<-summ.gfo[!duplicated(summ.gfo$Locus),]
glob.fst.out<-
	merge(summ.gfo[,2:3], glob.fst.out, by.x="Locus.ID", by.y="Locus")
glob.fst.out<-glob.fst.out[,-3]

glob.fst.out[glob.fst.out$Locus %in% no.dups.outliers$Locus,]
##USE LOSITAN ANALYSIS
##############################################################################
#****************************************************************************#
##########################FASTSTRUCTURE FST OUTLIERS##########################
#****************************************************************************#
##############################################################################
#FIRST: need to run genomewide_CIs C++ code
all.fst.files<-list.files("E://ubuntushare//stacks//populations_fstrugroups",
	pattern="fst_[0-9]-[0-9].txt")
all.summ.files<-list.files("E://ubuntushare//stacks//populations_fstrugroups",
	pattern="fst_[0-9]-[0-9]_summary.txt")
all.grp.outliers<-data.frame()
for(i in 1:length(all.fst.files))
{
	fst.file<-read.delim(paste(
		"E://ubuntushare//stacks//populations_fstrugroups//",
		all.fst.files[i], sep=""), header=T)
	summ.file<-read.delim(paste(
		"E://ubuntushare//stacks//populations_fstrugroups//",
		all.summ.files[i], sep=""), header=T)
	all.grp.outliers<-rbind(all.grp.outliers,
		fst.file[fst.file$SmoothFst >= summ.file$CI99smooth,])
}

grp.outliers<-all.grp.outliers[!duplicated(all.grp.outliers$Locus),]
length(levels(factor(grp.outliers$Chrom))) #number scaffolds
length(levels(as.factor(
	contigs[contigs$Scaffold %in% grp.outliers$Chrom,1]))) #num lgs

setwd("E://ubuntushare//stacks//populations_fstrugroups")
g1.2<-read.delim("fst_1-2.txt", header=T)
g1.3<-read.delim("fst_1-3.txt", header=T)
g1.4<-read.delim("fst_1-4.txt", header=T)
g1.5<-read.delim("fst_1-5.txt", header=T)
g3.2<-read.delim("fst_3-2.txt", header=T)
g3.4<-read.delim("fst_3-4.txt", header=T)
g4.2<-read.delim("fst_4-2.txt", header=T)
g5.2<-read.delim("fst_5-2.txt", header=T)
g5.3<-read.delim("fst_5-3.txt", header=T)
g5.4<-read.delim("fst_5-4.txt", header=T)

g1.2.ci<-read.delim("fst_1-2_summary.txt", header=T)
g1.3.ci<-read.delim("fst_1-3_summary.txt", header=T)
g1.4.ci<-read.delim("fst_1-4_summary.txt", header=T)
g1.5.ci<-read.delim("fst_1-5_summary.txt", header=T)
g3.2.ci<-read.delim("fst_3-2_summary.txt", header=T)
g3.4.ci<-read.delim("fst_3-4_summary.txt", header=T)
g4.2.ci<-read.delim("fst_4-2_summary.txt", header=T)
g5.2.ci<-read.delim("fst_5-2_summary.txt", header=T)
g5.3.ci<-read.delim("fst_5-3_summary.txt", header=T)
g5.4.ci<-read.delim("fst_5-4_summary.txt", header=T)

jpeg("E://Docs//PopGen//groups_fsts.jpeg", width=480*6, height=480*9, res=300)
par(mfrow=c(5,2), mar=c(0,12,1,0), oma=c(4,10,0,0))
plot.fsts.scaffs(g1.2, "Group 1-Group 2", g1.2.ci)
plot.fsts.scaffs(g1.3, "Group 1-Group 3", g1.3.ci)
plot.fsts.scaffs(g1.4, "Group 1-Group 4", g1.4.ci)
plot.fsts.scaffs(g1.5, "Group 1-Group 5", g1.5.ci)
plot.fsts.scaffs(g3.2, "Group 3-Group 2", g3.2.ci)
plot.fsts.scaffs(g3.4, "Group 3-Group 4", g3.4.ci)
plot.fsts.scaffs(g4.2, "Group 4-Group 2", g4.2.ci)
plot.fsts.scaffs(g5.2, "Group 5-Group 2", g5.2.ci)
plot.fsts.scaffs(g5.3, "Group 5-Group 3", g5.3.ci)
plot.fsts.scaffs(g5.4, "Group 5-Group 4", g5.4.ci)

mtext(side=1, "Genomic Location", outer = TRUE, line=1)
mtext(side=2, "Smoothed FST value", outer=TRUE, line=6)

dev.off()

#are any of these the same as the original pop gen analysis?
length(levels(as.factor(grp.outliers[grp.outliers$Locus 
	%in% no.dups.outliers$Locus,1])))

##############################################################################
#****************************************************************************#
########################LESS STRINGENT POPULATIONS FSTS#######################
#****************************************************************************#
##############################################################################
#FIRST: need to run genomewide_CIs C++ code
nolim.fst.files<-list.files("E://ubuntushare//stacks//populations_nopopslimit",
	pattern="fst_[A-Z].*[A-Z].txt")
nolim.summ.files<-list.files("E://ubuntushare//stacks//populations_nopopslimit",
	pattern="fst_[A-Z].*[A-Z]_summary.txt")
all.nolim.outliers<-data.frame()
for(i in 1:length(nolim.fst.files))
{
	fst.file<-read.delim(paste(
		"E://ubuntushare//stacks//populations_nopopslimit//",
		nolim.fst.files[i], sep=""), header=T)
	summ.file<-read.delim(paste(
		"E://ubuntushare//stacks//populations_nopopslimit//",
		nolim.summ.files[i], sep=""), header=T)
	all.nolim.outliers<-rbind(all.nolim.outliers,
		fst.file[fst.file$SmoothFst >= summ.file$CI99smooth,])
}

nolim.outliers<-all.nolim.outliers[!duplicated(all.nolim.outliers$Locus),]
length(levels(factor(nolim.outliers$Chrom))) #number scaffolds
length(levels(as.factor(
	contigs[contigs$Scaffold %in% nolim.outliers$Chrom,1]))) #num lgs
length(levels(as.factor(
	nolim.outliers[nolim.outliers$Locus %in% grp.outliers$Locus,1])))
length(levels(as.factor(
	nolim.outliers[nolim.outliers$Locus %in% no.dups.outliers$Locus,1])))
dim(nolim.outliers)

nolim.summ<-read.table("E://ubuntushare//stacks//populations_nopopslimit//batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
ss<-split(nolim.summ, nolim.summ$Pop.ID)#create a list of data frames for each pop
nolim.sum.stats<-as.data.frame(cbind(
	unlist(lapply(ss, function(x){length(levels(as.factor(x$Locus.ID)))})),#num RAD loci
	unlist(lapply(ss, nrow)), #num SNPS
	unlist(lapply(ss, function(x){nrow(x[x$Private==1,])})),#num private
	unlist(lapply(ss, function(x){(nrow(x[x$P!=1,])/nrow(x))*100})),#%SNPS polymorphic
	unlist(lapply(ss, function(x){ 
			den<-length(levels(as.factor(x$Locus.ID)))
			num<-length(levels(as.factor(x[x$P!=1,]$Locus.ID)))
			return(num/den*100)
	})),#%RAD loci polymorphic
	unlist(lapply(ss, function(x){mean(x$P)})),#average frequencuy of major allele
	unlist(lapply(ss, function(x){mean(x$Obs.Het)})),#average observed heterozygosity
	unlist(lapply(ss, function(x){mean(x$Pi)})),#average nucleotide diversity
	unlist(lapply(ss, function(x){mean(x$Fis)}))#average Fis
))
colnames(nolim.sum.stats)<-c("N.RAD","N.SNP","Private","%Poly.SNP","%Poly.RAD","P","Hobs","Pi","Fis")	
write.table(sum.stats, "E://Docs//PopGen//nolim_fst_summary_stats.txt", sep='\t', )



#############################################################################
############################MALE-FEMALE LOCI#################################	
#############################################################################


m.het1<-m.f.summ.dat[which(m.f.summ.dat$Pop.ID == "male" &
	m.f.summ.dat$Obs.Het == 1),]
m.het1.f<-m.f.summ.dat[m.f.summ.dat$Locus.ID %in% m.het1$Locus.ID,]#recover females
f.het1<-m.f.summ.dat[which(m.f.summ.dat$Pop.ID == "female" &
	m.f.summ.dat$Obs.Het == 1),]
f.het1.m<-m.f.summ.dat[m.f.summ.dat$Locus.ID %in% f.het1$Locus.ID,]

f1.scaff<-f.het1.m[f.het1.m$Chr %in% use.contigs$Scaffold,]
m1.scaff<-m.het1.f[m.het1.f$Chr %in% use.contigs$Scaffold,]
lg1<-NULL
lg2<-NULL
for(i in 1:nrow(f1.scaff)){
	id<-levels(as.factor(
		use.contigs[use.contigs$Scaffold %in% factor(f1.scaff$Chr[i]),1]))
	lg1<-c(lg1, id)
}
for(i in 1:nrow(m1.scaff)){
	id<-levels(as.factor(
		use.contigs[use.contigs$Scaffold %in% factor(m1.scaff$Chr[i]),1]))
	lg2<-c(lg2, id)
}
f1.scaff<-cbind(f1.scaff,lg1)
m1.scaff<-cbind(m1.scaff,lg2)
f1.scaff<-droplevels(f1.scaff)
m1.scaff<-droplevels(m1.scaff)

f1.m0.all<-f.het1.m[f.het1.m$Locus.ID %in% 
	f.het1.m[which(f.het1.m$Pop.ID == "male" &
	f.het1.m$Obs.Het == 0),]$Locus.ID,]
f1.m0<-f1.scaff[f1.scaff$Locus.ID %in% 
	f1.scaff[which(f1.scaff$Pop.ID == "male" &
	 f1.scaff$Obs.Het == 0),]$Locus.ID,]#problem is some snps don't reflect RAD locus
f1.m0.split<-split(f1.m0, f1.m0$Locus.ID)

f0.m1.all<-m.het1.f[m.het1.f$Locus.ID %in% 
	m.het1.f[which(m.het1.f$Pop.ID == "female" &
	 m.het1.f$Obs.Het == 0),]$Locus.ID,]
f0.m1<-m1.scaff[m1.scaff$Locus.ID %in% 
	m1.scaff[which(m1.scaff$Pop.ID == "female" &
	 m1.scaff$Obs.Het == 0),]$Locus.ID,]#problem is some snps don't reflect RAD locus
f0.m1.split<-split(f0.m1, f0.m1$Locus.ID)

write.table(f1.m0.all, "f.het1.m.het0.txt", 
	col.names=TRUE, row.names=FALSE, sep="\t")
write.table(f1.m0, "f.het1.m.het0.lg.txt",  
	col.names=TRUE, row.names=FALSE, sep="\t")

write.table(f0.m1.all, "m.het1.f.het0.txt",  
	col.names=TRUE, row.names=FALSE, sep="\t")

write.table(f0.m1, "m.het1.f.het0.lg.txt",  
	col.names=TRUE, row.names=FALSE, sep="\t")


#######Returned to this: 4 June 2015
m.f.summ<-read.table("E://ubuntushare//stacks//populations_sex//batch_1.sumstats.tsv",
	sep='\t', skip=2, header=T, comment.char="")
m.f.fst<-read.delim("E://ubuntushare//stacks//populations_sex//batch_1.fst_female-male.tsv")
mf.fst.sum<-read.delim("E://ubuntushare//stacks//populations_sex//fst_female-male.txt")
mf.fst.ci<-read.delim("E://ubuntushare//stacks//populations_sex//fst_female-male_summary.txt")

sex.outliers<-mf.fst.sum[mf.fst.sum$SmoothFst >= mf.fst.ci$CI99smooth,]


sex.outliers<-sex.outliers[!(sex.outliers$Locus %in% no.dups.outliers$Locus),]

length(levels(factor(sex.outliers$Locus)))
length(levels(factor(sex.outliers$Chrom)))
length(levels(as.factor(
	contigs[contigs$Scaffold %in% sex.outliers$Chrom,1]))) #num lgs

jpeg("E://Docs//PopGen//male-female.jpeg", res=300, height=480*2, width=480*4)
par(mar=c(2,4,2,1), oma=c(2,4,2,1), las=1)
plot.fsts.scaffs(mf.fst.sum, "", mf.fst.ci)
mtext(side=1,"Genomic Location", outer = FALSE, cex=0.75)
mtext("Smoothed FST value", side=2, outer=FALSE, las=0, line=4, cex=0.75)
dev.off()

#do any of the outliers have unusual heterozygosity patterns?
out.summ<-m.f.summ[m.f.summ$Locus.ID %in% sex.outliers$Locus,]
out.summ.split<-split(out.summ, out.summ$Pop.ID)

het1.loc<-m.f.summ[(m.f.summ$Obs.Het==1),"Locus.ID"]
het1.summ<-m.f.summ[m.f.summ$Locus.ID %in% het1.loc,]
put.het1.loc<-het1.summ[het1.summ$Obs.Het <= 0.5,"Locus.ID"]
put.het1<-het1.summ[het1.summ$Locus.ID %in% put.het1.loc,]
put.het1.split<-split(put.het1, put.het1$Locus.ID)

het0.loc<-m.f.summ[(m.f.summ$Obs.Het==0),"Locus.ID"]
het0.summ<-m.f.summ[m.f.summ$Locus.ID %in% het0.loc,]
put.het0.loc<-het0.summ[het0.summ$Obs.Het >= 0.5,"Locus.ID"]
put.het0<-het1.summ[het0.summ$Locus.ID %in% put.het0.loc,]
put.het0.split<-split(put.het0, put.het0$Locus.ID)

###########Using PLINK's assoc results (treat male/female as case/control)
sex.ped<-as.data.frame(read.delim(
	"E://ubuntushare//stacks//populations_sex//batch_1.plink.ped",
	stringsAsFactors=F, header=F,colClasses="character"))
ped.names<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', sex.ped[,2])
ped.names<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', ped.names)
ped.sex<-sub('sample_\\w{4}(\\w+).*[_.].*','\\1', sex.ped[,2])
ped.sex[substr(ped.sex,1,2)=="DP"]<-"M"
ped.sex[substr(ped.sex,1,2)=="DB"]<-"F"
ped.sex[substr(ped.sex,1,2)=="NP"]<-"M"
ped.sex[substr(ped.sex,1,1)=="P"]<-"M"
ped.sex<-substr(ped.sex,1,1)
ped.sex[ped.sex=="F"]<-2
ped.sex[ped.sex=="M"]<-1
ped.new<-sex.ped
ped.new[,2]<-ped.names
ped.new[,5]<-ped.sex
ped.new[,6]<-ped.sex
write.table(ped.new, "E://ubuntushare//stacks//populations_sex//sex.plink.ped",
	col.names=F, row.names=F, quote=F, eol='\n', sep='\t')
sex.assoc<-as.data.frame(read.table(
	"E://ubuntushare//stacks//populations_sex//sex.plink.assoc",
	header=T))

sex.assoc.sig<-sex.assoc[sex.assoc$P<=0.05,]
sex.map<-read.table("E://ubuntushare//stacks//populations_sex//batch_1.plink.map")

sex.map.sig<-sex.map[sex.map$V2 %in% sex.assoc.sig$SNP,]
sex.map.sig<-sex.map.sig[match(sex.map.sig$V2), sex.assoc.sig$SNP),]
sex.sig<-merge(sex.assoc.sig, sex.map.sig, by.x="SNP", by.y="V2")
sex.sig<-cbind(sex.sig, as.numeric(sub('_\\d+','',sex.sig$SNP)))
colnames(sex.sig)[11:14]<-c("SCAFFOLD","CHROM","BPN","RADLOC")
dim(sex.sig[sex.sig$RADLOC %in% sex.outliers,])
sex.out.assoc<-merge(sex.sig, sex.outliers, by.x="RADLOC", by.y="Locus")

dim(sex.out.assoc[sex.out.assoc$SCAFFOLD %in% use.contigs$Scaffold,])
summary(as.factor(
	merge(sex.out.assoc, use.contigs, by.x="SCAFFOLD", by.y="Scaffold")$LG))
#none of the RAD loci are on the linkage map, but some of the scaffolds map back

##############################################################################
#****************************************************************************#
######################MANIPULATING FILES FOR ASSOCIATIONS#####################
#****************************************************************************#
##############################################################################
#fix the weird ones--ALAL, merge lines where necessary. It'll make life easier

raw.pheno<-read.table("E://ubuntushare//popgen.pheno.txt", sep="\t", header=T)
fem.keep<-raw.pheno[!is.na(raw.pheno$BandNum) & substr(raw.pheno$ID,5,5)=="F",]
db.keep<-raw.pheno[substr(raw.pheno$ID,5,6)=="DB",]
db.keep<-db.keep[(db.keep$Side=="LEFT" || 
	db.keep$Side == "Left"),]
not.fem.keep<-raw.pheno[substr(raw.pheno$ID,5,5)!="F" & 
	substr(raw.pheno$ID,5,6)!="DB",]
not.fem.keep<-not.fem.keep[(not.fem.keep$Side=="LEFT" || 
	not.fem.keep$Side == "Left"),]

pheno.dat<-rbind(fem.keep, not.fem.keep, db.keep)
#use plink file to get list of individuals to keep
ped<-read.table("E://ubuntushare//stacks//populations//batch_1.plink.ped", 
	skip = 1, stringsAsFactors=F, colClasses="character")
ped.names<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', ped[,2])
ped.names<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', ped.names)
pheno.dat$ID<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', pheno.dat$ID)
ped.sex<-sub('sample_\\w{4}(\\w+).*[_.].*','\\1', ped[,2])
ped.sex[substr(ped.sex,1,2)=="DP"]<-"M"
ped.sex[substr(ped.sex,1,2)=="DB"]<-"F"
ped.sex[substr(ped.sex,1,2)=="NP"]<-"M"
ped.sex[substr(ped.sex,1,1)=="P"]<-"M"
ped.sex<-substr(ped.sex,1,1)
ped.sex[ped.sex=="F"]<-2
ped.sex[ped.sex=="M"]<-1
ped.new<-ped
ped.new[,2]<-ped.names
ped.new[,5]<-ped.sex
ped.new[,6]<-pops.pheno$std.length

#first add length to phenotype row in ped file, row 6
#order the phenotype data based on ped filename order
pops.pheno<-pheno.dat[pheno.dat$ID %in% ped.names,] #this does not have the juveniles
pops.pheno<-pops.pheno[match(ped.names,pops.pheno$ID),]
pops.pheno<-replace(pops.pheno, is.na(pops.pheno),as.numeric(-9))
pops.pheno<-pops.pheno[,-8]
pops.pheno<-cbind(substr(pops.pheno$ID,1,4), pops.pheno)
colnames(pops.pheno)[1]<-"PopID"
ped.new[,6]<-pops.pheno$std.length



write.table(ped.new, "E://ubuntushare//stacks//populations//plink.pheno.ped",
	row.names=F, col.names=F, quote=F, sep="\t", eol="\n")

#plink format: FamID, IndID, PhenA,PhenB,PhenC,PhenD,PhenE
#can have a header row
plink.pheno<-as.data.frame(cbind(ped.new[,1],pops.pheno$ID, pops.pheno$SVL,
	pops.pheno$depth,	(pops.pheno$SnoutLength/pops.pheno$HeadLength), 
	pops.pheno$SountDepth,	pops.pheno$MBandArea, pops.pheno$BandNum))
colnames(plink.pheno)<-c("FID", "IID", "SVL","BodyDepth", "SnoutHeadProp", 
	"SnoutDepth", "MeanBandArea", "BandNum")
plink.pheno<-plink.pheno[-which(plink.pheno$IID==-9),]
write.table(plink.pheno, "E://ubuntushare//stacks//populations//pheno.txt", 
	col.names=T, row.names=F, quote=F, sep='\t', eol='\n')
#run plink --ped plink.pheno.ped --map batch_1.plink.map --pheno pheno.txt \
	--assoc --all-pheno --qt-means --noweb  --out plink.pheno --allow-no-sex
band.num.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.pheno.BandNum.qassoc",
	header=T)
dim(band.num.assoc[band.num.assoc$P<0.05,])

band.area.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.pheno.MeanBandArea.qassoc",
	header=T)
snout.depth.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.pheno.SnoutDepth.qassoc",
	header=T)
snout.head.pr.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.pheno.SnoutHeadProp.qassoc",
	header=T)
body.depth.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.pheno.BodyDepth.qassoc",
	header=T)
body.length.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.length.qassoc",
	header=T)
svl.assoc<-read.table(
	"E://ubuntushare//stacks//populations//plink.pheno.qassoc",
	header=T)

pheno.dist<-dist(pops.pheno[,3:10])


num.bands<-pops.pheno[pops.pheno$BandNum!=-9,c("PopID","ID","BandNum")]
numbands.lmer<-lmer(BandNum~PopID+(1|ID), data=num.bands)
pheno.lmer<-lmer(as.matrix(pops.pheno[,3:10])~pops.pheno[,1]+(1|pops.pheno[,2]))



#create male and female files for pst analysis
fem.pheno<-pops.pheno[pops.pheno$BandNum!=-9,]
fem.pheno<-replace(fem.pheno, fem.pheno==-9,NA)
write.table(fem.pheno, 
	"E://ubuntushare//calculate_pairwise_pst//calculate_pairwise_pst//fem.pheno.txt",
	sep="\t", quote=F, col.names=T, row.names=F)
mal.pheno<-pops.pheno[pops.pheno$BandNum==-9,]
mal.pheno<-replace(mal.pheno, mal.pheno==-9, NA)
mal.pheno<-mal.pheno[,-9]
mal.pheno<-mal.pheno[,-9]
mal.pheno<-mal.pheno[complete.cases(mal.pheno),]
write.table(mal.pheno, 
	"E://ubuntushare//calculate_pairwise_pst//calculate_pairwise_pst//mal.pheno.txt",
	sep="\t", quote=F, col.names=T, row.names=F)
trait.names<-c("SVL", "Standard Length", "Body Depth", 
	"Snout Length", "Snout Depth", "Head Length")

jpeg("trait_summaries.jpeg", res=300, height=21, width=14, units="in")
par(mfrow=c(6,2), cex=0.5, mar=c(1,3,1,1), oma=c(3,2,3,1))
for(i in 3:7){
	plot(factor(fem.pheno[,1]),as.numeric(fem.pheno[,i]), 
		xaxt='n',xlab="", las=1, cex.axis=2)
	legend("top", trait.names[i-2], bty="n", cex=2)
	plot(factor(mal.pheno[,1]),as.numeric(mal.pheno[,i]), 
		xaxt='n',xlab="", las=1, cex.axis=2)
	legend("top", trait.names[i-2], bty="n", cex=2)
}
par(mar=c(5,3,0,1))
p1<-plot(factor(fem.pheno[,1]),fem.pheno[,8], xaxt='n',xlab="", las=1, 
	cex.axis=2)
axis(1, at=seq(1,12,1), labels=NA)
text(x=seq(1,12,1), y=9, srt = 45, adj = 1,cex=2,
     labels = p1$names, xpd = TRUE)
legend("top", trait.names[6], bty="n", cex=2)
p2<-plot(factor(mal.pheno[,1]),mal.pheno[,8], xaxt='n',xlab="", las=1,
	cex.axis=2)
axis(1, at=seq(1,12,1), labels=NA)
text(x=seq(1,12,1), y=6, srt = 45, adj = 1,cex=2,
     labels = p2$names, xpd = TRUE)
legend("top", trait.names[6], bty="n", cex=2)

par(fig = c(0, 1, 0, 1), oma=c(1,1,0,1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-.6, 1.0752, "Females", bty='n', cex=2)
legend(0.5, 1.0752, "Males", bty='n', cex=2)
mtext("Length (mm)", 2, outer=T, line=-.75, cex=1)
dev.off()

#bands
jpeg("female_trait_summaries.jpeg", res=300, height=7, width=14, units="in")
par(mfrow=c(1,2))
p1<-plot(factor(fem.pheno[,1]),fem.pheno$BandNum, xlab="",xaxt='n', las=1)
axis(1, at=seq(1,12,1), labels=NA)
text(x=seq(1,12,1), y=min(p1$out-1), srt = 45, adj = 1,
     labels = p1$names, xpd = TRUE)
mtext("Number of Bands", 2, outer=FALSE, line = 2)
legend("top", "Band Number", bty="n")
p2<-plot(factor(fem.pheno[,1]),fem.pheno$MBandArea, xaxt='n',xlab="", las=1)
axis(1, at=seq(1,12,1), labels=NA)
text(x=seq(1,12,1), y=-0.175, srt = 45, adj = 1,
     labels = p2$names, xpd = TRUE)
mtext(expression(paste("Area (mm"^"2)", sep="")), 2, outer=FALSE, line=2)
legend("top", "Mean Band Area", bty="n")
dev.off()


#########################################################################
#***********************************************************************#
##########RECALCULATE FSTS OF PCADAPT LOCI FOR PST-FST ANALYSIS##########
#***********************************************************************#
#########################################################################
#Fst=(variance_p)/(pbar*(1-pbar)) 
#pbar=avg allele freq across pops
#variance_p= variance in allele freqs among pops
#OR Fst=(HT-HS)/HT
#Let's go with Fst=1-((sum(HE in each pop))/(num.pops*HEtotal))
#and HE=1-(sum(p2))

pa.split<-split(pa5.out.dat, pa5.out.dat$Locus.ID)

calculate.fst<-function(dataframe){
	each.pop<-split(dataframe, dataframe$Locus.ID)	
	pop.he<-lapply(each.pop, function(x){ he<-x$P*x$P })
	Hs<-sum(dataframe$Exp.Het*dataframe$N)/(2*sum(dataframe$N))
	pbar<-sum(dataframe$P)/(2*sum(dataframe$N))
	qbar<-sum(1-dataframe$P)/
		(2*sum(dataframe$N))
	Ht<-1-sum((pbar*pbar)+(qbar*qbar))
	Fst<-(Ht-Hs)/Ht
	return(Fst)
}
pop.names<-levels(factor(pa.split[[1]]$Pop.ID))
pop.names<-pop.names[match(rownames(dist), pop.names)]
pairwise.fst<-list()
averages<-as.data.frame(setNames(replicate(length(pop.names), 
	numeric(0), simplify=F), pop.names))
for(j in 1:(length(pop.names)-1)){
		for(k in (j+1):length(pop.names)){
			averages[j,k]<-0
}}
names(averages)<-pop.names

for(i in 1:length(pa.split)){
	these.names<-levels(factor(pa.split[[i]]$Pop.ID))
	pairwise.fst[[i]]<-as.data.frame(setNames(replicate(length(these.names), 
		numeric(0), simplify=F), these.names))
	this.row<-names(pa.split)[i]
	pa.split[[i]]<-pa.split[[i]][match(
		rownames(dist), pa.split[[i]]$Pop.ID),]
	pop.split<-split(pa.split[[i]], pa.split[[i]]$Pop.ID)
	for(j in 1:(length(pop.split)-1)){
		for(k in (j+1):length(pop.split)){
			fsts<-calculate.fst(as.data.frame(
				rbind(pop.split[[j]], pop.split[[k]])))
			pairwise.fst[[i]][j,k]<-fsts
			averages[j,k]<-averages[j,k]+fsts
			this.row<-c(this.row, fsts)
		}
	}
	names(pairwise.fst)[[i]]<-names(pa.split)[[i]]
	pairwise.fst[[i]]<-rbind(pairwise.fst[[i]],
		rep(NA, ncol(pairwise.fst[[i]])))
	rownames(pairwise.fst[[i]])<-colnames(pairwise.fst[[i]])	
}
averages<-rbind(averages,rep(NA, ncol(averages)))
rownames(averages)<-colnames(averages)
averages<-averages/length(pa.split)#average fst values for pcadapt outliers

fem.pf.outliers<-list(
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "MBandArea")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "BandNum")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SVL")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "std.length")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "depth")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SnoutLength")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SountDepth")]))),
		as.dist(t(averages)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "HeadLength")]))),
		as.dist(t(averages)), nrepet=9999)
)#none of these are significant. bummer.

fem.ibd<-list(
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "MBandArea")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "BandNum")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SVL")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "std.length")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "depth")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SnoutLength")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SountDepth")]))),
		as.dist(dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "HeadLength")]))),
		as.dist(dist), nrepet=9999)
)
sub.dist<-dist[-7,-7]
sub.dist<-sub.dist[-5,-5]
fem.sub.ibd<-list(
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "MBandArea")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "BandNum")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "SVL")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "std.length")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "depth")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "SnoutLength")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "SountDepth")]))),
		as.dist(sub.dist), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "HeadLength")]))),
		as.dist(sub.dist), nrepet=9999)
)
sub.avg<-averages[-7,-7]
sub.avg<-sub.avg[-5,-5]
fem.sub.pf<-list(
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "MBandArea")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "BandNum")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "SVL")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "std.length")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "depth")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "SnoutLength")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "SountDepth")]))),
		as.dist(t(sub.avg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.sub[,c("PopID", "HeadLength")]))),
		as.dist(t(sub.avg)), nrepet=9999)
)



##################################CALCULATE PST###########################
fem.var<-data.frame(between=as.numeric(), within=as.numeric())
count<-1
for(i in 3:10){
	fem.var[count,]<-as.data.frame(VarCorr(
		lmer(fem.pheno[,i]~(1|fem.pheno[,1]),REML=FALSE)))[,"vcov"]
	count<-count+1
}
rownames(fem.var)<-colnames(fem.pheno)[3:10]
fem.pst<-fem.var$between/(fem.var$between+(2*fem.var$within))
fem.var<-cbind(fem.var, fem.pst)

#test with data from Sokal and Rohlf
#data from PA Thomas
test.dat<-data.frame(host=c(rep(1,8), rep(2,10),rep(3,13), rep(4,6)),
	width=c(380, 376,360,368,372,366,374,382,
	350,356,358,376,338,342,366,350,344,364,
	354,360,362,352,366,372,362,344,342,358,351,348,348,
	376,344,342,372,374,360))

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

pf.ba<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "MBandArea")]))),
	as.dist(fsts), nrepet=9999)
pf.bn<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "BandNum")]))),
	as.dist(fsts), nrepet=9999)
pf.sv<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SVL")]))),
	as.dist(fsts), nrepet=9999)
pf.ln<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "std.length")]))),
	as.dist(fsts), nrepet=9999)
pf.dp<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "depth")]))),
	as.dist(fsts), nrepet=9999)
pf.sl<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SnoutLength")]))),
	as.dist(fsts), nrepet=9999)
pf.sd<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SountDepth")]))),
	as.dist(fsts), nrepet=9999)
pf.hl<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "HeadLength")]))),
	as.dist(fsts), nrepet=9999)

pm.sv<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "SVL")]))),
	as.dist(fsts), nrepet=9999)
pm.ln<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "std.length")]))),
	as.dist(fsts), nrepet=9999)
pm.dp<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "depth")]))),
	as.dist(fsts), nrepet=9999)
pm.sl<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "SnoutLength")]))),
	as.dist(fsts), nrepet=9999)
pm.sd<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "SountDepth")]))),
	as.dist(fsts), nrepet=9999)
pm.hl<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "HeadLength")]))),
	as.dist(fsts), nrepet=9999)

####PCA to capture variation??
prcomp(fem.pheno)

##########USE DATA FROM POPULATIONS RUN WITH FASTSTRUCTURE GROUPS###########

fst.mat.fsg<-read.table("E://ubuntushare//stacks//populations_fstrugroups//fst_summary_fstrugroups.txt",
	 header=T, row.names=1, sep='\t')
fst.mat.fst<-t(fst.mat.fsg)

#use fast.groups.5 from population_structure.R
fem.pheno<-cbind(fem.pheno,fast.groups.5[
	sub('(sample_)([A-Z]{5,7}[0-9]?[0-9]?[0-9]?)(.fq)?(_align)',
	'\\2',fast.groups.5$inds) %in% fem.pheno$ID,2])

fem.fsg.pst<-list(
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "MBandArea")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "BandNum")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "SVL")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "std.length")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "depth")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "SnoutLength")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "SountDepth")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999),
	mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("V2", "HeadLength")]))),
		as.dist(t(fst.mat.fsg)), nrepet=9999)
)


#############################P-MATRIX#########################################
#standardize traits first
standardize.trait<-function(trait){
	trait/mean(trait)
}

#GLOBAL
std.fem<-as.data.frame(cbind(fem.pheno[,1:2], standardize.trait(fem.pheno[,3]), 
	standardize.trait(fem.pheno[,4]),standardize.trait(fem.pheno[,5]),
	standardize.trait(fem.pheno[,6]),standardize.trait(fem.pheno[,7]),
	standardize.trait(fem.pheno[,8]),standardize.trait(fem.pheno[,9]),
	standardize.trait(fem.pheno[,10])))
colnames(std.fem)<-colnames(fem.pheno)
pmat.fem.global<-cov(std.fem[,3:10])

#pca of traits
fem.phen.uns.rda<-rda(fem.pheno[,3:10])
fem.phen.std.rda<-rda(std.fem[,3:10])
#pca of p-matrices
fem.pmat.uns.rda<-rda(cov(fem.pheno[,3:10]))
fem.pmat.std.rda<-rda(pmat.fem.global)

fem.phen.uns.rda$CA$v #shows loadings per trait


std.mal<-as.data.frame(cbind(mal.pheno[,1:2], standardize.trait(mal.pheno[,3]), 
	standardize.trait(mal.pheno[,4]),standardize.trait(mal.pheno[,5]),
	standardize.trait(mal.pheno[,6]),standardize.trait(mal.pheno[,7]),
	standardize.trait(mal.pheno[,8])))
colnames(std.mal)<-colnames(mal.pheno)
pmat.mal.global<-cov(std.mal[,3:8])
pm.m.g.rda<-rda(pmat.mal.global)

#BETWEEN POPS
calc.pmat<-function(phen.dat.list, dim1, dim2){
	pmat<-lapply(phen.dat.list, function(x){
		cov(x[,dim1:dim2])
	})
	return(pmat)
}

#females
fem.pheno$PopID<-factor(fem.pheno$PopID)
pops.fem.unstd.dat<-split(fem.pheno, fem.pheno$PopID)
pops.fem.std.dat<-lapply(pops.fem.unstd.dat, function(x){
	y<-as.data.frame(cbind(x[,1:2], standardize.trait(x[,3]), 
	standardize.trait(x[,4]),standardize.trait(x[,5]),
	standardize.trait(x[,6]),standardize.trait(x[,7]),
	standardize.trait(x[,8]),standardize.trait(x[,9]),
	standardize.trait(x[,10])))}
)
pmat.fem.std.pops<-calc.pmat(pops.fem.std.dat,3,10)
pmat.fem.unst.pops<-calc.pmat(pops.fem.unstd.dat,3,10)

#males
mal.pheno$PopID<-factor(mal.pheno$PopID)
pops.mal.unstd.dat<-split(mal.pheno, mal.pheno$PopID)
pops.mal.std.dat<-lapply(pops.mal.unstd.dat, function(x){
	y<-as.data.frame(cbind(x[,1:2], standardize.trait(x[,3]), 
	standardize.trait(x[,4]),standardize.trait(x[,5]),
	standardize.trait(x[,6]),standardize.trait(x[,7]),
	standardize.trait(x[,8])))}
)
pmat.mal.std.pops<-calc.pmat(pops.mal.std.dat,3,8)
pmat.mal.unst.pops<-calc.pmat(pops.mal.unstd.dat,3,8)

#write pmatrices to file
#following format of presentation in Bertram et al 2011
#variance on diagonal, covariance lower, correlations upper, 
#with p1,p2,and p3 to right

for(i in 1:length(pmat.fem.unst.pops)){
	dat<-pmat.fem.unst.pops[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.fem.unst.pops[[i]])[
		upper.tri(cor(pmat.fem.unst.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.fem.unst.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(fem.pheno)[3:10],
		"p1","p2","p3")
	rownames(dat)<-colnames(fem.pheno)[3:10]
	write.table(dat, 
		paste("pop.",names(pmat.fem.unst.pops)[i], ".pmat_unstd.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}
for(i in 1:length(pmat.fem.std.pops)){
	dat<-pmat.fem.std[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.fem.std.pops[[i]])[
			upper.tri(cor(pmat.fem.std.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.fem.std.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(fem.pheno)[3:10], "p1","p2","p3")
	rownames(dat)<-colnames(fem.pheno)[3:10]
	write.table(dat, 
		paste("pop.",names(pmat.fem.std)[i], ".pmat.std.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}


for(i in 1:length(pmat.mal.unst.pops)){
	dat<-pmat.mal.unst.pops[[i]]
	dat[upper.tri(dat)]<-
		as.numeric(cor(pmat.mal.unst.pops[[i]])[
		upper.tri(cor(pmat.mal.unst.pops[[i]]))])
	dat<-as.matrix(
		cbind(dat,eigen(pmat.mal.unst.pops[[i]])$vectors[,1:3]))
	colnames(dat)<-c(colnames(mal.pheno)[3:8],
		"p1","p2","p3")
	rownames(dat)<-colnames(mal.pheno)[3:8]
	write.table(dat, 
		paste("pop.",names(pmat.mal.unst.pops)[i], 
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
	colnames(dat)<-c(colnames(mal.pheno)[3:8], "p1","p2","p3")
	rownames(dat)<-colnames(fem.pheno)[3:8]
	write.table(dat, 
		paste("pop.",names(pmat.mal.std.pops)[i],
			".male.pmat.std.txt", sep=""), 
		sep='\t',
		eol='\n', row.names=T, col.names=T, quote=F)

}

fem.pmat.rda<-lapply(pmat.fem.unst.pops, rda)

###ANGLES###
library(boot)
lead.ueigvec<-function(data){
	ev.1<-eigen(data)$vectors
	e1.norm <-ev.1[,1]/sqrt(sum(ev.1[,1]^2))
	return(e1.norm)
}
calculate.angle<-function(vec1, vec2){
	angle<-acos(vec1%*%vec2)*180/pi
	return(angle)
}
compare.pmatrix.angles<-function(df1, df2){
	ev1<-lead.ueigvec(df1)
	ev2<-lead.ueigvec(df2)
	angle<-calculate.angle(ev1,ev2)
}

angles<-matrix(nrow=length(pmat.fem.std.pops), ncol=length(pmat.fem.std.pops))
for(i in 1:(length(pmat.fem.std.pops)-1)){
	for(j in (i+1):(length(pmat.fem.std.pops))){
		angles[i,j]<-compare.pmatrix.angles(pmat.fem.std.pops[[i]], 
			pmat.fem.std.pops[[j]])
		
	}
}
row.names(angles)<-names(pmat.fem.std.pops)
colnames(angles)<-names(pmat.fem.std.pops)
write.table(angles, "angles.txt", 
	sep='\t', eol='\n', quote=F,
	row.names=T, col.names=T)

mal.angles<-matrix(nrow=length(pmat.mal.std.pops), 
	ncol=length(pmat.mal.std.pops))
for(i in 1:(length(pmat.mal.std.pops)-1)){
	for(j in (i+1):(length(pmat.mal.std.pops))){
		mal.angles[i,j]<-compare.pmatrix.angles(pmat.mal.std.pops[[i]], 
			pmat.mal.std.pops[[j]])
		
	}
}
row.names(mal.angles)<-names(pmat.mal.std.pops)
colnames(mal.angles)<-names(pmat.mal.std.pops)
write.table(mal.angles, "male.angles.txt", 
	sep='\t', eol='\n', quote=F,
	row.names=T, col.names=T)

###bootstrapping###
angle.for.boot<-function(df, i){
	dat<-df[i,-2]
	dat$PopID<-factor(dat$PopID)
	split.data<-split(dat, dat$PopID)
	std.data<-lapply(split.data, function(x){
		y<-as.data.frame(cbind(standardize.trait(x[,2]),
			standardize.trait(x[,3]), 
			standardize.trait(x[,4]),standardize.trait(x[,5]),
			standardize.trait(x[,6]),standardize.trait(x[,7]),
			standardize.trait(x[,8]),standardize.trait(x[,9])))
		return(y)
		}
	)
	pmat<-calc.pmat(std.data,1,8)
	angle<-compare.pmatrix.angles(pmat[[1]],pmat[[2]])	
	return(angle)
}
boot.results<-list()
boot.angles<-matrix(nrow=length(pmat.fem.std.pops), ncol=length(pmat.fem.std.pops))
pmat.names<-NULL
count<-1
for(i in 1:(length(pmat.fem.std.pops)-1)){
	for(j in (i+1):(length(pmat.fem.std.pops))){
		name1<-names(pmat.fem.std.pops)[i]
		name2<-names(pmat.fem.std.pops)[j]
		boot.data<-fem.pheno[fem.pheno$PopID==name1 | fem.pheno$PopID==name2,]
		rownames(boot.data)<-NULL
		boot.res<-boot(boot.data, angle.for.boot, R=999)
		ci<-boot.ci(boot.res, 0.95, type = "basic")
		boot.results[[count]]<-as.matrix(cbind(boot.res$t0,
			ci$basic[,4], ci$basic[,5]))
		pmat.names<-c(pmat.names, paste(name1,"-",name2,sep=""))
		boot.angles[i,j]<-boot.res$t0
		count<-count+1
	}
}
names(boot.results)<-pmat.names

boot.out<-data.frame(angle=numeric(length(boot.results)),
	low.ci=numeric(length(boot.results)),
	upp.ci=numeric(length(boot.results)))
for(i in 1:length(boot.results)){
	boot.out[i,]<-boot.results[[i]]
}
row.names(boot.out)<-pmat.names
write.table(boot.out, "angle.bootstrap.results.txt", 
	sep='\t', eol='\n', quote=F,
	row.names=T, col.names=T)



###############################BLOWS METHOD#################################
#A=first three eigenvectors of P for first species
#B=first three eigenvectors of P for second species
#S=t(A)Bt(B)A
#sim=sum(eigen(S))/3

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

#########TEST ON TWO TRAITS##########
two.traits<-fem.pheno[,c(1,3,4)]
two.traits.split<-split(two.traits, two.traits$PopID)
two.traits.split.std<-lapply(two.traits.split, function(x){
	y<-as.data.frame(cbind(x[,1],standardize.trait(x[,2]),
	standardize.trait(x[,3])))
})
two.traits.pmat<-calc.pmat(two.traits.split.std,2,3)
plot.pmatrix(two.traits.pmat[[1]])
dat.cov<-cov(eigen(two.traits.pmat[[1]])$vectors[,1:2])
plot(c(-1,1), c(-1,1) , asp=1, type="n", xlim=c(-1,1), ylim=c(-1,1), 
		cex.lab=1.5, xlab="",ylab="", las=1)

plot.eigenvectors(dat.cov, 1)
	plot.eigenvectors(dat.cov, 2)

#######################################
#unstandardized
fem.sim.unstd.mat<-generate.sim.mat(pmat.fem.unst.pops)
fem.sim.unstd.fst<-mantel.rtest(as.dist(fem.sim.unstd.mat), 
	as.dist(t(pwise.fst)),	nrepet=9999)
mal.sim.unstd.mat<-generate.sim.mat(pmat.mal.unst.pops)
mal.sim.unstd.fst<-mantel.rtest(as.dist(mal.sim.unstd.mat), 
	as.dist(t(pwise.fst)),	nrepet=9999)
#standardized
fem.sim.std.mat<-generate.sim.mat(pmat.fem.std)
fem.sim.std.fst<-mantel.rtest(as.dist(fem.sim.std.mat), as.dist(t(pwise.fst)),
	nrepet=9999)

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
all.pmax<-lapply(pmat.fem.unst.pops, find.pmax)

#standardized
pmax.std.mat<-vector.correlations(pmat.fem.std)
blows.sim.std.mat<-fem.sim.std.mat
blows.sim.std.mat[upper.tri(blows.sim.std.mat)]<-
	pmax.std.mat[upper.tri(pmax.std.mat)]
write.table(blows.sim.std.mat, "standardized_p_sim_matrix.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)
#unstandardized
pmax.unstd.mat<-vector.correlations(pmat.fem.unst.pops)
blows.sim.unstd.mat<-fem.sim.unstd.mat
blows.sim.unstd.mat[upper.tri(blows.sim.unstd.mat)]<-
	pmax.unstd.mat[upper.tri(pmax.unstd.mat)]
write.table(blows.sim.unstd.mat, "unstd_p_sim_matrix.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)

pmax.mal.unstd.mat<-vector.correlations(pmat.mal.unst.pops)
blows.mal.sim.unstd.mat<-mal.sim.unstd.mat
blows.mal.sim.unstd.mat[upper.tri(blows.mal.sim.unstd.mat)]<-
	pmax.mal.unstd.mat[upper.tri(pmax.mal.unstd.mat)]
write.table(blows.mal.sim.unstd.mat, "unstd_p_sim_matrix_male.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)


#####################TENSORS############################
#Adapted from Aguirre et al. 2013 supplemental material
library(gdata);library(matrixcalc);library(MCMCglmm)

pmat.mal.unst.pops<-calc.pmat(pops.mal.unstd.dat,3,10)
pmat.fem.unst.pops<-calc.pmat(pops.fem.unstd.dat,3,10)

covtensor<-function(matrix.list){
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

#variance along e1 and e2, 1 for each population of each replicate line
f.e11.proj <- lapply(pmat.fem.unst.pops, proj, b = f.e1.L)
f.e21.proj <- lapply(pmat.fem.unst.pops, proj, b = f.e2.L)
f.e31.proj <- lapply(pmat.fem.unst.pops, proj, b = f.e3.L)

m.e11.proj <- lapply(pmat.mal.unst.pops, proj, b = m.e1.L)
m.e21.proj <- lapply(pmat.mal.unst.pops, proj, b = m.e2.L)
m.e31.proj <- lapply(pmat.mal.unst.pops, proj, b = m.e3.L)

#summary plots
jpeg("blows.method4.summary.jpeg", width=14, height=14, units="in", res=300)
layout(matrix(c(1,1,2,3),2,2,byrow=F))
layout.show(3)
#plot eigenvalues of non-zero eigentensors for S
par(mar=c(5,4,2,0))
plot(fem.tensor$s.alpha[,1:fem.tensor$nonzero], ylab="alpha",
	xlab="", xaxt="n", las=1)#3 until it levels off.
axis(1, at=seq(1,fem.tensor$nonzero,1), 
	labels=paste("E",1:fem.tensor$nonzero, sep=""))
points(mal.tensor$s.alpha[,1:mal.tensor$nonzero],pch=6)
legend("topright", pch=c(1,6), c("Males", "Females"))

#plot the variance in each population in the direction of e11,e21, and e31
par(mar=c(3,5,2,1))
plot(x=seq(1,12,1),f.e11.proj, xaxt="n", las=1,ylim=c(-50,200),
	xlab="", ylab="lambda")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.names, par("usr")[1]-75,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),f.e11.proj, lty=2)
points(x=seq(1,12,1),f.e21.proj, pch=19)
lines(x=seq(1,12,1),f.e21.proj, lty=1)
points(x=seq(1,12,1),f.e31.proj, pch=15)
lines(x=seq(1,12,1),f.e31.proj,lty=4)
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4), c("e1","e2","e3"))
text(x=2,y=200, "Females")

par(mar=c(5,5,0,1))
plot(x=seq(1,12,1),m.e11.proj, xaxt="n", las=1,ylim=c(-50,200),
	xlab="Population", ylab="lambda")
axis(1, at=seq(1,12,1), labels=F)
text(x=seq(1,12,1), labels=pop.names, par("usr")[1]-75,
	srt=-45, xpd=TRUE)
lines(x=seq(1,12,1),m.e11.proj, lty=2)
points(x=seq(1,12,1),m.e21.proj, pch=19)
lines(x=seq(1,12,1),m.e21.proj, lty=1)
points(x=seq(1,12,1),m.e31.proj, pch=15)
lines(x=seq(1,12,1),m.e31.proj,lty=4)
legend("bottomright", pch=c(1,19,15), lty=c(2,1,4), c("e1","e2","e3"))
text(x=2,y=200, "Males")
dev.off()


