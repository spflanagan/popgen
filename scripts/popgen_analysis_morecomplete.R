#Author: Sarah P. Flanagan
#Date: 27 July 2015
#Purpose: Analyze Population genetics data
#made from used analyses from fst_cline_analysis and population_structure

rm(list=ls())

library(ade4)
library(lme4)
library(maps);library(gplots)
library(mapdata)
library(vegan)
library(boot);library(car)
library(adegenet)
library(scales)
library(gdata);library(matrixcalc)

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

dist<-read.table("E://Docs//PopGen//geographical_distances.txt", 
	header=T, row.names=1, sep='\t')
pwise.fst<-read.table("E:\\Docs\\PopGen\\fst_summary.txt",
	 header=T, row.names=1, sep='\t')
env.dat<-read.table("E://ubuntushare//bayenv2//env_data_bayenv_raw.txt")

setwd("E://Docs//PopGen")
#############################################################################
#***************************************************************************#
#################################FUNCTIONS###################################
#***************************************************************************#
#############################################################################

#***************************************************************************#
#PLOT ANY GENOME-WIDE STATISTIC
#***************************************************************************#
plot.genome.wide<-function(bp,var,y.max,x.max, rect.xs=NULL,y.min=0,x.min=0, 
	plot.new=FALSE, plot.axis=TRUE, rect.color="white",plot.rect=TRUE, 
	pt.cex=1, pt.col="black"){
	#********************************************
	#this function plots a variable without scaffold info. 
	#feed it the basepair (x) values and variable (y) values 
	#*********************************************
	if(plot.new==TRUE){ par(new=new) }
	plot(bp, var,xlab="",ylab="", new=plot.new,
		type="n", bg="transparent", axes=F, bty="n", 
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	if(plot.rect==TRUE){
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
	}
	if(plot.axis){
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
	points(bp, var, pch=19, cex=pt.cex,col=pt.col,
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

#***************************************************************************#
##AUTOMATED FST PLOTTING
#***************************************************************************#
plot.fsts.scaffs<-function(dat, dat.name, ci.dat=NULL, pt.cex=1,y.lab=NULL,
	col.pts=NULL, col.pt.col="dark green", col.pt.pch=8, col.pt.cex=2){
	#********************************************
	#this function plots every scaffold. 
	#dat is a fst file from stacks_genomewideCIs
	#dat should have a column named "Chr", one called "BP"
	#at least one other column is needed: dat.name
	#*********************************************
	byscaff<-split(dat, factor(dat$Chr))
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
				dat.name])
		}}
	}

	x.min<-min(addition.values)
	x.max<-max(addition.values)
	y.max<-max(dat[,dat.name])+0.5*max(dat[,dat.name])
	if(min(dat[,dat.name]) < 0) {
		y.min<-min(dat[,dat.name]) + 0.5*min(dat[,dat.name])
	} else {
		y.min<-0
	}

	plot(new.x[[1]], byscaff[[1]][,dat.name], 
		xlim=c(x.min,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(byscaff)){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		rect(rect.xs[j,1],y.min,rect.xs[j,2],y.max, 
			col=rect.color, border=NA)
	}
	for(j in 1:length(byscaff)){
		points(new.x[[j]],byscaff[[j]][,dat.name], pch=19, cex=pt.cex,
			xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	}
		#plot.genome.wide(new.x[[j]], 
		#	byscaff[[j]][,dat.name],rect.xs=rect.xs[j,],
		#	y.max,x.max, y.min=y.min,x.min=x.min, 
		#	plot.new=FALSE, plot.axis=FALSE, pt.cex=0.25)
		
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=10)),
		xlim=c(x.min,x.max),ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
	if(is.null(y.lab)){
		mtext(dat.name,2,las=1, line=4,cex=.75)}
	else{
		mtext(y.lab,2,las=1, line=4,cex=.75)}
	if(!is.null(ci.dat)){
		clip(x.min,x.max,y.min,y.max)
		abline(h=ci.dat$CI99smooth,col="red")
		if(!is.null(col.pts)){
			clip(x.min,x.max,ci.dat$CI99smooth, y.max)
			points(col.pt.x, col.pt.y,
				col=col.pt.col, pch=col.pt.pch, cex=col.pt.cex)
		}
	}
	return(as.data.frame(cbind(old.bp=dat$BP,new.bp=unlist(new.x), 
		locus=dat$Locus)))
}

#***************************************************************************#
##REORDER A DATAFRAME
#***************************************************************************#
reorder.df<-function(dat,order.list){
	#dat has to have the grouping IDs in row 1
	#those grouping ids must match the factors in order.list
	dat.sep<-split(dat, dat[,1])
	dat.new<-dat.sep[[order.list[1]]]
	for(i in 2:length(order.list)){
		dat.new<-rbind(dat.new, dat.sep[[order.list[i]]])
	}
	return(dat.new)
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
#uses ade4 package

t.pwise.fsts<-t(pwise.fsts)

ibd<-mantel.rtest(as.dist(t(dist)),as.dist(t(pwise.fst)))
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
pca1<-glPca(dat.plink, parallel=FALSE, nf=5)
scatter(pca1)

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")

#or my custom pca plotting skillzz
ind.names<-dimnames(pca1$scores)[[1]]
pop.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
	"FLAB","FLPB","FLHB","FLCC")
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

jpeg("adegenet.dapc.k3.jpeg", height=7, width=7, units="in", res=300)
par(mar=c(5,4,4,2)+0.1, oma=c(5,4,4,2)+0.1)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
mtext("Discriminant Axis 1", 1, line = 2, outer=TRUE)
mtext("Discriminant Axis 2",2, line = 2, outer=TRUE)
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

#########################################################################
#***********************************************************************#
################################STRUCTURE################################
#***********************************************************************#
#########################################################################
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

####RUN STRUCTURE ON ONLY ONE POP..TXSP is largest
subset.txsp<-inds[pop.id=="TXSP"]
write.table(subset.txsp,"E://ubuntushare//stacks//populations//txsp.indlist", 
	quote=FALSE, col.names=F, row.names=F, eol="\n", sep="\t")


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
ld.map<-read.table("stacks//populations//subset.map",	header=F)
pcadapt.outliers<-ld.map[rep.snps,]
colnames(pcadapt.outliers)<-c("chr","loc","dist","bp")
pa.out.radloc<-sub('(\\d+)_\\d+','\\1',pcadapt.outliers$loc)

summ.dat<-read.table("E://ubuntushare//stacks//populations//batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
pa.out.dat<-summ.dat[summ.dat$Locus.ID %in% pa.out.radloc,]
pa.out.dat$Chr<-factor(pa.out.dat$Chr)
#are any on the linkage map?
length(levels(as.factor(use.contigs[use.contigs$Scaffold %in% pa.out.dat$Chr,1])))

#plot individual scores #used to be height=width=12
jpeg("pcadapt.scores1.2.jpeg", height=7, width=7, units="in", res=300)
plot(as.numeric(scores[1,]),as.numeric(scores[2,]),pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("bottomright", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0117)", 1, line = 2)
mtext("PC2 (rho2: 0.0091)",2, line = 2)
dev.off()

jpeg("pcadapt.scores1.3.jpeg", height=7, width=7, units="in", res=300)
plot(as.numeric(scores[1,]),as.numeric(scores[3,]),pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("topright", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(12), 0.5), ncol=3)
mtext("PC1 (rho2: 0.0117)", 1, line = 2)
mtext("PC3 (rho2: 0.0018)",2, line = 2)
dev.off()


##############################################################################
#****************************************************************************#
#############################GLOBAL FST OUTLIERS##############################
#****************************************************************************#
##############################################################################
#Import global values generated by C++ program calculate_global_fsts
#pruned loci found in all 12 populations, 1833 of them.
global.fsts<-read.table(
	"E://ubuntushare//stacks//populations//ld.hwe.sub.globalstats.txt",
	header=T, sep='\t')

fst.ci.99<-mean(global.fsts$Fst) + 2.57583*sd(global.fsts$Fst)
#what are the outliers??
glob.fst.out<-global.fsts[global.fsts$Fst >= fst.ci.99,]

#compare to other anlyses
summ.gfo<-summ.dat[summ.dat$Locus.ID %in% glob.fst.out$Locus,]
summ.gfo<-summ.gfo[!duplicated(summ.gfo$Locus),]
glob.fst.out<-
	merge(summ.gfo[,2:3], glob.fst.out, by.x="Locus.ID", by.y="Locus")
glob.fst.out<-glob.fst.out[,-3]
glob.fst.out[glob.fst.out$Locus %in% no.dups.outliers$Locus,]

#non-pruned loci in all 12 populations
np.global.fsts<-read.table(
	"E://ubuntushare//stacks//populations//allpops.nonpruned.globalstats.txt",
	header=T, sep='\t')
np.fst.ci.99<-mean(np.global.fsts$Fst) + 2.57583*sd(np.global.fsts$Fst)
np.glob.fst.out<-np.global.fsts[np.global.fsts$Fst >= np.fst.ci.99,]

#non-pruned loci in any of the 12 populations
any.global.fsts<-read.table(
	"E://ubuntushare//stacks//populations_nopopslimit//anypop.nonpruned.globalstats.txt",
	header=T, sep='\t')
any.fst.ci.99<-mean(any.global.fsts$Fst) + 2.57583*sd(any.global.fsts$Fst)
any.glob.fst.out<-np.global.fsts[any.global.fsts$Fst >= any.fst.ci.99,]

#compare these
dim(glob.fst.out[glob.fst.out %in% any.glob.fst.out,])#0
np.glob.fst.out[np.glob.fst.out$Locus %in% any.glob.fst.out$Locus,]#230
any.glob.fst.out[!(any.glob.fst.out$Locus %in% np.glob.fst.out$Locus),]#1
dim(np.glob.fst.out[np.glob.fst.out$Locus %in% no.dups.outliers$Locus,])#141
#************************************PLOT*************************************#
#plot the anypops one (more background? stronger signal?)
all.scaff<-split(np.global.fsts, np.global.fsts$Chrom)
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
y.max<-max(np.global.fsts$Fst)+0.1*max(np.global.fsts$Fst)
y.min<-min(np.global.fsts$Fst)-0.1*min(np.global.fsts$Fst)
if(min(np.global.fsts$Fst) < 0) {
	y.min<-min(np.global.fsts$Fst) - 0.1*min(np.global.fsts$Fst)
} else {
	y.min<-0
}

#jpeg("nonpruned.global.fst.jpeg", width=24, height=12, units="in", res=300)
jpeg("nonpruned.global.fst.small.jpeg", width=14, height=4.5, units="in", res=300)
par(mar=c(0,3,1,0), oma=c(1,3,0,0),cex=2)
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
		all.scaff[[i]]$Fst,plot.rect=FALSE,
		y.max,x.max, rect.xs[i,],y.min=0,x.min=0, 
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$Fst >= np.fst.ci.99,]
	#temp.out.prune<-all.scaff[[i]][all.scaff[[i]]$Locus %in% 
	points(temp.sig$BP, temp.sig$Fst, col="red", pch=19, cex=0.5)
}
axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
	ylim = c(y.min, y.max), pos=0,
	las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
#axis(1, at=c(0,rect.xs[,2]), xlim=c(0,(rect.xs[nrow(rect.xs),2]+5)), 
#	ylim=c(0,y.max), las=1, tck=-0.01, labels=FALSE, pos=0, xpd=TRUE,
#	par("usr")[3]-0.00051,srt = 38, adj = 1, xpd = TRUE,cex=0.75)
legend("top",c("Outlier", "Neutral"), col=c("red","black"),bty="n",
	horiz=TRUE,pt.cex=0.5, pch=19,cex=0.75)
mtext(side=1, "Genomic Location", outer = TRUE, cex=1.5)
mtext(side=2, "Global Fst", outer=TRUE, line=1,cex=1.5)
#clip(x.min,x.max,y.min,y.max)
#abline(h=np.fst.ci.99)
dev.off()


#########################################################################
#***********************************************************************#
################################BAYENV2##################################
#***********************************************************************#
#########################################################################

#**************************STARTING WITH PLINK FILES********************#
ped<-read.table("E:\\ubuntushare\\stacks\\populations\\ld.hwe.sub.ped", 
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

write.table(ped,"E://ubuntushare//stacks//populations//bayenv.plink.ped", 
	row.names=F, col.names=F, quote=F, sep="\t",eol="\n")

clust.plink<-cbind(ped.pops, ped[,2],ped.pops)
write.table(clust.plink, 
	"E://ubuntushare//stacks//populations//plink.clust.txt",
	col.names=F, row.names=F, quote=F, sep="\t", eol="\n")


#plink.map -> numbers instead of scaffold_#
map<-read.table("E://ubuntushare//stacks//populations//ld.hwe.sub.map", 
	skip = 1)
chr.nums<-sub('scaffold_','',map[,1])
map[,1]<-chr.nums

write.table(map, "E://ubuntushare//stacks//populations//bayenv.plink.map", 
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
#and NCHROBS-MAC for every snp at every pop
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
matrix.files<-list.files("E://ubuntushare//bayenv2//",pattern="matrix")
matrices<-list()
for(i in 1:length(matrix.files))
{
	matrices[[i]]<-as.matrix(read.table(
		paste("E://ubuntushare//bayenv2//",matrix.files[i],sep=""), 
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

avg.matrix<-matrices[[1]]
for(i in 2:length(matrices)){
	avg.matrix<-avg.matrix+matrices[[i]]
}
avg.matrix<-avg.matrix/length(matrices)
write.table(avg.matrix, "E://ubuntushare//bayenv2//all.6348.matrix",
	sep="\t", eol="\n", quote=F, row.names=F, col.names=F)
#************************************SNPFILEs******************************#
#for SNPFILE, need just one file per SNP apparently.
#want to use all of the snps...need to get map with those inds.
all.snps.ped<-read.table("E://ubuntushare//stacks//populations//plink.ped", 
	header=F, stringsAsFactors=F)
ped.pop<-sub('sample_(\\w{4})\\w+.*[_.].*','\\1', all.snps.ped[,2])
all.snps.clust<-cbind(ped.pop,all.snps.ped[,2],ped.pop)
write.table(all.snps.clust, 
	"E://ubuntushare//stacks//populations//all.6348.clust.txt", 
	sep="\t", eol="\n", quote=F, row.names=F, col.names=F)
#then need to run plink --file batch_1.plink --freq --within all.6348.clust.txt \
	#--allow-no-sex --noweb --out all.bayenv.plink

#read in frequency per pop
all.snps.frq<-read.table("E://ubuntushare//stacks//populations//all.6348.frq.strat", 
	header=T, stringsAsFactors=F)
all.snps.frq<-cbind(all.snps.frq,all.snps.frq$NCHROBS-all.snps.frq$MAC)
colnames(all.snps.frq)[ncol(all.snps.frq)]<-"NAC"
pop.order<-levels(as.factor(all.snps.frq$CLST))
snp.names<-split(all.snps.frq$SNP,all.snps.frq$CLST)[[1]]

mac.by.pop<-as.data.frame(split(all.snps.frq$MAC,all.snps.frq$CLST))
rownames(mac.by.pop)<-snp.names
nac.by.pop<-as.data.frame(split(all.snps.frq$NAC,all.snps.frq$CLST))
rownames(nac.by.pop)<-snp.names
snps.dat<-interleave(mac.by.pop,nac.by.pop)

for(i in 1:nrow(mac.by.pop)){
	snps.dat<-interleave(mac.by.pop[i,],nac.by.pop[i,])
	write.table(snps.dat, 
		paste("E://ubuntushare//bayenv2//snpfiles//all.6348.",
			snp.names[i],sep=""), 
		col.names=F,row.names=F,quote=F,sep="\t",eol="\n")
}
#***************************PROCESS THE OUTPUT**************************#
####Environmental Associations####
all.6348.map<-read.table(
	"E://ubuntushare//stacks//populations//all.6348.plink.map",
	header=F)
bf.env<-read.table("E://ubuntushare//bayenv2//bf_environ.env_data_bayenv.txt")
colnames(bf.env)<-c("locus", "Temp_BF", "Temp_rho", "Temp_rs", 
	"Salinity_BF", "Salinity_rho", "Salinity_rs", "coll.temp_BF", 
	"coll.temp_rho", "coll.temp_rs", "coll.sal_BF", "coll.sal_rho", 
	"coll.sal_rs", "seagrass_BF", "seagrass_rho","seagrass_rs")
bf.env$locus<-sub('./snpfiles/all.6348.(\\d+)','\\1', bf.env$locus)

bf<-bf.env[,c(1,2,5,8,11,14)] #not normally distr
bf$rad.loc<-sub('(\\d+)(_\\d+)','\\1',bf$locus)
rho<-bf.env[,c(1,3,6,9,12,15)] #normal
rho$rad.loc<-sub('(\\d+)(_\\d+)','\\1',rho$locus)
rs<-bf.env[,c(1,4,7,10,13,16)] #normal
rs$rad.loc<-sub('(\\d+)(_\\d+)','\\1',rs$locus)


bf.scaff<-merge(all.6348.map, bf, by.x="V2", by.y="locus")
colnames(bf.scaff)[1:4]<-c("locus","scaffold","dist","BP")
#focus on Bayes Factors, because of Lotterhos & Whitlock (2015)
bf.co<-apply(bf[,-1],2,quantile,0.95)
temp.bf.sig<-bf[bf$Temp_BF>bf.co["Temp_BF"],c(1,2)]
sal.bf.sig<-bf[bf$Sal_BF>bf.co["Salinity_BF"],c(1,3)]
ctemp.bf.sig<-bf[bf$ctemp_bf>bf.co["coll.temp_bf"],c(1,4)]
csal.bf.sig<-bf[bf$csal_bf>bf.co["coll.sal_bf"],c(1,5)]
grass.bf.sig<-bf[bf$seagrass_bf>bf.co["seagrass_bf"],c(1,6)]


#compare to other analyses
dim(temp.bf.sig[temp.bf.sig$rad.loc %in% pa.out.radloc,])#0
dim(sal.bf.sig[sal.bf.sig$rad.loc %in% pa.out.radloc,])
dim(grass.bf.sig[grass.bf.sig$rad.loc %in% pa.out.radloc,])

temp.lositan<-temp.bf.sig[temp.bf.sig$rad.loc %in% lositan.sig$RAD.loc,]#15
dim(sal.bf.sig[sal.bf.sig$rad.loc %in% lositan.sig$RAD.loc,])
dim(grass.bf.sig[grass.bf.sig$rad.loc %in% lositan.sig$RAD.loc,])

temp.fst<-temp.bf.sig[temp.bf.sig$rad.loc %in% np.glob.fst.out$Locus,]#223
dim(sal.bf.sig[sal.bf.sig$rad.loc %in% np.glob.fst.out$Locus,])
dim(grass.bf.sig[grass.bf.sig$rad.loc %in% np.glob.fst.out$Locus,])

temp.xtx<-temp.bf.sig[temp.bf.sig$rad.loc %in% xtx.1$rad.loc,]#61
dim(sal.bf.sig[sal.bf.sig$rad.loc %in% xtx.1$rad.loc,])
dim(grass.bf.sig[grass.bf.sig$rad.loc %in% xtx.1$rad.loc,])

dim(temp.lositan[temp.lositan$rad.loc %in% temp.fst$rad.loc,])
dim(temp.lositan[temp.lositan$rad.loc %in% temp.xtx$rad.loc,])
dim(temp.fst[temp.fst$rad.loc %in% temp.xtx$rad.loc,])
#************************************PLOT*************************************#
bf.plot<-bf.scaff
bf.plot$Temp_BF<-log10(bf.plot$Temp_BF)
all.scaff<-split(bf.plot, bf.scaff$scaffold)
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

#jpeg("bf.env.assoc.jpeg", width=24, height=12, units="in", res=300)
jpeg("bf.env.assoc.small.jpeg", width=14, height=4.5, units="in", res=300)
par(mar=c(0,3,1,0), oma=c(1,3,0,0),cex=2)
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
		y.max,x.max, rect.xs[i,],y.min=0,x.min=0, pt.col="gray",
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in% temp.bf.sig$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col="gray32", pch=19, cex=0.5)
	#mark ones also fsts
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in%
		temp.fst$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col="red", pch=1, cex=0.75)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in%
		temp.lositan$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col="blue", pch=0, cex=0.75)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in%
		temp.xtx$locus,]
	points(temp.sig$BP, temp.sig$Temp_BF, col="purple", pch=6, cex=0.75)
	
}
axis(2, at = seq(y.min,y.max+0.00002,round((y.max-y.min)/2, digits=8)),
	ylim = c(y.min, y.max), pos=0,
	las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
#axis(1, at=c(0,rect.xs[,2]), xlim=c(0,(rect.xs[nrow(rect.xs),2]+5)), 
#	ylim=c(0,y.max), las=1, tck=-0.01, labels=FALSE, pos=0, xpd=TRUE,
#	par("usr")[3]-0.00051,srt = 38, adj = 1, xpd = TRUE,cex=0.75)
legend("top", c("Top 5% BF", "Fst outliers", "Lositan Outliers", "XtX 1%"), 
	col=c("gray32","red", "blue", "purple"),
	pch=c(19,1,0,6), pt.cex=c(0.5,0.75,0.75,0.75), bty='n',
	horiz=TRUE, cex=0.85)
mtext(side=1, "Genomic Location", outer = TRUE, cex=1.5)
mtext(side=2, "log10(BF)", line=1,outer=TRUE,cex=1.5)
#clip(x.min,x.max,y.min,y.max)
#abline(h=np.fst.ci.99)
dev.off()



#compare to each other
temp.sig<-temp.bf.sig[temp.bf.sig$locus %in% ctemp.bf.sig$locus,]#237
sal.sig<-sal.bf.sig[sal.bf.sig$locus %in% csal.bf.sig$locus,]#275
temp.sal<-temp.sig[temp.sig$locus %in% sal.sig$locus,]#127
all.sig<-temp.sal[temp.sal$locus %in% grass.bf.sig$locus,]#120

#remove ones significant in all of them
temp.sig<-temp.sig[!(temp.sig$locus %in% all.sig$locus),]#117
sal.sig<-sal.sig[!(sal.sig$locus %in% all.sig$locus),]#91
temp.sal.sig<-temp.sal[!(temp.sal$locus %in% all.sig$locus),]#7
grass.sig<-grass.bf.sig[!(grass.bf.sig$locus %in% all.sig$locus),]#494


########################Population differentiation##########################
#used all 6348 SNPs found in 12 populations
xtx<-read.table("E://ubuntushare//bayenv2//XtX_out.env_data_bayenv.txt")
colnames(xtx)<-c("locus", "XtX")
xtx$locus<-sub('./snpfiles/all.6348.(\\d+)','\\1', xtx$locus)
xtx$rad.loc<-sub('(\\d+)(_\\d+)','\\1',xtx$locus)
xtx.scaff<-merge(all.6348.map, xtx, by.x="V2", by.y="locus")
colnames(xtx.scaff)<-c("locus","scaffold","dist","BP","XtX", "rad.loc")
xtx.1<-xtx.scaff[xtx.scaff$XtX >= quantile(xtx.scaff$XtX,0.99),]
xtx.fst<-xtx.1[np.glob.fst.out$Locus %in% xtx.1$rad.loc,]

#************************************PLOT*************************************#
#plot the anypops one (more background? stronger signal?)
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

#jpeg("xtx.jpeg", width=24, height=12, units="in", res=300)
jpeg("xtx.small.jpeg", width=14, height=4.5, units="in", res=300)
par(mar=c(0,2,1,0), oma=c(1,2,0,0),cex=2)
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
		y.max,x.max, rect.xs[i,],y.min=0,x.min=0, 
		plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
	temp.sig<-all.scaff[[i]][all.scaff[[i]]$locus %in% xtx.1$locus,]
	points(temp.sig$BP, temp.sig$XtX, col="red", pch=19, cex=0.5)
	#mark ones also fsts
	temp.xtx.fst<-all.scaff[[i]][all.scaff[[i]]$locus %in%
		xtx.fst$locus,]
	points(temp.xtx.fst$BP, temp.xtx.fst$XtX, col="blue", pch=1, cex=0.75)
}
axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
	ylim = c(y.min, y.max), pos=0,
	las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
#axis(1, at=c(0,rect.xs[,2]), xlim=c(0,(rect.xs[nrow(rect.xs),2]+5)), 
#	ylim=c(0,y.max), las=1, tck=-0.01, labels=FALSE, pos=0, xpd=TRUE,
#	par("usr")[3]-0.00051,srt = 38, adj = 1, xpd = TRUE,cex=0.75)
legend("top", c("XtX outliers", "Fst outliers"), col=c("red", "blue"),
	pch=c(19,1), pt.cex=c(0.5,0.75), bty='n', horiz=TRUE, cex=0.85)
mtext(side=1, "Genomic Location", outer = TRUE, line=0,cex=1.5)
mtext(side=2, "XtX", outer=TRUE,line=0,cex=1.5)
#clip(x.min,x.max,y.min,y.max)
#abline(h=np.fst.ci.99)
dev.off()


#########################################################################
#***********************************************************************#
##############################LOSITAN####################################
#***********************************************************************#
#########################################################################
lositan.loci<-read.delim(
	"E://ubuntushare//pop_structure//lositan//ld.hwe.lositan.loci", header=T)
lositan.ci<-read.delim(
	"E://ubuntushare//pop_structure//lositan//ld.hwe.lositan.ci", 
	header=T, sep='\t')

plot_ci <- function(df.ci, bcolor, mcolor, tcolor) {
  lines(df.ci[,1],df.ci[,2], type='l', col=bcolor)#low
  lines(df.ci[,1],df.ci[,3], type='l', col=mcolor, lty=3)#med
  lines(df.ci[,1],df.ci[,4], type='l', col=tcolor)#upp
}

plot_loci <- function(df.loc, color) {
  points(df.loc[,2],df.loc[,3], col=color, pch=19)
}

lositan.loci<-cbind(lositan.loci, sub('(\\d+)(_\\d+)','\\1',lositan.loci[,1]))
colnames(lositan.loci)[5]<-"RAD.loc"

lositan.sig<-lositan.loci[(lositan.loci[,4] <= 0.05) | 
	(lositan.loci[,4]>=0.95),]
#comparisons
lositan.allfsts<-lositan.sig[lositan.sig$RAD.loc 
	%in% glob.fst.out$Locus,]
lositan.pcadapt<-lositan.sig[lositan.sig$RAD.loc %in% pa.out.radloc,]
pcadapt.allfsts<-glob.fst.out[glob.fst.out$Locus %in% pa.out.radloc,]
lositan.xtx<-lositan.sig[lositan.sig$RAD.loc %in% xtx.1$rad.loc,]
pcadapt.xtx<-xtx.1[xtx.1$rad.loc %in% pa.out.radloc,]
allfsts.xtx<-glob.fst.out[glob.fst.out$Locus %in% xtx.1$rad.loc,]
lositan.temp<-lositan.sig[lositan.sig$RAD.loc %in% temp.bf.sig$rad.loc,]


#jpeg("Fst_het_lositan.jpeg", res=300, height=7, width=14, units="in")
jpeg("Fst_het_lositan.small.jpeg", res=300, height=4.5, width=14, units="in")
plot(-10,ylim=c(0,0.25),xlim=c(0,1),axes=F, xlab="", ylab="", las=1,cex=2)
plot_ci(lositan.ci, 'black', 'gray', 'black')
plot_loci(lositan.loci, 'black')
points(lositan.sig$Het, lositan.sig$Fst, cex=1, col="red", pch=8)
points(lositan.allfsts$Het, lositan.allfsts$Fst, cex=3, col="green")
points(lositan.pcadapt$Het, lositan.pcadapt$Fst, cex=3, col="purple", pch=6)
points(lositan.temp$Het, lositan.temp$Fst, cex=3, col="blue", pch=0)
legend("top",bty="n",ncol=1,cex=1.5,
	c("99% Quantile", "Median"),lty=c(1,3), col=c("black", "gray"))
legend("topright", bty="n",ncol=2,cex=1.5,
	c("Lositan", "Lositan & Fst Out", "Lositan & PCAdapt", "Lositan & Temp"),
	pch=c(8,1,6,0), col=c("red", "green", "purple", "blue"))
axis(1, las=1,tck = -0.05, cex.axis=1.5,pos=0)
axis(2, las=1,tck = -0.05,cex.axis=1.5,pos=0)
mtext("He",1,outer=T,cex=1.5,line=-2)
mtext("Fst",2,outer=T,cex=1.5,line=-2)
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




#############################################################################
#****************************************************************************#
############################MALE-FEMALE LOCI#################################	
#****************************************************************************#
#############################################################################
#use populations and calc_fst_cis program output
m.f.summ<-read.table("E://ubuntushare//stacks//populations_sex//batch_1.sumstats.tsv",
	sep='\t', skip=2, header=T, comment.char="")
m.f.fst<-read.delim("E://ubuntushare//stacks//populations_sex//batch_1.fst_female-male.tsv")
mf.fst.sum<-read.delim("E://ubuntushare//stacks//populations_sex//fst_female-male.txt")
mf.fst.ci<-read.delim("E://ubuntushare//stacks//populations_sex//fst_female-male_summary.txt")

sex.outliers<-mf.fst.sum[mf.fst.sum$SmoothFst >= mf.fst.ci$CI99smooth,]
dim(sex.outliers[sex.outliers$Locus %in% np.glob.fst.out$Locus,])

sex.outliers<-sex.outliers[!(sex.outliers$Locus %in% np.glob.fst.out$Locus),]

length(levels(factor(sex.outliers$Locus)))
length(levels(factor(sex.outliers$Chrom)))
length(levels(as.factor(
	contigs[contigs$Scaffold %in% sex.outliers$Chrom,1]))) #num lgs

jpeg("E://Docs//PopGen//male-female.jpeg", res=300, height=480*2, width=480*4)
par(mar=c(2,4,2,1), oma=c(2,4,2,1), las=1)
plot.fsts.scaffs(mf.fst.sum, "SmoothFst", mf.fst.ci)
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

#***********************Using PLINK's assoc results**************************#
#####***REWRITE FILES***#####
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
#### RUN PLINK!! (treat male/female as case/control) ####

sex.assoc<-as.data.frame(read.table(
	"E://ubuntushare//stacks//populations_sex//sex.plink.assoc",
	header=T))

sex.assoc.sig<-sex.assoc[sex.assoc$P<=0.05,]
sex.map<-read.table("E://ubuntushare//stacks//populations_sex//batch_1.plink.map")

sex.map.sig<-sex.map[sex.map$V2 %in% sex.assoc.sig$SNP,]
sex.map.sig<-sex.map.sig[match(sex.map.sig$V2, sex.assoc.sig$SNP),]
sex.sig<-merge(sex.assoc.sig, sex.map.sig, by.x="SNP", by.y="V2")
sex.sig<-cbind(sex.sig, as.numeric(sub('_\\d+','',sex.sig$SNP)))
colnames(sex.sig)[11:14]<-c("SCAFFOLD","CHROM","BPN","RADLOC")

sex.out.both<-sex.sig[sex.sig$RADLOC %in% sex.outliers$Locus,]#104
dim(sex.sig[sex.sig$RADLOC %in% np.glob.fst.out$Locus,])#91
dim(sex.out.both[sex.out.both$RADLOC %in% np.glob.fst.out$Locus,])#4

#sex.out.assoc<-merge(sex.sig, sex.outliers, by.x="RADLOC", by.y="Locus")

#linkage map comparison (unnecessary)
dim(sex.out.assoc[sex.out.assoc$SCAFFOLD %in% use.contigs$Scaffold,])
summary(as.factor(
	merge(sex.out.assoc, use.contigs, by.x="SCAFFOLD", by.y="Scaffold")$LG))
#none of the RAD loci are on the linkage map, but some of the scaffolds map back

####PLOT####
jpeg("male.female.fst.jpeg",height=4,width=8, units="in",res=300)
par(mar=c(0,3,1,1), oma=c(0,3,1,1), las=1)
mf.plot<-plot.fsts.scaffs(mf.fst.sum, "SmoothFst", pt.cex=1,y.lab="")
points(mf.plot[mf.plot$locus %in% sex.outliers$Locus,"new.bp"],
	sex.outliers$SmoothFst, col="red", pch=19)
points(mf.plot[mf.plot$locus %in% sex.out.both$RADLOC,"new.bp"],
	mf.fst.sum[mf.fst.sum$Locus %in% sex.out.both$RADLOC,"SmoothFst"], 
	col="blue", cex=1.5)
legend("topright", pch=c(19,19,1), col=c("black","red", "blue"),
	c("non-significant", "Fst outliers", "PLINK association"),cex=0.75,
	bty="n", horiz=TRUE)
mtext(side=1,"Genomic Location", outer = FALSE, cex=0.75)
mtext("Smoothed FST value", side=2, outer=FALSE, las=0, line=4, cex=0.75)
dev.off()


##############################################################################
#****************************************************************************#
####################################PST-FST###################################
#****************************************************************************#
##############################################################################
#************create male and female files for pst analysis******************#
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
ped<-read.table("E://ubuntushare//stacks//populations//batch_1.plink.ped", 
	skip = 1, stringsAsFactors=F, colClasses="character")
ped.names<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', ped[,2])
ped.names<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', ped.names)
pheno.dat$ID<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', pheno.dat$ID)
pops.pheno<-pheno.dat[pheno.dat$ID %in% ped.names,] #this does not have the juveniles
pops.pheno<-pops.pheno[match(ped.names,pops.pheno$ID),]
pops.pheno<-replace(pops.pheno, is.na(pops.pheno),as.numeric(-9))
pops.pheno<-pops.pheno[,-8]
pops.pheno<-cbind(substr(pops.pheno$ID,1,4), pops.pheno)
colnames(pops.pheno)[1]<-"PopID"
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

#******PLOT TRAIT SUMMARIES********#
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



##################################CALCULATE PST###########################
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

##########read in the files I wrote########
mal.pheno<-read.table(header=T,
	"E://ubuntushare//calculate_pairwise_pst//calculate_pairwise_pst//mal.pheno.txt")
fem.pheno<-read.table(header=T,
	"E://ubuntushare//calculate_pairwise_pst//calculate_pairwise_pst//fem.pheno.txt")

#reorder them to match dist file
fem.pheno.sep<-split(fem.pheno, fem.pheno$PopID)
fem.pheno.new<-rbind(fem.pheno.sep$TXSP,fem.pheno.sep$TXCC,fem.pheno.sep$TXCB,
	fem.pheno.sep$ALST,fem.pheno.sep$FLSG,fem.pheno.sep$FLKB,
	fem.pheno.sep$FLFD,fem.pheno.sep$FLSI,fem.pheno.sep$FLAB,
	fem.pheno.sep$FLPB,fem.pheno.sep$FLHB,fem.pheno.sep$FLCC)
mal.pheno.sep<-split(mal.pheno, mal.pheno$PopID)
mal.pheno.new<-rbind(mal.pheno.sep$TXSP,mal.pheno.sep$TXCC,mal.pheno.sep$TXCB,
	mal.pheno.sep$ALST,mal.pheno.sep$FLSG,mal.pheno.sep$FLKB,
	mal.pheno.sep$FLFD,mal.pheno.sep$FLSI,mal.pheno.sep$FLAB,
	mal.pheno.sep$FLPB,mal.pheno.sep$FLHB,mal.pheno.sep$FLCC)

pf.ba<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "MBandArea")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.bn<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "BandNum")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.sv<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "SVL")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.ln<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "std.length")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.dp<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "depth")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.sl<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "SnoutLength")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.sd<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "SountDepth")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pf.hl<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno.new[,c("PopID", "HeadLength")]))),
	as.dist(t(pwise.fst)), nrepet=9999)

pm.sv<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno.new[,c("PopID", "SVL")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pm.ln<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno.new[,c("PopID", "std.length")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pm.dp<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno.new[,c("PopID", "depth")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pm.sl<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno.new[,c("PopID", "SnoutLength")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pm.sd<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno.new[,c("PopID", "SountDepth")]))),
	as.dist(t(pwise.fst)), nrepet=9999)
pm.hl<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno.new[,c("PopID", "HeadLength")]))),
	as.dist(t(pwise.fst)), nrepet=9999)

########USE DATA FROM LESS STRINGENT POPULATIONS RUN (NO POPS LIMIT)##########
pwise.fst.more<-read.table("E:\\ubuntushare\\stacks\\populations_nopopslimit\\batch_1.fst_summary.tsv",
	header=T, row.names=1, sep='\t')
	pwise.fst.more<-rbind(pwise.fst.more,rep(NA, ncol(pwise.fst.more)))
	rownames(pwise.fst.more)<-colnames(pwise.fst.more)

f.more.ba<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "MBandArea")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.bn<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "BandNum")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.sv<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SVL")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.ln<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "std.length")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.dp<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "depth")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.sl<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SnoutLength")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.sd<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "SountDepth")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
f.more.hl<-mantel.rtest(as.dist(t(pairwise.pst(fem.pheno[,c("PopID", "HeadLength")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)

m.more.sv<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "SVL")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
m.more.ln<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "std.length")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
m.more.dp<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "depth")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
m.more.sl<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "SnoutLength")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
m.more.sd<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "SountDepth")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)
m.more.hl<-mantel.rtest(as.dist(t(pairwise.pst(mal.pheno[,c("PopID", "HeadLength")]))),
	as.dist(t(pwise.fst.more)), nrepet=9999)


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

#################PCA, PLOTTING#############################
bands.split<-split(fem.pheno[,c(1,2,9,10)],fem.pheno$PopID)
fem.pheno.split<-split(fem.pheno[,1:8],fem.pheno$PopID)
mal.pheno.split<-split(mal.pheno,mal.pheno$PopID)
#pca per pop
band.pca<-lapply(bands.split, function(x){ rda(x[,3:4]) })
fem.pheno.pca<-lapply(fem.pheno.split, function(x){ rda(x[,3:8]) })
mal.pheno.pca<-lapply(mal.pheno.split, function(x){ rda(x[,3:8]) })

#extract eigenvalue
band.eig<-NULL
band.u<-NULL
for(i in 1:length(band.pca)){
	band.eig<-rbind(band.eig,band.pca[[i]]$CA$eig)
	rownames(band.eig)[i]<-names(band.pca)[i]
	band.u<-as.data.frame(rbind(band.u,as.data.frame(
		cbind(as.character(bands.split[[i]]$PopID),band.pca[[i]]$CA$u[,1]))))
}
band.dist<-as.matrix(as.dist(band.eig[,1]))#ignore warnings
band.u.sep<-split(band.u, band.u[,1])
band.u.new<-rbind(band.u.sep$TXSP,band.u.sep$TXCC,band.u.sep$TXCB,
	band.u.sep$ALST,band.u.sep$FLSG,band.u.sep$FLKB,
	band.u.sep$FLFD,band.u.sep$FLSI,band.u.sep$FLAB,
	band.u.sep$FLPB,band.u.sep$FLHB,band.u.sep$FLCC)
band.u.new$V2<-as.numeric(band.u.new$V2)

fem.pheno.eig<-NULL
fem.pheno.u<-NULL
for(i in 1:length(fem.pheno.pca)){
	fem.pheno.eig<-rbind(fem.pheno.eig,fem.pheno.pca[[i]]$CA$eig)
	rownames(fem.pheno.eig)[i]<-names(fem.pheno.pca)[i]
	fem.pheno.u<-as.data.frame(rbind(fem.pheno.u,as.data.frame(
		cbind(as.character(fem.pheno.split[[i]]$PopID),
		fem.pheno.pca[[i]]$CA$u[,1]))))

}
fem.dist<-as.matrix(as.dist(fem.pheno.eig[,1]))#ignore warnings
fem.u.sep<-split(fem.pheno.u, fem.pheno.u[,1])
fem.u.new<-rbind(fem.u.sep$TXSP,fem.u.sep$TXCC,fem.u.sep$TXCB,
	fem.u.sep$ALST,fem.u.sep$FLSG,fem.u.sep$FLKB,
	fem.u.sep$FLFD,fem.u.sep$FLSI,fem.u.sep$FLAB,
	fem.u.sep$FLPB,fem.u.sep$FLHB,fem.u.sep$FLCC)
fem.u.new$V2<-as.numeric(fem.u.new$V2)

mal.pheno.eig<-NULL
mal.pheno.u<-NULL
for(i in 1:length(mal.pheno.pca)){
	mal.pheno.eig<-rbind(mal.pheno.eig,mal.pheno.pca[[i]]$CA$eig)
	rownames(mal.pheno.eig)[i]<-names(mal.pheno.pca)[i]
	mal.pheno.u<-as.data.frame(rbind(mal.pheno.u,as.data.frame(
		cbind(as.character(mal.pheno.split[[i]]$PopID),
		mal.pheno.pca[[i]]$CA$u[,1]))))
}
mal.dist<-as.matrix(as.dist(mal.pheno.eig[,1]))#ignore warnings
mal.u.sep<-split(mal.pheno.u, mal.pheno.u[,1])
mal.u<-rbind(mal.u.sep$TXSP,mal.u.sep$TXCC,mal.u.sep$TXCB,
	mal.u.sep$ALST,mal.u.sep$FLSG,mal.u.sep$FLKB,
	mal.u.sep$FLFD,mal.u.sep$FLSI,mal.u.sep$FLAB,
	mal.u.sep$FLPB,mal.u.sep$FLHB,mal.u.sep$FLCC)
mal.u$V2<-as.numeric(mal.u$V2)

mantel.rtest(as.dist(t(pairwise.pst(data.frame(cbind(rownames(band.eig),
	as.numeric(band.eig[,1])))))),
	as.dist(t(dist)), nrepet=9999)

mantel.rtest(band.dist,as.dist(t(dist)))

env.u<-rda(t(env.dat))$CA$u
env.u.new<-env.u[match(pop.list,rownames(env.u)),1]
env.dist<-dist(env.u.new)

band.pst<-pairwise.pst(band.u.new)
fem.pst<-pairwise.pst(fem.u.new)
mal.pst<-pairwise.pst(mal.u)

###PLOT###
jpeg("pst.fst.dist.jpeg",height=7,width=7, units="in", res=300)
par(las=1, oma=c(1,1,1,1), mar=c(3,3,1,3))
plot(dist[upper.tri(dist)], pwise.fst[upper.tri(pwise.fst)], pch=19,
	ylim=c(0,1),xlab="",ylab="")
points(dist[upper.tri(dist)],band.pst[upper.tri(band.pst)], pch=6,col="darkgreen")
points(dist[upper.tri(dist)],fem.pst[upper.tri(fem.pst)],pch=4,col="red")
points(dist[upper.tri(dist)],mal.pst[upper.tri(mal.pst)],pch=5,col="blue")
#points(dist[upper.tri(dist)],env.dist[lower.tri(env.dist)],pch=15,col="purple")
legend("topright", col=c("black","blue","red","darkgreen"),pch=c(19,5,4,6),
	c("Fst","Male PCA Pst", "Female PCA Pst", "Female Bands PCA Pst"))
axis(4)
mtext("Distance (miles)",1, outer=T, line=-0.5)
mtext("Smoothed Pairwise Fst",2, las=0, outer=T, line=-0.5)
mtext("Pairwise Pst",4, outer=T, las=0,line=-0.5)
dev.off()

##############################################################################
#****************************************************************************#
################################P-MATRIX######################################
#****************************************************************************#
##############################################################################

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
fem.phen.uns.rda<-rda(fem.pheno.new[,3:10])
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
fem.pheno.new$PopID<-factor(fem.pheno.new$PopID)
pops.fem.unstd.dat<-split(fem.pheno.new, fem.pheno.new$PopID)
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
mal.pheno.new$PopID<-factor(mal.pheno.new$PopID)
pops.mal.unstd.dat<-split(mal.pheno.new, mal.pheno.new$PopID)
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
#females
pmat.f.u.p<-pmat.fem.unst.pops[match(pop.list, names(pmat.fem.unst.pops))]
fem.sim.unstd.mat<-generate.sim.mat(pmat.f.u.p)
fem.sim.unstd.fst<-mantel.rtest(as.dist(fem.sim.unstd.mat), 
	as.dist(t(pwise.fst)),	nrepet=9999)
pmax.unstd.mat<-vector.correlations(pmat.fem.unst.pops)
blows.sim.unstd.mat<-fem.sim.unstd.mat
blows.sim.unstd.mat[upper.tri(blows.sim.unstd.mat)]<-
	pmax.unstd.mat[upper.tri(pmax.unstd.mat)]
write.table(blows.sim.unstd.mat, "unstd_p_sim_matrix.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)
#males
pmat.m.u.p<-pmat.mal.unst.pops[match(pop.list, names(pmat.mal.unst.pops))]
mal.sim.unstd.mat<-generate.sim.mat(pmat.m.u.p)
mal.sim.unstd.fst<-mantel.rtest(as.dist(mal.sim.unstd.mat), 
	as.dist(t(pwise.fst)),	nrepet=9999)
pmax.mal.unstd.mat<-vector.correlations(pmat.mal.unst.pops)
blows.mal.sim.unstd.mat<-mal.sim.unstd.mat
blows.mal.sim.unstd.mat[upper.tri(blows.mal.sim.unstd.mat)]<-
	pmax.mal.unstd.mat[upper.tri(pmax.mal.unstd.mat)]
write.table(blows.mal.sim.unstd.mat, "unstd_p_sim_matrix_male.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)
#standardized
fem.sim.std.mat<-generate.sim.mat(pmat.fem.std)
fem.sim.std.fst<-mantel.rtest(as.dist(fem.sim.std.mat), as.dist(t(pwise.fst)),
	nrepet=9999)#unstandardized
pmax.std.mat<-vector.correlations(pmat.fem.std)
blows.sim.std.mat<-fem.sim.std.mat
blows.sim.std.mat[upper.tri(blows.sim.std.mat)]<-
	pmax.std.mat[upper.tri(pmax.std.mat)]
write.table(blows.sim.std.mat, "standardized_p_sim_matrix.txt", col.names=T,
	row.names=T, eol='\n', sep='\t', quote=F)

#######PLOT##########
make a heatmap, males above the diagonal, females below
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



