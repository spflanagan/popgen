#Author: Sarah P. Flanagan
#Date: 2 April 2015
#Purpose: Analyze fst cline data
#Look for signs of differentiation between north and south
#but not south-south or north-north.

rm(list=ls())

library(ade4)
library(lme4)
library(maps)
library(mapdata)
library(vegan)
library(boot)
library(car)
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
*****************************************************************************#
################################CLINE ANALYSIS################################
#****************************************************************************#
##############################################################################
#using data from the C++ program fst_geographical_clines

#clines: txsp.alst,flsg.flsi, flab.flcc 
#north-north: alst.flcc, alst.flsg,flcc.flsg
#south-south: txsp.flab, txsp.flsi,flab.flsi
sig.in.all<-sig.diff[complete.cases(sig.diff),]
sig.all.fst<-cbind(sig.in.all[,1:3], 
	sig.in.all[,grep("_fst", names(sig.in.all))])
sig.all.smf<-cbind(sig.in.all[,1:3], 
	sig.in.all[,grep("_smooth", names(sig.in.all))])

sig.all.cline<-sig.all.fst[apply(#outer apply selects 'all' that are true
	apply(sig.all.fst[7:24], 1, ">", 0),#every cell is T/F >0
	2, all),]#from the columns (apply transposes df)
write.table(sig.all.cline$Locus, row.names=F, col.names=F,
	"E://ubuntushare//fst_geographical_clines//fst_geographical_clines//sig_loci.txt", sep='\t')
write.table(sig.all.cline, row.names=F,col.names=T,
	"E://ubuntushare//fst_geographical_clines//fst_geographical_clines//sig_fsts_out.txt", sep='\t')

	
##############################################################################
*****************************************************************************#
######################CLINE ANALYSIS: POPULATIONS OUTPUT######################
#****************************************************************************#
##############################################################################
#USE DATA FROM CLINE-CLASSIFIED POPULATIONS
setwd("E://ubuntushare//stacks//populations_clines")
clines.files<-list.files(pattern=paste("fst_[[:alpha:]]+", "-",
	"[[:alpha:]].*.txt", sep=""))
clines.summ<-clines.files[grep("_summary",clines.files)]
clines.files<-clines.files[grep("_summary",clines.files, invert=TRUE)]

northmid.north<-read.delim("fst_north-mid-north.txt", header=T)
northmid.southmid<-read.delim("fst_north-mid-south-mid.txt", header=T)
northmid.south<-read.delim("fst_north-mid-south.txt", header=T)
southmid.north<-read.delim("fst_south-mid-north.txt", header=T)
southmid.south<-read.delim("fst_south-mid-south.txt", header=T)
south.north<-read.delim("fst_south-north.txt", header=T)

nm.n.ci<-read.delim("fst_north-mid-north_summary.txt", header=T)
nm.sm.ci<-read.delim("fst_north-mid-south-mid_summary.txt", header=T)
nm.s.ci<-read.delim("fst_north-mid-south_summary.txt", header=T)
sm.n.ci<-read.delim("fst_south-mid-north_summary.txt", header=T)
sm.s.ci<-read.delim("fst_south-mid-south_summary.txt", header=T)
s.n.ci<-read.delim("fst_south-north_summary.txt", header=T)
#UNSURE WHAT TO DO WITH THESE

cline.outliers<-data.frame()
for(i in 1:length(clines.files))
{
	fst.file<-read.delim(paste(
		"E://ubuntushare//stacks//populations_clines//",
		clines.files[i], sep=""), header=T)
	summ.file<-read.delim(paste(
		"E://ubuntushare//stacks//populations_clines//",
		clines.summ[i], sep=""), header=T)
	cline.outliers<-rbind(cline.outliers,
		fst.file[fst.file$SmoothFst >= summ.file$CI99smooth,])
}

cline.no.dups.outliers<-cline.outliers[!duplicated(cline.outliers$Locus),]
length(levels(factor(cline.no.dups.outliers$Chrom))) #number scaffolds
length(levels(as.factor(
	contigs[contigs$Scaffold %in% cline.no.dups.outliers$Chrom,1]))) #num lgs





##############################################################################
*****************************************************************************#
###########################CLINE FST OUTLIER ANALYSIS#########################
#****************************************************************************#
##############################################################################
setwd("E://ubuntushare//stacks//populations")
ALST.FLCC<-read.delim("fst_ALST-FLCC.txt", header=T)#north
ALST.FLSG<-read.delim("fst_ALST-FLSG.txt", header=T)#north
FLAB.FLCC<-read.delim("fst_FLAB-FLCC.txt", header=T)#cline
FLAB.FLSI<-read.delim("fst_FLAB-FLSI.txt", header=T)#south
FLCC.FLSG<-read.delim("fst_FLCC-FLSG.txt", header=T)#north
FLSG.FLSI<-read.delim("fst_FLSG-FLSI.txt", header=T)#cline
TXSP.ALST<-read.delim("fst_TXSP-ALST.txt", header=T)#cline
TXSP.FLAB<-read.delim("fst_TXSP-FLAB.txt", header=T)#south
TXSP.FLSI<-read.delim("fst_TXSP-FLSI.txt", header=T)#south

ALST.FLCC.ci<-read.delim("fst_ALST-FLCC_summary.txt", header=T)#north
ALST.FLSG.ci<-read.delim("fst_ALST-FLSG_summary.txt", header=T)#north
FLAB.FLCC.ci<-read.delim("fst_FLAB-FLCC_summary.txt", header=T)#cline
FLAB.FLSI.ci<-read.delim("fst_FLAB-FLSI_summary.txt", header=T)#south
FLCC.FLSG.ci<-read.delim("fst_FLCC-FLSG_summary.txt", header=T)#north
FLSG.FLSI.ci<-read.delim("fst_FLSG-FLSI_summary.txt", header=T)#cline
TXSP.ALST.ci<-read.delim("fst_TXSP-ALST_summary.txt", header=T)#cline
TXSP.FLAB.ci<-read.delim("fst_TXSP-FLAB_summary.txt", header=T)#south
TXSP.FLSI.ci<-read.delim("fst_TXSP-FLSI_summary.txt", header=T)#south


jpeg("E://Docs//PopGen//clines_fst.jpeg", width=480*9, height=480*6, res=300)
par(mfrow=c(3,3), mar=c(0,10,1,0), oma=c(4,10,0,0))
plot.fsts.scaffs(TXSP.ALST, "Texas Gulf \nCline", TXSP.ALST.ci)
plot.fsts.scaffs(ALST.FLCC, "Texas Gulf & \nFlorida Atlantic \nNorth", ALST.FLCC.ci)
plot.fsts.scaffs(TXSP.FLAB, "Texas Gulf & \nFlorida Atlantic \nSouth", TXSP.FLAB.ci)
plot.fsts.scaffs(FLSG.FLSI, "Florida Gulf \nCline", FLSG.FLSI.ci)
plot.fsts.scaffs(ALST.FLSG, "Texas Gulf & \nFlorida Gulf \nNorth", ALST.FLSG.ci)
plot.fsts.scaffs(TXSP.FLSI, "Texas Gulf & \nFlorida Gulf \nSouth", TXSP.FLSI.ci)
plot.fsts.scaffs(FLAB.FLCC, "Florida Atlantic \nCline", FLAB.FLCC.ci)
plot.fsts.scaffs(FLCC.FLSG, "Florida Gulf & \nFlorida Atlantic \nNorth", FLCC.FLSG.ci)
plot.fsts.scaffs(FLAB.FLSI, "Florida Gulf & \nFlorida Atlantic \nSouth", FLAB.FLSI.ci)

mtext(side=1, "Genomic Location", outer = TRUE, line=1)
mtext(side=2, "Smoothed FST value", outer=TRUE, line=6)

dev.off()
#only problem with this: the scaffolds don't end up being the same length
#for each pop.
setwd("E://Docs//PopGen")
##IDENTIFY ALL OUTLIERS
#list.cline.99<-list(
txsp.alst.99<-	TXSP.ALST[TXSP.ALST$SmoothFst > TXSP.ALST.ci$CI99smooth,]#,
alst.flcc.99<-	ALST.FLCC[ALST.FLCC$SmoothFst > ALST.FLCC.ci$CI99smooth,]#,
txsp.flab.99<-	TXSP.FLAB[TXSP.FLAB$SmoothFst > TXSP.FLAB.ci$CI99smooth,]#,
flsg.flsi.99<-	FLSG.FLSI[FLSG.FLSI$SmoothFst > FLSG.FLSI.ci$CI99smooth,]#,
alst.flsg.99<-	ALST.FLSG[ALST.FLSG$SmoothFst > ALST.FLSG.ci$CI99smooth,]#,
txsp.flsi.99<-	TXSP.FLSI[TXSP.FLSI$SmoothFst > TXSP.FLSI.ci$CI99smooth,]#,
flab.flcc.99<-	FLAB.FLCC[FLAB.FLCC$SmoothFst > FLAB.FLCC.ci$CI99smooth,]#,
flcc.flsg.99<-	FLCC.FLSG[FLCC.FLSG$SmoothFst > FLCC.FLSG.ci$CI99smooth,]#,
flab.flsi.99<-	FLAB.FLSI[FLAB.FLSI$SmoothFst > FLAB.FLSI.ci$CI99smooth,]#)

all.loci<-rbind(txsp.alst.99[,1:3],alst.flcc.99[,1:3],txsp.flab.99[,1:3],
	flsg.flsi.99[,1:3],alst.flsg.99[,1:3],txsp.flsi.99[,1:3],
	flab.flcc.99[,1:3],flcc.flsg.99[,1:3],flab.flsi.99[,1:3])

all.loci<-all.loci[!duplicated(all.loci[,1:3]),]

all.cline.99<-merge(all.loci, txsp.alst.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"TXSP.ALST"
all.cline.99<-merge(all.cline.99, alst.flcc.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"ALST.FLCC"
all.cline.99<-merge(all.cline.99, txsp.flab.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"TXSP.FLAB"
all.cline.99<-merge(all.cline.99, flsg.flsi.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"FLSG.FLSI"
all.cline.99<-merge(all.cline.99, alst.flsg.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"ALST.FLSG"
all.cline.99<-merge(all.cline.99, txsp.flsi.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"TXSP.FLSI"
all.cline.99<-merge(all.cline.99, flab.flcc.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"FLAB.FLCC"
all.cline.99<-merge(all.cline.99, flcc.flsg.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"FLCC.FLSG"
all.cline.99<-merge(all.cline.99, flab.flsi.99[,-4], by=c("Locus","Chrom","BP"),all=T)
names(all.cline.99)[length(all.cline.99)]<-"FLAB.FLSI"

cline.outliers<-all.cline.99[complete.cases(all.cline.99[,c(4,7,10)]),]
north.outliers<-all.cline.99[complete.cases(all.cline.99[,c(5,8,11)]),]
south.outliers<-all.cline.99[complete.cases(all.cline.99[,c(6,9,12)]),]

########################20 MAY 2015 RETURN TO THE ANALYSIS####################
#none of the clines have the same exact SNPs, or even the same scaffolds
#but let's see if they have and small values in north-north
#c1=cline1=flab.flcc
#c2=cline2=flsg.flsi
#c3=cline3=txsp.alst
#n1=north1=alst.flcc
#n2=north2=alst.flsg
#n3=north3=flcc.flsg
#s1=south1=flab.flsi
#s2=south2=txsp.flab
#s3=south3=txsp.flsi

n1.c1<-ALST.FLCC[ALST.FLCC$Locus %in% flab.flcc.99$Locus,]
n1.c1.low<-n1.c1[n1.c1$SmoothFst < ALST.FLCC.ci$MeanSmoothFst,]
n1.c2<-ALST.FLCC[ALST.FLCC$Locus %in% flsg.flsi.99$Locus,]
n1.c2.low<-n1.c2[n1.c2$SmoothFst < ALST.FLCC.ci$MeanSmoothFst,]
n1.c3<-ALST.FLCC[ALST.FLCC$Locus %in% txsp.alst.99$Locus,]
n1.c3.low<-n1.c3[n1.c3$SmoothFst < ALST.FLCC.ci$MeanSmoothFst,]

n2.c1<-ALST.FLSG[ALST.FLSG$Locus %in% flab.flcc.99$Locus,]
n2.c1.low<-n2.c1[n2.c1$SmoothFst < ALST.FLSG.ci$MeanSmoothFst,]
n2.c2<-ALST.FLSG[ALST.FLSG$Locus %in% flsg.flsi.99$Locus,]
n2.c2.low<-n2.c2[n2.c2$SmoothFst < ALST.FLSG.ci$MeanSmoothFst,]
n2.c3<-ALST.FLSG[ALST.FLSG$Locus %in% txsp.alst.99$Locus,]
n2.c3.low<-n2.c3[n2.c3$SmoothFst < ALST.FLSG.ci$MeanSmoothFst,]

n3.c1<-FLCC.FLSG[FLCC.FLSG$Locus %in% flab.flcc.99$Locus,]
n3.c1.low<-n3.c1[n3.c1$SmoothFst < FLCC.FLSG.ci$MeanSmoothFst,]
n3.c2<-FLCC.FLSG[FLCC.FLSG$Locus %in% flsg.flsi.99$Locus,]
n3.c2.low<-n3.c2[n3.c2$SmoothFst < FLCC.FLSG.ci$MeanSmoothFst,]
n3.c3<-FLCC.FLSG[FLCC.FLSG$Locus %in% txsp.alst.99$Locus,]
n3.c3.low<-n3.c3[n3.c3$SmoothFst < FLCC.FLSG.ci$MeanSmoothFst,]

s1.c1<-FLAB.FLSI[FLAB.FLSI$Locus %in% flab.flcc.99$Locus,]
s1.c1.low<-s1.c1[s1.c1$SmoothFst < FLAB.FLSI.ci$MeanSmoothFst,]
s1.c2<-FLAB.FLSI[FLAB.FLSI$Locus %in% flsg.flsi.99$Locus,]
s1.c2.low<-s1.c2[s1.c2$SmoothFst < FLAB.FLSI.ci$MeanSmoothFst,]
s1.c3<-FLAB.FLSI[FLAB.FLSI$Locus %in% txsp.alst.99$Locus,]
s1.c3.low<-s1.c3[s1.c3$SmoothFst < FLAB.FLSI.ci$MeanSmoothFst,]

s2.c1<-TXSP.FLAB[TXSP.FLAB$Locus %in% flab.flcc.99$Locus,]
s2.c1.low<-s2.c1[s2.c1$SmoothFst < TXSP.FLAB.ci$MeanSmoothFst,]
s2.c2<-TXSP.FLAB[TXSP.FLAB$Locus %in% flsg.flsi.99$Locus,]
s2.c2.low<-s2.c2[s2.c2$SmoothFst < TXSP.FLAB.ci$MeanSmoothFst,]
s2.c3<-TXSP.FLAB[TXSP.FLAB$Locus %in% txsp.alst.99$Locus,]
s2.c3.low<-s2.c3[s2.c3$SmoothFst < TXSP.FLAB.ci$MeanSmoothFst,]

s3.c1<-TXSP.FLSI[TXSP.FLSI$Locus %in% flab.flcc.99$Locus,]
s3.c1.low<-s3.c1[s3.c1$SmoothFst < TXSP.FLSI.ci$MeanSmoothFst,]
s3.c2<-TXSP.FLSI[TXSP.FLSI$Locus %in% flsg.flsi.99$Locus,]
s3.c2.low<-s3.c2[s3.c2$SmoothFst < TXSP.FLSI.ci$MeanSmoothFst,]
s3.c3<-TXSP.FLSI[TXSP.FLSI$Locus %in% txsp.alst.99$Locus,]
s3.c3.low<-s3.c3[s3.c3$SmoothFst < TXSP.FLSI.ci$MeanSmoothFst,]

north.low.outliers<-rbind(n1.c2.low, n1.c2.low,n1.c3.low,n2.c1.low,
	n2.c2.low,n2.c3.low,n3.c1.low,n3.c2.low,n3.c3.low)
south.low.outliers<-rbind(s1.c2.low, s1.c2.low,s1.c3.low,s2.c1.low,
	s2.c2.low,s2.c3.low,s3.c1.low,s3.c2.low,s3.c3.low)
all.low.outliers<-rbind(north.low.outliers, south.low.outliers)

length(levels(as.factor(
	contigs[contigs$Scaffold %in% all.low.outliers$Chrom,1]))) #num lg
summary(as.factor(contigs[contigs$Scaffold %in% all.low.outliers$Chrom,1]))

#select those also found in non-outlier analysis
length(levels(as.factor(all.low.outliers[all.low.outliers$Locus 
	%in% pa.out.dat$Locus.ID,1]))) #29715 is a match
length(levels(as.factor(all.low.outliers[all.low.outliers$Chrom 
	%in% pa.out.dat$Chr,1])))
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

#############################################################################
###############EXAMINE ALLELE FREQ PER POP OF SIG. CLINE LOCI################
#############################################################################
sig.loci<-levels(as.factor(sig.all.cline$Locus))
sig.snp<-levels(as.factor(sig.all.cline$BP))
sig.dat<-summ.dat[summ.dat$Locus.ID %in% sig.loci,]
sig.dat<-sig.dat[sig.dat$BP %in% sig.snp,]
sig.split<-split(sig.dat, sig.dat$Locus.ID)
mar.coor$lon<-(-1*mar.coor$lon)
#adjust coordinates for ease of mapping
mar.coor[5,3]<-mar.coor[5,3]-1
mar.coor[6,2]<-mar.coor[6,2]-2
mar.coor[6,3]<-mar.coor[6,3]-1
mar.coor[7,2]<-mar.coor[7,2]-1.5#moves down
mar.coor[7,3]<-mar.coor[7,3]-2#moves left
mar.coor[8,2]<-mar.coor[8,2]-2
mar.coor[8,3]<-mar.coor[8,3]-1.5
mar.coor[9,2]<-mar.coor[9,2]-1.5
mar.coor[9,3]<-mar.coor[9,3]-1.5
mar.coor[10,2]<-mar.coor[10,2]-2
mar.coor[11,2]<-mar.coor[11,2]-1


plot.map.bg<-function(fill.col="gray98",border="grey"){
	map("worldHires", "usa",xlim=c(-100,-78), ylim=c(23,35), 
		col="gray98", fill=TRUE, border="grey")
	map("worldHires", "mexico",xlim=c(-100,-78), ylim=c(23,35), 
		col="gray98", fill=TRUE, add=TRUE, border="grey")
	map("worldHires", "cuba",xlim=c(-100,-78), ylim=c(23,35), 
		col="gray98", fill=TRUE, add=TRUE, border="grey")
}
setwd("E://Docs//PopGen")
jpeg("cline_allele_freqs_15.jpg", width=480*16, height=480*16, res=300)
par(mfrow=c(4,4),  mar = c(0,0,0,0), oma=c(0,0,0,0))
for(i in 1:length(sig.split)){
	
	plot.map.bg()
	if(length(split(sig.split[[i]],sig.split[[i]]$BP))>1)
	{
		plots<-split(sig.split[[i]],sig.split[[i]]$BP)
		colors<-seq(1,length(plots),1)
		for(j in 1:length(plots)){
			plots[[j]]<-plots[[j]][match(mar.coor$site, 
				plots[[j]]$Pop.ID),]
			text(mar.coor$lon+1, mar.coor$lat+0.5*j, labels=plots[[j]]$P, 
				col=colors[j], cex=0.75)
			text(mar.coor$lon-0.5, mar.coor$lat+0.5*j, labels=
				paste(plots[[j]]$P.Nuc, ":"), col=colors[j], cex=0.75)
			text(x=-89+2.1*j, 34, paste(".",plots[[j]]$BP), 
				col=colors[j], cex=1.5)
		}
		text(x=-89-1, 34, paste("Locus", plots[[1]]$Locus.ID),
			 col="black",cex=1.5)

	}else{
		sig.split[[i]]<-sig.split[[i]][match(mar.coor$site, 
			sig.split[[i]]$Pop.ID),]
		text(mar.coor$lon+1, mar.coor$lat+0.5, labels=sig.split[[i]]$P, 
			col="black", cex=1)
		text(mar.coor$lon-0.5, mar.coor$lat+0.5, labels=
			paste(sig.split[[i]]$P.Nuc, ":"), col="black", cex=1)
		text(x=-89, 34, paste("Locus", sig.split[[1]]$Locus.ID, ".", 
			sig.split[[i]]$BP), col="black", cex=1.5)
	}
}
dev.off()


sig.split[[2]]<-sig.split[[2]][match(mar.coor$site, sig.split[[2]]$Pop.ID),]
plot.map.bg()
text(mar.coor$lon+1, mar.coor$lat, labels=sig.split$`10031`$P, 
	col="dark green", cex=1)
text(mar.coor$lon, mar.coor$lat, labels=paste(sig.split$`10031`$P.Nuc, ":"), 
	col="dark green", cex=1)
text(x=-89, 34, paste("Locus", sig.split[[2]]$Locus.ID, ".", sig.split[[2]]$BP),
	col="dark green", cex=1.5)

sig.split[[3]]<-sig.split[[3]][match(mar.coor$site, sig.split[[3]]$Pop.ID),]
map("worldHires", "usa",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, border="grey")
map("worldHires", "mexico",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, add=TRUE, border="grey")
map("worldHires", "cuba",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, add=TRUE, border="grey")
text(mar.coor$lon+1, mar.coor$lat, labels=sig.split[[3]]$P, 
	col="purple", cex=1)
text(mar.coor$lon, mar.coor$lat, labels=paste(sig.split[[3]]$P.Nuc, ":"), 
	col="purple", cex=1)
text(x=-89, 34, paste("Locus", sig.split[[3]]$Locus.ID, ".", sig.split[[3]]$BP),
	col="purple", cex=1.5)

sig.split[[4]]<-sig.split[[4]][match(mar.coor$site, sig.split[[4]]$Pop.ID),]
map("worldHires", "usa",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, border="grey")
map("worldHires", "mexico",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, add=TRUE, border="grey")
map("worldHires", "cuba",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, add=TRUE, border="grey")
text(mar.coor$lon+1, mar.coor$lat, labels=sig.split[[4]]$P, 
	col="dark blue", cex=1)
text(mar.coor$lon, mar.coor$lat, labels=paste(sig.split[[4]]$P.Nuc, ":"), 
	col="dark blue", cex=1)
text(x=-89, 34, paste("Locus", sig.split[[4]]$Locus.ID, ".", sig.split[[4]]$BP),
	col="dark blue", cex=1.5)

sig.split[[5]]<-sig.split[[5]][match(mar.coor$site, sig.split[[5]]$Pop.ID),]
map("worldHires", "usa",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, border="grey")
map("worldHires", "mexico",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, add=TRUE, border="grey")
map("worldHires", "cuba",xlim=c(-100,-78), ylim=c(23,35), 
	col="gray98", fill=TRUE, add=TRUE, border="grey")
text(mar.coor$lon+1, mar.coor$lat, labels=sig.split[[5]]$P, 
	col="dark orange", cex=1)
text(mar.coor$lon, mar.coor$lat, labels=paste(sig.split[[5]]$P.Nuc, ":"), 
	col="dark orange", cex=1)
text(x=-89, 34, paste("Locus", sig.split[[5]]$Locus.ID, ".", sig.split[[5]]$BP),
	col="dark orange", cex=1.5)
dev.off()

##############################################################################
#****************************************************************************#
######################PLOT FST, FIS, PI USING LINKAGE MAP#####################
#****************************************************************************#
##############################################################################

summ.sub<-summ.dat[summ.dat$Chr %in% use.contigs$Scaffold,]
#Want to sort by chromosome in the order of use.contigs
chroms<-as.vector(summ.sub$Chr)
for(i in 1:nrow(use.contigs)){
	chroms[chroms %in% use.contigs[i,3]]<-use.contigs[i,1]
}

summ.sub<-cbind(summ.sub, chroms)#this should work!
summ.sub$chroms<-as.numeric(summ.sub$chroms)
summ.sub<-summ.sub[order(summ.sub$chroms),]

#####PREP DATA FOR PLOTTING#####
sub.chr<-split(summ.sub, summ.sub$chroms)
plot.chroms<-lapply(split(summ.sub, summ.sub$Pop.ID), function(x){split(x, x$chroms)})

last.max<-0
rect.xs<-NULL
addition.values<-0
for(i in 1:length(sub.chr)){
	new.max<-last.max+round(max(sub.chr[[i]]$BP), -2)
	rect.xs<-rbind(rect.xs,c(last.max, new.max))
	addition.values<-c(addition.values, new.max)
	last.max<-new.max
}

for(i in 1:length(plot.chroms)){
	for(j in 1:length(plot.chroms[[i]])){
		plot.chroms[[i]][[j]]$BP<-plot.chroms[[i]][[j]]$BP+addition.values[j]
	}
}

x.max<-max(addition.values)
y.max<-max(summ.sub$Smoothed.Pi)+0.1*max(summ.sub$Smoothed.Pi)#I don't like the way this looks on graph
y.max<-0.04
if(min(summ.sub$Smoothed.Pi) < 0) {
	y.min<-min(summ.sub$Smoothed.Pi) - 0.1*min(summ.sub$Smoothed.Pi)
} else {
	y.min<-0
}
pop.names<-c("Alabama Saltwater", "Florida Anne's Beach", 
	"Florida Cape Canaveral", "Florida Fort Desot","Florida Harbor Branch",
	"Florida Keaton Beach","Florida Palatka Bay", "Florida St. George",
	"Florida Sanibel Island", "Texas Christmas Bay", "Texas Corpus Christi",
	"Texas South Padre")
x.lab<-((rect.xs[,2]-rect.xs[,1])/2)+addition.values[1:21]

#PI
jpeg("marine_all_pi_lgs.jpg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,10,1,0), oma=c(4,10,0,0))

for(i in 1:length(plot.chroms)){
	plot(plot.chroms[[i]][[1]]$BP, plot.chroms[[i]][[1]]$Smoothed.Pi, 
		xlim=c(0,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(plot.chroms[[i]])){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(plot.chroms[[i]][[j]]$BP, 
			plot.chroms[[i]][[j]]$Smoothed.Pi,
			y.max,x.max, rect.xs[j,],y.min=y.min,x.min=0, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.75)
	mtext(pop.names[i],2,las=1, line=2)
}
par(mgp=c(3,0.1,0))
axis(1, at=c(0, x.max), labels=c("",""), lwd.ticks=0, pos=0)
axis(1, at=x.lab, labels=seq(1,21,1), ylab="", xlab="",  
	xlim=c(0,x.max), tck=-0.01, cex.axis=0.75, pos=0)
mtext(side=1, "Position on Linkage Group", outer = TRUE, line=1.25)
mtext(side=2, "Smoothed pi value for each population", outer=TRUE, line=6)
dev.off()

#UNSMOOTHED PI
y.max<-max(summ.sub$Pi)+0.1*max(summ.sub$Pi)#I don't like the way this looks on graph
if(min(summ.sub$Pi) < 0) {
	y.min<-min(summ.sub$Pi) - 0.1*min(summ.sub$Pi)
} else {
	y.min<-0
}
pop.names<-c("Alabama Saltwater", "Florida Anne's Beach", 
	"Florida Cape Canaveral", "Florida Fort Desot","Florida Harbor Branch",
	"Florida Keaton Beach","Florida Palatka Bay", "Florida St. George",
	"Florida Sanibel Island", "Texas Christmas Bay", "Texas Corpus Christi",
	"Texas South Padre")
x.lab<-((rect.xs[,2]-rect.xs[,1])/2)+addition.values[1:21]

#PI
jpeg("marine_all_pi_lgs_unsmooth.jpg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,10,1,0), oma=c(4,10,0,0))

for(i in 1:length(plot.chroms)){
	plot(plot.chroms[[i]][[1]]$BP, plot.chroms[[i]][[1]]$Pi, 
		xlim=c(0,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(plot.chroms[[i]])){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(plot.chroms[[i]][[j]]$BP, 
			plot.chroms[[i]][[j]]$Pi,
			y.max,x.max, rect.xs[j,],y.min=y.min,x.min=0, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.75)
	mtext(pop.names[i],2,las=1, line=2)
}
par(mgp=c(3,0.1,0))
axis(1, at=c(0, x.max), labels=c("",""), lwd.ticks=0, pos=0)
axis(1, at=x.lab, labels=seq(1,21,1), ylab="", xlab="",  
	xlim=c(0,x.max), tck=-0.01, cex.axis=0.75, pos=0)
mtext(side=1, "Position on Linkage Group", outer = TRUE, line=1.25)
mtext(side=2, "Pi value for each population", outer=TRUE, line=6)
dev.off()


#FIS
y.max<-max(summ.sub$Smoothed.Fis)+0.1*max(summ.sub$Smoothed.Fis)
if(min(summ.sub$Smoothed.Fis) < 0) {
	y.min<-min(summ.sub$Smoothed.Fis) + 0.1*min(summ.sub$Smoothed.Fis)
} else {
	y.min<-0
}

jpeg("marine_all_fis_lgs.jpg", width=480*4, height=480*8, res=300)
par(mfrow=c(12,1), mar=c(0,10,1,0), oma=c(4,10,0,0))

for(i in 1:length(plot.chroms)){
	plot(plot.chroms[[i]][[1]]$BP, plot.chroms[[i]][[1]]$Smoothed.Fis, 
		xlim=c(0,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(plot.chroms[[i]])){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(plot.chroms[[i]][[j]]$BP, 
			plot.chroms[[i]][[j]]$Smoothed.Fis,
			y.max,x.max, rect.xs[j,],y.min=y.min,x.min=0, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.75)
	mtext(pop.names[i],2,las=1, line=2)
}
par(mgp=c(3,0.1,0))
axis(1, at=c(0, x.max), labels=c("",""), lwd.ticks=0, pos=y.min)
axis(1, at=x.lab, labels=seq(1,21,1), ylab="", xlab="",  
	xlim=c(0,x.max), tck=-0.01, cex.axis=0.75, pos=y.min)
mtext(side=1, "Position on Linkage Group", outer = TRUE, line=1.25)
mtext(side=2, "Smoothed Fis value for each population", outer=TRUE, line=6)
dev.off()

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

#############################################################################
############################FIS = 1 WEIRDOS##################################	
#############################################################################
fis.outlier<-summ.dat[summ.dat$Smoothed.Fis > -6.68,]#3rd quantile = -6.679
fis.1<-summ.dat[summ.dat$Fis ==1,]
fis.1.loc<-levels(as.factor(fis.1$Locus.ID))
summary(fis.1$Pop.ID) #show distribution among 

fis.1.private<-fis.1[fis.1$Private==1,]
fis.1.s3<-fis.1[fis.1$Chr=="scaffold_3",]

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

scaff<-fis.1[fis.1$Chr %in% use.contigs$Scaffold,]
lg<-NULL
for(i in 1:nrow(scaff)){
	id<-levels(as.factor(
		use.contigs[use.contigs$Scaffold %in% factor(scaff$Chr[i]),1]))
	lg<-c(lg, id)
}
scaff<-cbind(scaff,lg)
fis.1.lg.4<-scaff[scaff$lg==4,]
#are individuals with fis=1 on lg 4 all males or all females??
#use male-female fst comparison
m.f.summ.dat<-read.table("G://Popgen//stacks//populations_sex//batch_1.sumstats.tsv",
	sep='\t', skip=2, header=T, comment.char="")

m.f.matches<-m.f.summ.dat[m.f.summ.dat$Locus.ID %in% scaff$Locus.ID,]
#so not m/f

#filter scaff
filt.fis1.scaff<-scaff[scaff$N>=35,]
filt.fis1.scaff<-filt.fis1.scaff[filt.fis1.scaff$P<=0.95,]
filt.fis1.scaff<-droplevels(filt.fis1.scaff)
#not sure what to do from here.

#############################################################################
############################P = 0.5 WEIRDOS##################################	
#############################################################################
p.5<-summ.dat[summ.dat$P == 0.5,]
p.5.f<-p.5[p.5$N>=35,]
p.5.loc<-levels(as.factor(p.5$Locus.ID))
summary(fis.1$Pop.ID) #show distribution among 

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

p.scaff<-p.5.f[p.5.f$Chr %in% use.contigs$Scaffold,]
lg<-NULL
for(i in 1:nrow(p.scaff)){
	id<-levels(as.factor(
		use.contigs[use.contigs$Scaffold %in% factor(p.scaff$Chr[i]),1]))
	lg<-c(lg, id)
}
p.scaff<-cbind(p.scaff,lg)
p.scaff<-droplevels(p.scaff)
summary(p.scaff)

#could THESE be sex linked??

p.mf.matches<-m.f.summ.dat[m.f.summ.dat$Locus.ID %in% p.scaff$Locus.ID,]
p.mf.matches.priv<-p.mf.matches[p.mf.matches$Private ==1,]
p.f.p.5<-p.mf.matches[p.mf.matches$Pop.ID=="female",]
p.f.p.5<-p.f.p.5[p.f.p.5$P==0.5,]
p.m.p.5<-p.mf.matches[p.mf.matches$Pop.ID=="male",]
p.m.p.5<-p.m.p.5[p.m.p.5$P==0.5,]

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

#############################################################################
############################CRAPPY OLD CLINE ANALYSIS########################
#############################################################################
sig.in.all<-sig.diff[complete.cases(sig.diff),]
sig.all.fst<-cbind(sig.in.all[,1:3], 
	sig.in.all[,grep("_fst", names(sig.in.all))])
sig.all.smf<-cbind(sig.in.all[,1:3], 
	sig.in.all[,grep("_smooth", names(sig.in.all))])

clines.all<-sig.all.fst[,1:6]
north.all<-sig.all.fst[,7:15]
north.all<-cbind(sig.all.fst[,1:3], north.all)
south.all<-sig.all.fst[,16:24]
south.all<-cbind(sig.all.fst[,1:3], south.all)
clines.diff<-sig.all.fst[abs(clines.all[,4:6])> abs(north.all[,4:12])
	& abs(clines.all[,4:6])> abs(south.all[,4:12]),]

clines.diff.1<-sig.all.fst[abs(sig.all.fst[,4])< abs(sig.all.fst[,7])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,8])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,9])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,10])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,11])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,12])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,12])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,13])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,14])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,15])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,16])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,17])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,18])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,19])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,21])
	& abs(sig.all.fst[,4])> abs(sig.all.fst[,20]),]

clines.diff.2<-sig.all.fst[abs(sig.all.fst[,5])< abs(sig.all.fst[,7])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,8])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,9])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,10])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,11])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,12])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,12])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,13])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,14])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,15])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,16])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,17])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,18])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,19])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,21])
	& abs(sig.all.fst[,5])> abs(sig.all.fst[,20]),]

clines.diff.3<-sig.all.fst[abs(sig.all.fst[,6])< abs(sig.all.fst[,7])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,8])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,9])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,10])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,11])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,12])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,12])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,13])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,14])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,15])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,16])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,17])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,18])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,19])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,21])
	& abs(sig.all.fst[,6])> abs(sig.all.fst[,20]),]#has 0

clines.diff.1.2<-clines.diff.2[clines.diff.2$Locus %in% clines.diff.1$Locus
	& clines.diff.2$Scaffold %in% clines.diff.1$Scaffold
	& clines.diff.2$BP %in% clines.diff.1$BP,]

scaff1189<-sig.all.fst[sig.all.fst$Scaffold=="scaffold_1189",]
scaff468<-sig.all.fst[sig.all.fst$Scaffold=="scaffold_468",]
scaff7<-sig.all.fst[sig.all.fst$Scaffold=="scaffold_7",]
scaff854<-sig.all.fst[sig.all.fst$Scaffold=="scaffold_854",]

write.table(clines.diff.1.2$Locus, row.names=F, col.names=F,
	"E://ubuntushare//fasta_from_stacks_catalog//fasta_from_stacks_catalog//sig_loci.txt", sep='\t')

#these are the sig ones from clines
scaff1189<-summ.dat[summ.dat$Chr=="scaffold_1189",]
scaff468<-summ.dat[summ.dat$Chr=="scaffold_468",]
scaff7<-summ.dat[summ.dat$Chr=="scaffold_7",]
scaff854<-summ.dat[summ.dat$Chr=="scaffold_854",]

scaff468$BP<-scaff468$BP+50000
scaff7$BP<-scaff7$BP+134000
scaff854$BP<-scaff854$BP+1370000
sig.scaff<-rbind(scaff1189, scaff468, scaff7, scaff854)
scaff.ss<-split(sig.scaff, sig.scaff$Pop.ID)

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

##########PLOT#######
library(car)

plot.eigenvectors<-function(dat.cov, col.num){
	#Plots the eigenvectors (scaled by the eigenvalues) 
	j = 1.96*sqrt(eigen(dat.cov)$values[col.num])
	delx1 = 0 + j*eigen(dat.cov)$vectors[1,col.num]
	delx2 = 0 + j*eigen(dat.cov)$vectors[2,col.num]
	delx11 = 0 - j*eigen(dat.cov)$vectors[1,col.num]
	delx22 = 0 - j*eigen(dat.cov)$vectors[2,col.num]
	x.values=c(0 , delx1, delx11)
	y.values=c(0 , delx2, delx22)
	lines(x.values, y.values, lwd=2, col="black")
}
plot.pmatrix<-function(pmatrix){
	#get leading eigenvector for bands
	#and leading eigenvector for length
	#Plots the P-matrix 95% confidence ellipse and eigenvectors.
	plot(c(-1,1), c(-1,1) , asp=1, type="n", xlim=c(-1,1), ylim=c(-1,1), 
		cex.lab=1.5, xlab="",ylab="", las=1)
	dat.cov<-cov(eigen(pmatrix)$vectors[,1:2])
	ellipse(center=c(0, 0), shape=dat.cov, 
		radius=1.96, center.cex=0.1, 
		lwd=2, col="black", add=TRUE)
	plot.eigenvectors(dat.cov, 1)
	plot.eigenvectors(dat.cov, 2)
}
jpeg("P-matrices.jpeg", res=300, height=14, width=14, units="in")
par(mfrow=c(3,4), mar=c(4,4,0.25,0.25))
plot.pmatrix(pmat.fem.std[["TXSP"]])
legend("top", "TXSP", bty="n")
plot.pmatrix(pmat.fem.std[["TXCC"]])
legend("top", "TXCC", bty="n")
plot.pmatrix(pmat.fem.std[["TXCB"]])
legend("top", "TXCB", bty="n")
plot.pmatrix(pmat.fem.std[["ALST"]])
legend("top", "ALST", bty="n")
plot.pmatrix(pmat.fem.std[["FLSI"]])
legend("top", "FLSI", bty="n")
plot.pmatrix(pmat.fem.std[["FLFD"]])
legend("top", "FLFD", bty="n")
mtext(expression("p"[2]),2, outer=T, line=-1.5)
plot.pmatrix(pmat.fem.std[["FLKB"]])
legend("top", "FLKB", bty="n")
plot.pmatrix(pmat.fem.std[["FLSG"]])
legend("top", "FLSG", bty="n")
plot.pmatrix(pmat.fem.std[["FLAB"]])
legend("top", "FLAB", bty="n")
plot.pmatrix(pmat.fem.std[["FLPB"]])
legend("top", "FLPB", bty="n")
plot.pmatrix(pmat.fem.std[["FLHB"]])
legend("top", "FLHB", bty="n")
plot.pmatrix(pmat.fem.std[["FLCC"]])
legend("top", "FLCC", bty="n")
mtext(expression("p"[1]),1, outer=T, line=-1)
dev.off()

jpeg("P-matrices.male.jpeg", res=300, height=14, width=14, units="in")
par(mfrow=c(3,4), mar=c(4,4,0.25,0.25))
plot.pmatrix(pmat.mal.std.pops[["TXSP"]])
legend("top", "TXSP", bty="n")
plot.pmatrix(pmat.mal.std.pops[["TXCC"]])
legend("top", "TXCC", bty="n")
plot.pmatrix(pmat.mal.std.pops[["TXCB"]])
legend("top", "TXCB", bty="n")
plot.pmatrix(pmat.mal.std.pops[["ALST"]])
legend("top", "ALST", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLSI"]])
legend("top", "FLSI", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLFD"]])
legend("top", "FLFD", bty="n")
mtext(expression("p"[2]),2, outer=T, line=-1.5)
plot.pmatrix(pmat.mal.std.pops[["FLKB"]])
legend("top", "FLKB", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLSG"]])
legend("top", "FLSG", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLAB"]])
legend("top", "FLAB", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLPB"]])
legend("top", "FLPB", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLHB"]])
legend("top", "FLHB", bty="n")
plot.pmatrix(pmat.mal.std.pops[["FLCC"]])
legend("top", "FLCC", bty="n")
mtext(expression("p"[1]),1, outer=T, line=-1)
dev.off()


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

########################RANDOM SKEWERS#######################################
#don't use this analysis.
library(evolqg)
rand.skewers<-RandomSkewers(pmat.fem.std.pops)
rs.out<-matrix(nrow=nrow(rand.skewers$probabilities), 
	ncol=ncol(rand.skewers$probabilities))
rs.out[upper.tri(rs.out)]<-rand.skewers$correlations[
	lower.tri(rand.skewers$correlations)]
rs.out[lower.tri(rs.out)]<-rand.skewers$probabilities[
	lower.tri(rand.skewers$probabilities)]
colnames(rs.out)<-colnames(rand.skewers$probabilities)
rownames(rs.out)<-rownames(rand.skewers$probabilities)
write.table(rs.out, "random.skewers.txt", sep="\t", eol='\n', quote=F,
	row.names=T, col.names=T)

mal.rand.skewers<-RandomSkewers(pmat.mal.std.pops)
mrs.out<-matrix(nrow=nrow(mal.rand.skewers$probabilities), 
	ncol=ncol(mal.rand.skewers$probabilities))
mrs.out[upper.tri(mrs.out)]<-mal.rand.skewers$correlations[
	lower.tri(mal.rand.skewers$correlations)]
mrs.out[lower.tri(mrs.out)]<-mal.rand.skewers$probabilities[
	lower.tri(mal.rand.skewers$probabilities)]
colnames(mrs.out)<-colnames(mal.rand.skewers$probabilities)
rownames(mrs.out)<-rownames(mal.rand.skewers$probabilities)
write.table(mrs.out, "mal.random.skewers.txt", sep="\t", eol='\n', quote=F,
	row.names=T, col.names=T)
#############################################################################
#IGNORE----------------->ATTEMPTS TO USE MIXED MODELS<-----------------IGNORE
#############################################################################
#Yij = u + Ii + Rj(i)
#re-create data matrix w/ columns PopID, ID, Trait, and Value
new.fem<-as.data.frame(rbind(as.matrix(unname(cbind(std.fem[,c(1,2,3)], 
	rep("SVL",nrow(std.fem))))),
as.matrix(unname(cbind(std.fem[,c(1,2,4)],rep("Length",nrow(std.fem))))), 
as.matrix(unname(cbind(std.fem[,c(1,2,5)],rep("Depth",nrow(std.fem))))),
as.matrix(unname(cbind(std.fem[,c(1,2,6)],rep("SnoutLength",nrow(std.fem))))),
as.matrix(unname(cbind(std.fem[,c(1,2,7)], rep("SnoutDepth",nrow(std.fem))))),
as.matrix(unname(cbind(std.fem[,c(1,2,8)],rep("HeadLength",nrow(std.fem))))),
as.matrix(unname(cbind(std.fem[,c(1,2,9)],rep("BandArea",nrow(std.fem))))),
as.matrix(unname(cbind(std.fem[,c(1,2,10)],rep("BandNum",nrow(std.fem)))))
))
colnames(new.fem)<-c("PopID", "IndID", "TraitValue", "Trait")
new.fem$PopID<-factor(new.fem$PopID)

pmatrix.fem<-lmer(TraitValue~Trait+PopID+(1|IndID), data=new.fem, REML=TRUE)
pmatrix.cov<-summary(pmatrix.fem)$vcov

new.pops.fem<-as.data.frame(rbind(as.matrix(unname(cbind(pops.std.fem[,c(1,2,3)], 
	rep("SVL",nrow(pops.std.fem))))),
as.matrix(unname(cbind(pops.std.fem[,c(1,2,4)],rep("Length",nrow(pops.std.fem))))), 
as.matrix(unname(cbind(pops.std.fem[,c(1,2,5)],rep("Depth",nrow(pops.std.fem))))),
as.matrix(unname(cbind(pops.std.fem[,c(1,2,6)],rep("SnoutLength",nrow(pops.std.fem))))),
as.matrix(unname(cbind(pops.std.fem[,c(1,2,7)], rep("SnoutDepth",nrow(pops.std.fem))))),
as.matrix(unname(cbind(pops.std.fem[,c(1,2,8)],rep("HeadLength",nrow(pops.std.fem))))),
as.matrix(unname(cbind(pops.std.fem[,c(1,2,9)],rep("BandArea",nrow(pops.std.fem))))),
as.matrix(unname(cbind(pops.std.fem[,c(1,2,10)],rep("BandNum",nrow(pops.std.fem)))))
))
colnames(new.pops.fem)<-c("PopID", "IndID", "TraitValue", "Trait")
ALST.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="ALST",], REML=TRUE)
FLAB.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLAB",], REML=TRUE)
FLCC.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLCC",], REML=TRUE)
FLFD.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLFD",], REML=TRUE)
FLHB.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLHB",], REML=TRUE)
FLKB.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLKB",], REML=TRUE)
FLPB.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLPB",], REML=TRUE)
FLSG.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLSG",], REML=TRUE)
FLSI.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="FLSI",], REML=TRUE)
TXCB.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="TXCB",], REML=TRUE)
TXCC.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="TXCC",], REML=TRUE)
TXSP.lmer<-lmer(TraitValue~Trait+(1|IndID), data=new.fem[new.fem$PopID=="TXSP",], REML=TRUE)




pops.std.fem<-NULL
for(i in 1:length(list.pops.std.fem)){
	pops.std.fem<-rbind(pops.std.fem, list.pops.std.fem[[i]])
}
colnames(pops.std.fem)<-colnames(fem.pheno)

pops.stdf.split<-split(pops.std.fem, fem.pheno$PopID)




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


