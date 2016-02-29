#Author: Sarah P. Flanagan
#Date: 27 July 2015
#Purpose: Un-used analyses of population genetics data
#made from un-used analyses from fst_cline_analysis and population_structure

rm(list=ls())

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

#########################################################################
#***********************************************************************#
###########COMPARE GROUP ASSIGNMENTS FROM ALL OF THE PROGRAMS############
#***********************************************************************#
#########################################################################
fast.groups.3 #faststructure
adegenet.groups #adegenet

########
#generate ld.hwe whitelist for Stacks populations
ld.hwe.map<-read.table("E://ubuntushare//stacks//populations//ld.hwe.sub.map",
	sep='\t')
whitelist<-sub('(\\d+)(_)(\\d+)','\\1',ld.hwe.map[,2])
write.table(whitelist[!duplicated(whitelist)], 
	"E://ubuntushare//stacks//populations//ld.hwe.whitelist.txt",
	sep='\t', eol='\n', quote=F, col.names=F, row.names=F)
write.table(ld.hwe.map[,2], 
	"E://ubuntushare//stacks//populations//ld.hwe.snplist.txt",
	sep='\t', eol='\n', quote=F, col.names=F, row.names=F)



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





