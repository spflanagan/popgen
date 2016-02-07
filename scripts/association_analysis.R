#Author: Sarah P. Flanagan
#Date: 23 May 2015
#Purpose: Look for associations with phenotypes

rm(list=ls())

##############################################################################
#****************************************************************************#
######################MANIPULATING FILES FOR ASSOCIATIONS#####################
#****************************************************************************#
##############################################################################
#fix the weird ones--ALAL, merge lines where necessary. It'll make life easier

raw.pheno<-read.table("B://ubuntushare//current_analysis//popgen.pheno.txt", sep="\t", header=T)
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
ped<-read.table("B://ubuntushare//current_analysis//batch_1.plink.ped", 
	skip = 1, stringsAsFactors=F, colClasses="character")
ped.names<-sub('sample_(\\w{4}\\w+).*[_.].*','\\1', ped[,2])
ped.names<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', ped.names)
pheno.dat$ID<-sub('([[:alpha:]]{5,7})([[:digit:]]{1})$', '\\10\\2', pheno.dat$ID)
pops.pheno<-pheno.dat[pheno.dat$ID %in% ped.names,] #this does not have the juveniles
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


write.table(ped.new, "B://ubuntushare//current_analysis//plink.pheno.ped",
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
##############################################################################
#****************************************************************************#
###########################ANALYZING PLINK OUTPUT#############################
#****************************************************************************#
##############################################################################

band.num.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.pheno.BandNum.qassoc",
	header=T)
dim(band.num.assoc[band.num.assoc$P<0.05,])

band.area.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.pheno.MeanBandArea.qassoc",
	header=T)
snout.depth.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.pheno.SnoutDepth.qassoc",
	header=T)
snout.head.pr.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.pheno.SnoutHeadProp.qassoc",
	header=T)
body.depth.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.pheno.BodyDepth.qassoc",
	header=T)
body.length.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.length.qassoc",
	header=T)
svl.assoc<-read.table(
	"B://ubuntushare//current_analysis//plink.pheno.qassoc",
	header=T)

num.05<-band.num.assoc[band.num.assoc$P<0.05,"SNP"]
area.05<-band.area.assoc[band.area.assoc$P < 0.05, "SNP"]

num.005<-band.num.assoc[band.num.assoc$P<0.005,"SNP"]
area.005<-band.area.assoc[band.area.assoc$P < 0.005, "SNP"]


library(lme4)
library(vegan)
#create male and female files for pst analysis
fem.pheno<-pops.pheno[pops.pheno$BandNum!=-9,]
fem.pheno<-replace(fem.pheno, fem.pheno==-9,NA)
write.table(fem.pheno, 
	"B://ubuntushare//current_analysis//fem.pheno.txt",
	sep="\t", quote=F, col.names=T, row.names=F)
mal.pheno<-pops.pheno[pops.pheno$BandNum==-9,]
mal.pheno<-replace(mal.pheno, mal.pheno==-9, NA)
mal.pheno<-mal.pheno[,-9]
mal.pheno<-mal.pheno[,-9]
mal.pheno<-mal.pheno[complete.cases(mal.pheno),]
write.table(mal.pheno, 
	"B://ubuntushare//current_analysis//mal.pheno.txt",
	sep="\t", quote=F, col.names=T, row.names=F)
trait.names<-c("SVL", "Standard Length", "Body Depth", 
	"Snout Length", "Snout Depth", "Head Length")
colors<-as.numeric(fem.pheno$PopID)
plot(trait.pca, type="n", display="sites", ann=F)
text(scores(trait.pca, display="sites"), fem.pheno$ID, cex=0.75, col=colors)





