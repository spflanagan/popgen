#Author: Sarah P. Flanagan
#Date: May 2015
#Purpose: Perform population structure analyses using several methods
#including those that correct for isolation by distance



#########################################################################
#***********************************************************************#
################################ADEGENET#################################
#***********************************************************************#
#########################################################################
library(adegenet)
library(scales)

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
################################STRUCTURE################################
#***********************************************************************#
#########################################################################
#may not need to use R for this.

#########################################################################
#***********************************************************************#
################################BAYENV2##################################
#***********************************************************************#
#########################################################################

#**************************STARTING WITH PLINK FILES********************#
ped<-read.table("E://ubuntushare//stacks//populations//batch_1.plink.ped", 
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

write.table(ped,"E://ubuntushare//stacks//populations//plink.ped", 
	row.names=F, col.names=F, quote=F, sep="\t",eol="\n")

clust.plink<-cbind(ped.pops, ped[,2],ped.pops)
write.table(clust.plink, "E://ubuntushare//stacks//populations//plink.clust.txt",
	col.names=F, row.names=F, quote=F, sep="\t", eol="\n")


#plink.map -> numbers instead of scaffold_#
map<-read.table("E://ubuntushare//stacks//populations//batch_1.plink.map", 
	skip = 1)
chr.nums<-sub('scaffold_','',map[,1])
map[,1]<-chr.nums

write.table(map, "E://ubuntushare//stacks//populations//plinkn.map", 
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
freq<-read.table("E://ubuntushare//stacks//populations//ld.subset.frq.strat", 
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
library(gdata)
snpsfile<-interleave(mac.by.pop,nac.by.pop)

write.table(mac.by.pop, "E://ubuntushare//bayenv2//SNPFILE", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n")

write.table(snpsfile, "E://ubuntushare//bayenv2//SNPSFILE", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n")

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
image(matrices[[10]])#these all look pretty similar

fsts<-read.table("E://ubuntushare//stacks//populations//batch_1.fst_summary.tsv",
	header=T,row.names=1,sep='\t')
fsts<-rbind(fsts, rep(NA,12))
rownames(fsts)[12]<-"TXCC"


fst.t<-t(fsts)
fsts[lower.tri(fsts)]<-fst.t[lower.tri(fst.t)]
fsts[is.na(fsts)]<-1 #fill in the diagonal
#order pops to be same as matrices
fsts<-fsts[order(rownames(fsts),pop.order),order(colnames(fsts),pop.order)]

cor(fsts,matrices[[1]])
#these correlations are pretty low

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