summ.dat<-read.table("results/stacks/populations/batch_1.sumstats.tsv",
	sep='\t', skip=12, header=T, comment.char="")
ld.hwe<-read.table("results/stacks/populations/ld.subset.hwe.snplist")

lh.list<-do.call(rbind,lapply(as.list(as.character(ld.hwe$V1)), function(x) {
	do.call(rbind,strsplit(x,"_")) }))
colnames(lh.list)<-c("Locus","Col")
summ.prune<-summ.dat[summ.dat$Locus.ID %in% ld.hwe$V1,]
#I'm just checking to see loci and scaffolds etc.

summ.prune$SNPID<-paste(summ.prune$Chr,".",summ.prune$BP,sep="")
summ.prune<-summ.prune[!duplicated(summ.prune$SNPID),]

summ.dat$SNPLoc<-paste(summ.dat$Locus.ID,"_",summ.dat$Col,sep="")
summ.prune<-summ.dat[summ.dat$SNPLoc %in% ld.hwe$V1,]
summ.prune<-summ.prune[!duplicated(summ.prune$SNPLoc),]