#Author: Sarah P. Flanagan
#Last updated: 8 May 2016
#Date: 8 May 2016
#Purpose: Analyze Nerophis ophidion Population genetics data 

rm(list=ls())

library(ade4)
library(adegenet)
library(pcadapt)
library(scales)
library(gdata)

setwd("B:/ubuntushare/popgen/nerophis/")
source("../../SCA/scripts/plotting_functions.R")

pop.list<-c("SEW","LEM","GEL","STR","GTL","FIN")
pop.list<-c("LEM","GEL","STR","GTL","FIN")
pop.labels<-c("DKW","DKE","DEN","SEE","FIS")
#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################

pwise.fst.sub<-read.table("stacks/fst_summary_subset.txt",
	 header=T, row.names=1, sep='\t')
global.fst<-read.delim("stacks/pruned.globalstats.txt")
nerophis.coor<-read.csv("collecting_sites.csv")
sumstats<-read.delim("stacks/batch_3.sumstats.tsv",skip=5,header=T)
sumstats$Locus<-paste(sumstats$Locus.ID,sumstats$Col,sep=".")
geo.dist<-as.matrix(read.delim("nerophis_distances.txt",
	header=T,row.names=1,sep='\t'))

#########################################################################
#***********************************************************************#
########################POPULATION STRUCTURE#############################
#***********************************************************************#
#########################################################################

#******************************ADEGENET*********************************#
dat.plink<-read.PLINK("stacks/subset.raw",parallel=FALSE)
#look at alleles
png("Missingness_noSEW.png",height=7, width=7,units="in",res=300)
glPlot(dat.plink, posi="topleft")
dev.off()
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

pop<-substr(ind.names, 1,3)
colors<-pop
colors[colors=="FIN"]<-rainbow(5)[1]
colors[colors=="GEL"]<-rainbow(5)[2]
colors[colors=="GTL"]<-rainbow(5)[3]
colors[colors=="LEM"]<-rainbow(5)[4]
#colors[colors=="SEW"]<-rainbow(6)[5]
colors[colors=="STR"]<-rainbow(5)[5]
pop.list<-levels(as.factor(pop))

jpeg("subset.adegenet.pca1.2.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,1], pca1$scores[,2], pch=16, cex=2,lwd=1.3,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("topleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(5), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

jpeg("subset.adegenet.pca1.3.jpeg",res=300,height=7,width=7,units="in")
plot(pca1$scores[,2], pca1$scores[,3], pch=16, cex=2,
	col=alpha(colors, 0.5), ylab="", xlab="")
legend("topleft", pop.list, pch=19, pt.cex=2,
	col=alpha(rainbow(5), 0.5), ncol=3)
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2)
mtext(paste("PC3: ", round(pca1$eig[3]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 2)
dev.off()

#discriminant analysis of principal components (DAPC)
dat.clust<-find.clusters(dat.plink, parallel=FALSE, n.pca=20, n.clust=NULL,
	choose.n.clust=FALSE, max.n.clust=6)#
dapc1<-dapc(dat.plink, dat.clust$grp, n.pca=20,n.da=3, parallel=F)
png("E:/Docs/PopGen/adegenet.dapc.png",height=7,width=7,units="in",res=300)
scatter(dapc1, scree.da=FALSE, bg="white", posi.pca="topleft", legend=TRUE)
dev.off()
compoplot(dapc1)

#output k=5 clusters
adegenet.groups<-as.data.frame(cbind(names(dat.clust$grp), dat.clust$grp))
#get the discriminant analysis loadings
adegenet.da<-merge(adegenet.groups,dapc1$ind.coord,by=0)
adegenet.da$pop<-substr(adegenet.da$V1, 1,3)
adegenet.da$colors[adegenet.da$pop=="FIN"]<-rainbow(5)[1]
adegenet.da$colors[adegenet.da$pop=="GEL"]<-rainbow(5)[2]
adegenet.da$colors[adegenet.da$pop=="GTL"]<-rainbow(5)[3]
adegenet.da$colors[adegenet.da$pop=="LEM"]<-rainbow(5)[4]
#adegenet.da$colors[adegenet.da$pop=="SEW"]<-rainbow(5)[5]
adegenet.da$colors[adegenet.da$pop=="STR"]<-rainbow(5)[5]

adegenet.da$shape<-as.numeric(adegenet.da$V2)
adegenet.da$shape[adegenet.da$V2=="1"]<-21
adegenet.da$shape[adegenet.da$V2=="2"]<-22
adegenet.da$shape[adegenet.da$V2=="3"]<-23
adegenet.da$shape[adegenet.da$V2=="4"]<-24
adegenet.da$shape[adegenet.da$V2=="5"]<-25

plot(adegenet.da$LD1,adegenet.da$LD2,pch=as.numeric(adegenet.da$shape),
	bg=alpha(adegenet.da$colors,0.5),col=alpha("black",0.5),ylab="",xlab="",cex=2)
legend("topright",pch=c(21,22,23,24,25),pt.cex=2,
	c("Group 1","Group 2","Group 3", "Group 4", "Group 5"),
	col=alpha("black",0.5),ncol=3)
mtext("Discriminant Axis 1",1,line=2,cex=1.3)
mtext("Discriminant Axis 2",2,line=1.5,las=0,cex=1.3)
text(x=-8,y=4.5,"Adegenet,\nBest K = 5")

par(fig=c(0.5,1,0,0.45),new=T)
pca1$pops<-substr(rownames(pca1$scores),1,3)
pca1$colors[pca1$pops=="FIN"]<-rainbow(6)[1]
pca1$colors[pca1$pops=="GEL"]<-rainbow(6)[2]
pca1$colors[pca1$pops=="GTL"]<-rainbow(6)[3]
pca1$colors[pca1$pops=="LEM"]<-rainbow(6)[4]
pca1$colors[pca1$pops=="SEW"]<-rainbow(6)[5]
pca1$colors[pca1$pops=="STR"]<-rainbow(6)[6]
plot(pca1$scores[,1], pca1$scores[,2], pch=21, cex=2,lwd=1.3,
	bg=alpha(pca1$colors, 0.5),col=alpha("black",0.5), ylab="", xlab="")
mtext(paste("PC1: ", round(pca1$eig[1]/sum(pca1$eig)*100, 2), "%", sep=""), 
	1, line = 2,cex=1.3)
mtext(paste("PC2: ", round(pca1$eig[2]/sum(pca1$eig)*100, 2), "%", sep=""), 
	2, line = 1.5,las=0,cex=1.3)
text(x=1,y=2.8,"Adegenet,\nBest K = 3")

#*****************************STRUCTURE***********************************#
str.path<-"structure/nop_str/admix/Results/"
structure.k2<-read.table(paste(str.path,"admix_run_2_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k2$V1<-sub('([A-Z]{3}).*','\\1', structure.k2$V1)
structure.k3<-read.table(paste(str.path,"admix_run_3_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k3$V1<-sub('([A-Z]{3}).*','\\1', structure.k3$V1)
structure.k4<-read.table(paste(str.path,"admix_run_4_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k4$V1<-sub('([A-Z]{3}).*','\\1', structure.k4$V1)
structure.k5<-read.table(paste(str.path,"admix_run_5_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k5$V1<-sub('([A-Z]{3}).*','\\1', structure.k5$V1)
structure.k6<-read.table(paste(str.path,"admix_run_6_f_clusters.txt",sep=""),
	sep='\t', header=F)
structure.k6$V1<-sub('([A-Z]{3}).*','\\1', structure.k6$V1)


all.colors<-c("palegreen","goldenrod1","orchid3","tomato","darkblue", "forest green")

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

tapply(structure.k6$V2,structure.k6$V1,max)#V2=TX,V3=FLKB,V4=TXCB,V5=FLAB,V6=FLCC
str6<-data.frame(structure.k6$V1,structure.k6$V2, structure.k6$V4,structure.k6$V3,
	structure.k6$V5,structure.k6$V6)

png("structure_k2-6.png",height=10,width=7,units="in",res=300)
par(mfrow=c(5,length(pop.list)),mar=c(0.5,0,1,0),oma=c(1,3,1,0))
plotting.structure(str2,2,pop.list, make.file=FALSE, 
	colors=all.colors[c(1,5)],xlabel=F,ylabel="STRUCTURE\nK=2")
plotting.structure(str3,3,pop.list, make.file=FALSE, 
	colors=all.colors[c(1,3,5)],xlabel=F,ylabel="STRUCTURE\nK=3")
plotting.structure(str4,4,pop.list, make.file=FALSE,
	colors=all.colors[c(1,2,3,5)],xlabel=F,ylabel="STRUCTURE\nK=4")
plotting.structure(str5,5,pop.list, make.file=FALSE,
	colors=all.colors[c(1,2,3,5,6)],xlabel=F,ylabel="STRUCTURE\nK=5")
plotting.structure(str6,6,pop.list, make.file=FALSE, colors=all.colors,
	xlabel=F,ylabel="STRUCTURE\nK=6")
dev.off()
#It looks like k=2 is best...they are kind of panmictic.



#########################################################################
#***********************************************************************#
#####################POPULATION DIFFERENTIATION##########################
#***********************************************************************#
#########################################################################

#**************************GLOBAL FSTS********************************#
png("global_fsts.png",height=4,width=10,units="in",res=300)
par(mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(global.fst$Fst,pch=19,ylab=expression(italic(F)[italic(ST)]),
	xlab="SNP Index")
dev.off()

#**************************STACKS FSTS********************************#
setwd("stacks")
stacks.files<-data.frame(file=list.files(pattern="batch_3.fst_\\w{3}-\\w{3}"))
stacks.files$pop1<-gsub("batch_3.fst_(\\w{3})-\\w{3}.tsv","\\1",stacks.files[,1])
stacks.files$pop2<-gsub("batch_3.fst_\\w{3}-(\\w{3}).tsv","\\1",stacks.files[,1])

png("stacks_AMOVAfsts.png",height=10,width=10,units="in",res=300)
par(mfrow=c(length(pop.list),length(pop.list)),mar=c(0.5,1,0.5,0.5),
	oma=c(2,2,2,2))
for(i in 1:length(pop.list)){
	for(j in 1:length(pop.list)){
		
		file<-as.character(stacks.files[stacks.files$pop1==pop.list[i] &
			stacks.files$pop2==pop.list[j] | 
			stacks.files$pop1==pop.list[j] &
			stacks.files$pop2==pop.list[i] ,"file"])
		if(length(file)>0 & i <=j){
			print(pop.list[i])
			dat<-read.delim(file)
			plot(dat$AMOVA.Fst,pch=19,axes=F,xlab="",ylab="")
			axis(2,las=1,pos=0)
			abline(h=mean(dat$AMOVA.Fst),lty=2,col="dodgerblue")
		} else{
			plot(x=c(0,1),y=c(0,1),axes=F,type='n')
		}
	}
}
mtext(expression(italic(F)[italic(ST)]),2,outer=T)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#column labels
text(x=-0.85,y=1,pop.list[1])
text(x=-0.5,y=1,pop.list[2])
text(x=-0.15,y=1,pop.list[3])
text(x=0.15,y=1,pop.list[4])
text(x=0.5,y=1,pop.list[5])
text(x=0.85,y=1,pop.list[6])
#row labels
text(x=-1,y=0.85,pop.list[1])
text(x=-1,y=0.5,pop.list[2])
text(x=-1,y=0.15,pop.list[3])
text(x=-1,y=-0.15,pop.list[4])
text(x=-1,y=-0.5,pop.list[5])
text(x=-1,y=-0.85,pop.list[6])
dev.off()


#########################################################################
##########################BAYENV-FILTER WOD##############################
#########################################################################
setwd("./wod_files")
wod.files<-list.files(pattern="_out.txt")

wod.dat<-data.frame()

for(i in 1:length(wod.files)){
  wod.dat<-rbind(wod.dat,read.delim(wod.files[i], sep="\t"))
}
#keep only those from 2004-2014 (2014 is newest)
wod.filt<-wod.dat[wod.dat$Year >= 2004,]


mar.coor<-read.csv("../collecting_sites.csv")
adj.coor<-as.data.frame(cbind(site=as.character(mar.coor$Location), 
      lat.l=as.numeric(mar.coor$Latitude-0.5), 
      lat.r=as.numeric(mar.coor$Latitude),
	lat.h=as.numeric(mar.coor$Latitude+0.5), 
      lon.l=as.numeric(mar.coor$Longitute-0.5),
      lon.r=as.numeric(mar.coor$Longitute),
	lon.h=as.numeric(mar.coor$Longitute+0.5)), stringsAsFactors = FALSE)
adj.coor[,2:7]<-as.data.frame(sapply(adj.coor[,2:7],as.numeric))

#restrict to the actual coordinates
wod.rest<-list()
for(i in 1:nrow(adj.coor)){
  wod.temp<-as.data.frame(subset(wod.filt,
                                 round(Long,1) <= adj.coor$lon.h[i] & round(Long,1) >= adj.coor$lon.l[i]))
  wod.rest[[i]]<-as.data.frame(subset(wod.temp,  
                                      Lat <= adj.coor$lat.h[i] & Lat >= adj.coor$lat.l[i]))
}
names(wod.rest)<-paste(mar.coor$Location)
averages<-lapply(wod.rest, function(x){
  avgs<-apply(x,2,mean,na.rm=T)
  n<-nrow(x)
  return(list("avgs"=avgs, "n"=n))
})

environ.file<-NULL
for(i in 1:length(averages)){
  temp.avg<-append(averages[[i]]$avgs,averages[[i]]$n)
  environ.file<-cbind(environ.file,temp.avg)
}
colnames(environ.file)<-mar.coor$Location
rownames(environ.file)<-c(names(averages[[1]]$avgs), "n")
environ.file<-environ.file[-9,]#remove pH
environ.file<-environ.file[-9,]#remove Oxygen


##include temperature variance
temp.var<-lapply(wod.rest,function(x){
  tempvar<-var(x$Temp)
  return(tempvar)
})
temp.var<-t(do.call("rbind",temp.var))
rownames(temp.var)<-c("tempvar")
std.tempvar<-(temp.var-mean(temp.var))/sd(temp.var)
environ.file<-rbind(environ.file,temp.var)
#write environmental data to file

write.table(environ.file, "bayenv/wod_data_nop_bayenv.txt",
	sep='\t',eol="\t\n", quote=F)

#**********************************BAYENV2***********************************#
#####STARTING WITH PLINK FILES
ped<-read.table("stacks/subset.ped", 
	stringsAsFactors=F, colClasses="character")
ped.pops<-substr(ped[,2],1,3)
ped.sex<-sub('(\\w{3})(\\w)(\\d+)','\\2', ped[,2])
ped.sex[ped.sex=="F"]<-2
ped.sex[ped.sex=="M"]<-1
ped.sex[ped.sex=="I"]<-0
ped[,1]<-ped.pops
ped[,5]<-ped.sex
write.table(ped,"bayenv/bayenv.plink.ped", 
	row.names=F, col.names=F, quote=F, sep=" ",eol="\n")

clust.plink<-cbind(ped.pops, ped[,2],ped.pops)
write.table(clust.plink, 
	"stacks/plink.clust.txt",
	col.names=F, row.names=F, quote=F, sep="\t", eol="\n")
#Then plink --ped bayenv.plink.ped --map subset.map \
#--extract plink.snplist --out bayenv --noweb --allow-no-sex --recode \
#--freq --within plink.clust.txt 

#####CONVERT PLINK TO BAYENV2
freq<-read.table("bayenv/bayenv.frq.strat", 
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

write.table(snpsfile, "bayenv/nop.snpsfile", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n") #bayenv SNPSFILE

#NOW RUN MATRIX ESTIMATION: run_bayenv2_matrix_general.sh

#####check Bayenv2 matrix
matrix.files<-list.files("bayenv/",pattern="matrix")
matrices<-list()
for(i in 1:length(matrix.files))
{
	matrices[[i]]<-as.matrix(read.table(
		paste("bayenv/",matrix.files[i],sep=""), 
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
all.snps.ped<-read.table("stacks/batch_3.plink.ped", header=F, stringsAsFactors=F)
ped.pop<-sub('(\\w{3})\\w+\\d+','\\1', all.snps.ped[,2])
all.snps.clust<-cbind(all.snps.ped[,1],all.snps.ped[,2],all.snps.ped[,1])
write.table(all.snps.clust, "bayenv/all.clust.txt", sep="\t", eol="\n", quote=F,
	row.names=F, col.names=F)
#then need to run plink --file stacks/batch_3.plink --freq --within bayenv/all.clust.txt \
	--allow-no-sex --noweb --out bayenv/all.bayenv.plink

#read in frequency per pop
all.snps.frq<-read.table("bayenv/all.bayenv.plink.frq.strat", 
	header=T, stringsAsFactors=F)
freq<-cbind(all.snps.frq,all.snps.frq$NCHROBS-all.snps.frq$MAC)
colnames(freq)[ncol(freq)]<-"NAC"
pop.order<-levels(as.factor(freq$CLST))
snp.names<-split(freq$SNP,freq$CLST)[[1]]

mac.by.pop<-as.data.frame(split(freq$MAC,freq$CLST))
rownames(mac.by.pop)<-snp.names
nac.by.pop<-as.data.frame(split(freq$NAC,freq$CLST))
rownames(nac.by.pop)<-snp.names
snpsfile<-interleave(mac.by.pop,nac.by.pop)

write.table(snpsfile, "bayenv/nop.all.snpsfile", 
	col.names=F,row.names=F,quote=F,sep="\t",eol="\n")
for(i in seq(1,(nrow(snpsfile)/2),2)){
	write.table(snpsfile[i:(i+1),],
		paste("bayenv/snpfiles/",rownames(snpsfile)[i],sep=""),
		col.names=F,row.names=F,quote=F,sep='\t',eol='\n')
}
#####ENVFILE
env.raw<-read.table("bayenv/wod_data_nop_bayenv.txt")
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
	"bayenv/env_data_bayenv_std.txt",
	sep='\t',quote=F,col.names=F,row.names=F,eol='\n')

##Are they correlated with distance?
colnames(env.std)<-colnames(env.raw)
rownames(env.std)<-rownames(env.raw)
env.dist<-as.matrix(vegdist(t(env.raw)))
env.dist<-env.dist[rownames(geo.dist),colnames(geo.dist)]
mantel.rtest(as.dist(t(geo.dist)),as.dist(env.dist),999)
#Monte-Carlo test
#Observation: 0.01898832 
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
#Based on 999 replicates
#Simulated p-value: 0.438 





