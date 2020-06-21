---
title: Supplemental Material for "Multiple colonizations of freshwater by the Gulf pipefish reveal a shared genomic signature of adaptation"
preprint: false
author: 
  - name: Sarah P. Flanagan
    affilnum: 1
    corresponding: true
    email: spflanagan.phd@gmail.com
  - name: Emily Rose
    affilnum: 2
  - name: Adam Jones
    affilnum: 3
affiliation:
  - affilnum: 1
    affil: School of Biological Sciences, University of Canterbury, 4800 Private Bag, Christchurch 8140 New Zealand
  - affilnum: 2
    affil: Department of Biology, The University of Tampa, Tampa, FL 33606 USA
  - affilnum: 3
    affil: Department of Biological Sciences, University of Idaho, Moscow, ID 83844 USA
abstract: >
  This document includes supplementary material for the paper. 
  In this document, we walk through the generation of each of the figures in the manuscript, including all of the additional analyses that were used to decide on the final presentation of results.
  This document shows the population structure analyses (Fsts, PCAdapt, Admixture, and Treemix), which were used to generate Figure 1.
  The second major section of this document deals with the identification of outliers and the interpretation of outliers from multiple analyses.
  These analyses combined to create Figure 2 and Figure 3. 
header-includes: >
  \usepackage{lipsum}
  \usepackage{float}
  \floatplacement{figure}{H}
bibliography: programs.bib
output:
  bookdown::pdf_document2:
    toc: true
    toc_depth: 2
    fig_caption: yes
    keep_tex: yes
    number_sections: no
    template: manuscript.latex
  html_document: null
  word_document: null
fontsize: 11pt
capsize: normalsize
csl: molecular-ecology.csl
documentclass: article
spacing: singlespacing
---

# Overview of the study {-}

The initial analyses are in `200_fwsw_analysis.Rmd` and conducted the analyses on a dataset generated from comparing lumped 'freshwater' and 'saltwater' populations, containing SNPs found in 50% of individuals and with a minor allele frequency of at least 5%. The revised paper will instead focus on two datasets:

1. A dataset with all 16 populations, generated from all pairwise comparisons of populations, containing SNPs found in every population, in 75% of individuals, and with a minor allele frequency of at least 5%.

2. A dataset containing only the 4 freshwater populations (TXFW, LAFW, ALFW, FLFW) and their nearest saltwater populations (TXCC, ALST, FLCC -- note ALST is the nearest neighbor to both ALFW and LAFW). This dataset also contains SNPs found in 75% of individuals with a minor allele frequency of at least 5%. 




```r
source("../../gwscaR/R/gwscaR.R")
source("../../gwscaR/R/gwscaR_plot.R")
source("../../gwscaR/R/gwscaR_utility.R")
source("../../gwscaR/R/gwscaR_fsts.R")
source("../../gwscaR/R/gwscaR_popgen.R")
source("../../gwscaR/R/vcf2dadi.R")
source("../R/203_treemix_plotting_funcs.R")#I've modified these functions
library(knitr)
library(scales)
library(kableExtra)
library(vegan)
library(RColorBrewer)
```


```r
pop.list<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
	"FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLLG")
pop.labs<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLFW")
fw.list<-c("TXFW","LAFW","ALFW","FLLG")
sw.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB",
	"FLFD","FLSI","FLAB","FLPB","FLHB","FLCC")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
lgn<-seq(1,22)
all.colors<-c(rep("black",2),"#2166ac","black","#2166ac","black","#2166ac",
        rep("black",8),"#2166ac")
#grp.colors<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33','#f781bf')
grp.colors<-c('#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837')
col_vector<-c(red='#e6194b', green='#3cb44b', blue='#4363d8',yellow='#ffe119', 
              cyan='#46f0f0',orange='#f58231', teal='#008080', purple='#911eb4',  
              magenta='#f032e6', lime='#bcf60c', pink='#fabebe',  lavendar='#e6beff', 
              brown='#9a6324', olive='#808000', apricot='#ffd8b1',maroon='#800000', 
              mint='#aaffc3', navy='#000075', beige='#fffac8', grey='#808080', 
              white='#ffffff', black='#000000')

col_vector<-c('#762a83','#762a83',"#2166ac",'#762a83',"#2166ac",'#af8dc3',
              "#2166ac",'#e7d4e8','#e7d4e8','#e7d4e8','#e7d4e8','#7fbf7b',
              '#1b7837','#1b7837','#1b7837',"#2166ac")

ppi<-data.frame(Pop=pop.labs,cols = col_vector,
                pch=rep(c(15,16,17,18),4))
ppi$pch[grep("FW",ppi$Pop)]<-c(15,16,17,18)
```



## Map of the populations

First, we'll plot these populations on a map


```r
library(maps);library(gplots);library(mapdata)
mar.coor<-read.csv("marine_coordinates_revised.csv", header=T)
fw.coor<-read.csv("fw_coordinates.csv", header=T)
```

```r
jpeg("all_sites_map.jpeg", res=300, height=7,width=14, units="in")
#pdf("all_sites_map.pdf",height=7,width=14)
par(oma=c(0,0,0,0),mar=c(0,0,0,0),pin=c(7,7))
map("worldHires", "usa",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray90", mar=c(0,0,0,0),fill=TRUE, res=300,myborder=0)
map("worldHires", "mexico",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray95", fill=TRUE, add=TRUE)
points(mar.coor$lon, mar.coor$lat,  col="black", cex=2, pch=19)
points(-1*fw.coor$lon, fw.coor$lat,  col="cornflowerblue", cex=2, pch=18)
abline(h=c(25,30,35),lty=3)
abline(v=c(-80,-85,-90,-95),lty=3)
text(x=c(-99.5,-99.5),y=c(25,30),c("25N","30N"),cex=1.75)
text(x=c(-80,-85,-90,-95),y=rep(31.8,4),c("80W","85W","90W","95W"),cex=1.75)
text(y=26,x=-90,"Gulf of Mexico",cex=1.75)
text(y=25.5,x=-98.5,"Mexico",cex=1.75)
text(x=-91,y=31,"USA",cex=1.75)
text(x=-78,y=29.5,"Atlantic Ocean",cex=1.75)
text(x=-96.4,y=26,"TXSP",font=2,cex=1.75)
text(x=-96.6,y=27.2,"TXCC",font=2,cex=1.75)
text(x=-95.6,y=28.3,"TXFW",font=2,col="cornflowerblue",cex=1.75)
text(x=-94.4,y=29,"TXCB",font=2,cex=1.75)
text(x=-90.5,y=29.8,"LAFW",font=2,col="cornflowerblue",cex=1.75)
text(x=-88,y=30,"ALST",font=2,cex=1.75)
text(x=-87,y=30.75,"ALFW",font=2,col="cornflowerblue",cex=1.75)
text(x=-85,y=29.4,"FLSG",font=2,cex=1.75)
text(x=-83.7,y=29,"FLKB",font=2,cex=1.75)
text(x=-83.4,y=27.6,"FLFD",font=2,cex=1.75)
text(x=-82.4,y=26,"FLSI",font=2,cex=1.75)
text(x=-79.9,y=24.8,"FLAB",font=2,cex=1.75)
text(x=-79.2,y=26.8,"FLPB",font=2,cex=1.75)
text(x=-79.4,y=27.2,"FLHB",font=2,cex=1.75)
text(x=-79.9,y=28.5,"FLCC",font=2,cex=1.75)
text(x=-80.9,y=29.5,"FLFW",font=2,col="cornflowerblue",cex=1.75)
dev.off()
```

## Phenotypic variation



```r
raw.pheno<-read.table("../sw_results/popgen.pheno.txt", sep="\t", header=T)
	raw.pheno$PopID<-gsub("(\\w{4})\\w+","\\1",raw.pheno$ID)
	raw.pheno<-raw.pheno[raw.pheno$PopID %in% pop.list,]
	raw.pheno$sex<-gsub("\\w{4}(\\w)\\w+","\\1",raw.pheno$ID)
	raw.pheno$TailLength<-raw.pheno$std.length-raw.pheno$SVL
	raw.pheno$HeadLength<-raw.pheno$HeadLength-raw.pheno$SnoutLength

fem.pheno<-raw.pheno[raw.pheno$sex %in% c("F","D"),-8]
	fem.pheno<-fem.pheno[,c(11,1,10,2,12,4,5,6,7,8,9)]
	fem.pheno<-fem.pheno[order(match(fem.pheno$PopID,pop.list)),]
	write.table(fem.pheno,"fem.pheno.txt",sep='\t',row.names=F,col.names=T,
		quote=F)
	
mal.pheno<-raw.pheno[raw.pheno$sex %in% c("P","N"),-8]
	mal.pheno<-mal.pheno[,c(11,1,10,2,12,4,5,6,7)]
	mal.pheno<-mal.pheno[order(match(mal.pheno$PopID,pop.list)),]
	write.table(mal.pheno,"mal.pheno.txt",sep='\t',row.names=F,col.names=T,
		quote=F)
```

```r
fem.pheno<-read.table("fem.pheno.txt",header=T)
	fem.pheno<-fem.pheno[!is.na(fem.pheno$BandNum),]
mal.pheno<-read.table("mal.pheno.txt",header=T)
```


```r
fem.pheno$PopID<-factor(fem.pheno$PopID)
fem.pheno<-fem.pheno[!is.na(fem.pheno$BandNum),]
mal.pheno$PopID<-factor(mal.pheno$PopID)
bands.pcdat<-fem.pheno[!is.na(fem.pheno$BandNum),
	c("PopID","ID","MBandArea","BandNum")]
# run pcas
band.pca<-rda(bands.pcdat[,3:4])
fem.pheno.pca<-rda(fem.pheno[,4:9])
mal.pheno.pca<-rda(mal.pheno[,4:9])
```

```r
####extract eigenvalue
band.eig<-band.pca$CA$eig
band.pc<-band.eig/sum(band.eig)*100

#extract PC scores
band.u<-data.frame(bands.pcdat[,1:2],
                   "BandPC1"=band.pca$CA$u[,1],stringsAsFactors=F)
band.u.sep<-split(band.u, band.u[,1])
band.u.new<-rbind(band.u.sep$TXSP,band.u.sep$TXCC,band.u.sep$TXCB,
	band.u.sep$ALST,band.u.sep$FLSG,band.u.sep$FLKB,
	band.u.sep$FLFD,band.u.sep$FLSI,band.u.sep$FLAB,
	band.u.sep$FLPB,band.u.sep$FLHB,band.u.sep$FLCC)

fem.pheno.eig<-fem.pheno.pca$CA$eig
fem.pheno.pc<-fem.pheno.eig/sum(fem.pheno.eig)*100

#extract PC scores
fem.pheno.u<-data.frame(fem.pheno[,1:2],
	"FemBodyPC1"=fem.pheno.pca$CA$u[,1],stringsAsFactors=F)
fem.u.sep<-split(fem.pheno.u, fem.pheno.u[,1])
fem.u.new<-rbind(fem.u.sep$TXSP,fem.u.sep$TXCC,fem.u.sep$TXCB,
	fem.u.sep$ALST,fem.u.sep$FLSG,fem.u.sep$FLKB,
	fem.u.sep$FLFD,fem.u.sep$FLSI,fem.u.sep$FLAB,
	fem.u.sep$FLPB,fem.u.sep$FLHB,fem.u.sep$FLCC)

mal.pheno.eig<-mal.pheno.pca$CA$eig
mal.pheno.pc<-mal.pheno.eig/sum(mal.pheno.eig)*100

#extract PC scores
mal.u<-data.frame(mal.pheno[,1:2],"MalBodyPC1"=mal.pheno.pca$CA$u[,1],
	stringsAsFactors=F)
mal.u.sep<-split(mal.u, mal.u[,1])
mal.u.new<-rbind(mal.u.sep$TXSP,mal.u.sep$TXCC,mal.u.sep$TXCB,
	mal.u.sep$ALST,mal.u.sep$FLSG,mal.u.sep$FLKB,
	mal.u.sep$FLFD,mal.u.sep$FLSI,mal.u.sep$FLAB,
	mal.u.sep$FLPB,mal.u.sep$FLHB,mal.u.sep$FLCC)
```


```r
# females
fem.pop<-as.character(bands.pcdat$PopID)
fem.pop[fem.pop=="FLLG"]<-"FLFW"
fem.colors<-as.character(fem.pop)
fem.pch<-as.character(fem.pop)
fw.fem.col<-as.character(fem.pop[fem.pop %in% fw.list])
for(i in 1:length(fem.pop)){
  fem.colors[i]<-as.character(ppi[ppi$Pop %in% fem.pop[i],"cols"])
  fem.pch[i]<-as.numeric(as.character(ppi[ppi$Pop %in% fem.pop[i],"pch"]))
}
fem.pch<-as.numeric(fem.pch)

# males
mal.pop<-as.character(mal.pheno$PopID)
mal.pop[mal.pop=="FLLG"]<-"FLFW"
mal.colors<-as.character(mal.pop)
mal.pch<-as.character(mal.pop)
fw.mal.col<-as.character(mal.pop[mal.pop %in% fw.list])
for(i in 1:length(mal.pop)){
  mal.colors[i]<-as.character(ppi[ppi$Pop %in% mal.pop[i],"cols"])
  mal.pch[i]<-as.character(ppi[ppi$Pop %in% mal.pop[i],"pch"])
}
mal.pch<-as.numeric(mal.pch)


fw.fem.rows<-which(fem.pheno$PopID %in% fw.list)
fw.mal.rows<-which(mal.pheno$PopID %in% fw.list)
```



```r
ptCex<-2

par(mfrow=c(2,3),oma=c(2,2,2,2),mar=c(2,2,2,2),lwd=1.3)
mp<-plot(mal.pheno.pca,type="n",xlim=c(-3,3),ylim=c(-8.2,4)
	,xlab="",ylab="",las=1,cex.axis=1.5)
points(mal.pheno.pca,col=alpha(mal.colors,0.5),cex=ptCex,pch=mal.pch)
mtext(paste0("PC1 (",round(mal.pheno.pc[1],2),"%)"),1,line=2)
mtext(paste0("PC2 (",round(mal.pheno.pc[2],2),"%)"),2,line=2.5)
legend("top",bty='n',c("Male Body Traits"),cex=1.5)

fp<-plot(fem.pheno.pca,type="n",xlab="",ylab="",las=1,cex.axis=1.5,ylim=c(-4,12),
	xlim=c(-3,3))
points(fem.pheno.pca,col=alpha(fem.colors,0.5),cex=ptCex,pch=fem.pch)
mtext(paste0("PC1 (",round(fem.pheno.pc[1],2),"%)"),1,line=2)
mtext(paste0("PC2 (",round(fem.pheno.pc[2],2),"%)"),2,line=2.5)
legend("top",bty='n',c("Female Body Traits"),cex=1.5)

bp<-plot(band.pca,type="n",xlab="",ylab="",las=1,
         cex.axis=1.5,xlim=c(-2,2),ylim=c(-3,1))
points(band.pca,pch=fem.pch,col=alpha(fem.colors,0.5),cex=ptCex)
mtext(paste0("PC1 (",round(band.pc[1],2),"%)"),1,line=2)
mtext(paste0("PC2 (",round(band.pc[2],2),"%)"),2,line=2.5)
legend("top",bty='n',c("Female Band Traits"),cex=1.5)


plot(mp$sites[fw.mal.rows,],type="n",xlab="",ylab="",las=1,cex.axis=1.5)
abline(h=0,lty=3)
abline(v=0,lty=3)
points(mp$sites[fw.mal.rows,],xlim=c(-0.1,0.1),ylim=c(-.2,.2),
	col=alpha(mal.colors[fw.mal.rows],0.5),cex=ptCex,pch=mal.pch[fw.mal.rows])
mtext(paste0("PC1 (",round(mal.pheno.pc[1],2),"%)"),1,line=2)
mtext(paste0("PC2 (",round(mal.pheno.pc[2],2),"%)"),2,line=2.5)


plot(fp$sites[fw.fem.rows,],type="n",xlab="",ylab="",las=1,
	cex.axis=1.5,ylim=c(-4,12),xlim=c(-3,3))
abline(h=0,lty=3)
abline(v=0,lty=3)
points(fp$sites[fw.fem.rows,],
	col=alpha(fem.colors[fw.fem.rows],0.5),cex=ptCex,
	pch=fem.pch[fw.fem.rows])
mtext(paste0("PC1 (",round(fem.pheno.pc[1],2),"%)"),1,line=2)
mtext(paste0("PC2 (",round(fem.pheno.pc[2],2),"%)"),2,line=2.5)


plot(bp$sites[fw.fem.rows,],type="n",xlab="",ylab="",las=1,
	cex.axis=1.5,xlim=c(-2,2),ylim=c(-3,1))
points(bp$sites[fw.fem.rows,],
	pch=fem.pch[fw.fem.rows],col=alpha(fem.colors[fw.fem.rows],0.5),
	cex=ptCex)
mtext(paste0("PC1 (",round(band.pc[1],2),"%)"),1,line=2)
mtext(paste0("PC2 (",round(band.pc[2],2),"%)"),2,line=2.5)
abline(h=0,lty=3)
abline(v=0,lty=3)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), 
    mar = c(0, 0, 0, 0), new = TRUE, cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", legend=ppi$Pop, 
	col=as.character(ppi$col),
	pt.cex=ptCex,bty='n',pch=ppi$pch, ncol=8)
```

\begin{figure}[H]
\includegraphics{../figs/plotPCA-1} \caption{Principal components analysis of morphological traits in S. scovelli reveals that phenotypic variation among populations is not based on habitat type. The top set of panels show the results of the PCA with all 16 populations, color-coded by populations and point shape. The bottom set of panels show the same PCA results, but with different x- and y-axis scaling and without the saltwater populations plotted, to facilitate visualizing the differences among saltwater populations. The left panels show male body traits (SVL, tail length, trunk depth, head length, snout length, and snout depth), the middle panels show those same traits in females, and the right panels show the female band traits (band number and band area).}(\#fig:plotPCA)
\end{figure}

## Table of summary statistics



```r
pop_map<-read.delim("../fwsw_pops_map.txt",header = FALSE,stringsAsFactors = FALSE)
ful_vcf<-parse.vcf(
  "filter_rad_20191014@1654/14_filtered/radiator_data_20191014@1710.vcf")
colnames(ful_vcf)<-gsub("\\-","_",colnames(ful_vcf))
```

```r
sub_vcf<-parse.vcf("stacks/populations_subset75/batch_2.pruned.vcf")
sub_all_vcf<-parse.vcf("stacks/populations_subset75/all_pops_subset75/batch_2.vcf")
# combine the two
sub<-merge(sub_vcf,sub_all_vcf,by="ID",all = TRUE)
rmv<-grep(".y",colnames(sub))
sub<-sub[,-rmv]
colnames(sub)<-gsub(".x","",colnames(sub))
```


```r
pop_summaries<-do.call(rbind,lapply(pop.list,function(pop,ful,sub){
  ful_dat<-ful[,c(1:9,grep(pop,colnames(ful)))]
  sub_dat<-sub[,c(1:9,grep(pop,colnames(sub)))]
  # calc observed het values
  ful_ho<-apply(ful_dat,1,calc.het)
  sub_ho<-apply(sub_dat,1,calc.het)
  # estimate allele freqs
  ful_afs<-do.call(rbind,apply(ful_dat,1,calc.afs.vcf))
  sub_afs<-do.call(rbind,apply(sub_dat,1,calc.afs.vcf))
  # save data frame
  dat<-data.frame(pop=pop,
                  preg_full=length(grep(paste0(pop,"P"),colnames(ful))),
                  preg_sub=length(grep(paste0(pop,"P"),colnames(sub))),
                  nonp_full=length(grep(paste0(pop,"NP"),colnames(ful))),
                  nonp_sub=length(grep(paste0(pop,"NP"),colnames(sub))),
                  nfem_full=length(grep(paste0(pop,"F"),colnames(ful))),
                  nfem_sub=length(grep(paste0(pop,"F"),colnames(sub))),
                  njuv_full=length(grep(paste0(pop,"J"),colnames(ful)))+
                    length(grep(paste0(pop,"DB"),colnames(ful))),
                  njuv_sub=length(grep(paste0(pop,"J"),colnames(sub))) + 
                    length(grep(paste0(pop,"DB"),colnames(sub))),
                  ho_full = mean(ful_ho,na.rm = TRUE),
                  hov_full = var(ful_ho,na.rm = TRUE),
                  ho_sub = mean(sub_ho,na.rm = TRUE),
                  hov_sub = var(sub_ho,na.rm = TRUE),
                  poly_ful = nrow(ful_afs[ful_afs$RefFreq<1,])/nrow(ful_afs)*100,
                  poly_sub = nrow(sub_afs[sub_afs$RefFreq<1,])/nrow(sub_afs)*100,
                  p_full = mean(ful_afs$RefFreq,na.rm=TRUE),
                  pv_full = var(ful_afs$RefFreq,na.rm=TRUE),
                  p_sub = mean(sub_afs$RefFreq,na.rm=TRUE),
                  pv_sub = var(sub_afs$RefFreq,na.rm=TRUE))
  return(dat)
},ful=ful_vcf,sub=sub))
write.csv(pop_summaries,"population_summaries.csv",row.names = FALSE,quote=FALSE)
```

We can merge the genetic info with the environmental info. 

```r
pop_summaries<-read.csv("population_summaries.csv")
env.data<-data.frame(t(read.csv("bayenv/env_data_raw.csv",row.names = 1)))
env.data$pop<-rownames(env.data)
pop_summaries<-merge(env.data,pop_summaries,by="pop")
```


```r
pretty_sum<-data.frame(Population = pop_summaries$pop,
                       Temperature = pop_summaries$temp,
                       Salinity = pop_summaries$salinity,
                       SeagrassDensity=pop_summaries$seagrass,
                       N_Pregnant = paste0(pop_summaries$preg_full,
                                           " (",pop_summaries$preg_sub,")"),
                       N_NonPregnant = paste0(pop_summaries$nonp_full,
                                              " (",pop_summaries$nonp_sub,")"),
                       N_Female = paste0(pop_summaries$nfem_full,
                                         " (",pop_summaries$nfem_sub,")"),
                       N_Juvenile = paste0(pop_summaries$njuv_full,
                                           " (",pop_summaries$njuv_sub,")"),
                       H_o = paste0(round(pop_summaries$ho_full,digits=3),
                                    "\u00B1",round(pop_summaries$hov_full,digits=3),
                                   " (",round(pop_summaries$ho_sub,digits=3),
                                   "\u00B1",round(pop_summaries$hov_sub,digits=3),")"),
                       MinorAlleleFrequency = paste0(
                         round(1-pop_summaries$p_full,digits=3),
                         "\u00B1",round(pop_summaries$pv_full,digits=3),
                         " (",round(1-pop_summaries$p_sub,digits=3),
                         "\u00B1",round(pop_summaries$pv_sub,digits = 3),")"),
                       PercentPolymorphicLoci = paste0(
                         round(pop_summaries$poly_ful,digits=1),
                         " (",round(pop_summaries$poly_sub,digits=1),")"),
                       stringsAsFactors = FALSE)
pretty_sum$Population[pretty_sum$Population=="FLLG"]<-"FLFW"
write.table(pretty_sum,"Table1_populationSummaries.txt",
            sep='\t',col.names = TRUE,quote=FALSE,row.names = FALSE)
```


## Minor allele frequencies


```r
locus.info<-colnames(ful_vcf[1:9])
fw.afs<-lapply(fw.list,function(pop,vcf){
  this.vcf<-cbind(vcf[,locus.info],vcf[,grep(pop,colnames(vcf))])
  this.afs<-do.call(rbind,apply(this.vcf,1,calc.afs.vcf))
}, vcf=ful_vcf)
names(fw.afs)<-c("TXFW","LAFW","ALFW","FLFW")
sw.afs<-lapply(sw.list,function(pop,vcf){
  this.vcf<-cbind(vcf[,locus.info],vcf[,grep(pop,colnames(vcf))])
  this.afs<-do.call(rbind,apply(this.vcf,1,calc.afs.vcf))
},vcf=ful_vcf)
names(sw.afs)<-sw.list
all.afs<-c(fw.afs,sw.afs)
minAF<-lapply(all.afs,function(x){
  mins<-apply(x,1,function(row){ return(min(row[4],row[6])) } )
  return(as.numeric(mins))
})
names(minAF)<-names(all.afs)
```


```r
par(mfrow=c(4,4),mar=c(2,2,1,0),oma=c(2,3,0.5,0.5))
for(i in 1:length(pop.labs)){
  if(pop.labs[i] %in% names(fw.afs)){
    color<-"cornflowerblue"
  }else{
    color<-"black"
  }
  hist(minAF[[pop.labs[i]]],ylab="",xlab="",main="",
       xlim=c(0,0.5),ylim=c(0,10000),axes=F,col=color,
       breaks=seq(0,0.5,0.05))
  axis(1,pos=0,cex.axis=2)
  if(i %in% c(1,5,9,13)){
    axis(2,pos=0,las=1,cex.axis=2,labels = seq(0,10,2),at=seq(0,10000,2000))
  }else{
    axis(2,pos=0,las=1,labels = FALSE,cex.axis=2)
    }
  mtext(pop.labs[i],3,col=color,cex=2*0.75,line=-1)
}
mtext("Minor Allele Frequency",1,outer=TRUE,cex=1.75*0.75)
mtext("Number of SNPs (x 1000)",2,outer = TRUE,cex=1.75*0.75,line=1)
```

![Minor allele frequency distributions of the full dataset (7433 SNPs) for each population. Freshwater populations are plotted in blue. The histograms show the number of SNPs with various frequencies of the reference alleles. All populations are skewed towards having small minor allele frequencies, but the TXFW and FLFW have additional reductions in genetic variation.](../figs/minorAFplot-1.png)


# Population structure of all 16 populations  {-}

## FSTs

This will calculate the pairwise $F_{ST}$s (but it's slow)


```r
pop_map<-read.delim("../fwsw_pops_map.txt",header = FALSE,stringsAsFactors = FALSE)
fst_mat<-matrix(NA,
                nrow=length(unique(pop_map$V2)),ncol=length(unique(pop_map$V2)),
                dimnames = list(pop.list,pop.list))
vcf<-parse.vcf("filter_rad_20191014@1654/14_filtered/radiator_data_20191014@1710.vcf")
colnames(vcf)<-gsub("\\-","_",colnames(vcf)) # fix individual names
# loop through each population pair
for(i in 1:(nrow(fst_mat)-1)){
  for(j in (i+1):ncol(fst_mat)){
    
    if(rownames(fst_mat)[i] != colnames(fst_mat)[j]){ # sanity check
      map1<-pop_map[pop_map[,2]==rownames(fst_mat)[i],]
      map2<-pop_map[pop_map[,2]==colnames(fst_mat)[j],]
      pwfsts<-gwsca(vcf,colnames(vcf)[1:9],map1[,1],map2[,1])
      fst_mat[i,j]<-mean(pwfsts$Fst,na.rm = TRUE)
    }
  }
}
colnames(fst_mat)<-rownames(fst_mat)<-pop.labs
write.table(fst_mat,"pairwise_fsts_full.txt",sep='\t',
            col.names = TRUE,row.names=TRUE,quote=FALSE)
```


```r
full_fsts<-read.delim("stacks/populations_whitelist/batch_2.fst_summary.tsv",
                      row.names = 1)
full_fsts<-rbind(full_fsts,TXSP=rep(NA,ncol(full_fsts))) #add the final row
Tfull_fsts<-t(full_fsts)
full_fsts[lower.tri(full_fsts)]<-Tfull_fsts[lower.tri(Tfull_fsts)] # now it's symmetric
full_fsts<-full_fsts[pop.list,pop.list]
colnames(full_fsts)<-rownames(full_fsts)<-pop.labs
fst_mat<-as.matrix(full_fsts)
```



```r
library(lattice); library(grid); library(RColorBrewer)
poporder<-read.delim("treemix/poporder")
colors<-poporder$colors
poporder<-poporder$poporder

nr<-treemix.cov.plot("treemix/fwsw_k100b",poporder)
dimnames(nr)[[1]]<-dimnames(nr)[[2]]<-pop.labs
# colors
colors<-c("black","darkgrey","grey","lightgrey","cornflowerblue")
pal<-colorRampPalette(colors)
ncol=80
cols<-pal(ncol)
rev.colors<-c("cornflowerblue","lightgrey","grey","darkgrey","black")
rev.pal<-colorRampPalette(rev.colors)
rev.cols<-rev.pal(ncol)

hm.height<-list(x=2,units="in")#2.2/3.8
hm.width<-list(x=2.4,units="in")#2.4 in RStudio/3.9

heatmaps.name<-"../figs/fst_heatmaps.png"

png(heatmaps.name,height=5,width=8,units="in",res=300)

fst.lv<-levelplot(as.matrix(fst_mat),col.regions=cols,alpha.regions=0.7,
                  scales = list(x=list(rot=90),tck = 0),xlab="",ylab="")
print(fst.lv,split=c(1,1,2,1),more=TRUE,panel.width=hm.width,
      panel.height=hm.height,cex=2)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(italic(F)[ST]), 0.2, 0, hjust=0.5, vjust=1.2,gp=gpar(cex=0.75))
trellis.unfocus()

nr.lv<-levelplot(nr,col.regions=cols,alpha.regions=0.7,
                 scales = list(x=list(rot=90),tck = 0),xlab="",ylab="")
print(nr.lv,split=c(2,1,2,1),more=FALSE,newpage=FALSE,panel.width=hm.width,
      panel.height=hm.height,cex=2)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text("Treemix", 0.2, 0, hjust=0.5, vjust=1.2,gp=gpar(cex=0.75))
trellis.unfocus()

dev.off()
```

![Heatmaps depicting population structure. In all graphs, dark colors depict similarity between populations and light grey and blue depict populations with high differentiation. The left panel shows pairwise FST values calculated by the populations module in Stacks (Catchen et al. 2013). The right panel shows covariances between populations as calculated by TreeMix (Pickrell and Pritchard 2012).](../figs/fst_heatmaps.png)


## PCAdapt {-}



```r
library(pcadapt)
```


```r
filename<-read.pcadapt("filter_rad_20191014@1654/14_filtered/radiator_data_20191014@1710.vcf",
                       type="vcf")
```

```
## No variant got discarded.
## Summary:
## 
## 	- input file:				filter_rad_20191014@1654/14_filtered/radiator_data_20191014@1710.vcf
## 	- output file:				/tmp/RtmpWYtABi/file21ac60c77de8.pcadapt
## 
## 	- number of individuals detected:	605
## 	- number of loci detected:		7433
## 
## 7433 lines detected.
## 605 columns detected.
```

```r
x<-pcadapt(filename, K=20)
plot(x,option="screeplot") #K=7
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/pcadaptChoose-1} \caption{Scree plot from PCAdapt, specifying keeping 20 PC axes.}(\#fig:pcadaptChoose)
\end{figure}


```r
pa<-pcadapt(filename,K=7)
saveRDS(pa,"fwsw_all_pcadapt.RDS")
pa.props<-round((pa$singular.values/sum(pa$singular.values))*100,2)
kable(pa.props,caption="Proportion of variation explained by all 7 of the retained PC axes in PCAdapt")
```

\begin{table}

\caption{(\#tab:pcadaptAnalyze)Proportion of variation explained by all 7 of the retained PC axes in PCAdapt}
\centering
\begin{tabular}[t]{r}
\hline
x\\
\hline
33.26\\
\hline
20.25\\
\hline
16.71\\
\hline
10.09\\
\hline
7.58\\
\hline
6.68\\
\hline
5.43\\
\hline
\end{tabular}
\end{table}

```r
ind_dat<-read.table(
  "filter_rad_20191014@1654/14_filtered/individuals.qc.stats_20191014@1654.tsv",
  header=T, stringsAsFactors = F)
pops<-ind_dat$STRATA	
grp<-pops
grp[grp=="TXFW" | grp=="LAFW" | grp=="ALFW" | grp=="FLLG"]<-"freshwater"
grp[grp!="freshwater"]<-"saltwater"

#colors
pap<-data.frame(Pop=pops,cols=pops,pch=pops,grp=grp,stringsAsFactors = F)
pap$Pop[pap$Pop == "FLLG"]<-"FLFW"
for(i in 1:nrow(pap)){
  pap[i,"cols"]<-as.character(ppi[ppi$Pop %in% pap[i,"Pop"],"cols"])
}
for(i in 1:nrow(pap)){
  pap[i,"pch"]<-as.numeric(ppi[ppi$Pop %in% pap[i,"Pop"],"pch"])
}
write.table(pap,"pcadapt_colp.txt",col.names=TRUE,sep='\t',quote=F)
```



```r
#plot
par(mfrow=c(2,3),oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(pa$scores[,1],pa$scores[,2],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),
     pch=as.numeric(pap$pch),	cex=1.5)
mtext(paste("PC1 (",pa.props[1],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC2 (",pa.props[2],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[,3],pa$scores[,4],col=alpha(pap$cols,0.5),
     bg=alpha(pap$cols,0.75),pch=as.numeric(pap$pch),cex=1.5)
mtext(paste("PC3 (",pa.props[3],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC4 (",pa.props[4],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[,5],pa$scores[,6],col=alpha(pap$cols,0.5),
     bg=alpha(pap$cols,0.75),pch=as.numeric(pap$pch),cex=1.5)
mtext(paste("PC5 (",pa.props[5],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC6 (",pa.props[6],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[grp=="freshwater",1],pa$scores[grp=="freshwater",2],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),
     pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
	cex=1.5)
mtext(paste("PC1 (",pa.props[1],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC2 (",pa.props[2],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[grp=="freshwater",3],pa$scores[grp=="freshwater",4],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),
     pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
	cex=1.5)
mtext(paste("PC3 (",pa.props[3],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC4 (",pa.props[4],"%)",sep=""),2,line = 2,cex=0.75)
plot(pa$scores[grp=="freshwater",5],pa$scores[grp=="freshwater",6],
     col=alpha(pap$cols[pap$grp=="freshwater"],0.5),
     bg=alpha(pap$cols[pap$grp=="freshwater"],0.75),
     pch=as.numeric(pap$pch[pap$grp=="freshwater"]),
	cex=1.5)
mtext(paste("PC5 (",pa.props[5],"%)",sep=""),1,line = 2,cex=0.75)
mtext(paste("PC6 (",pa.props[6],"%)",sep=""),2,line = 2,cex=0.75)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("top", legend=ppi$Pop, pch=as.numeric(ppi$pch), pt.cex=1.5,cex=0.85,
       col=alpha(ppi$cols, 0.5),pt.bg=alpha(ppi$cols,0.25), ncol=8,bty='n')
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/plotPcadaptInitial-1} \caption{Principal components analysis of genotypes in S. scovelli reveals population structure due to geographic distance and habitat type. The top set of panels show the results of the PCA with all 16 populations, color-coded by populations and point shape. The bottom set of panels show the same PCA results, but with different x- and y-axis scaling and without the saltwater populations plotted, to facilitate visualizing the differences among saltwater populations. The left panels show the first and second PC axes, whcih together account for 53.5% of the varation, the middle panels show the 3rd and 4th PC axes (another 26.8% of the variation), and the right panels the fifth and sixth axes (another % of variation).}(\#fig:plotPcadaptInitial)
\end{figure}




## Admixture scree plots {-}



```r
admixK<-read.delim("admixture/K_CVs.txt",header = FALSE)
admixK$K<-as.numeric(gsub(".*\\(K=(\\d+)\\).*","\\1",admixK$V1))
admixK$CV<-as.numeric(gsub("^.*\\: (\\d+\\.\\d+)$","\\1",admixK$V1))

admixK<-admixK[order(admixK$K),]

plot(admixK$K,admixK$CV,pch=19,type = "b",lty=1,xlab = "K",ylab="CV",las=1,lwd=2)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixScree-1} \caption{Admixture screeplot for K=1 through K=16. The coefficient of variation (CV) is shown on the y-axis.}(\#fig:admixScree)
\end{figure}

Looks like $K=5$ or $K=7$ are the best, let's look at all of the $K=2$ through $K=7$.


```r
library(RColorBrewer)
famfile<-"admixture/fwsw_all_filt.fam"

poporderFile<-"treemix/poplist"
poporderDF<-read.table(poporderFile,col.names = c("Pop"),stringsAsFactors = F)
poporderDF$orderNum<-1:nrow(poporderDF)
```


```r
admixPlotting<-function(qfile,K,famfile="admixture/fwsw_all_filt.fam",
                        poporder=poporderDF){
  # read files in 
  famTable<- read.table(famfile, col.names = 
                        c("Pop","Ind","Father","Mother","Sex","phenotype"),
                      stringsAsFactors = F)[1:2]

  qtbl<-read.table(qfile,stringsAsFactors = F)
  
  # create useful tables
  mergedAdmixtureTable <- cbind(qtbl, famTable)
  mergedAdmixTabOrderNs <- merge(mergedAdmixtureTable,poporder,by="Pop")
  ordered <- mergedAdmixTabOrderNs[order(mergedAdmixTabOrderNs$orderNum),]
  
  plotting.structure(ordered[,1:(ncol(ordered)-2)],k = K,
                     pop.order = poporder$Pop,make.file = FALSE)
  admix<-ordered[,1:(ncol(ordered)-2)]
  return(admix)
}
```


```r
admixK2<-admixPlotting("admixture/fwsw_all_filt.2.Q",2)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixK2-1} \caption{Admixture plot for K=2. The colors represent different genetic populations.}(\#fig:admixK2)
\end{figure}


```r
admixK3<-admixPlotting("admixture/fwsw_all_filt.3.Q",3)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixK3-1} \caption{Admixture plot for K=3. The colors represent different genetic populations.}(\#fig:admixK3)
\end{figure}

```r
admixK4<-admixPlotting("admixture/fwsw_all_filt.4.Q",4)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixK4-1} \caption{Admixture plot for K=4. The colors represent different genetic populations.}(\#fig:admixK4)
\end{figure}


```r
admixK5<-admixPlotting("admixture/fwsw_all_filt.5.Q",5)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixk5-1} \caption{Admixture plot for K=5. The colors represent different genetic populations.}(\#fig:admixk5)
\end{figure}

```r
write.table(admixK5,"admixture/admixK5.txt",sep = '\t',
            quote = FALSE,col.names = TRUE,row.names = FALSE)
```


```r
admixK6<-admixPlotting("admixture/fwsw_all_filt.6.Q",6)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixk6-1} \caption{Admixture plot for K=6. The colors represent different genetic populations.}(\#fig:admixk6)
\end{figure}


```r
admixK7<-admixPlotting("admixture/fwsw_all_filt.7.Q",7)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/admixk7-1} \caption{Admixture plot for K=7. The colors represent different genetic populations.}(\#fig:admixk7)
\end{figure}

```r
write.table(admixK7,"admixture/admixK7.txt",sep = '\t',quote = FALSE,col.names = TRUE,row.names = FALSE)
```


## Treemix analysis {-}



```r
poporder<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST",
            "ALFW","FLSG","FLKB","FLFD","FLSI","FLAB",
            "FLPB","FLHB","FLCC","FLLG")
colors<-poporder
colors[colors %in% "FLLG"]<-grp.colors[6]
colors[colors %in% c("FLPB","FLHB","FLCC")]<-grp.colors[6]
colors[colors %in% c("FLAB")]<-grp.colors[5]
colors[colors %in% c("FLSI","FLFD","FLKB","FLSG")]<-grp.colors[3]
colors[colors %in% c("ALST","ALFW","LAFW")]<-grp.colors[2]
colors[colors %in% c("TXSP","TXCC","TXFW","TXCB")]<-grp.colors[1]
write.table(cbind(poporder,colors),"poporder",quote=F,sep='\t')
```

### Tree with no migration edges and no root

It's informative to plot the treemix tree that does not have any migration edges added for comparison later.


```r
# unrooted tree
library(ape)
tre<-read.tree("treemix/unrooted_consensus.newick")
png("../figs/treemix_unrooted_consense.png",height=8,width=8,units="in",res=300)
plot(tre)
dev.off()
write.tree(tre,'treemix/unrooted_consensus.newick')
```


```r
knitr::include_graphics('../figs/treemix_unrooted_consense.png')
```

\begin{figure}[H]
\includegraphics[width=8in,]{../figs/treemix_unrooted_consense} \caption{The consensus tree from running Treemix without any migration edges  and no root.}(\#fig:treemixUnrooted)
\end{figure}


### Tree with FLAB as root

We ran Treemix assuming the FLAB was the root with 100 bootstrap replicates. We then used PHYLIP's consense program to build a consensus tree, assuming that it was a rooted tree. We also re-saved it to file so it would be in one line and thus compatible with treemix. 


```r
# rooted tree
tre<-read.tree("treemix/rooted_consensus.newick")
png("../../figs/treemix_rooted_consense.png",height=8,width=8,units="in",res=300)
plot(tre)
dev.off()
write.tree(tre,'treemix/rooted_consensus.newick')
```


```r
knitr::include_graphics('../figs/treemix_unrooted_consense.png')
```

\begin{figure}[H]
\includegraphics[width=8in,]{../figs/treemix_unrooted_consense} \caption{The consensus tree from running Treemix without any migration edges and FLAB as root.}(\#fig:treemixRooted)
\end{figure}


### Choosing the optimal number of migration edges

To choose the optimal number of migration edges, we used the R package optM [@fitak_optm:_2019].



```r
library(OptM)
tmOpt<-optM("treemix/migrations/")
evanno_treemix(tmOpt)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/plotOptM-1} \caption{Plot showing the comparison of Treemix number of migration edges using the Evanno method.}(\#fig:plotOptM)
\end{figure}


```r
lliks<-list.files(path="treemix/migrations/",pattern=".llik",full.names = TRUE)
likes<-data.frame(do.call(rbind,lapply(lliks,function(file){
  likdat<-read.delim(file,header=FALSE,row.names=1,sep=':')
  migs<-as.numeric(gsub(".*m(\\d).*","\\1",file))
  return(cbind(migs=migs,loglikelihood=likdat[2,]))
})))
rownames(likes)<-lliks

llikMean<-tapply(likes$loglikelihood,likes$migs,mean)
llikSEM<-tapply(likes$loglikelihood,likes$migs,function(x){
  return(sqrt(var(x)/length(x)))
})

plot(0:5,as.numeric(llikMean),pch=19,cex=2,
     xlab="Number of migration edges",
     ylab="Log likelihood",ylim=c(-8000,1200))
points(0:5,llikMean,lwd=2,type='l')
arrows(x0 = 0:5,y0=c(llikMean-llikSEM),
       x1=0:5,y1=c(llikMean+llikSEM),code=3,angle=0,lwd=2)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/treemixLlik-1} \caption{Average log likelihoods of treemix bootstrap replicates with 0 through 5 migration edges. Shown are the means (of 100 bootstraps) and the standard error of the mean.}(\#fig:treemixLlik)
\end{figure}


```r
m1s<-list.files(pattern="m1.*treeout",path="treemix/migrations/",full.names = TRUE)
edges<-do.call(rbind,lapply(m1s,function(file){
  treeout<-scan(file,what="character",sep='\n')
  edge<-do.call(rbind,strsplit(treeout[2],' '))
  return(edge)
}))
rownames(edges)<-m1s
```


```r
starts<-gsub("[[:digit:]]","\\1",edges[,5])
starts<-gsub(":","",starts)
starts<-gsub("\\.","",starts)
starts<-gsub("e-","",starts)
starts<-gsub(",",", ",starts)
starTab<-data.frame(tree=names(summary(as.factor(starts))),
                    counts=summary(as.factor(starts)))
rownames(starTab)<-NULL
kable(starTab,"latex",booktabs=TRUE,row.names=FALSE,longtable=TRUE,
      col.names = c("tree location","number of bootstraps"),
      caption="Number of bootstraps with the one migration edge beginning at each of these points on the population trees.") %>%
  kable_styling(latex_options=c("HOLD_position","repeat_header"))  %>%
  column_spec(1,width = "30em")
```


\begin{longtable}[t]{>{\raggedright\arraybackslash}p{30em}r}
\caption{(\#tab:M1Starts)Number of bootstraps with the one migration edge beginning at each of these points on the population trees.}\\
\toprule
tree location & number of bootstraps\\
\midrule
\endfirsthead
\caption[]{(\#tab:M1Starts)Number of bootstraps with the one migration edge beginning at each of these points on the population trees. \textit{(continued)}}\\
\toprule
tree location & number of bootstraps\\
\midrule
\endhead
\
\endfoot
\bottomrule
\endlastfoot
(((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLHB, ((FLCC, FLLG), FLPB))) & 1\\
(((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLHB, (FLCC, (FLPB, FLLG)))) & 1\\
(((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, ((FLHB, FLCC), FLPB))) & 11\\
(((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, (FLCC, (FLHB, FLPB)))) & 9\\
(((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, (FLHB, (FLCC, FLPB)))) & 28\\
\addlinespace
(((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, (FLPB, (FLCC, FLHB)))) & 7\\
((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), ((FLKB, FLSG), (FLSI, FLFD))), (FLLG, ((FLHB, FLCC), FLPB))) & 1\\
((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), ((FLKB, FLSG), (FLSI, FLFD))), (FLLG, (FLCC, (FLHB, FLPB)))) & 1\\
((((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST), ((FLKB, FLSG), (FLSI, FLFD))), (FLLG, (FLHB, (FLCC, FLPB)))) & 8\\
(((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, (FLCC, (FLHB, FLPB)))) & 6\\
\addlinespace
(((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, (FLCC, (FLPB, FLHB)))) & 1\\
(((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), (FLKB, FLSG)), (FLSI, FLFD)), (FLLG, (FLHB, (FLCC, FLPB)))) & 7\\
(((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), (FLSG, FLKB)), (FLSI, FLFD)), (FLLG, (FLCC, (FLHB, FLPB)))) & 1\\
((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), ((FLKB, FLSG), (FLSI, FLFD))), (FLLG, (FLCC, (FLHB, FLPB)))) & 1\\
((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), ((FLKB, FLSG), (FLSI, FLFD))), (FLLG, (FLHB, (FLCC, FLPB)))) & 1\\
\addlinespace
((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), ((FLKB, FLSG), (FLSI, FLFD))), (FLLG, (FLPB, (FLCC, FLHB)))) & 1\\
((((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST), (FLSG, (FLKB, (FLSI, FLFD)))), (FLLG, (FLHB, (FLCC, FLPB)))) & 1\\
(((FLKB, (FLSG, ((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST))), (FLSI, FLFD)), (FLLG, (FLCC, (FLHB, FLPB)))) & 1\\
(((FLKB, (FLSG, ((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST))), (FLSI, FLFD)), (FLLG, (FLHB, (FLCC, FLPB)))) & 5\\
(((FLKB, (FLSG, ((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST))), (FLSI, FLFD)), (FLLG, (FLHB, (FLPB, FLCC)))) & 1\\
\addlinespace
(((FLKB, (FLSG, ((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST))), (FLSI, FLFD)), (FLLG, (FLCC, (FLHB, FLPB)))) & 1\\
(((FLKB, (FLSG, ((LAFW, (ALFW, ((TXCB, (TXCC, TXSP)), TXFW))), ALST))), (FLSI, FLFD)), (FLLG, (FLHB, (FLCC, FLPB)))) & 2\\
(((FLKB, (FLSI, FLFD)), (FLSG, ((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST))), (FLLG, (FLCC, (FLHB, FLPB)))) & 1\\
(((FLSG, ((((TXCB, (TXCC, TXSP)), TXFW), (ALFW, LAFW)), ALST)), (FLKB, (FLSI, FLFD))), (FLLG, (FLHB, (FLCC, FLPB)))) & 1\\
(FLLG, ((FLHB, FLCC), FLPB)) & 1\\
\addlinespace
FLLG & 1\\
TXFW & 1\\*
\end{longtable}


```r
stops<-gsub("[[:digit:]]","\\1",edges[,6])
stops<-gsub(":","",stops)
stops<-gsub("\\.","",stops)
stops<-gsub("e-","",stops)
kable(summary(as.factor(stops)),"latex",booktabs=TRUE,
      col.names = c("number of bootstraps"),
      caption="Number of bootstraps with the one migration edge ending at each of these points on the population trees.")%>%
  kable_styling(latex_options=c("HOLD_position"))
```

\begin{table}[H]

\caption{(\#tab:M1Stops)Number of bootstraps with the one migration edge ending at each of these points on the population trees.}
\centering
\begin{tabular}[t]{lr}
\toprule
  & number of bootstraps\\
\midrule
(TXCB,(TXCC,TXSP)) & 1\\
FLCC & 1\\
FLPB & 99\\
\bottomrule
\end{tabular}
\end{table}

The overwhelming majority of bootstrapped trees with one migration edge have that edge leading from the main branch of the tree to FLPB (see the above tables for the common start and end points).




```r
bestM1<-gsub("\\.llik","",
             rownames(likes[likes$migs==1,])[
               which.max(likes$loglikelihood[likes$migs==1])])
png("../figs/treemix_comparison.png",height = 4,width=8,units="in",res=300)
par(mfrow=c(1,3),mar=c(1,2,1,2),oma=c(1,6,1,2),xpd=TRUE)
t0<-plot_tree("treemix/migrations/fwsw_FLAB_m0",scale=T,mbar=F,cex = 1.5,
              lwd=2,disp=0.002,scadj=0.05)
t1<-plot_tree("treemix/migrations/fwsw_FLAB_m1",scale=T,mbar=T,cex = 1.5,
              lwd=2,mig_left=FALSE,disp=0.0002,scadj=0.05)
t2<-plot_tree(bestM1,scale=T,mbar=T,cex = 1.5,
              lwd=2,mig_left=FALSE,disp=0.0002,scadj=0.05)
dev.off()
```


```r
knitr::include_graphics('../figs/treemix_comparison.png')
```

\begin{figure}[H]
\includegraphics[width=0.9\linewidth,]{../figs/treemix_comparison} \caption{The plot of the tree with FLAB as root but no migration edges (left) compared to the tree with FLAB as root and two migration edges. The drift parameter is plotted on the x-axis, and migration edges are colored based on the migration weight.}(\#fig:treemixCompare)
\end{figure}




### Threepop and fourpop analysis
 
We investigated the two migration edges that are in the best Treemix model using threepop and fourpop analyses. In the threepop analysis, significantly negative f3 statistics mean that the first pop in the list (A in A;B,C) is admixed [@pickrellInferencePopulationSplits2012; @reichReconstructingIndianPopulation2009]. Therefore, with the threepop analysis we want to look for the tree (A;B,C) where A corresponds to the end of an arrow and (B,C) corresponds to where an arrow begins. In the fourpop analysis, a significantly non-zero value indicates gene flow in the tree [@pickrellInferencePopulationSplits2012; @reichReconstructingIndianPopulation2009]. 


```r
threepop<-data.frame(do.call(rbind,
                             strsplit(grep(
                               ";",readLines("treemix/fwsw_threepop.txt"),
                               value = TRUE),' ')),
                     stringsAsFactors = FALSE)
colnames(threepop)<-c("pops","f3_stat","f3_se","f3_z")

fourpop<-data.frame(do.call(rbind,
                            strsplit(grep(";",
                                          readLines("treemix/fwsw_fourpop.txt"),
                                          value = TRUE),' ')),
                     stringsAsFactors = FALSE)
colnames(fourpop)<-c("pops","f4_stat","f4_se","f4_z")
```


The migration edge we investigated indicated potential migration from the ancestral Florida branch to the FLPB branch. We first investigated all of the three-population trees with FLPB and FLAB and found that the majority of these trees are positive, which suggests that FLPB is not admixed with the other Florida populations. Several trees are negative, but their standard errors overlap with zero, which suggests that FLPB may experience some admixtre with those populations. Unsurprisingly, given the admixture and pcadapt results, these trees are those with other Atlantic Florida populations (FLHB, FLCC, and FLLG/FLFW). 


```r
FLPB_edge<-threepop[grep("FLPB;.*FLAB.*",threepop$pops),]
FLPB_edge$f3_stat<-as.numeric(FLPB_edge$f3_stat)
FLPB_edge$f3_se<-as.numeric(FLPB_edge$f3_se)
FLPB_edge$lowSE<-FLPB_edge$f3_stat-FLPB_edge$f3_se
FLPB_edge$uppSE<-FLPB_edge$f3_stat+FLPB_edge$f3_se
FLPB_edge$f3_z<-as.numeric(FLPB_edge$f3_z)

ymax<-max(abs(c(FLPB_edge$lowSE,FLPB_edge$uppSE)))+0.001

# all the relevant trees
plot(1:nrow(FLPB_edge),as.numeric(FLPB_edge$f3_stat),
     ylim=c(ymax*-1,ymax),axes=FALSE,
     xlab="",ylab="f3 statistic +/- SE")
abline(h=0,lty=2)
arrows(1:nrow(FLPB_edge),FLPB_edge$lowSE,
       1:nrow(FLPB_edge),FLPB_edge$uppSE,code=3,angle=0)
axis(2,at=seq(-0.05,0.05,0.01),las=1)
axis(1,pos=-0.01,at=1:nrow(FLPB_edge),
     labels = FLPB_edge$pops,las=2,cex.axis=0.75)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/f3Edge1-1} \caption{Plot of the f3 statistic for three-population trees including both FLPB and FLAB. Error bars show the standard errors generated from block jackknifes. Trees indicative of admixture in FLPB are those with significantly negative f3 statistics.}(\#fig:f3Edge1)
\end{figure}

To investigate the results of the fourpop analysis for this migration edge, we focused on four-population trees including both FLAB and FLPB and other Florida populations. All of these four population trees show evidence of admixture, which is unsurprising given that these trees include populations that are in the same population clusters in the admixture and pcadapt results.


```r
f4s<-fourpop[c(grep("FLPB,FL.*;FLAB,FL.*",fourpop$pops),
               grep("FL.*,FLPB;FLAB,FL.*",fourpop$pops),
               grep("FLPB,FL.*;FL.*,FLAB",fourpop$pops),
               grep("FL.*,FLPB;FL.*,FLAB",fourpop$pops)),]
kable(f4s,"latex",booktabs=TRUE,longtable=TRUE,
      caption = "Four-population trees with FLPB and FLAB and their f4 statistic, standard error, z-score, and p-value." )%>%
  kable_styling(latex_options=c("HOLD_position","repeat_header"))
```


\begin{longtable}[t]{lllll}
\caption{(\#tab:f4Edge1)Four-population trees with FLPB and FLAB and their f4 statistic, standard error, z-score, and p-value.}\\
\toprule
  & pops & f4\_stat & f4\_se & f4\_z\\
\midrule
\endfirsthead
\caption[]{(\#tab:f4Edge1)Four-population trees with FLPB and FLAB and their f4 statistic, standard error, z-score, and p-value. \textit{(continued)}}\\
\toprule
  & pops & f4\_stat & f4\_se & f4\_z\\
\midrule
\endhead
\
\endfoot
\bottomrule
\endlastfoot
7569 & FLSG,FLPB;FLAB,FLHB & 0.0108426 & 0.0006659 & 16.2826\\
7572 & FLSG,FLPB;FLAB,FLCC & 0.0108142 & 0.000676126 & 15.9943\\
7575 & FLSG,FLPB;FLAB,FLLG & 0.00958893 & 0.000746125 & 12.8516\\
7674 & FLKB,FLPB;FLAB,FLHB & 0.0108759 & 0.000673324 & 16.1526\\
7677 & FLKB,FLPB;FLAB,FLCC & 0.0108676 & 0.000684369 & 15.8797\\
\addlinespace
7680 & FLKB,FLPB;FLAB,FLLG & 0.00966079 & 0.000769578 & 12.5534\\
7734 & FLFD,FLPB;FLAB,FLHB & 0.0111537 & 0.000685198 & 16.2781\\
7737 & FLFD,FLPB;FLAB,FLCC & 0.0111212 & 0.000694793 & 16.0064\\
7740 & FLFD,FLPB;FLAB,FLLG & 0.00984951 & 0.000772256 & 12.7542\\
7764 & FLSI,FLPB;FLAB,FLHB & 0.0114181 & 0.000690251 & 16.5419\\
\addlinespace
7767 & FLSI,FLPB;FLAB,FLCC & 0.01141 & 0.00070008 & 16.2981\\
7770 & FLSI,FLPB;FLAB,FLLG & 0.0100674 & 0.000766841 & 13.1284\\
13029 & FLSG,FLPB;FLAB,FLHB & 0.0108426 & 0.0006659 & 16.2826\\
13032 & FLSG,FLPB;FLAB,FLCC & 0.0108142 & 0.000676126 & 15.9943\\
13035 & FLSG,FLPB;FLAB,FLLG & 0.00958893 & 0.000746125 & 12.8516\\
\addlinespace
13134 & FLKB,FLPB;FLAB,FLHB & 0.0108759 & 0.000673324 & 16.1526\\
13137 & FLKB,FLPB;FLAB,FLCC & 0.0108676 & 0.000684369 & 15.8797\\
13140 & FLKB,FLPB;FLAB,FLLG & 0.00966079 & 0.000769578 & 12.5534\\
13194 & FLFD,FLPB;FLAB,FLHB & 0.0111537 & 0.000685198 & 16.2781\\
13197 & FLFD,FLPB;FLAB,FLCC & 0.0111212 & 0.000694793 & 16.0064\\
\addlinespace
13200 & FLFD,FLPB;FLAB,FLLG & 0.00984951 & 0.000772256 & 12.7542\\
13224 & FLSI,FLPB;FLAB,FLHB & 0.0114181 & 0.000690251 & 16.5419\\
13227 & FLSI,FLPB;FLAB,FLCC & 0.01141 & 0.00070008 & 16.2981\\
13230 & FLSI,FLPB;FLAB,FLLG & 0.0100674 & 0.000766841 & 13.1284\\
7464 & FLSG,FLPB;FLKB,FLAB & 0.00230133 & 0.000211957 & 10.8575\\
\addlinespace
7509 & FLSG,FLPB;FLFD,FLAB & 0.00177997 & 0.000177424 & 10.0323\\
7538 & FLSG,FLPB;FLSI,FLAB & 0.00179893 & 0.000166984 & 10.7731\\
7614 & FLKB,FLPB;FLFD,FLAB & 0.00185955 & 0.000176447 & 10.5389\\
7643 & FLKB,FLPB;FLSI,FLAB & 0.00188169 & 0.000170658 & 11.0261\\
7703 & FLFD,FLPB;FLSI,FLAB & 0.00222764 & 0.000196885 & 11.3145\\
\addlinespace
12924 & FLSG,FLPB;FLKB,FLAB & 0.00230133 & 0.000211957 & 10.8575\\
12969 & FLSG,FLPB;FLFD,FLAB & 0.00177997 & 0.000177424 & 10.0323\\
12998 & FLSG,FLPB;FLSI,FLAB & 0.00179893 & 0.000166984 & 10.7731\\
13074 & FLKB,FLPB;FLFD,FLAB & 0.00185955 & 0.000176447 & 10.5389\\
13103 & FLKB,FLPB;FLSI,FLAB & 0.00188169 & 0.000170658 & 11.0261\\
\addlinespace
13163 & FLFD,FLPB;FLSI,FLAB & 0.00222764 & 0.000196885 & 11.3145\\*
\end{longtable}


```r
f4s$lowSE<-as.numeric(f4s$f4_stat)-as.numeric(f4s$f4_se)
f4s$uppSE<-as.numeric(f4s$f4_stat)+as.numeric(f4s$f4_se)
ymax<-max(abs(c(f4s$lowSE,FLPB_edge$uppSE)))+0.001

# all the relevant trees
plot(1:nrow(f4s),as.numeric(f4s$f4_stat),
     ylim=c(ymax*-1,ymax),axes=FALSE,
     xlab="",ylab="f4 statistic +/- SE")
abline(h=0,lty=2)
arrows(1:nrow(f4s),f4s$lowSE,
       1:nrow(f4s),f4s$uppSE,code=3,angle=0)
axis(2,at=seq(-0.05,0.05,0.01),las=1)
axis(1,pos=-0.01,at=1:nrow(f4s),
     labels = f4s$pops,las=2,cex.axis=0.75)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/f4Edge1Fig-1} \caption{Plot of the f4 statistic for four-population trees including both FLPB and FLAB and other Florida populations. Error bars show the standard errors generated from block jackknifes. Trees indicative of admixture in the tree are those with significantly non-zero f4 statistics.}(\#fig:f4Edge1Fig)
\end{figure}


## Make figure


```r
grp.colors<-c('#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837') 
grp7colors<-c('#762a83','#9970ab','#c2a5cf','#d9f0d3','#a6dba0','#5aae61','#1b7837')
grp5colors<-c('#762a83','#c2a5cf','#a6dba0','#5aae61','#1b7837')
```


```r
#admixture 
admixK5<-read.delim("admixture/admixK5.txt",header = TRUE)
admixK7<-read.delim("admixture/admixK7.txt",header = TRUE)
#pcadapt
pa<-readRDS("fwsw_all_pcadapt.RDS")
pa.props<-round((pa$singular.values/sum(pa$singular.values))*100,2)
pap<-read.delim("pcadapt_colp.txt",sep='\t')
# map
library(jpeg)
img<-readJPEG("all_sites_map.jpeg")

# stuff for treemix
source("../R/203_treemix_plotting_funcs.R") #I've modified the functions from treemix
library(lattice); library(grid); library(RColorBrewer)
poporder<-read.delim("treemix/poporder")
colors<-poporder$colors
poporder<-poporder$poporder

d <- read.table("treemix/migrations/fwsw_FLAB_m1.vertices.gz", 
                as.is  = T, comment.char = "", quote = "")
branch.cols<-rep("black",nrow(d))
branch.cols[d[,2] %in% c("TXFW","ALFW","LAFW","FLLG")]<-"cornflowerblue"

tip.names<-as.vector(d[d[,5] == "TIP",2])
tip.names<-data.frame(Original=tip.names,Replacement=tip.names,stringsAsFactors = FALSE)
tip.names$Replacement[tip.names$Replacement=="FLLG"]<-"FLFW"


col_vector<-c(TXSP='#762a83',TXCC='#762a83',TXFW="#2166ac",TXCB='#762a83',LAFW="#2166ac",ALST='#af8dc3',
              ALFW="#2166ac",FLSG='#e7d4e8',FLKB='#e7d4e8',FLFD='#e7d4e8',FLSI='#e7d4e8',FLAB='#7fbf7b',
              FLPB='#1b7837',FLHB='#1b7837',FLCC='#1b7837',FLFW="#2166ac")
```



```r
npop<-length(pop.list)
pseq<-1:npop
m<-matrix(c(rep(1,16),rep(2,6),
            3:18,rep(2,6),
            19:34,rep(2,6),
            rep(35,8),rep(36,8),rep(37,6)),
          nrow=4,ncol=npop+6,byrow = T)
jpeg("../figs/NewPopStructure_v1.jpeg",res=300,height=8,width=10,units="in")
#set the layout
layout(mat=m,heights=c(6,1,1,6))
#MAP
#open an empty plot window with coordinates
par(oma=c(1.5,3.5,1,2),mar=c(0,0,0,0),xpd=NA)
plot(1:14,ty="n",axes=FALSE,xlab="",ylab="",xpd=TRUE)
#specify the position of the image through bottom-left and top-right coords
rasterImage(img,1,1,14,14,xpd=TRUE)
text(x = 1.5,y=11.5,"A",cex = 2,font =2)
text(x = 0.5,y=1.2,"C",cex = 2,font =2)

#Treemix
par(mar=c(0,0,0,1))
t2<-plot_tree("treemix/migrations/fwsw_FLAB_m1",scale=T,mbar=T,cex = 1.5,
              lwd=2,mig_left=FALSE,disp=0.0002,scadj=0.05,
              tip.order=tip.names,branch.cols = branch.cols,
              tip.col=col_vector[tip.names$Replacement])
text(x = 0.035,y=0.9,"B",cex = 2,font =2)
#STRUCTURE
par(mar=c(1,0,0,0))
plotting.structure(admixK5,5,pop.order = poporder, pop.list, make.file=FALSE, 
                   xlabcol = all.colors,plot.new=FALSE,
                   colors=grp5colors[c(3,5,4,1,2)],xlabel=FALSE,
                   ylabel=expression(atop(italic(K)==5)),lab.cex=0.85)
plotting.structure(admixK7,7,pop.order = poporder,pop.labs, make.file=FALSE,
                   plot.new=FALSE,colors=grp7colors[c(6,7,1,2,5,3,4)],xlabel=TRUE,
                   xlabcol = all.colors,
                   ylabel=expression(atop(italic(K)==7)),lab.cex=0.85)

#PCADAPT
par(mar=c(2,2,2,2))
plot(pa$scores[,1],pa$scores[,2],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),
     pch=as.numeric(pap$pch),	cex=3,bty="L",xlab="",ylab="",cex.axis=1.5)

mtext(paste("PC1 (",pa.props[1],"%)",sep=""),1,line = 2.5,cex=1)
mtext(paste("PC2 (",pa.props[2],"%)",sep=""),2,line = 2.5,cex=1)
text(x = -0.07,y=0.1,"D",cex = 2,font =2)

plot(pa$scores[,3],pa$scores[,4],col=alpha(pap$cols,0.5),bg=alpha(pap$cols,0.75),
     pch=as.numeric(pap$pch),
     cex=3,	bty="L",xlab="",ylab="",cex.axis=1.5)

mtext(paste("PC3 (",pa.props[3],"%)",sep=""),1,line = 2.5,cex=1)
mtext(paste("PC4 (",pa.props[4],"%)",sep=""),2,line = 2.5,cex=1)

plot(1:10,ty="n",axes=FALSE,xlab="",ylab="",xpd=TRUE)
legend("bottom", legend=ppi$Pop, pch=as.numeric(ppi$pch), pt.cex=3,cex=1.5,
       col=alpha(ppi$cols, 0.5),pt.bg=alpha(ppi$cols,0.25), ncol=2,bty='n')
dev.off()
```

![Population Structure summary plot (Figure 1 in main document)](../figs/NewPopStructure_v1.jpeg)


## Make table



```r
full_fsts<-read.delim("stacks/populations_whitelist/batch_2.fst_summary.tsv",
                      row.names = 1)
full_fsts<-rbind(full_fsts,TXSP=rep(NA,ncol(full_fsts))) #add the final row
Tfull_fsts<-t(full_fsts)
full_fsts[lower.tri(full_fsts)]<-Tfull_fsts[lower.tri(Tfull_fsts)] # now it's symmetric
full_fsts<-full_fsts[pop.list,pop.list]
colnames(full_fsts)<-rownames(full_fsts)<-pop.labs
fst_mat<-as.matrix(full_fsts)
```


```r
cov<-read.table(gzfile("treemix/unrooted/fwsw_ML_consensus.cov.gz"), as.is = T, head = T, quote = "", 
                comment.char = "")
#reorder
covplot <- data.frame(matrix(nrow = nrow(cov), ncol = ncol(cov)))
for(i in 1:length(pop.list)){
  for( j in 1:length(pop.list)){
    
    covplot[i, j] = cov[which(names(cov)==pop.list[i]), which(names(cov)==pop.list[j])]
    rownames(covplot)[i]<-pop.list[i]
    colnames(covplot)[j]<-pop.list[j]
  }
}
covplot<-as.matrix(covplot)
```

```r
library(xlsx); library(RColorBrewer); library(scales)

table2<-fst_mat
table2[lower.tri(table2)]<-covplot[lower.tri(covplot)]
diag(table2)<-diag(covplot)
table2<-round(table2,digits = 4)

# first export the data
sheetName <- "FstCov"
file<-"table2_fst_cov.xlsx"
write.xlsx(table2,file,sheetName = sheetName)

wb<-loadWorkbook(file)
sheets <- getSheets(wb)               
sheet <- sheets[[sheetName]]          
rows <- getRows(sheet, rowIndex=2:(nrow(table2)+1)) # 1st row is headers
cells <- getCells(rows, colIndex = 2:(ncol(table2)+1)) # 1st col is rownames         
values<-lapply(cells,getCellValue)

# set the colors
pal<-colorRampPalette(c("#deebf7","#3182bd"))
cols<-matrix(nrow = nrow(table2),ncol=nrow(table2))
cols[upper.tri(cols)]<-pal(10)[as.numeric(cut(table2[upper.tri(table2)],breaks = 10))]
pal<-colorRampPalette(c("#78c679","#f7fcb9"))
cols[lower.tri(cols)]<-pal(10)[as.numeric(cut(table2[lower.tri(table2)],breaks = 10))]
pal<-colorRampPalette(c("#f7f7f7","#969696"))
diag(cols)<-pal(10)[as.numeric(cut(diag(table2),breaks = 10))]

for(i in 1:nrow(table2)){
  for(j in 1:ncol(table2)){
    csij<-CellStyle(wb,fill=Fill(foregroundColor = alpha(cols[i,j])))
    setCellStyle(cells[[paste(i+1,j+1,sep=".")]],cellStyle = csij)
  }
}

saveWorkbook(wb, file)
```



# Fst Outliers  {-}

The alignments were done with a preliminary genome assembly that is different from the published, updated one. The data that I need to convert are the vcf and the stacks files. 




```r
convert.agp<-function(locus=NULL,old.agp,old.scf,new.agp,scf.agp,
                      chr=NULL,bp=NULL,id=NULL){
  
  if(!is.null(locus)){
    chr<-locus$`#CHROM`
    bp<-locus$POS
    id<-locus$ID
  }else{
    bp<-as.numeric(unlist(bp))
    chr<-as.character(chr)
    id<-as.character(id)
  }
  component<-as.data.frame(old.agp[old.agp$object == chr & 
                                     old.agp$object_beg <= bp & 
                                     old.agp$object_end >= bp,],
                           stringsAsFactors=FALSE)
  if(nrow(component)>0){
    # it's found on one of the LGs
    comp.id<-component$component_id
    if(comp.id != 100){
      #make sure it's an actual scaffold as a component
      comp.bp<-as.numeric(as.character(component$component_beg))+
        (bp-as.numeric(as.character(component$object_beg)))-1
      #sanity check - is it a reasonable size?
      if(comp.bp<as.numeric(as.character(component$component_end))){ 
        updated<-new.agp[new.agp$component_id%in%comp.id & 
                  as.numeric(as.character(new.agp$component_beg)) <=comp.bp & 
                  as.numeric(as.character(new.agp$component_end)) >= comp.bp,]
        if(nrow(updated)==0){ #if you didn't find it, check scaffold
          updated<-scf.agp[scf.agp$object%in%comp.id & 
                  as.numeric(as.character(scf.agp$object_beg)) <=comp.bp & 
                  as.numeric(as.character(scf.agp$object_end)) >= comp.bp,]
          updated.bp<-comp.bp
          updated.chr<-as.character(comp.id)
        } else{
          updated.bp<-updated$object_beg+comp.bp
          updated.chr<-as.character(updated$object)  
        }
      }else {
        print("WARNING: position in component larger than component")
        updated.bp<-comp.id
        updated.chr<-as.character(comp.id)
      }
    }else{
      print(paste("WARNING: locus ",id, " is not on a scaffold",sep=""))
      updated.bp<-bp
      updated.chr<-NA
    }
    out<-data.frame(Locus=id,OrigChr=chr,OrigBP=bp,
                    NewChr=updated.chr,NewBP=updated.bp,
                    stringsAsFactors = FALSE)
  }else{
    #it's not on an LG - let's check the scaffolds
    component<-as.data.frame(old.scf[old.scf$object == chr & 
                                       old.scf$object_beg <= bp & 
                                       old.scf$object_end >= bp,],
                             stringsAsFactors=FALSE)
    if(nrow(component)>0){
      #then we found it
      #check to make sure my bp makes sense
      if(bp < max(old.scf[old.scf$object==chr,"object_end"])){
        comp.bp<-bp
        comp.id<-as.character(chr)
        #look for it in the new assembly
        updated<-new.agp[new.agp$component_id%in%comp.id & 
                  as.numeric(as.character(new.agp$component_beg)) <=comp.bp &   
                  as.numeric(as.character(new.agp$component_end)) >= comp.bp,]
        if(nrow(updated)==0){ #if you didn't find it, check scaffold
          updated<-scf.agp[scf.agp$object%in%comp.id & 
                  as.numeric(as.character(scf.agp$object_beg)) <=comp.bp & 
                  as.numeric(as.character(scf.agp$object_end)) >= comp.bp,]
          updated.bp<-comp.bp
          updated.chr<-as.character(comp.id)
        } else{
          updated.bp<-updated$object_beg+comp.bp
          updated.chr<-as.character(updated$object)  
        }
      } else {
          print("WARNING: position in scaffold larger than scaffold")
          updated.bp<-NA
          updated.chr<-NA
      }
      out<-data.frame(Locus=id,OrigChr=chr,OrigBP=bp,
                      NewChr=updated.chr,NewBP=updated.bp,
                      stringsAsFactors = FALSE)
    }else{
      out<-data.frame(Locus=id,OrigChr=chr,OrigBP=bp,
                      NewChr=NA,NewBP=NA,
                      stringsAsFactors = FALSE)
      print(paste("WARNING: locus ", id, " not found",sep=""))
    }
  }
  
  return(out)
}
convert.stacks<-function(stacks.fst,outname,lgs,ssc.agp,sscf.agp,chr.agp,scf.agp){
  stacks.fst$Chr<-as.character(stacks.fst$Chr)
  for(i in 1:nrow(stacks.fst)){
     convert<-convert.agp(old.agp=ssc.agp,old.scf=sscf.agp,
                          new.agp=chr.agp[chr.agp$W=="W",],scf.agp = scf.agp,
                          chr=as.character(stacks.fst$Chr[i]),
                          bp=stacks.fst$BP[i],id=as.character(stacks.fst$Locus.ID[i]))
    stacks.fst[i,"Chr"]<-as.character(convert["NewChr"])
    stacks.fst[i,"BP"]<-convert["NewBP"]
  }
  # reorder by chrom
  scaffs<-levels(as.factor(stacks.fst$Chr))
  scaffs[1:22]<-lgs
  upd.fst<-do.call(rbind,lapply(scaffs,function(lg){
    this.chr<-stacks.fst[stacks.fst$Chr==lg,]
    this.chr<-this.chr[order(as.numeric(this.chr$BP)),]
    return(this.chr)
  }))
  write.table(upd.fst,outname,col.names = TRUE,row.names = FALSE,quote=FALSE,sep='\t')
  print(by(upd.fst,upd.fst$Chr,function(chr){ return(max(chr$BP)/1000000) })[lgs])
  return(upd.fst)
}
```


```r
# old agps
ssc.agp<-read.delim("../../scovelli_genome/SSC_genome.agp",
                    comment.char="#",header=FALSE)
colnames(ssc.agp)<-c("object","object_beg","object_end","part_number","W",
                     "component_id","component_beg","component_end","orientation")
sscf.agp<-read.delim("../../scovelli_genome/SSC_scaffolds.agp",
                     comment.char="#",header=FALSE)
colnames(sscf.agp)<-c("object","object_beg","object_end","part_number","W",
                      "component_id","component_beg","component_end","orientation")
# new scaffold and chrom level agps
scf.agp<-read.delim(gzfile("../../scovelli_genome/ssc_2016_12_20_scafflevel.agp.gz"),
                    comment.char="#",header=FALSE)
chr.agp<-read.delim(gzfile("../../scovelli_genome/ssc_2016_12_20_chromlevel.agp.gz"),
                    comment.char="#",header=FALSE)
colnames(scf.agp)<-c("object","object_beg","object_end","part_number","W",
                     "component_id","component_beg","component_end","orientation")
colnames(chr.agp)<-c("object","object_beg","object_end","part_number","W",
                     "component_id","component_beg","component_end","orientation")
```


```r
vcf<-parse.vcf("stacks/populations_subset75/batch_2.pruned.vcf")
converted<-data.frame(Locus=integer(),OrigChr=character(),OrigBP=integer(),
                      NewChr=character(),NewBP=integer(),
                      stringsAsFactors = FALSE)
for(i in 1:nrow(vcf)){
  converted[i,]<-convert.agp(locus=vcf[i,],old.agp=ssc.agp,old.scf=sscf.agp,
                             new.agp=chr.agp[chr.agp$W=="W",],scf.agp = scf.agp)
}
```


```r
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
```



```r
new.vcf<-as.data.frame(vcf,stringsAsFactor=FALSE)
for(i in 1:nrow(vcf)){
  new.vcf$POS[i]<-converted$NewBP[i]
  new.vcf$`#CHROM`[i]<-as.character(converted$NewChr[i])
}
write.table(new.vcf,"converted_subset.vcf",sep='\t',
            quote=FALSE,col.names = TRUE,row.names = FALSE)
```

```r
snps<-read.delim("stacks/populations_subset75/pruned_snps.txt")

fwsw.fl<-read.delim("stacks/populations_subset75/batch_2.fst_FLCC-FLFW.tsv")
fwsw.tx<-read.delim("stacks/populations_subset75/batch_2.fst_TXCC-TXFW.tsv")
fwsw.al<-read.delim("stacks/populations_subset75/batch_2.fst_ALFW-ALST.tsv")
fwsw.la<-read.delim("stacks/populations_subset75/batch_2.fst_ALST-LAFW.tsv")


fwsw.tx<-fwsw.tx[paste(fwsw.tx$Locus.ID,fwsw.tx$BP) %in% paste(snps$V1,snps$V2),]
fwsw.fl<-fwsw.fl[paste(fwsw.fl$Locus.ID,fwsw.fl$BP) %in% paste(snps$V1,snps$V2),]
fwsw.al<-fwsw.al[paste(fwsw.al$Locus.ID,fwsw.al$BP) %in% paste(snps$V1,snps$V2),]
fwsw.la<-fwsw.la[paste(fwsw.la$Locus.ID,fwsw.la$BP) %in% paste(snps$V1,snps$V2),]


upd.tx<-convert.stacks(fwsw.tx,"stacks/populations_subset75/converted.fst_TXCC-TXFW.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
upd.al<-convert.stacks(fwsw.al,"stacks/populations_subset75/converted.fst_ALFW-ALST.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
upd.la<-convert.stacks(fwsw.al,"stacks/populations_subset75/converted.fst_ALST-LAFW.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
upd.fl<-convert.stacks(fwsw.fl,"stacks/populations_subset75/converted.fst_FLCC-FLFW.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
```


```r
swsw.fl<-read.delim("stacks/populations_whitelist/batch_2.fst_FLCC-FLHB.tsv")
swsw.tx<-read.delim("stacks/populations_whitelist/batch_2.fst_TXCB-TXCC.tsv")
swsw.al<-read.delim("stacks/populations_whitelist/batch_2.fst_ALST-FLSG.tsv")

upd.st<-convert.stacks(swsw.tx,"stacks/converted.fst_TXCB-TXCC.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
upd.sa<-convert.stacks(swsw.al,"stacks/converted.fst_ALST-FLSG.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
upd.sf<-convert.stacks(swsw.fl,"stacks/converted.fst_FLCC-FLHB.txt",
                       lgs,ssc.agp,sscf.agp,chr.agp,scf.agp)
```


## Stacks {-}


```r
fwsw.al<-read.delim("stacks/populations_subset75/converted.fst_ALFW-ALST.txt")
fwsw.la<-read.delim("stacks/populations_subset75/converted.fst_ALST-LAFW.txt")
fwsw.tx<-read.delim("stacks/populations_subset75/converted.fst_TXCC-TXFW.txt")
fwsw.fl<-read.delim("stacks/populations_subset75/converted.fst_FLCC-FLFW.txt")
kable(cbind(nrow(fwsw.al),nrow(fwsw.la),nrow(fwsw.tx),nrow(fwsw.fl)),
      col.names=c("Alabama","Louisiana","Texas","Florida"),
      caption="The number of SNPs in each dataset")
```

\begin{table}

\caption{(\#tab:readStacksFsts)The number of SNPs in each dataset}
\centering
\begin{tabular}[t]{r|r|r|r}
\hline
Alabama & Louisiana & Texas & Florida\\
\hline
6918 & 6918 & 3687 & 2049\\
\hline
\end{tabular}
\end{table}

These datasets do not contain the full 12103 SNPs (\@ref(tab:readStacksFsts)) because some of those SNPs are fixed in the populations or do not pass coverage or minimum allele frequency thresholds.


```r
source("../R/205_popgenPlotting.R")
fst_dat<-list(fwsw.al,fwsw.la,fwsw.tx,fwsw.fl)
fsts<-plot_multiple_LGs(list_fsts = fst_dat,fst_name = "Corrected.AMOVA.Fst",
                        bp_name="BP",chr_name="Chr",lgs=lgs,
                        plot_labs=list("ALFW vs ALST","ALST vs LAFW",
                                       "TXFW vs TXCC","FLFW vs FLCC"),
                        pt_cols = list(c(grp.colors[3],grp.colors[2]),
                                       c(grp.colors[2],grp.colors[3]),
                                       c(grp.colors[1],grp.colors[2]),
                                       c(grp.colors[6],grp.colors[5])),
                        ncol=1,addSmooth = FALSE,pch=19,y.lim = c(0,1),
                        pt.cex=1,axis.size = 1)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/plotStacksFsts-1} \caption{Manhattan plot of pairwise AMOVA-corrected Fst values from Stacks for each freshwater - nearest saltwater population pair. The x-axis corresponds to genomic locations, with chromosomes labelled. To the right are loci that mapped to unanchored scaffolds.}(\#fig:plotStacksFsts)
\end{figure}

So this generated a plot for each pairwise comparison. We could look for shared outliers and see if we can find anything.


```r
tx.sig<-fwsw.tx[fwsw.tx$Fisher.s.P<0.01,"Locus.ID"]
la.sig<-fwsw.la[fwsw.la$Fisher.s.P<0.01,"Locus.ID"]
al.sig<-fwsw.al[fwsw.al$Fisher.s.P<0.01,"Locus.ID"]
fl.sig<-fwsw.fl[fwsw.fl$Fisher.s.P<0.01,"Locus.ID"]
length(tx.sig[(tx.sig %in% c(la.sig,al.sig,fl.sig))])
```

```
## [1] 84
```

```r
length(la.sig[(la.sig %in% c(tx.sig,al.sig,fl.sig))])
```

```
## [1] 189
```

```r
length(al.sig[(al.sig %in% c(la.sig,tx.sig,fl.sig))])
```

```
## [1] 189
```

```r
all.shared<-fl.sig[fl.sig %in% la.sig & fl.sig %in% al.sig & fl.sig %in% tx.sig]

gulf_shared<-tx.sig[tx.sig %in% la.sig & tx.sig %in% al.sig]
```

There are 10 outliers (as determined by Fisher's P from stacks < 0.01). But because of the large pairwise Fsts between the Florida populations, we'll focus on the shared SNPs in the Texas, Alabama, and Louisiana analyeses, in which we have 46

As a point of comparison, we can repeat this with the similar pairwise saltwater-saltwater $F_{ST}$ comparisons.


```r
swsw.fl<-read.delim("stacks/converted.fst_TXCB-TXCC.txt")
swsw.tx<-read.delim("stacks/converted.fst_ALST-FLSG.txt")
swsw.al<-read.delim("stacks/converted.fst_FLCC-FLHB.txt")
fst_dat<-list(swsw.fl,swsw.tx,swsw.al)
fsts<-plot_multiple_LGs(list_fsts = fst_dat,fst_name = "Corrected.AMOVA.Fst",
                        bp_name="BP",chr_name="Chr",lgs=lgs,
                        plot_labs=list("TXSP vs TXCC","ALST vs FLSG","FLHB vs FLCC"),
                        pt_cols = list(c(grp.colors[1],grp.colors[2]),
                                       c(grp.colors[2],grp.colors[3]),
                                       c(grp.colors[6],grp.colors[5])),
                        ncol=1,addSmooth = FALSE,pch=19,y.lim = c(0,1),
                        pt.cex=1,axis.size = 1)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/plotStacksFstsSWSW-1} \caption{Manhattan plot of pairwise AMOVA-corrected Fst values from Stacks for the saltwater populations nearest to freshwater populations compared to their nearest saltwater neightbor. The x-axis corresponds to genomic locations, with chromosomes labelled. To the right are loci that mapped to unanchored scaffolds.}(\#fig:plotStacksFstsSWSW)
\end{figure}

So this generated a plot for each pairwise comparison. We could look for shared outliers and see if we can find anything.


```r
tx.sig<-swsw.tx[swsw.tx$Fisher.s.P<0.01,"Locus.ID"]
al.sig<-swsw.al[swsw.al$Fisher.s.P<0.01,"Locus.ID"]
fl.sig<-swsw.fl[swsw.fl$Fisher.s.P<0.01,"Locus.ID"]
all.shared<-fl.sig[fl.sig %in% al.sig & fl.sig %in% tx.sig]
```

There are 1 outliers (as determined by Fisher's P from stacks < 0.01).

## Permutations {-}


```r
permute.gwsca<-function(vcf,map1,nperms,z=1.96, maf.cutoff = 0.05,cov.thresh=0.2){
  # calculate the actuals
  actual_fsts<-gwsca(vcf,colnames(vcf)[1:9],
                     map1[map1[,2] %in% unique(map1[,2])[1],1],
                     map1[map1[,2] %in% unique(map1[,2])[2],1],
                     maf.cutoff=maf.cutoff,prop.ind.thresh=cov.thresh)
  # do the permutations
  perm_fsts<-lapply(1:nperms,function(i,vcf,map1){
    perm_map<-map1
    perm_map[,2]<-perm_map[,2][permute::shuffle(perm_map[,2])]
    perm_dat<-gwsca(vcf,colnames(vcf)[1:9],
                     perm_map[perm_map[,2] %in% unique(perm_map[,2])[1],1],
                     perm_map[perm_map[,2] %in% unique(perm_map[,2])[2],1],
                     maf.cutoff,cov.thresh)
   
    return(perm_dat)
  },vcf=vcf,map1=map1)
  
  # calculate stats
  fsts<-t(do.call(rbind,lapply(perm_fsts,'[[',"Fst"))) #extract permuted fsts
  perm_fst_mu<-rowMeans(fsts)
  perm_fst_in<-NULL
  for(i in 1:nrow(actual_fsts)){
    pmax<-max(fsts[i,] )
    pmin<-min(fsts[i,] )
    if(actual_fsts[i,"Fst"] > pmax | actual_fsts[i,"Fst"] < pmin ){
      perm_fst_in[i]<-1
    }else{
      perm_fst_in[i]<-0
    }
  }
  
  fst_dat<-data.frame(cbind(actual_fsts,
                            n_perms=nperms,
                            mean_perm=perm_fst_mu,
                            act_in_perm=perm_fst_in))
  return(fst_dat)
}
```

```r
vcf<-parse.vcf("converted_subset.vcf")
popmap<-data.frame(inds=colnames(vcf)[10:ncol(vcf)],
                   pops=gsub("sample_(\\w{4}).*","\\1",colnames(vcf)[10:ncol(vcf)]),
                   stringsAsFactors = FALSE)
pwise_maps<-list(popmap[popmap$pops %in% c("TXFW","TXCC"),],
                 popmap[popmap$pops %in% c("FLLG","FLCC"),],
                 popmap[popmap$pops %in% c("ALFW","ALST"),],
                 popmap[popmap$pops %in% c("LAFW","ALST"),])

permuted_fsts<-lapply(pwise_maps,permute.gwsca,vcf=vcf,nperms=1000, maf.cutoff=0)
saveRDS(permuted_fsts,"permuted_fsts.RDS")
```

Now let's visualize it.


```r
plot_fst_hists<-function(perms,plot_lab=NULL,cols=NULL,permlab="mean_perm",
                         reallab="Fst",baseplot=TRUE,inset=NULL){
  require(scales)
  if(is.null(plot_lab)){
    plot_lab<-""
  }
  if(is.null(cols)){
    cols<-c("grey","black")
  } else if(length(cols)==1){
    cols<-c("dark grey",cols)
  }
  #inset<-par()$fig
  #browser()
  if(isTRUE(baseplot)){
    hist(perms[,permlab],col=alpha(cols[1],0.5),border = alpha(cols[1],0.5),
         xlim=c(0,1),breaks = seq(0,1,0.1),main = plot_lab,
         xlab=expression(italic(F)[ST]),
         ylab="Number of SNPs")
    hist(perms[,reallab],col=alpha(cols[2],0.5),border = alpha(cols[2],0.5),
         xlim=c(0,1),breaks = seq(0,1,0.1),main = "",
         xlab=expression(italic(F)[ST]),
         ylab="Number of SNPs",add=TRUE)
  }
  if(!is.null(inset)){ # add an inset
    # adjust the fig coordinates
    ifig<-c(inset[1]+0.25*(inset[2]-inset[1]),inset[2], 
            inset[3]+0.25*(inset[4]-inset[3]), inset[4])
    par(fig = ifig,new=TRUE) # start x, end x, start y, end y (percent plotting space)
    hist(perms[,permlab][perms[,reallab]>0],col=alpha(cols[1],0.5),
         border = alpha(cols[1],0.5),
         xlim=c(0,1),breaks = seq(0,1,0.1),main = "",xlab="",
         ylab="")
    box() #give it a box
    hist(perms[,reallab,][perms[,reallab]>0],col=alpha(cols[2],0.5),
         border = alpha(cols[2],0.5),
         xlim=c(0,1),breaks = seq(0,1,0.1),main = "",xlab="",
         ylab="",add=TRUE)
  }
  invisible(par()$fig)
}
```

```r
plot_labs<-list("TXFW vs TXCC","FLFW vs FLCC","ALFW vs ALST","ALST vs LAFW")
pt_cols<-list(TXTX=grp.colors[1],FLFL=grp.colors[6],
              ALAL=grp.colors[3],ALLA=grp.colors[2])
png("../figs/permuted_fsts.png",pointsize = 16,height=7,width=8,units="in",res=300)
par(mfrow=c(2,2),new=FALSE,mar=c(4,4,3,1))
#plot the base
pars<-mapply(plot_fst_hists,perms=permuted_fsts,plot_lab=plot_labs,
             cols=pt_cols,SIMPLIFY = FALSE)
# add the insets
ipars<-mapply(plot_fst_hists,perms=permuted_fsts,plot_lab=plot_labs,cols=pt_cols,
              inset=pars,
              MoreArgs = list(baseplot=FALSE))
dev.off()
```

![Histograms of pairwise Fst values generated by the permutations. The inset in each shows the same data but with a smaller range on the y-axis to provide an imporved visualisation.](../figs/permuted_fsts.png)

Now let's start to aggregate everything.


```r
permuted_fsts<-readRDS("permuted_fsts.RDS")

perm_dat<-do.call(cbind,permuted_fsts)
colnames(perm_dat)[colnames(perm_dat) %in% "act_in_perm"]<-c("perm_TX","perm_FL","perm_AL","perm_LA")
perm_dat<-perm_dat[,c("Chrom","Pos","perm_TX","perm_FL","perm_AL","perm_LA")]
perm_dat$Pos<-as.numeric(perm_dat$Pos)

vcf<-parse.vcf("converted_subset.vcf")
fw_SNPinfo<-data.frame(ID=vcf$ID,Chrom=vcf$`#CHROM`,Pos=vcf$POS,BP=vcf$POS-1,
                       REF=vcf$REF,ALT=vcf$ALT,
                       perm_TX=permuted_fsts[[1]]$act_in_perm,
                       perm_FL=permuted_fsts[[2]]$act_in_perm,
                       perm_AL=permuted_fsts[[3]]$act_in_perm,
                       perm_LA=permuted_fsts[[4]]$act_in_perm,
                       stringsAsFactors = FALSE)
```




```r
# add Fsts
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.al,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),
                                               "Corrected.AMOVA.Fst")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_AL"
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.la,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),
                                               "Corrected.AMOVA.Fst")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_LA"
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.tx,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),
                                               "Corrected.AMOVA.Fst")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_TX"
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.fl,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),
                                               "Corrected.AMOVA.Fst")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_FL"
# add p-values
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.al,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),"Fisher.s.P")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_AL_P"
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.la,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),"Fisher.s.P")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_LA_P"
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.tx,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),"Fisher.s.P")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_TX_P"
fw_SNPinfo<-merge(fw_SNPinfo,fwsw.fl,by.x=c("Chrom","BP"),by.y=c("Chr","BP"),
                  all.x=TRUE,all.y = FALSE)[,c(colnames(fw_SNPinfo),"Fisher.s.P")] 
colnames(fw_SNPinfo)[ncol(fw_SNPinfo)]<-"stacks_FL_P"
saveRDS(fw_SNPinfo,"fw_SNPinfo.RDS")
```

## Fst plot


```r
library(UpSetR);library(scales);library(ggplot2)
library(grid);library(gwscaR);library(gridGraphics)
source("../R/upset_hacked.R")
source("../R/205_popgenPlotting.R")
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
```


Then we'll plot the Stacks Fst outliers.


```r
pop.list<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
	"FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLLG")
pop.labs<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLFW")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
lgn<-seq(1,22)
cols<-c(perm=alpha('#e41a1c',0.75),sal=alpha('#377eb8',0.75),pc=alpha('#a65628',0.75),
        stacks=alpha('#f781bf',0.75),xtx=alpha('#ff7f00',0.75))
grp7colors<-c('#762a83','#9970ab','#c2a5cf','#d9f0d3','#a6dba0','#5aae61','#1b7837')
png("../figs/FstOutliers_NoBayenv.png",height=8,width=8.5,units="in",res=300,pointsize=20)
par(mfrow=c(4,1),oma=c(1,1,0.5,1),mar=c(2,2,1,1),xpd=TRUE)
# plot TX
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_TX",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,
                   pch=19,pt.cols = c(grp7colors[1],grp7colors[2]),pt.cex = 2)
#points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#       plot_dat$stacks_TX[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#       col=cols["sal"],cex=1,pch=2)
#points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#       plot_dat$stacks_TX[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#       col=cols["xtx"],cex=1,pch=3,lwd=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       plot_dat$stacks_TX[
         rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       col=cols["perm"],cex=2,pch=5,lwd=3)
#points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
#       plot_dat$stacks_TX[plot_dat$pcadaptQ<0.01],
#       col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 )],
       plot_dat$stacks_TX[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05)],
       col=cols["stacks"],cex=2,pch=8,lwd=3)
axis(2,las=1,pos=-1500000,cex.axis=1)
mtext("Texas",2,line=1.2,cex=1)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15,cex=1)

# plot AL
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_AL",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,
                   pch=19,pt.cols = c(grp7colors[3],"lightgrey"),pt.cex = 2)
# points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#        plot_dat$stacks_AL[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#        col=cols["xtx"],cex=1,pch=3,lwd=2)
# points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#        plot_dat$stacks_AL[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#        col=cols["sal"],cex=1,pch=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       plot_dat$stacks_AL[
         rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       col=cols["perm"],cex=2,pch=5,lwd=3)
# points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
#        plot_dat$stacks_AL[plot_dat$pcadaptQ<0.01],
#        col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05)],
       plot_dat$stacks_AL[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05)],
       col=cols["stacks"],cex=2,pch=8,lwd=3)
axis(2,las=1,pos=-1500000,cex.axis=1)
mtext("Alabama",2,line=1.2,cex=1)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15,cex=1)

# plot LA
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_LA",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19,
                   pt.cols = c("lightgrey",grp7colors[3]),pt.cex = 2)
# points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#        plot_dat$stacks_LA[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#        col=cols["sal"],cex=1,pch=2)
# points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#        plot_dat$stacks_LA[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#        col=cols["xtx"],cex=1,pch=3,lwd=2)
# points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
#        plot_dat$stacks_LA[plot_dat$pcadaptQ<0.01],
#        col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       plot_dat$stacks_LA[
         rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       col=cols["perm"],cex=2,pch=5,lwd=3)

points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05)],
       plot_dat$stacks_LA[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 )],
       col=cols["stacks"],cex=2,pch=8,lwd=3)
axis(2,las=1,pos=-1500000,cex.axis=1)
mtext("Louisiana",2,line=1.2,cex=1)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15,cex=1)


# FL
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_FL",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19,
                   pt.cols = c(grp7colors[6],grp7colors[7]),pt.cex = 2)
# points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#        plot_dat$stacks_FL[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
#        col=cols["sal"],cex=1,pch=2)
# points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#        plot_dat$stacks_FL[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
#        col=cols["xtx"],cex=1,pch=3,lwd=2)
# points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
#        plot_dat$stacks_FL[plot_dat$pcadaptQ<0.01],
#        col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       plot_dat$stacks_FL[
         rowSums(plot_dat[,c("perm_TX","perm_AL","perm_LA")])==3],
       col=cols["perm"],cex=2,pch=5,lwd=3)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05)],
       plot_dat$stacks_FL[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05)],
       col=cols["stacks"],cex=2,pch=8,lwd=3)
axis(2,las=1,pos=-1500000,cex.axis=1)
mtext("Florida",2,line=1.2,cex=1)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15,cex=1)

# add outside legend

opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
            mar=c(0, 0, 0, 0), new=TRUE)
on.exit(par(opar))
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("top",
       c(expression("Permutation"~italic("F")["ST"]),
         expression("Stacks"~italic("F")["ST"])
         #"PCAdapt",expression(italic("X")^T~italic("X")),"Salinity BF"
         ),
       #xjust = 0.5,x.intersp = 0.5,
       col = cols[c("perm","stacks")],
       pt.bg=cols[c("perm","stacks")],
       #pch=c(4,5,0,3,2),
       pch=c(5,8),pt.lwd =2,pt.cex=2,cex=1,
       bty='n',ncol=2)
dev.off()
```


```r
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
outliers<-list( permutations=fw_SNPinfo$ID[
                 rowSums(fw_SNPinfo[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
               
               Alabama=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05)], 
               Louisiana=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_LA_P < 0.05)],
               Texas=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_TX_P < 0.05)],
               Florida=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_FL_P < 0.05)],
               sharedStacks=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05 & 
                                                  fw_SNPinfo$stacks_LA_P < 0.05 &
                 fw_SNPinfo$stacks_TX_P < 0.05)])
cols<-c(permutations='#e41a1c',salBF='#377eb8',pcadapt='#a65628',
        xtx='#ff7f00',sharedStacks='#f781bf',
        Alabama='#af8dc3',Louisiana='#e7d4e8',Texas='#762a83',Florida='#1b7837')

png("../figs/sharedStacks.png",res=300,height=8,width=7.5,units="in",pointsize = 20)
upset(fromList(outliers),sets=c("Texas","Alabama","Louisiana","Florida"),
      point.size=7,line.size=2,mainbar.y.label = "Shared Stacks Outliers",
      sets.x.label = "# Outliers",text.scale=rep(3,6),
      sets.bar.color =cols[c("Texas","Alabama","Louisiana","Florida")],
      margin1scale = 0.2,
      sets.pt.color=cols[c("Texas","Alabama","Louisiana","Florida")],
      keep.order = TRUE)
dev.off()
```

```
## png 
##   2
```


```r
outliers<-list(permutations=fw_SNPinfo$ID[
                 rowSums(fw_SNPinfo[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
               Alabama=fw_SNPinfo$ID[which(fw_SNPinfo$perm_AL==1)], 
               Louisiana=fw_SNPinfo$ID[which(fw_SNPinfo$perm_LA ==1)],
               Texas=fw_SNPinfo$ID[which(fw_SNPinfo$perm_TX==1)],
               Florida=fw_SNPinfo$ID[which(fw_SNPinfo$perm_FL==1)],
               sharedStacks=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05 & 
                                                  fw_SNPinfo$stacks_LA_P < 0.05 &
                                                  fw_SNPinfo$stacks_TX_P < 0.05)])
cols<-c(permutations='#e41a1c',salBF='#377eb8',pcadapt='#a65628',
        xtx='#ff7f00',sharedStacks='#f781bf',
        Alabama='#af8dc3',Louisiana='#e7d4e8',Texas='#762a83',Florida='#1b7837')
png("../figs/sharedPerms.png",res=300,height=8,width=8.5,units="in",pointsize = 20)
upset(fromList(outliers),sets=c("Texas","Alabama","Louisiana","Florida","sharedStacks"),
      point.size=7,line.size=2,mainbar.y.label = "Shared Permutation Outliers",
      sets.x.label = "# Outliers",text.scale=rep(3,6),
      sets.bar.color =cols[c("Texas","Alabama","Louisiana","Florida","sharedStacks")],
      margin1scale = 0.2,
      sets.pt.color=cols[c("Texas","Alabama","Louisiana","Florida","sharedStacks")],
      keep.order = TRUE)
dev.off()
```

```
## png 
##   2
```



```r
library(magick);library(multipanelfigure)

image_files <- c("../figs/FstOutliers_NoBayenv.png",
                 "../figs/sharedStacks.png",
                 "../figs/sharedPerms.png")


png("../figs/fstPlots.png",height=4,width=6.5,units="in",res=300,pointsize = 20)
figure <- multi_panel_figure(
  width = c(4, 2),
  height = c(1.8,1.8),
  unit = "inches",
  row_spacing = 0.0,
  column_spacing = 0,
    
)
(figure %<>% fill_panel(image_files[1],row = 1:2,scaling="fit",
                        allow_panel_overwriting = TRUE) )
#(figure %<>% fill_panel(image_files[2],column=2,row=1:2,scaling="fit"))
(figure %<>% fill_panel(image_files[2], column=2, row=1,scaling="fit",
                        allow_panel_overwriting = TRUE) )
(figure %<>% fill_panel(image_files[3],column=2,row=2,scaling="fit",
                        allow_panel_overwriting = TRUE) )
dev.off()
```

## PCAdapt {-}

For this analysis, to maintain consistency with the other outlier analyses, I'm using the subestted dataset. So I need to run PCAdapt [@luu_pcadapt:_2017] another time.


```r
vcf<-parse.vcf("converted_subset.vcf")
write.table("##fileformat=VCFv","pcadapt_fw/fwsw.pruned.vcf",quote=FALSE,
            col.names = FALSE,row.names = FALSE)
suppressWarnings(write.table(vcf,"pcadapt_fw/fwsw.pruned.vcf",
                             quote=FALSE,append = TRUE,
                             row.names = FALSE,col.names = TRUE,sep='\t'))
```




```r
library(pcadapt)
#need to remove the first line with a # 
filename<-read.pcadapt("pcadapt_fw/fwsw.pruned.vcf",type="vcf") 
```

```
## No variant got discarded.
## Summary:
## 
## 	- input file:				pcadapt_fw/fwsw.pruned.vcf
## 	- output file:				/tmp/RtmpWYtABi/file21ac69a1c01c.pcadapt
## 
## 	- number of individuals detected:	303
## 	- number of loci detected:		12103
## 
## 12103 lines detected.
## 303 columns detected.
```

```r
x<-pcadapt(filename, K=20)
plot(x,option="screeplot")
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/pcadaptOutliers-1} \caption{Screeplot for the subsetted dataset with 20 PC axes retained in the analysis.}(\#fig:pcadaptOutliers)
\end{figure}

$K=4$ seems like the best choice here to keep values to the left of the straight line (or could be $K=6$). 


```r
# Organize pop info
pops<-gsub("sample_(\\w{4}).*","\\1",colnames(vcf)[10:ncol(vcf)])	
grp<-pops
grp[grp=="TXFW" | grp=="LAFW" | grp=="ALFW" | grp=="FLLG"]<-"freshwater"
grp[grp!="freshwater"]<-"saltwater"
```



```r
res<-pcadapt(filename,K=4)
plot(res, option="manhattan")
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/PcadaptK4-1-1} \caption{Manhattan plot for the PCAdapt outlier analysis with the subsetted dataset and 4 PC axes retained in the analysis.}(\#fig:PcadaptK4-1)
\end{figure}

```r
plot(res, option="qqplot")
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/PcadaptK4-2-1} \caption{Q-Q plot for the PCAdapt outlier analysis with the subsetted dataset and 4 PC axes retained in the analysis.}(\#fig:PcadaptK4-2)
\end{figure}

```r
plot(res, option="stat.distribution")
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/PcadaptK4-3-1} \caption{Distribution of the  for the PCAdapt outlier analysis with the subsetted dataset and 4 PC axes retained in the analysis.}(\#fig:PcadaptK4-3)
\end{figure}

```r
plot(res, option="scores",pop=pops)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/PcadaptK4-4-1} \caption{Plot of the scores from the PCAdapt outlier analysis with the subsetted dataset and 4 PC axes retained in the analysis.}(\#fig:PcadaptK4-4)
\end{figure}

The PCAdapt vignette recommends displaying the loadings and evaluate if loadings are clustered in single or several genomic regions


```r
par(mfrow = c(2, 2))
for (i in 1:4)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/PcadaptLoadings-1} \caption{Plots of the loadings for the four PC axes according to genomic position (on the x-axis).}(\#fig:PcadaptLoadings)
\end{figure}

This suggests that loadings are not clustered (assuming these are grouped by space), so we don't need to worry about LD thinning. Now let's look chromosome by chromosome:


```r
par(mfrow=c(6,4),mar=c(3,3,2,1.5))
l<-lapply(lgs, function(lg,vcf){
  plot(res$loadings[which(vcf$`#CHROM` %in% lg), 1], pch = 19, cex = .3, 
       xlab = paste0("Position on ", lg), ylab = "Loadings PC 1")
  mtext(lg,3,outer=FALSE)
},vcf=vcf)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/PcadaptLoadingsChr-1} \caption{Inspecting the loadings for the four PC axes according to genomic position (on the x-axis) for each chromosome individually.}(\#fig:PcadaptLoadingsChr)
\end{figure}

None of the LGs seem to have huge clusters of outliers so we can move on, lumping them all together.

We need to choose a cutoff for outlier detection. I'll use the qvalue approach, which identifies outliers with a false discovery rate of $\alpha$, which I'm setting here to 0.05.


```r
library(qvalue)
qval <- qvalue(res$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
snp_pc<-get.pc(res,outliers) # Get the PCs associated with outliers
```

We identified 181 outliers with this analysis, which are associated with 4 of the 4 clusters. If we look at the distribution of these, though, we see that most are associated with PC 1


```r
table(snp_pc$PC)
```

```
## 
##   1   2   3   4 
## 115  36  10  20
```

Now we can add the qvalues to the fw_SNPinfo dataframe


```r
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
fw_SNPinfo$pcadaptQ<-qval
fw_SNPinfo$pcadaptPC<-get.pc(res,1:length(qval))$PC
saveRDS(fw_SNPinfo,"fw_SNPinfo.RDS")
```

I should note that for some of these PCAdapt gives "NA" -- not sure what causes this behaviour but there it is. It's true for the residuals and everything else. It looks to be due to low allele frequencies -- though stacks should have been run with a minimum allele frequency cutoff, so this is perplexing. 

# Bayenv {-}

Investigate the environmental data


```r
env.data<-read.csv("bayenv/env_data_raw.csv",row.names = 1)
env.data<-rbind(env.data,pop=c(rep("SW",12),rep("FW",4)))
env.data<-as.data.frame(t(env.data))
wilcox.test(as.numeric(env.data$temp)~env.data$pop) #ties, but p=0.539
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  as.numeric(env.data$temp) by env.data$pop
## W = 18.5, p-value = 0.5392
## alternative hypothesis: true location shift is not equal to 0
```

```r
wilcox.test(as.numeric(env.data$seagrass)~env.data$pop) #ties, but p=0.897
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  as.numeric(env.data$seagrass) by env.data$pop
## W = 25.5, p-value = 0.8968
## alternative hypothesis: true location shift is not equal to 0
```


This analysis is of just the freshwater and saltwater populations. First I ran bayenv using the script `run_bayenv2_matrix_general.sh`.


```bash
# set up correct file formats
../../scripts/run_bayenv2_matrix_general.sh FILEMANIP \
bayenv/sub75.pruned.clust stacks/populations_subset75/batch_2.pruned.ped \
stacks/populations_subset75/batch_2.pruned.map bayenv ~/Programs/bayenv 7
# estimate matrices
../../scripts/run_bayenv2_matrix_general.sh MATRIX bayenv ~/Programs/bayenv/ 7
```

At this point I looked at the matrices

![Heatmap plots of each of the ten replicate Bayenv matrices, run without the Florida populations. Colors represent (co)variances, with larger values shown in darker colors](../fwsw_results/bayenv/fwsw75_pruned.png)

They all looked fine so I chose the randomly selected a matrix to use for the remainder of the analyses.


```bash
# Create the SNPFILES
../../scripts/run_bayenv2_matrix_general.sh SNPFILES SNPSFILE SNPFILES
# Run bayenv
nohup ../../scripts/run_bayenv2_matrix_general.sh BAYENV \
~/Programs/bayenv/ matrix env_data_sub75.txt 7 3 SNPFILES > bayenv.log &
```

### Analyze Bayenv output

Once Bayenv was finished running, I first aggregated all of the output. 


```r
get_bayenv_results<-function(dir,env_vars){
  # process the variable names
  var_names<-unlist(lapply(env_vars,function(var){
    nms<-c(paste0(var,"_BF"),paste0(var,"_rho"),paste0(var,"_rs"))
    return(nms)
  }))
  # list all the files
  bf.files<-list.files(pattern="bf",path = dir,full.names = TRUE)
  xtx.files<-list.files(pattern="xtx",path=dir,full.names = TRUE)
  # get the bf files
  bf.dat<-do.call(rbind,lapply(bf.files,function(filename){
    bf<-read.table(filename,header = FALSE)
  }))
  colnames(bf.dat)<-c("locus", var_names)
  
  # xtx files
  xtx.files<-list.files(pattern="xtx",path=dir,full.names = TRUE)
  xtx.dat<-do.call(rbind,lapply(xtx.files,function(filename){
    xtx<-read.table(filename,header = FALSE,stringsAsFactors = FALSE)
  }))
  colnames(xtx.dat)<-c("locus","XtX")
  
  # combine the two
  bayenv.dat<-merge(xtx.dat,bf.dat,by="locus")
  return(bayenv.dat)
}
```

The SNP names are uninformative, just the row number the SNP was in. We can make these better using the freq info


```r
bayenv_dat<-get_bayenv_results(dir="bayenv/SNPFILES",
                               env_vars=c("temp","salinity","seagrass"))

freq<-read.table("bayenv/bayenv.frq.strat",header=T, stringsAsFactors=F)
#want to get $MAC for every snp at every pop 
#and NCHROBS-MAC for every stnp at every pop
freq<-cbind(freq,freq$NCHROBS-freq$MAC)
colnames(freq)[ncol(freq)]<-"NAC"
snp.names<-split(freq$SNP,freq$CLST)[[1]]
snp.names<-gsub("(\\d+)_\\d+","\\1",snp.names)

snp_dat<-data.frame(ID=snp.names,loc=seq(1,length(snp.names)*2,2))
bayenv_dat$locus<-as.numeric(gsub("SNPFILES\\/(\\d+)","\\1",bayenv_dat$locus))
bayenv_dat<-bayenv_dat[order(bayenv_dat$locus),]
bayenv_dat<-merge(snp_dat,bayenv_dat,by.x="loc",by.y="locus")
colnames(bayenv_dat)[1:2]<-c("index","SNPID")

pmap<-read.delim("stacks/populations_subset75/batch_2.pruned.map",header = FALSE) 
# the position location is wrong in this ^ file
pmap$locus<-gsub("(\\d+)_\\d+","\\1",pmap[,2])

bayenv_dat<-merge(pmap[,c(1,4,5)],bayenv_dat,by.x="locus",by.y="SNPID")
colnames(bayenv_dat)[2:3]<-c("Chrom","SNPID")
write.table(bayenv_dat,"bayenv/bayenv_output.txt",sep="\t",
            col.names = TRUE,row.names = FALSE,quote = FALSE)
```

Merge it with other snp info

```r
bayenv_dat<-read.delim("bayenv/bayenv_output.txt",header = TRUE)
bayenv_dat$logSalBF<-log(bayenv_dat$salinity_BF)
bayenv_dat$logTemBF<-log(bayenv_dat$temp_BF)
bayenv_dat$logSegBF<-log(bayenv_dat$seagrass_BF)
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
fw_SNPinfo<-merge(fw_SNPinfo,
                  bayenv_dat[,c("SNPID","XtX","logSalBF","logTemBF","logSegBF")],
                  by.x="ID",by.y="SNPID",all.x=TRUE)
saveRDS(fw_SNPinfo,"fw_SNPinfo.RDS")
```

Now we can investigate the outliers etc.


```r
#taken directly from fwsw_analysis.R
bayenv.dat<-read.delim("bayenv/bayenv_output.txt",header=T)
# calculate quantiles for bayes factors
#focus on Bayes Factors, because of Lotterhos & Whitlock (2015)
bf.co<-apply(bayenv.dat[,grep("BF",colnames(bayenv.dat))],2,quantile,0.99) 
temp.bf.sig<-bayenv.dat[bayenv.dat$temp_BF>bf.co["temp_BF"],]
sal.bf.sig<-bayenv.dat[bayenv.dat$salinity_BF>bf.co["salinity_BF"],]
grass.bf.sig<-bayenv.dat[bayenv.dat$seagrass_BF>bf.co["seagrass_BF"],]
#get the log transformed Bayes Factors
bayenv.dat$logSal<-log(bayenv.dat$salinity_BF)
bayenv.dat$logTemp<-log(bayenv.dat$temp_BF)
bayenv.dat$logSeagrass<-log(bayenv.dat$seagrass_BF)

# xtx
xtx.sig<-bayenv.dat[bayenv.dat$XtX > quantile(bayenv.dat$XtX,0.99),]
```

There are 0 overlapping outliers between temperature-, salinity-, and seagrass-associated loci.

But if we only care about salinity ones, there are 122 outliers. Are any of those XtX outliers too? 0 overlap -- not bad!




```r
library(scales)
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
cols<-c(perm=alpha('#e41a1c',0.75),sal=alpha('#377eb8',0.75),pc=alpha('#a65628',0.75),
        stacks=alpha('#f781bf',0.75),xtx=alpha('#ff7f00',0.75))
#png("../figs/BayenvOutliers.png",height=5,width=8.5,units="in",res=300,pointsize=12)
par(mfrow=c(2,1),oma=c(1,2,1,1),mar=c(2,2,1,1),xpd=TRUE)
# plot XtX
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "XtX",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19)
points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       plot_dat$XtX[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       col=cols["xtx"],cex=0.75,pch=3)
points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       plot_dat$XtX[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       col=cols["sal"],cex=0.85,pch=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       plot_dat$XtX[
         rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       col=cols["perm"],cex=1,pch=4,lwd=2)
points(plot_dat$plot.pos[plot_dat$stacks_AL < 0.05 & plot_dat$stacks_LA < 0.05 &
                           plot_dat$stacks_TX < 0.05 & plot_dat$stacks_FL < 0.05],
       plot_dat$XtX[plot_dat$stacks_AL < 0.05 & plot_dat$stacks_LA < 0.05 &
                           plot_dat$stacks_TX < 0.05 & plot_dat$stacks_FL < 0.05],
       col=cols["stacks"],cex=1,pch=5,lwd=2)
points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
       plot_dat$XtX[plot_dat$pcadaptQ<0.01],
       col=cols["pc"],cex=1,pch=0,lwd=2)
axis(2,las=1)
mtext(expression(italic("X")^"T"~italic("X")),2,line=2)

# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=0)

# plot Bayes Factors
plot_dat<-fst.plot(plot_dat,scaffs.to.plot = lgs,fst.name = "logSalBF",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19)
points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       plot_dat$logSalBF[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       col=cols["sal"],cex=0.85,pch=2)
points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       plot_dat$logSalBF[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       col=cols["xtx"],cex=0.75,pch=3)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       plot_dat$logSalBF[
         rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       col=cols["perm"],cex=1,pch=4,lwd=2)
points(plot_dat$plot.pos[plot_dat$stacks_AL < 0.05 & plot_dat$stacks_LA < 0.05 &
                           plot_dat$stacks_TX < 0.05 & plot_dat$stacks_FL < 0.05],
       plot_dat$logSalBF[plot_dat$stacks_AL < 0.05 & plot_dat$stacks_LA < 0.05 &
                           plot_dat$stacks_TX < 0.05 & plot_dat$stacks_FL < 0.05],
       col=cols["stacks"],cex=1,pch=5,lwd=2)
points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
       plot_dat$logSalBF[plot_dat$pcadaptQ<0.01],
       col=cols["pc"],cex=1,pch=0,lwd=2)
axis(2,las=1)
mtext("log(Salinity Bayes Factors)",2,line=2)

# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-5)

# add outside legend

opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
            mar=c(0, 0, 0, 0), new=TRUE)
on.exit(par(opar))
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("top",c(expression("Permutation"~italic("F")["ST"]),
         expression("Stacks"~italic("F")["ST"]),
         "PCAdapt",expression(italic("X")^T~italic("X")),"Salinity BF"),
       xjust = 0.5,x.intersp = 0.5,
       col = cols[c("perm","stacks","pc","xtx","sal")],
       pt.bg=cols[c("perm","stacks","pc","xtx","sal")],
       pch=c(4,5,0,3,2),bty='n',ncol=5)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/plotBayenv-1} \caption{Manhattan plots of the Bayenv XTX statistic and the log of the Bayes factors associated with salinity. Shown are only loci that mapped to regions on chromosomes. Grey points represent loci that are not outlier in any analyses. Colored points represent outliers in both the Bayenv analyses and the other outlier analyses.}(\#fig:plotBayenv)
\end{figure}

```r
#dev.off()
```



### Bayenv without Florida

First I need to subset the ped file so that it doesn't contain the Florida samples.


```r
ped<-read.delim("stacks/populations_subset75/batch_2.pruned.ped",header=FALSE,sep=' ')
keep_ped<-ped[grep("FL",ped$V2,invert=TRUE),]
keep_ped$V1<-as.numeric(as.factor(gsub("sample_(\\w{4}).*","\\1",keep_ped$V2)))
write.table(keep_ped,"stacks/populations_subset75/noFL_subset75.ped",col.names=FALSE,row.names=FALSE,quote=FALSE)
```

And create a cluster file (i.e., pop map)


```r
clust<-keep_ped[,1:2]
clust$clust<-gsub("sample_(\\w{4}).*","\\1",clust$V2)
write.table(clust,"bayenv/noFL/sub75.noFL.clust",col.names = FALSE,row.names = FALSE,quote=FALSE)
```

And I want to use the correctly specified locations


Then I can use the `run_bayenv2_matrix_general.sh` script to run the matrices.


```bash
# set up correct file formats creating a clust file - run from bayenv/ dir
../../scripts/run_bayenv2_matrix_general.sh FILEMANIP \
noFL/sub75.noFL.clust ../stacks/populations_subset75/noFL_subset75.ped \
../stacks/populations_subset75/batch_2.pruned.map noFL/ ~/Programs/bayenv 5
# estimate matrices
../../scripts/run_bayenv2_matrix_general.sh MATRIX noFL ~/Programs/bayenv/ 5
```

They're all really similar-looking so I used a randomly-selected one (it was #9). 

![Heatmap plots of each of the ten replicate Bayenv matrices. Colors represent (co)variances, with larger values shown in darker colors](../fwsw_results/bayenv/noFL/fwsw75_noFL.png)



```bash
# Create the SNPFILES
../../scripts/run_bayenv2_matrix_general.sh SNPFILES noFL/SNPSFILE noFL/SNPFILES
```

Before running the analysis, I need to subset and standardize the environmental data too.


```r
env_dat<-read.csv("bayenv/env_data_sub75.csv",row.names=1)

# make sure it's in the same order as the plinkfile
freq<-read.table("bayenv/noFL/bayenv.frq.strat",header=T, stringsAsFactors=F)
freq<-cbind(freq,freq$NCHROBS-freq$MAC)
colnames(freq)[ncol(freq)]<-"NAC"
pop.order<-levels(as.factor(freq$CLST))
snp.names<-split(freq$SNP,freq$CLST)[[1]]

env_dat<-env_dat[,pop.order]
std_dat<-t(apply(env_dat,1,function(x){
  stands<-(x-mean(x))/sd(x)
  return(stands)
}))

write.table(std_dat[,pop.order],
            "bayenv/noFL/env_data_noFL.txt",sep='\t',col.names = FALSE,row.names = FALSE)
```



```bash
# Run bayenv
nohup ../../scripts/run_bayenv2_matrix_general.sh BAYENV \
~/Programs/bayenv/ noFL/matrix noFL/env_data_noFL.txt 5 3 noFL/SNPFILES > ../../logs/bayenv_noFL.log 2>&1 &
```


```r
bayenv_noFL<-get_bayenv_results(dir="bayenv/noFL/SNPFILES/",
                               env_vars=c("temp","salinity","seagrass"))
# match up SNP names and IDS
snp_names<-read.delim("bayenv/noFL/SNPFILES_names.txt",header=FALSE)
snpsfile<-read.delim("bayenv/noFL/SNPSFILE",header = FALSE)

snp_names$BayenvID<-as.numeric(rownames(snpsfile)[seq(1,(nrow(snpsfile)-1),2)])
colnames(snp_names)[1]<-"SNP_name"
snp_names$ID<-as.numeric(gsub("(\\d+)_\\d+","\\1",snp_names$SNP_name))
snp_names$Column<-as.numeric(gsub("(\\d+)_(\\d+)","\\2",snp_names$SNP_name))

bayenv_noFL$BayenvID<-as.numeric(gsub(".*/(\\d+)","\\1",bayenv_noFL$locus))
bayenv_noFL<-merge(bayenv_noFL,snp_names)

# get positional info from the master information file
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
bayenv_noFL<-merge(bayenv_noFL,fw_SNPinfo[,c(1:4,21:24)],by="ID")
colnames(bayenv_noFL)[4:13]<-paste0(colnames(bayenv_noFL[,4:13]),"_noFL")
colnames(bayenv_noFL)[19:22]<-paste0(colnames(bayenv_noFL[,19:22]),"_FL")
colnames(bayenv_noFL)<-gsub("\\.\\w","",colnames(bayenv_noFL))

# calculate logs of Bayes factors
bayenv_noFL$logSalBF_noFL<-log(bayenv_noFL$salinity_BF_noFL)
bayenv_noFL$logTemBF_noFL<-log(bayenv_noFL$temp_BF_noFL)
bayenv_noFL$logSegBF_noFL<-log(bayenv_noFL$seagrass_BF_noFL)

# save to file
write.table(bayenv_noFL,"bayenv/bayenv_output_noFL.txt",sep="\t",
            col.names = TRUE,row.names = FALSE,quote = FALSE)
```



```r
library(scales)
bayenv_noFL<-read.delim("bayenv/bayenv_output_noFL.txt",header=TRUE)

cols<-c(perm=alpha('#e41a1c',0.75),sal=alpha('#377eb8',0.75),pc=alpha('#a65628',0.75),
        stacks=alpha('#f781bf',0.75),xtx=alpha('#ff7f00',0.75))
#png("../figs/BayenvOutliers.png",height=5,width=8.5,units="in",res=300,pointsize=12)
par(mfrow=c(2,1),oma=c(1,2,1,1),mar=c(2,2,1,1),xpd=TRUE)
# plot XtX
plot_dat<-fst.plot(bayenv_noFL,scaffs.to.plot = lgs,fst.name = "XtX_noFL",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19)
points(plot_dat$plot.pos[plot_dat$XtX_noFL>=quantile(plot_dat$XtX_noFL,0.99)],
       plot_dat$XtX_noFL[plot_dat$XtX_noFL>=quantile(plot_dat$XtX_noFL,0.99)],
       col=cols["xtx"],cex=0.75,pch=3)
points(plot_dat$plot.pos[plot_dat$logSalBF_noFL>=quantile(plot_dat$logSalBF_noFL,0.99)],
       plot_dat$XtX_noFL[plot_dat$logSalBF_noFL>=quantile(plot_dat$logSalBF_noFL,0.99)],
       col=cols["sal"],cex=0.85,pch=2)
axis(2,las=1)
mtext(expression(italic("X")^"T"~italic("X")),2,line=2)

# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=0)

# plot Bayes Factors
plot_dat<-fst.plot(plot_dat,scaffs.to.plot = lgs,fst.name = "logSalBF_noFL",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19)
points(plot_dat$plot.pos[plot_dat$logSalBF_noFL>=quantile(plot_dat$logSalBF_noFL,0.99)],
       plot_dat$logSalBF_noFL[plot_dat$logSalBF_noFL>=quantile(plot_dat$logSalBF_noFL,0.99)],
       col=cols["sal"],cex=0.85,pch=2)
points(plot_dat$plot.pos[plot_dat$XtX_noFL>=quantile(plot_dat$XtX_noFL,0.99)],
       plot_dat$logSalBF_noFL[plot_dat$XtX_noFL>=quantile(plot_dat$XtX_noFL,0.99)],
       col=cols["xtx"],cex=0.75,pch=3)
axis(2,las=1)
mtext("log(Salinity Bayes Factors)",2,line=2)

# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-5)

# add outside legend

opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
            mar=c(0, 0, 0, 0), new=TRUE)
on.exit(par(opar))
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("top",c(expression(italic("X")^T~italic("X")),"Salinity BF"),
       xjust = 0.5,x.intersp = 0.5,
       col = cols[c("xtx","sal")],
       pt.bg=cols[c("xtx","sal")],
       pch=c(3,2),bty='n',ncol=5)
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/plotBayenvNoFL-1} \caption{Manhattan plots of the Bayenv XTX statistic and the log of the Bayes factors associated with salinity from the analysis without the Florida populations. Shown are only loci that mapped to regions on chromosomes. Grey points represent loci that are not outlier in any analyses. Colored points represent outliers in both the Bayenv analyses and the other outlier analyses.}(\#fig:plotBayenvNoFL)
\end{figure}

```r
#dev.off()
```


Let's compare the results of the analyses with and without the Florida populations.


```r
bayenv_cor<-cor(bayenv_noFL[,c("XtX_noFL",
                   "logSalBF_noFL",
                   "logTemBF_noFL",
                   "logSegBF_noFL",
                   "XtX_FL","logSalBF_FL","logTemBF_FL","logSegBF_FL")])
kable(bayenv_cor[grep("_FL",rownames(bayenv_cor)),grep("_noFL",colnames(bayenv_cor))],
      "latex",booktabs=TRUE,
      caption="Correlations between the Bayenv analyses without the Florida populations ('noFL') and with the Florida popualtions ('FL'). The comparison of each statistic with itself in each analysis is on the diagonal.")
```

\begin{table}

\caption{(\#tab:BayenvCompareFL)Correlations between the Bayenv analyses without the Florida populations ('noFL') and with the Florida popualtions ('FL'). The comparison of each statistic with itself in each analysis is on the diagonal.}
\centering
\begin{tabular}[t]{lrrrr}
\toprule
  & XtX\_noFL & logSalBF\_noFL & logTemBF\_noFL & logSegBF\_noFL\\
\midrule
XtX\_FL & 0.5482418 & 0.3697998 & 0.4264008 & 0.4234426\\
logSalBF\_FL & 0.4926238 & 0.5913697 & 0.3864321 & 0.4268952\\
logTemBF\_FL & 0.3977142 & 0.4089924 & 0.4567702 & 0.3327506\\
logSegBF\_FL & 0.3772125 & 0.3057937 & 0.2178636 & 0.3283923\\
\bottomrule
\end{tabular}
\end{table}

```r
par(mfrow=c(2,2))
plot(bayenv_noFL$XtX_FL, bayenv_noFL$XtX_noFL,
     xlab="XtX with Florida populations",ylab="XtX without Florida Populations",
     pch=19,cex=1.5,col=alpha("dark grey",0.75))
plot(bayenv_noFL$logSalBF_FL, bayenv_noFL$logSalBF_noFL,
     xlab="log of Salinity Bayes Factors with Florida populations",
     ylab="log of Salinity Bayes Factors without Florida Populations",
     pch=19,cex=1.5,col=alpha("dark grey",0.75))
plot(bayenv_noFL$logTemBF_FL, bayenv_noFL$logTemBF_noFL,
     xlab="log of Temperature Bayes Factors with Florida populations",
     ylab="log of Temperature Bayes Factors without Florida Populations",
     pch=19,cex=1.5,col=alpha("dark grey",0.75))
plot(bayenv_noFL$logSegBF_FL, bayenv_noFL$logSegBF_noFL,
     xlab="log of Seagrass cover Bayes Factors with Florida populations",
     ylab="log of Seagrass cover Bayes Factors without Florida Populations",
     pch=19,cex=1.5,col=alpha("dark grey",0.75))
```

\begin{figure}[H]
\includegraphics{202_fwsw_reanalysis_files/figure-latex/BayenvCompareFLPlots-1} \caption{Plots of the Bayenv statistics in the analysis with the Florida populations (x-axis) vs without the Florida populations (y-axis).}(\#fig:BayenvCompareFLPlots)
\end{figure}

Are any of the outliers the same in the two analyses?


```r
bothXtX<-bayenv_noFL$ID[bayenv_noFL$XtX_FL>=quantile(plot_dat$XtX_FL,0.99) &
                 bayenv_noFL$XtX_noFL>=quantile(plot_dat$XtX_noFL,0.99)]
bothSal<-bayenv_noFL$ID[bayenv_noFL$logSalBF_FL>=quantile(plot_dat$logSalBF_FL,0.99) &
                 bayenv_noFL$logSalBF_noFL>=quantile(plot_dat$logSalBF_noFL,0.99)]
```

There are 33 SNPs that are XtX outliers in both Bayenv analyses and 58 salinity-associated outliers in both.


## Bayenv Plots

First we'll get set up with the libraries, code, and data.


```r
library(UpSetR);library(scales);library(ggplot2)
library(grid);library(gwscaR);library(gridGraphics)
source("../R/upset_hacked.R")
source("../R/205_popgenPlotting.R")
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
outliers<-list(xtx=fw_SNPinfo$ID[fw_SNPinfo$XtX >= quantile(fw_SNPinfo$XtX,0.99,na.rm=TRUE)],
               salBF=fw_SNPinfo$ID[fw_SNPinfo$logSalBF>=
                                     quantile(fw_SNPinfo$logSalBF,0.99,na.rm=TRUE)],
               permutations=fw_SNPinfo$ID[
                 rowSums(fw_SNPinfo[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
               pcadapt=fw_SNPinfo$ID[which(fw_SNPinfo$pcadaptQ<0.01)],
               Alabama=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05)], 
               Louisiana=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_LA_P < 0.05)],
               Texas=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_TX_P < 0.05)],
               Florida=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_FL_P < 0.05)],
               sharedStacks=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05 & 
                                                  fw_SNPinfo$stacks_LA_P < 0.05 &
                 fw_SNPinfo$stacks_TX_P < 0.05)])
```

Then we'll plot the Stacks Fst outliers.


```r
pop.list<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
	"FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLLG")
pop.labs<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLFW")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
lgn<-seq(1,22)
cols<-c(perm=alpha('#e41a1c',0.75),sal=alpha('#377eb8',0.75),pc=alpha('#a65628',0.75),
        stacks=alpha('#f781bf',0.75),xtx=alpha('#ff7f00',0.75))
grp7colors<-c('#762a83','#9970ab','#c2a5cf','#d9f0d3','#a6dba0','#5aae61','#1b7837')
png("../figs/FstOutliers.png",height=8,width=8.5,units="in",res=300,pointsize=20)
par(mfrow=c(4,1),oma=c(1,1,0.5,1),mar=c(2,2,1,1),xpd=TRUE)
# plot TX
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_TX",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,
                   pch=19,pt.cols = c(grp7colors[1],grp7colors[2]),pt.cex = 1)
points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       plot_dat$stacks_TX[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       col=cols["sal"],cex=1,pch=2)
points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       plot_dat$stacks_TX[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       col=cols["xtx"],cex=1,pch=3,lwd=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       plot_dat$stacks_TX[
         rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       col=cols["perm"],cex=1,pch=4,lwd=2)
points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
       plot_dat$stacks_TX[plot_dat$pcadaptQ<0.01],
       col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       plot_dat$stacks_TX[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       col=cols["stacks"],cex=1,pch=5,lwd=2)
axis(2,las=1,pos=-1500000)
mtext("TXFW vs. TXCC",2,line=1,cex=0.65)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15)

# plot AL
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_AL",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,
                   pch=19,pt.cols = c(grp7colors[3],"lightgrey"),pt.cex = 1)
points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       plot_dat$stacks_AL[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       col=cols["xtx"],cex=1,pch=3,lwd=2)
points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       plot_dat$stacks_AL[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       col=cols["sal"],cex=1,pch=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       plot_dat$stacks_AL[
         rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       col=cols["perm"],cex=1,pch=4,lwd=2)
points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
       plot_dat$stacks_AL[plot_dat$pcadaptQ<0.01],
       col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       plot_dat$stacks_AL[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       col=cols["stacks"],cex=1,pch=5,lwd=2)
axis(2,las=1,pos=-1500000)
mtext("ALFW vs. ALST",2,line=1,cex=0.65)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15)

# plot LA
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_LA",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19,
                   pt.cols = c("lightgrey",grp7colors[3]),pt.cex = 1)
points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       plot_dat$stacks_LA[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       col=cols["sal"],cex=1,pch=2)
points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       plot_dat$stacks_LA[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       col=cols["xtx"],cex=1,pch=3,lwd=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       plot_dat$stacks_LA[
         rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       col=cols["perm"],cex=1,pch=4,lwd=2)
points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
       plot_dat$stacks_LA[plot_dat$pcadaptQ<0.01],
       col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       plot_dat$stacks_LA[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       col=cols["stacks"],cex=1,pch=5,lwd=2)
axis(2,las=1,pos=-1500000)
mtext("LAFW vs. ALST",2,line=1,cex=0.65)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15)


# FL
plot_dat<-fst.plot(fw_SNPinfo,scaffs.to.plot = lgs,fst.name = "stacks_FL",
                   chrom.name = "Chrom",bp.name = "Pos",axis.size = 0,pch=19,
                   pt.cols = c(grp7colors[6],grp7colors[7]),pt.cex = 1)
points(plot_dat$plot.pos[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       plot_dat$stacks_FL[plot_dat$logSalBF>=quantile(plot_dat$logSalBF,0.99)],
       col=cols["sal"],cex=1,pch=2)
points(plot_dat$plot.pos[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       plot_dat$stacks_FL[plot_dat$XtX>=quantile(plot_dat$XtX,0.99)],
       col=cols["xtx"],cex=1,pch=3,lwd=2)
points(plot_dat$plot.pos[
  rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       plot_dat$stacks_FL[
         rowSums(plot_dat[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
       col=cols["perm"],cex=1,pch=4,lwd=2)
points(plot_dat$plot.pos[plot_dat$pcadaptQ<0.01],
       plot_dat$stacks_FL[plot_dat$pcadaptQ<0.01],
       col=cols["pc"],cex=1,pch=0,lwd=2)
points(plot_dat$plot.pos[which(plot_dat$stacks_AL_P < 0.05 & 
                                 plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       plot_dat$stacks_FL[which(plot_dat$stacks_AL_P < 0.05 & 
                                  plot_dat$stacks_LA_P < 0.05 &
                           plot_dat$stacks_TX_P < 0.05 & 
                             plot_dat$stacks_FL_P < 0.05)],
       col=cols["stacks"],cex=1,pch=5,lwd=2)
axis(2,las=1,pos=-1500000)
mtext("FLFW vs. FLCC",2,line=1,cex=0.65)
# add the LG labels
midpts<-tapply(plot_dat$plot.pos,plot_dat$Chrom,median)
text(x=midpts[lgs],y=-0.15)

# add outside legend

opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
            mar=c(0, 0, 0, 0), new=TRUE)
on.exit(par(opar))
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("top",c(expression("Permutation"~italic("F")["ST"]),
         expression("Stacks"~italic("F")["ST"]),
         "PCAdapt",expression(italic("X")^T~italic("X")),"Salinity BF"),
       xjust = 0.5,x.intersp = 0.5,
       col = cols[c("perm","stacks","pc","xtx","sal")],
       pt.bg=cols[c("perm","stacks","pc","xtx","sal")],
       pch=c(4,5,0,3,2),bty='n',ncol=5)
dev.off()
```

And create upset plots to show overlap between pairwise Stacks comparisons and between all outlier methods.


```r
cols<-c(permutations='#e41a1c',salBF='#377eb8',pcadapt='#a65628',
        xtx='#ff7f00',sharedStacks='#f781bf',
        Alabama='#af8dc3',Louisiana='#e7d4e8',Texas='#762a83',Florida='#1b7837')
png("../figs/upsetOutliers.png",res=300,height=4,width=7.5,units="in",pointsize = 20)
upset(fromList(outliers),sets=c("permutations","salBF","xtx","pcadapt","sharedStacks"),
      point.size=3.5,line.size=2,mainbar.y.label = "Number of Shared Outliers",
      sets.x.label = " Number of Outliers",text.scale=c(1.5,1.5,1.5,1.5,1.5,1.5),
      sets.bar.color =cols[c("permutations","salBF","xtx","pcadapt","sharedStacks")],
      margin1scale = 0.2,
      sets.pt.color=cols[c("permutations","salBF","xtx","pcadapt","sharedStacks")])
dev.off()
```

And finally we'll stitch them together into a final image.


```r
library(magick);library(multipanelfigure)

image_files <- c("../figs/FstOutliers.png",
                 "../figs/upsetOutliers.png",
                 "../figs/sharedStacks.png")


png("../figs/fstPlots.png",height=4,width=8,units="in",res=300,pointsize = 11)
figure <- multi_panel_figure(
  width = c(4.5, 3),
  height = c(1.8,1.8),
  unit = "inches",
  row_spacing = 0.0,column_spacing = 0
)
(figure %<>% fill_panel(image_files[1],row = 1:2,scaling="fit",
                        allow_panel_overwriting = TRUE) )
(figure %<>% fill_panel(image_files[3], column=2, row=1,scaling="fit") )
(figure %<>% fill_panel(image_files[2],column=2,row=2,scaling="fit"))
dev.off()
```

![Multipanel figure showing the locations of outlier loci in the genome (A) as well as overlap between the outlier analyses (B and C). This figure is Figure 2 in the main text.](../figs/fstPlots.png)

## Annotations {-}



```r
annotate_snps<-function(snpDF,gff,genome.blast,ID="Locus.ID",
                        chrom="Chr",bp="BP",pos="Column")
{
  fw.sig.reg<-do.call(rbind,apply(snpDF,1,function(sig){
    this.gff<-gff[as.character(gff$seqname) %in% 
                    as.character(unlist(sig[chrom])),]
    description<-NA
    SSCID<-NA
    if(nrow(this.gff)>0){
      
      this.reg<-this.gff[which(this.gff$start <= as.numeric(sig[bp]) & 
                           this.gff$end >= as.numeric(sig[bp])),]
      if(nrow(this.reg) == 0){
        if(as.numeric(sig[bp]) > max(as.numeric(this.gff$end))){
          region<-"beyond.last.contig"
        }else{
          region<-NA
        }
      }else{
        if(length(grep("SSCG\\d+",this.reg$attribute))>0){
          geneID<-unique(gsub(".*(SSCG\\d+).*",
                              "\\1",
                              this.reg$attribute[grep("SSCG\\d+",
                                                      this.reg$attribute)]))
          gene<-genome.blast[genome.blast$sscv4_gene_ID %in% geneID,
                             "blastp_hit_description"]
        }else{
          geneID<-NA
          gene<-NA
        }
        # if there are multiples they'll be in separated by a semi-colon
        region<-paste(this.reg$feature,collapse = ";")
        description<-paste(gene,collapse=";")
        SSCID<-paste(geneID,collapse=";")
      }
    }else{
      region<-"scaffNotFound"
    }
    return(data.frame(Locus=sig[[ID]],Chr=sig[chrom],BP=sig[bp],SNPCol=sig[pos],
                      region=region, description=description,SSCID=SSCID,
                      row.names=NULL))
  }))
}
```


```r
gff.name<-"ssc_2016_12_20_chromlevel.gff.gz"
if(length(grep("gz",gff.name))>0){
  gff<-read.delim(gzfile(paste("../../scovelli_genome/",gff.name,sep="")),header=F)
} else{
  gff<-read.delim(paste("../../scovelli_genome/",gff.name,sep=""),header=F)
}
colnames(gff)<-c("seqname","source","feature","start","end","score",
                 "strand","frame","attribute")
genome.blast<-read.csv("../../scovelli_genome/ssc_2016_12_20_cds_nr_blast_results.csv",
                       skip=1,header=T)#I saved it as a csv
```

```r
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")

snp_annotate<-annotate_snps(fw_SNPinfo,gff,genome.blast,ID="ID",
                            chrom="Chrom",bp="BP",pos = "Pos")
snp_annotate$Locus<-as.character(snp_annotate$Locus)
fw_SNPinfo$ID<-as.character(fw_SNPinfo$ID)
fw_SNPinfo<-merge(fw_SNPinfo,snp_annotate[,-c(2,3)],
                  by.x="ID",by.y="Locus",all.x=TRUE)

saveRDS(fw_SNPinfo,"fw_SNPinfo.RDS")
```


```r
fw_SNPinfo<-readRDS("fw_SNPinfo.RDS")
outliers<-list(xtx=fw_SNPinfo$ID[fw_SNPinfo$XtX >= quantile(fw_SNPinfo$XtX,0.99,na.rm = TRUE)],
               salBF=fw_SNPinfo$ID[fw_SNPinfo$logSalBF>=
                                     quantile(fw_SNPinfo$logSalBF,0.99,na.rm = TRUE)],
               permutations=fw_SNPinfo$ID[
                 rowSums(fw_SNPinfo[,c("perm_TX","perm_FL","perm_AL","perm_LA")])==4],
               pcadapt=fw_SNPinfo$ID[which(fw_SNPinfo$pcadaptQ<0.01)],
               Alabama=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05)], 
               Louisiana=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_LA_P < 0.05)],
               Texas=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_TX_P < 0.05)],
               Florida=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_FL_P < 0.05)],
               sharedStacks=fw_SNPinfo$ID[which(fw_SNPinfo$stacks_AL_P < 0.05 & 
                                                  fw_SNPinfo$stacks_LA_P < 0.05 &
                 fw_SNPinfo$stacks_TX_P < 0.05 & fw_SNPinfo$stacks_FL_P < 0.05)])
out_snps<-fw_SNPinfo[fw_SNPinfo$ID %in% unlist(outliers),]
```


```r
annInfo<-fw_SNPinfo[,c("ID","region")]
annInfo$region<-as.character(annInfo$region)

annInfo$region[grep("UTR",annInfo$region)]<-"regulatory"
annInfo$region[grep("gene",annInfo$region)]<-"coding"
annInfo$region[annInfo$region %in% "contig"]<-"non-coding"
annInfo$region[annInfo$region %in% "scaffNotFound"]<-"unkown"

annInfo$outlier<-"not-outlier"
annInfo$outlier[annInfo$ID %in% unlist(outliers)]<-"outlier"
kable(table(annInfo$region,annInfo$outlier),booktabs=TRUE,
      caption="The number of SNPs in the subsetted dataset that were outliers in coding and non-coding regions of the genome.")
```

\begin{table}

\caption{(\#tab:summarizeAnnotations)The number of SNPs in the subsetted dataset that were outliers in coding and non-coding regions of the genome.}
\centering
\begin{tabular}[t]{lrr}
\toprule
  & not-outlier & outlier\\
\midrule
coding & 4469 & 764\\
non-coding & 3535 & 620\\
regulatory & 363 & 59\\
unkown & 1 & 0\\
\bottomrule
\end{tabular}
\end{table}

If our null hypothesis is that we have randomly selected outliers with equal probability from coding and non-coding regions, we can use a Fisher's exact test to calculate the probability of our observed distribution of outliers.


```r
fisher.test(table(annInfo$region,annInfo$outlier))
```

```
## 
## 	Fisher's Exact Test for Count Data
## 
## data:  table(annInfo$region, annInfo$outlier)
## p-value = 0.8654
## alternative hypothesis: two.sided
```

Based on this, we could conclude that our table of outliers in coding regions is more extreme than expected under the null. 

Now let's look at a more specific set of genes -- those from our salinity-associated gene set.


```r
put_genes<-read.delim("putative_genes.txt",stringsAsFactors = FALSE)
s<-strsplit(put_genes$Scovelli_geneID,",")
genes<-data.frame(Gene=rep(put_genes$Gene,sapply(s,length)),
                  SSCID=unlist(s),stringsAsFactors = FALSE)
dat<-do.call(rbind,apply(genes,1,function(gene,snpinfo){
  if(length(grep(gene["SSCID"],snpinfo$SSCID))>0){
    ID<-cbind(snpinfo$ID[grep(gene["SSCID"],snpinfo$SSCID)])
    out<-data.frame(ID,rbind(gene),row.names = NULL)
    return(out)
  }
},snpinfo=fw_SNPinfo))
dat<-unique(dat)

# combine any duplicate annotations
dati<-do.call(rbind,lapply(unique(dat$ID),function(id,dat){
  td<-dat[dat$ID %in% id,]
  return(data.frame(ID=id,
                    Gene=paste(td$Gene,collapse = ";"),
                    SSCID=paste(td$SSCID,collapse = ";"),
                    stringsAsFactors = FALSE))
},dat=dat))

dati$ID<-as.character(dati$ID)
fw_SNPinfo<-merge(fw_SNPinfo,dati,by="ID",all.x=TRUE)

colnames(fw_SNPinfo)[colnames(fw_SNPinfo)=="SSCID.x"]<-"SSCID"
saveRDS(fw_SNPinfo,"fw_SNPinfo.RDS")
```

```r
annInfo$salgene<-"not-putative"
annInfo$salgene[annInfo$ID %in% fw_SNPinfo$ID[!is.na(fw_SNPinfo$Gene)]]<-"putative"
kable(table(annInfo$outlier,annInfo$salgene),booktabs=TRUE,
      caption="The number of SNPs in the subsetted dataset that were outliers in putative salinity genes or not.")
```

\begin{table}

\caption{(\#tab:putativeGeneAnns)The number of SNPs in the subsetted dataset that were outliers in putative salinity genes or not.}
\centering
\begin{tabular}[t]{lrr}
\toprule
  & not-putative & putative\\
\midrule
not-outlier & 10248 & 117\\
outlier & 1714 & 24\\
\bottomrule
\end{tabular}
\end{table}

```r
fisher.test(table(annInfo$outlier,annInfo$salgene))
```

```
## 
## 	Fisher's Exact Test for Count Data
## 
## data:  table(annInfo$outlier, annInfo$salgene)
## p-value = 0.3967
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  0.7533815 1.9217177
## sample estimates:
## odds ratio 
##   1.226425
```

Putative salinity genes do not appear to be enriched with outliers.



# References {-}

