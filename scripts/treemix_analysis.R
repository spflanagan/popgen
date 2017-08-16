#' Title: Analyze the treemix output
#' Author: Sarah P. Flanagan
#' Date: 15 June 2017
setwd("~/sf_ubuntushare/popgen/fwsw_results/treemix/")
poporder<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST",
            "ALFW","FLSG","FLKB","FLFD","FLSI","FLAB",
            "FLPB","FLHB","FLCC","FLLG")
grp.colors<-c('#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837')
colors<-poporder
colors[colors %in% "FLLG"]<-grp.colors[6]
colors[colors %in% c("FLPB","FLHB","FLCC")]<-grp.colors[6]
colors[colors %in% c("FLAB")]<-grp.colors[5]
colors[colors %in% c("FLSI","FLFD","FLKB","FLSG")]<-grp.colors[3]
colors[colors %in% c("ALST","ALFW","LAFW")]<-grp.colors[2]
colors[colors %in% c("TXSP","TXCC","TXFW","TXCB")]<-grp.colors[1]
write.table(cbind(poporder,colors),"poporder",quote=F,sep='\t')
treemix.dir<-"~/Programs/treemix-1.13/"
source(paste(treemix.dir,"src/plotting_funcs.R",sep=""))#I've modified these functions
library(lattice); library(grid)
source("../../../gwscaR/R/gwscaR.R")
source("../../../gwscaR/R/gwscaR_plot.R")
source("../../../gwscaR/R/gwscaR_utility.R")
source("../../../gwscaR/R/gwscaR_fsts.R")
source("../../../gwscaR/R/gwscaR_popgen.R")


## Basic tree
png("FWSW_treemix.png",height=7,width=11,units="in",res=300)
par(mfrow=c(1,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
tree<-plot_tree("fwsw.basic",plotmig=F,scale=F,mbar=F,plus=0.05)
mtext("Drift parameter",1,line=2)
resid<-plot_resid("fwsw.basic","poporder",wcols="rb")


### FLLG as outgroup
m0<-treemix.cov.plot("fwsw.k100bFLLGr",poporder,split=c(1,1,3,2),more=TRUE)
m1<-treemix.cov.plot("fwsw.k100bFLLGrm1",poporder,split=c(2,1,3,2),more=TRUE)
m2<-treemix.cov.plot("fwsw.k100bFLLGrm2",poporder,split=c(3,1,3,2),more=TRUE)
m3<-treemix.cov.plot("fwsw.k100bFLLGrm3",poporder,split=c(1,2,3,2),more=TRUE)
m4<-treemix.cov.plot("fwsw.k100bFLLGrm4",poporder,split=c(2,2,3,2),more=TRUE)
m5<-treemix.cov.plot("fwsw.k100bFLLGrm5",poporder,split=c(3,2,3,2),more=FALSE)
par(mfrow=c(2,3))
r0<-plot_resid("fwsw.k100bFLLGr","poporder")
r1<-plot_resid("fwsw.k100bFLLGrm1","poporder")
r2<-plot_resid("fwsw.k100bFLLGrm2","poporder")
r3<-plot_resid("fwsw.k100bFLLGrm3","poporder")
r4<-plot_resid("fwsw.k100bFLLGrm4","poporder")
r5<-plot_resid("fwsw.k100bFLLGrm5","poporder")

png("migration_trees_treemix.png",height=6,width=11,units="in",res=300)
par(mfrow=c(2,3),mar=c(1,1,1,1),oma=c(1,1,1,1))
t0<-plot_tree("fwsw.k100bFLLGr",plotmig = F,plus=0.05,scale=F,mbar=F)
t1<-plot_tree("fwsw.k100bFLLGrm1",plus=0.05,scale=F,mbar=F)
t2<-plot_tree("fwsw.k100bFLLGrm2",plus=0.05,scale=F,mbar=F)
t3<-plot_tree("fwsw.k100bFLLGrm3",plus=0.05,scale=F,mbar=F)
t4<-plot_tree("fwsw.k100bFLLGrm4",plus=0.05,scale=F,mbar=F)
t5<-plot_tree("fwsw.k100bFLLGrm5",plus=0.05,scale=F,mbar=F)
dev.off()

#' Evaluate migration p-values
tree0<-read.table(gzfile("fwsw.k100bFLLGr.treeout.gz"), as.is  = T, comment.char = "", quote = "")
tree1<-read.table(gzfile("fwsw.k100bFLLGrm1.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree2<-read.table(gzfile("fwsw.k100bFLLGrm2.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree3<-read.table(gzfile("fwsw.k100bFLLGrm3.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree4<-read.table(gzfile("fwsw.k100bFLLGrm4.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree5<-read.table(gzfile("fwsw.k100bFLLGrm5.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)

png("FWSW_treemix_m3.png",height=7,width=7,units="in",res=300)
t3<-plot_tree("fwsw.k100bFLLGrm3",plus=0.05,scale=F,mbar=T)
dev.off()

##### FLPB as outgroup #####
m0<-treemix.cov.plot("fwsw.k100bFLPBr",poporder,split=c(1,1,3,2),more=TRUE)
m1<-treemix.cov.plot("fwsw.k100bFLPBrm1",poporder,split=c(2,1,3,2),more=TRUE)
m2<-treemix.cov.plot("fwsw.k100bFLPBrm2",poporder,split=c(3,1,3,2),more=TRUE)
m3<-treemix.cov.plot("fwsw.k100bFLPBrm3",poporder,split=c(1,2,3,2),more=TRUE)
m4<-treemix.cov.plot("fwsw.k100bFLPBrm4",poporder,split=c(2,2,3,2),more=TRUE)
m5<-treemix.cov.plot("fwsw.k100bFLPBrm5",poporder,split=c(3,2,3,2),more=FALSE)
par(mfrow=c(2,3))
r0<-plot_resid("fwsw.k100bFLPBr","poporder")
r1<-plot_resid("fwsw.k100bFLPBrm1","poporder")
r2<-plot_resid("fwsw.k100bFLPBrm2","poporder")
r3<-plot_resid("fwsw.k100bFLPBrm3","poporder")
r4<-plot_resid("fwsw.k100bFLPBrm4","poporder")
r5<-plot_resid("fwsw.k100bFLPBrm5","poporder")

png("migration_trees_treemix_FLPB.png",height=6,width=11,units="in",res=300)
par(mfrow=c(2,3),mar=c(1,1,1,1),oma=c(1,1,1,1))
t0<-plot_tree("fwsw.k100bFLPBr",plotmig = F,plus=0.05,scale=F,mbar=F)
t1<-plot_tree("fwsw.k100bFLPBrm1",plus=0.05,scale=F,mbar=F)
t2<-plot_tree("fwsw.k100bFLPBrm2",plus=0.05,scale=F,mbar=F)
t3<-plot_tree("fwsw.k100bFLPBrm3",plus=0.05,scale=F,mbar=F)
t4<-plot_tree("fwsw.k100bFLPBrm4",plus=0.05,scale=F,mbar=F)
t5<-plot_tree("fwsw.k100bFLPBrm5",plus=0.05,scale=F,mbar=F)
dev.off()

#' Evaluate migration p-values
tree0<-read.table(gzfile("fwsw.k100bFLPBr.treeout.gz"), as.is  = T, comment.char = "", quote = "")
tree1<-read.table(gzfile("fwsw.k100bFLPBrm1.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree2<-read.table(gzfile("fwsw.k100bFLPBrm2.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree3<-read.table(gzfile("fwsw.k100bFLPBrm3.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree4<-read.table(gzfile("fwsw.k100bFLPBrm4.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)
tree5<-read.table(gzfile("fwsw.k100bFLPBrm5.treeout.gz"), as.is  = T, comment.char = "", quote = "",skip=1)

d <- read.table("fwsw.k100bFLPBrm3.vertices.gz", as.is  = T, comment.char = "", quote = "")
branch.cols<-rep("black",nrow(d))
branch.cols[d[,2] %in% c("TXFW","ALFW","LAFW","FLLG")]<-"cornflowerblue"

tip.names<-as.vector(d[d[,5] == "TIP",2])
tip.names<-data.frame(Original=tip.names,Replacement=tip.names,stringsAsFactors = FALSE)
tip.names$Replacement[tip.names$Replacement=="FLLG"]<-"FLFW"

png("FWSW_treemix_m3_FLPB.png",height=7,width=7,units="in",res=300)
t3<-plot_tree("fwsw.k100bFLPBrm3","poporder",plus=0.05,scale=F,mbar=F,arrow=0.1,tip.order = tip.names)
ybar<-0.01
mcols = rev( heat.colors(150) )
mcols = mcols[50:length(mcols)]
ymi = ybar+0.15
yma = ybar+0.35
l = 0.2
w = l/100
xma = max(t3$d$x/20)
rect( rep(0.15, 100), ymi+(0:99)*w, rep(0.15+xma, 100), ymi+(1:100)*w, col = mcols, border = mcols)
text(0.15+xma+0.001, ymi, lab = "0", adj = 0, cex = 0.7)
text(0.15+xma+0.001, yma, lab = "0.5", adj = 0, cex =0.7)
text(0.15, yma+0.06, lab = "Migration", adj = 0 , cex = 0.6)
text(0.15, yma+0.03, lab = "weight", adj = 0 , cex = 0.6)
dev.off()
