#' Title: Analyze the treemix output
#' Author: Sarah P. Flanagan
#' Date: 15 June 2017
setwd("~/sf_ubuntushare/popgen/fwsw_results/treemix/")
poporder<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST",
            "ALFW","FLSG","FLKB","FLFD","FLSI","FLAB",
            "FLPB","FLHB","FLCC","FLLG")
write.table(poporder,"poporder",quote=F)
treemix.dir<-"~/Programs/treemix-1.13/"
source(paste(treemix.dir,"src/plotting_funcs.R",sep=""))#I've modified these functions
png("FWSW_treemix.png",height=7,width=11,units="in",res=300)
par(mfrow=c(1,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
tree<-plot_tree("fwsw.basic",plotmig=F,scale=F,mbar=F,plus=0.05)
mtext("Drift parameter",1,line=2)
resid<-plot_resid("fwsw.basic","poporder",wcols="rb")


library(lattice);library(RColorBrewer); library(grid)


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
