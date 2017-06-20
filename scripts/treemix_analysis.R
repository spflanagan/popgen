#' Title: Analyze the treemix output
#' Author: Sarah P. Flanagan
#' Date: 15 June 2017
setwd("~/sf_ubuntushare/popgen/fwsw_results/treemix/")
poporder<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST",
            "ALFW","FLSG","FLKB","FLFD","FLSI","FLAB",
            "FLPB","FLHB","FLCC","FLLG")
write.table(poporder,"poporder",quote=F)
treemix.dir<-"~/Programs/treemix-1.13/"
source(paste(treemix.dir,"src/plotting_funcs.R",sep=""))
png("FWSW_treemix.png",height=7,width=11,units="in",res=300)
par(mfrow=c(1,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
tree<-plot_tree("fwsw.basic",plotmig=F,scale=F,mbar=F,plus=0.05)
mtext("Drift parameter",1,line=2)
resid<-plot_resid("fwsw.basic","poporder",wcols="rb")

#I've modified these functions
library(lattice);library(RColorBrewer); library(grid)
colors<-c("blue","yellow","red")
pal<-colorRampPalette(colors)
ncol=80
cols<-pal(ncol)
#' create a function to plot covariances
#' @param stem The filename basic stem for that run
#' @param poporder The list of populations in the order to plot them.
#' @return cp The covariance matrix in the plotting order.
treemix.cov.plot<-function(stem,poporder){
  cov<-read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  #reorder
  covplot = data.frame(matrix(nrow = nrow(cov), ncol = ncol(cov)))
  for(i in 1:length(poporder)){
    for( j in 1:length(poporder)){
      
      covplot[i, j] = cov[which(names(cov)==poporder[i]), which(names(cov)==poporder[j])]
      rownames(covplot)[i]<-poporder[i]
      colnames(covplot)[j]<-poporder[j]
    }
  }
  cp<-as.matrix(covplot)
  cp[lower.tri(cp)]<-NA
  cp[upper.tri(cp)]<-covplot[upper.tri(covplot)]
  levelplot(cp,col.regions=cols,alpha.regions=0.7,
            scales = list(x=list(rot=90),tck = 0),xlab="",ylab="")
  trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
  grid.text("covariance", 0.2, 0, hjust=0.5, vjust=1.2)
  trellis.unfocus()
  return(cp)
}

par(mfrow=c(2,3))
m0<-treemix.cov.plot("fwsw.k100bFLLGr",poporder)
m1<-treemix.cov.plot("fwsw.k100bFLLGrm1",poporder)
m2<-treemix.cov.plot("fwsw.k100bFLLGrm2",poporder)
m3<-treemix.cov.plot("fwsw.k100bFLLGrm3",poporder)
m4<-treemix.cov.plot("fwsw.k100bFLLGrm4",poporder)
m5<-treemix.cov.plot("fwsw.k100bFLLGrm5",poporder)
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
