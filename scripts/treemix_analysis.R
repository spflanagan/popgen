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
cov<-read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
#reorder
covplot = data.frame(matrix(nrow = nrow(c), ncol = ncol(c)))
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