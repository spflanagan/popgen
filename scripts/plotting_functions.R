
#############################################################################
#***************************************************************************#
#################################FUNCTIONS###################################
#***************************************************************************#
#############################################################################

#***************************************************************************#
#PLOT ANY GENOME-WIDE STATISTIC
#***************************************************************************#
plotting.genome.wide<-function(bp,var,y.max,x.max, rect.xs,y.min=0,x.min=0, 
	plot.new=FALSE, plot.axis=TRUE, rect.color="white", pt.cex=1){
	par(new=new)
	plot(bp, var,xlab="",ylab="", 
		type="n", bg="transparent", axes=F, bty="n", 
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	num.rect<-nrow(rect.xs)
	if(is.null(num.rect)) {
		rect(rect.xs[1],y.min,rect.xs[2],y.max, 
				col=rect.color, border=NA)
	} else {
		for(i in 1:nrow(rect.xs)){
			rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
				col=rect.color, border=NA)
		}
	}
	if(plot.axis){
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
	points(bp, var, pch=19, cex=pt.cex,
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

#***************************************************************************#
##AUTOMATED FST PLOTTING
#***************************************************************************#
plotting.fsts.scaffs<-function(dat, dat.name, ci.dat=NULL, 
	col.pts=NULL, col.pt.col="dark green", col.pt.pch=8, col.pt.cex=2){
	byscaff<-split(dat, factor(dat$Chr))#dat is a fst file from stacks_genomewideCIs
	if(!is.null(col.pts)){
		col.pt.byscaff<-split(col.pts, factor(col.pts$Chr))}
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	for(i in 1:length(byscaff)){
		new.max<-last.max+round(max(byscaff[[i]]$BP), -2)
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	
	new.x<-list()
	col.pt.x<-NULL
	col.pt.y<-NULL
	for(i in 1:length(byscaff)){
		new.x[[i]]<-byscaff[[i]]$BP+addition.values[i]
		if(!is.null(col.pts)){
		j<-match(names(byscaff[i]),names(col.pt.byscaff))
		if(!is.na(j)){
			col.pt.x<-c(col.pt.x, byscaff[[i]][(
				byscaff[[i]]$Locus %in% 
				col.pt.byscaff[[j]]$Locus) &
				(byscaff[[i]]$BP %in% col.pt.byscaff[[j]]$BP)
				,"BP"]+addition.values[i])
			col.pt.y<-c(col.pt.y, byscaff[[i]][(
				byscaff[[i]]$Locus
				 %in% col.pt.byscaff[[j]]$Locus) &
				(byscaff[[i]]$BP %in% col.pt.byscaff[[j]]$BP),
				"SmoothFst"])
		}}
	}

	x.min<-min(addition.values)
	x.max<-max(addition.values)
	y.max<-max(dat$SmoothFst)+0.5*max(dat$SmoothFst)
	if(min(dat$SmoothFst) < 0) {
		y.min<-min(dat$SmoothFst) + 0.5*min(dat$SmoothFst)
	} else {
		y.min<-0
	}

	plot(new.x[[1]], byscaff[[1]]$SmoothFst, 
		xlim=c(x.min,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(byscaff)){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		plot.genome.wide(new.x[[j]], 
			byscaff[[j]]$SmoothFst,
			y.max,x.max, rect.xs[j,],y.min=y.min,x.min=x.min, 
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.25)
	}
		
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=10)),
		xlim=c(x.min,x.max),ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
	mtext(dat.name,2,las=1, line=4,cex=.75)
	if(!is.null(ci.dat)){
		clip(x.min,x.max,y.min,y.max)
		abline(h=ci.dat$CI99smooth,col="red")
		if(!is.null(col.pts)){
			clip(x.min,x.max,ci.dat$CI99smooth, y.max)
			points(col.pt.x, col.pt.y,
				col=col.pt.col, pch=col.pt.pch, cex=col.pt.cex)
		}
	}
}




#***************************************************************************#
#PLOT A STRUCTURE BARPLOT
#***************************************************************************#

plotting.structure<-function(structure.out, k, pop.order, 
	filename=paste("str.k",k,".jpeg",sep=""),make.file=TRUE,
	plot.new=TRUE,colors=NULL,xlabel=TRUE,ylabel=NULL){
	str.split<-split(structure.out,structure.out[,1])
	if(is.null(colors)){
		bar.colors<-rainbow(k,s=0.5)
	} else {
		bar.colors<-colors
	}
	if(make.file==TRUE){
		jpeg(filename,width=7, height=1.25, units="in", res=300)
		par(mfrow=c(1,length(str.split)))
	} 
	#par(mfrow=c(1,length(str.split)),mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
	for(i in 1:length(str.split)){
		pop.index<-pop.order[i]
		barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
			beside=FALSE, space=0,	border=NA, col=bar.colors,
			xlab="", ylab="", xaxt='n', yaxt='n')#, new=plot.new)
		if(xlabel==TRUE){
			mtext(pop.index, 1, line=0.5, cex=1, outer=F)}
		if(!is.null(ylabel)){
			if(i == 1) { mtext(ylabel,2,cex=1) }
		}	
	}
	if(make.file==TRUE) {dev.off()}
}

