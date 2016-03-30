
#############################################################################
#***************************************************************************#
#################################FUNCTIONS###################################
#***************************************************************************#
#############################################################################

#***************************************************************************#
#PLOT ANY GENOME-WIDE STATISTIC
#***************************************************************************#
plotting.genome.wide<-function(bp,var,y.max,x.max, rect.xs=NULL,y.min=0,x.min=0, 
	plot.new=FALSE, plot.axis=TRUE, rect.color="white",plot.rect=TRUE, 
	pt.cex=1, pt.col="black"){
	#********************************************
	#this function plots a variable without scaffold info. 
	#feed it the basepair (x) values and variable (y) values 
	#*********************************************
	if(plot.new==TRUE){ par(new=new) }
	plot(bp, var,xlab="",ylab="", new=plot.new,
		type="n", bg="transparent", axes=F, bty="n", 
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	if(plot.rect==TRUE){
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
	}
	if(plot.axis){
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
	points(bp, var, pch=19, cex=pt.cex,col=pt.col,
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

#***************************************************************************#
##AUTOMATED FST PLOTTING
#***************************************************************************#
plotting.fsts.scaffs<-function(dat, dat.name, ci.dat=NULL, pt.cex=1,y.lab=NULL,
	col.pts=NULL, col.pt.col="dark green", col.pt.pch=8, col.pt.cex=2){
	#********************************************
	#this function plots every scaffold. 
	#dat is a fst file from stacks_genomewideCIs
	#dat should have a column named "Chr", one called "BP"
	#at least one other column is needed: dat.name
	#*********************************************
	byscaff<-split(dat, factor(dat$Chr))
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
				dat.name])
		}}
	}

	x.min<-min(addition.values)
	x.max<-max(addition.values)
	y.max<-max(dat[,dat.name])+0.5*max(dat[,dat.name])
	if(min(dat[,dat.name]) < 0) {
		y.min<-min(dat[,dat.name]) + 0.5*min(dat[,dat.name])
	} else {
		y.min<-0
	}

	plot(new.x[[1]], byscaff[[1]][,dat.name], 
		xlim=c(x.min,x.max), ylim=c(y.min, y.max), bty="n",type="n",
		axes=F, xlab="", ylab="")
	for(j in 1:length(byscaff)){
		if(j%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		rect(rect.xs[j,1],y.min,rect.xs[j,2],y.max, 
			col=rect.color, border=NA)
	}
	for(j in 1:length(byscaff)){
		points(new.x[[j]],byscaff[[j]][,dat.name], pch=19, cex=pt.cex,
			xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	}
		#plot.genome.wide(new.x[[j]], 
		#	byscaff[[j]][,dat.name],rect.xs=rect.xs[j,],
		#	y.max,x.max, y.min=y.min,x.min=x.min, 
		#	plot.new=FALSE, plot.axis=FALSE, pt.cex=0.25)
		
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=10)),
		xlim=c(x.min,x.max),ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)
	if(is.null(y.lab)){
		mtext(dat.name,2,las=1, line=4,cex=.75)}
	else{
		mtext(y.lab,2,las=1, line=4,cex=.75)}
	if(!is.null(ci.dat)){
		clip(x.min,x.max,y.min,y.max)
		abline(h=ci.dat$CI99smooth,col="red")
		if(!is.null(col.pts)){
			clip(x.min,x.max,ci.dat$CI99smooth, y.max)
			points(col.pt.x, col.pt.y,
				col=col.pt.col, pch=col.pt.pch, cex=col.pt.cex)
		}
	}
	return(as.data.frame(cbind(old.bp=dat$BP,new.bp=unlist(new.x), 
		locus=dat$Locus)))
}

#***************************************************************************#
##REORDER A DATAFRAME
#***************************************************************************#
reorder.df<-function(dat,order.list){
	#dat has to have the grouping IDs in row 1
	#those grouping ids must match the factors in order.list
	dat.sep<-split(dat, dat[,1])
	dat.new<-dat.sep[[order.list[1]]]
	for(i in 2:length(order.list)){
		dat.new<-rbind(dat.new, dat.sep[[order.list[i]]])
	}
	return(dat.new)
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

