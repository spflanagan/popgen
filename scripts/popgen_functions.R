#Author: Sarah P. Flanagan
#Last Updated: 13 July 2016
#Functions for use in population genomics analyses

#############################################################################
#***************************************************************************#
#################################FUNCTIONS###################################
#***************************************************************************#
#############################################################################

#############################################################################
########################******CALCULATIONS*****##############################
#############################################################################
#***************************************************************************#
#CALCULATE PAIRWISE FSTS PER LOCUS
#***************************************************************************#
pairwise.fst<-function(ped,allele1,allele2,pop.order){
	#V1 of ped should be pop index
	ped.split<-split(ped[,c(allele1,allele2)], factor(ped[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(pop.order),numeric(0), simplify = F), pop.order))
	for(i in 1:(length(pop.order)-1)){
	  for(j in (i+1):length(pop.order)){
		pop1<-factor(ped.split[[pop.order[i]]][ped.split[[pop.order[i]]]!="0"])
		pop2<-factor(ped.split[[pop.order[j]]][ped.split[[pop.order[j]]]!="0"])
		freq1<-summary(pop1)/sum(summary(pop1))	
		freq2<-summary(pop2)/sum(summary(pop2))	
		freqall<-summary(as.factor(c(pop1,pop2)))/
			sum(summary(as.factor(c(pop1,pop2))))
		if(length(freq1)>1){ hs1<-2*freq1[1]*freq1[2] 
		} else {
			hs1<-0
		}
		if(length(freq2)>1){ hs2<-2*freq2[1]*freq2[2] 
		} else {
			hs2<-0
		}
		if(length(freqall)>1){
			hs<-mean(c(hs1,hs2))
			ht<-2*freqall[1]*freqall[2]
			fst<-(ht-hs)/ht
		}
		if(length(freqall)<=1){ fst<-1 }
		dat.var[pop.order[i],pop.order[j]]<-fst
	  }
	}
	dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
	rownames(dat.var)<-colnames(dat.var)
	return(as.matrix(dat.var))
}

#***************************************************************************#
#CALCULATE ISOLATION BY DISTANCE PER LOCUS
#***************************************************************************#
fst.ibd.byloc<-function(ped.file,dist.mat,pop.order){
	results.mantel<-data.frame()
	for(i in seq(7,ncol(ped.file),2)){
		res<-mantel.rtest(
			as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
			as.dist(t(dist.mat)), nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	colnames(results.mantel)<-c("Obs","P")
	return(results.mantel)
}

#***************************************************************************#
#CALCULATE PAIRWISE PST BETWEEN POPULATION PAIRS
#***************************************************************************#
pairwise.pst<-function(dat, pop.order){
	#first column must be pop id/grouping factor
	library(nlme)
	dat.split<-split(dat, factor(dat[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(pop.order),numeric(0), simplify = F), pop.order))
	for(i in 1:(length(pop.order)-1)){
	  for(j in (i+1):length(pop.order)){
		temp.data<-rbind(as.data.frame(dat.split[[pop.order[i]]]),
			as.data.frame(dat.split[[pop.order[j]]]))
		colnames(temp.data)<-c("PopID","Var")
		temp.data$PopID<-factor(temp.data$PopID)
		anv <- lme(fixed=Var ~ 1, random=~1|PopID,data=temp.data)
		varcomp <- VarCorr(anv)
		v.btwn<- as.numeric(varcomp[1])
		v.wthn <- as.numeric(varcomp[2])
		pst <- v.btwn/(v.btwn+2*v.wthn)
		dat.var[pop.order[i],pop.order[j]]<-pst
		#aov.var<-summary.aov(
		#	aov(temp.data[,2]~temp.data[,1]))[[1]]$`Sum Sq`
		#aov.df<-summary.aov(
		#	aov(temp.data[,2]~temp.data[,1]))[[1]]$`Df`
		#dat.var[pop.order[i],pop.order[j]]<-aov.var[2]/(aov.var[2]+
		#	(2*(aov.var[1]/(aov.df[2]-1))))	
	  }
	}
	dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
	rownames(dat.var)<-colnames(dat.var)
	return(dat.var)
}

#***************************************************************************#
#CALCULATE PAIRWISE PSTS FOR ALL TRAITS AND TEST FOR IBD
#***************************************************************************#
all.traits.pst.mantel<-function(trait.df,comp.df,id.index){
	results.mantel<-data.frame()
	for(i in 3:ncol(trait.df)){
		res<-mantel.rtest(
			as.dist(t(pairwise.pst(trait.df[,c(id.index,i)],pop.order))),
			as.dist(t(comp.df)), nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	rownames(results.mantel)<-colnames(trait.df)[3:ncol(trait.df)]
	colnames(results.mantel)<-c("Obs","P")
	return(results.mantel)
}


#***************************************************************************#
#COMPARE FST AND PST PER LOCUS
#***************************************************************************#
fst.pst.byloc<-function(ped.file,trait.df,pop.order,trait.ind){
	results.list<-list()
	for(j in 3:ncol(trait.df)){
	results.mantel<-data.frame()
	for(i in seq(7,ncol(ped.file),2)){
		res<-mantel.rtest(
			as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
			as.dist(t(pairwise.pst(trait.df[,c(trait.ind,j)],pop.order))),
			 nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	colnames(results.mantel)<-c("Obs","P")
	results.list<-append(results.list,data.frame(results.mantel))
	}
	#names(results.list)<-colnames(trait.df)[3:ncol(trait.df)]
	return(results.list)
}

#***************************************************************************#
#CALCULATE STANDARD ERROR OF THE MEAN
#***************************************************************************#
sem<-function(x){
	sem<-sd(x)/sqrt(length(x))
	return(sem)
}

#############################################################################
########################********PLOTTING*******##############################
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
plotting.fsts.scaffs<-function(fst.dat, fst.name="Fst",chrom.name="Chr",
	bp.name="BP", pt.lty=0,pt.col="grey7",new=T,
	ci.dat=NULL, pt.cex=1,y.lab=NULL,axis.size=0.5,scaffold.order=NULL,
	groups=NULL,print.names=FALSE,pt.pch=19,
	sig.col="dark green", col.pt.pch=8, col.pt.cex=2){
		if(!is.null(scaffold.order)){
		scaff.ord<-scaffold.order$component_id
		lgs<-scaffold.order$object
	} else{
		scaff.ord<-levels(factor(fst.dat[,chrom.name]))
		lgs<-scaff.ord
	}
	if(!is.null(groups)){
		lgs<-groups
		scaff.ord<-groups
	}
	all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	xlist<-NULL
	xs<-NULL
	for(i in 1:length(scaff.ord)){
		all.scaff[[scaff.ord[i]]]<-
			all.scaff[[scaff.ord[i]]][order(all.scaff[[scaff.ord[i]]][,bp.name]),]	
		all.scaff[[scaff.ord[i]]][,bp.name]<-
			seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1)
		xs<-c(xs, seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1))
		new.max<-max(xs)
		#scaffold.order[i,"new_start"]<-last.max
		#scaffold.order[i,"new_end"]<-new.max
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		rownames(rect.xs)[i]<-scaff.ord[i]
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	#change BP to plot
	x.max<-max(xs)
	x.min<-min(xs)
	y.max<-max(fst.dat[,fst.name])+0.1*max(fst.dat[,fst.name])
	y.min<-min(fst.dat[,fst.name])-0.1*min(fst.dat[,fst.name])
	if(min(fst.dat[,fst.name]) < 0) {
		y.min<-min(fst.dat[,fst.name]) - 0.1*min(fst.dat[,fst.name])
	} else {
		y.min<-0
	}
	displacement<-y.min-((y.max-y.min)/30)
	if(new==T){
		plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), 
			ylim=c(y.min, y.max), 
			bty="n",type="n",	axes=F, xlab="", ylab="")
		for(i in 1:nrow(rect.xs)){
			if(i%%2 == 0) {
				rect.color<-"white"
			} else {
				rect.color<-"gray75"
			}
			rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
				col=rect.color, border=NA)
			if(print.names==T){
				text(x=mean(all.scaff[[scaff.ord[i]]][
					all.scaff[[scaff.ord[i]]]$Chrom==rownames(rect.xs)[i],
					bp.name]),
					y=displacement,labels=rownames(rect.xs)[i],
					adj=1,xpd=T,srt=45)
			}
		}
	}
	for(i in 1:length(scaff.ord)){
		if(pt.lty==0){
			points(all.scaff[[scaff.ord[i]]][,bp.name], 
				all.scaff[[scaff.ord[i]]][,fst.name], 
				pch=pt.pch, cex=0.5,col=pt.col,
				xlim=c(x.min,x.max),ylim=c(y.min, y.max))
		} else {
			lines(all.scaff[[scaff.ord[i]]][,bp.name], 
				all.scaff[[scaff.ord[i]]][,fst.name], 
				lty=pt.lty, cex=0.5,col=pt.col,
				xlim=c(x.min,x.max),ylim=c(y.min, y.max))

		}
		temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] >= ci.dat[1],]
		points(temp.sig[,bp.name], temp.sig[,fst.name], 
			col=sig.col[1], pch=col.pt.pch, cex=col.pt.pch)
		temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] <= ci.dat[2],]
		points(temp.sig[,bp.name], temp.sig[,fst.name], 
			col=sig.col[2], pch=col.pt.pch, cex=col.pt.pch)
	}
	if(axis.size>0){
		axis(2, at = seq(round(y.min,2),round(y.max,2),
			round((y.max-y.min)/2, digits=2)),
			ylim = c(y.min, y.max), pos=0,
			labels=seq(round(y.min,2),round(y.max,2),
				round((y.max-y.min)/2, digits=2)),
			las=1,tck = -0.01, xlab="", ylab="", cex.axis=axis.size)
	}
	xes<-do.call("rbind",all.scaff)
	return(xes)
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
	filename=paste("str.k",k,".jpeg",sep=""),make.file=TRUE,lab.cex=1,
	plot.new=TRUE,colors=NULL,xlabel=TRUE,ylabel=NULL,xlabcol=NULL){
	str.split<-split(structure.out,structure.out[,1])
	if(is.null(colors)){
		bar.colors<-rainbow(k,s=0.5)
	} else {
		bar.colors<-colors
	}
	if(is.null(xlabcol)){
	  xlabcol<-rep("black",length(str.split))
	}
	if(make.file==TRUE){
		jpeg(filename,width=7, height=1.25, units="in", res=300)
		par(mfrow=c(1,length(str.split)),mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
	} else { 
	  if(plot.new==TRUE){
  	  par(mfrow=c(1,length(str.split)),mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
	  }
  }
	for(i in 1:length(str.split)){
		pop.index<-pop.order[i]
		barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
			beside=FALSE, space=0,	border=NA, col=bar.colors,
			xlab="", ylab="", xaxt='n', yaxt='n')#, new=plot.new)
		if(xlabel==TRUE){
			mtext(pop.index, 1, line=0.5, cex=lab.cex, outer=F,col=xlabcol[i])}
		if(!is.null(ylabel)){
			if(i == 1) { mtext(ylabel,2,cex=lab.cex) }
		}	
	}
	if(make.file==TRUE) {dev.off()}
}

