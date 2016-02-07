#Author: Sarah P. Flanagan
#Date: 21 August 2014
#Purpose: Plot Fsts from stacks populations
#first run the data through stacks_genomeCIs (C++ program) to generate 
#genome-wide confidence intervals.

rm(list = ls())
setwd("C://Users//Sarah//Documents//Research//Scovelli11//RAD")
#setwd("C://Users//Sarah//Documents//Research//PopGen//ddRAD//Plate2//analysis")

fst.dat<-read.delim("fst_Female-Male.txt", sep='\t', header=T)
#fst.dat<-read.delim("fst_FLCC-FLFD.txt", sep='\t', header=T)
ci.dat<-read.delim("fst_Female-Male_summary.txt", sep='\t', header=T)
#ci.dat<-read.delim("fst_FLCC-FLFD_summary.txt", sep='\t', header=T)
sig.scaffolds<-NULL

plot.fsts<-function(data, sub, ci, sig.only){
	name<-paste(sub,".png", sep="")
	x<-data$BP[data$Chrom==sub]
	y<-data$SmoothFst[data$Chrom==sub]
	if(sig.only == TRUE && max(y) >= ci[1,6])
	{
		if(max(y)+0.075*max(y) < ci[1,7]){
			ymax<-ci[1,7]+0.075*ci[1,7]}
		else{
			ymax<-max(y)+0.075*max(y)}
		png(name)
		par(mar=c(4,4,1,1),oma=c(2,1,1,1), bty="l", mgp=c(2,0.5,0))
		plot(x, y, xlab = "", ylab="",pch=19,
			xlim=c(0,max(x)), ylim=c(0,ymax))
		clip(0,max(x),0,ymax)
		abline(h=ci[1,6], lwd=4, lty=2, col="blue")#95
		abline(h=ci[1,7], lwd=4, lty=3, col="red")#99
		axis(1, xlab = "Location on Chromosome")
		axis(2, ylab = "Smoothed Fst")
		title(xlab="Location on Chromosome", 
			ylab = "Smoothed Fst")
		legend("topright", c("99% genome-wide CI", "95% genome-wide CI"),
			col=c("red", "blue"), lty=c(2,3), box.lwd=0.5, 
			cex=0.75, ncol=2)
		mtext(text=sub, side = 1, outer=TRUE)
		dev.off()
		
	}
	if(sig.only == FALSE)
	{
		if(max(y)+0.075*max(y) < ci[1,7]){
			ymax<-ci[1,7]+0.075*ci[1,7]}
		else{
			ymax<-max(y)+0.075*max(y)}
		png(name)
		par(mar=c(4,4,1,1),oma=c(2,1,1,1), bty="l", mgp=c(2,0.5,0))
		plot(x, y, xlab = "", ylab="",
			pch=19,xlim=c(0,max(x)), ylim=c(0,ymax))
		clip(0,max(x),0,ymax)
		abline(h=ci[1,6], lwd=4, lty=2, col="blue")#95
		abline(h=ci[1,7], lwd=4, lty=3, col="red")#99
		title(xlab="Location on Chromosome", 
			ylab = "Smoothed Fst")
		legend("topright", c("99% genome-wide CI", "95% genome-wide CI"),
			col=c("red", "blue"), lty=c(2,3), box.lwd=0.5, 
			cex=0.75, ncol=2)
		mtext(text=sub, side = 1, outer=TRUE)
		
		dev.off()
	}
}
#subtest<-levels(fst.dat$Chrom)[1]
#plot.fsts(fst.dat, subtest, ci.dat,FALSE)#for testing purposes

#produce a graph for each chromosome
for(x in 1:length(levels(fst.dat$Chrom)))
{
	sub<-levels(fst.dat$Chrom)[x]
	#plot.fsts(fst.dat, sub, ci.dat, sig.only=TRUE)
	if(max(fst.dat$SmoothFst[fst.dat$Chrom==sub]) >= ci.dat[1,6])
	{
		sig.scaffolds<-c(sig.scaffolds, sub)
	}
}
write.table(sig.scaffolds, "Female-Male_sig_scaffolds.txt", sep='\t', 
	col.names=FALSE, row.names=FALSE)

###Figure out how to pass variables (like pathname and function choices) 
#to r script, and also how to set defaults for the function parameters
