rm(list=ls())
library(ggplot2)
setwd("C:/Users/Sarah/OneDrive/blastresults/popgen/bayenv")


go.plot<-function(file.list, file.name, analyses=NULL){
	dat<-read.table(file.list[1],skip=1,sep='\t')
	if(is.null(analyses)){
		analysis.names<-gsub("(\\w+)_\\w+.*","\\1",file.list)
	} else {
		analysis.names<-analyses
	}
	colnames(dat)<-c("GO","Number")
	dat$Analysis<-analysis.names[1]
	for(i in 2:length(file.list)){
		d<-read.table(file.list[i],skip=1,sep='\t')
		colnames(d)<-c("GO","Number")
		d$Analysis<-analysis.names[i]
		dat<-rbind(dat,d)
	}

	for(i in 1:length(levels(dat$GO))){
		t<-dat[dat$GO %in% levels(dat$GO)[i],]
		if(nrow(t)<length(analysis.names)){
			a.new<-analysis.names[!(analysis.names %in% t$Analysis)]
			g.new<-rep(t$GO[1],length(a.new))
			n.new<-rep(0,length(a.new))
			t.add<-data.frame(GO=g.new,Number=n.new,Analysis=a.new)
			dat<-rbind(dat,t.add)
		}
	}
	dat<-dat[order(dat$GO),]
	library(ggplot2)
	file.name<-paste(file.name,'.jpeg',sep="")
	jpeg(file.name,height=9,width=9,units="in",res=300)
	par(mar=c(2,2,2,2),oma=c(2,2,2,2),cex=2,lwd=1.3)
	p<-ggplot(dat,aes(factor(GO),Number,fill = factor(Analysis))) + 
		geom_bar(stat="identity",position="dodge") + 
		scale_fill_brewer(palette="Set1",name="Analysis") +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		coord_flip() +
		xlab("Gene Ontology") + ylab("Number")
	print(p)
	dev.off()
	return(dat)
}	


bio.files<-list.files(pattern="bio.txt")
bio2.files<-list.files(pattern="bio2.txt")
cell.files<-list.files(pattern="cell.txt")
cell2.files<-list.files(pattern="cell2.txt")
mol.files<-list.files(pattern="mol.txt")
mol2.files<-list.files(pattern="mol2.txt")

analysis.list<-c("Collection Salinity","Collection Temperature",
	"Seagrass Cover", "10-yr Salinity","10-yr Temperature",
	"10-yr Temperature Variance")

bio.dat<-go.plot(bio.files,"BF_Biology",analysis.list)
bio2.dat<-go.plot(bio2.files,"BF_Biology2",analysis.list)
cell.dat<-go.plot(cell.files,"BF_Cell",analysis.list)
cell2.dat<-go.plot(cell2.files,"BF_Cell2",analysis.list)
mol.dat<-go.plot(mol.files,"BF_Molecular",analysis.list)
mol2.dat<-go.plot(mol2.files,"BF_Molecular2",analysis.list)