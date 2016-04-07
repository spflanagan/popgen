
library(maps);library(gplots)
library(mapdata)

setwd("B:/ubuntushare/popgen/sw_results/migrate")
mar.coor<-read.csv("../marine_coordinates_revised.csv", header=T)

migrate<-read.table("migrate_out_all.txt")
colnames(migrate)<-c("Locus","Parameter","2.5","25","Mode","75","97.2",
	"Median","Mean")
migrate.order<-c("FLKB","FLPB","TXSP","ALST","FLAB","FLCC","FLFD","FLHB",
	"FLSG","FLSI","TXCB","TXCC")
migrate$ParameterType<-substr(migrate$Parameter,1,1)
thetas<-migrate[migrate$ParameterType=="T",]
thetas$Pop<-migrate.order
thetas<-thetas[match(mar.coor$site,thetas$Pop),]
theta.cex<-thetas$Mean*1000
m<-migrate[migrate$ParameterType=="M",]
m$start<-migrate.order[as.numeric(gsub("M_(\\d+)->\\d+","\\1",m$Parameter))]
m$end<-migrate.order[as.numeric(gsub("M_\\d+->(\\d+)","\\1",m$Parameter))]
pairwise.pops<-rbind(c("FLCC","FLHB"),c("FLHB","FLCC"),c("FLHB","FLPB"),
	c("FLPB","FLHB"),c("FLPB","FLAB"),c("FLAB","FLPB"),c("FLAB","FLSI"),
	c("FLSI","FLAB"),c("FLSI","FLFD"),c("FLFD","FLSI"),c("FLFD","FLKB"),
	c("FLKB","FLFD"),c("FLKB","FLSG"),c("FLSG","FLKB"),c("FLSG","ALST"),
	c("ALST","FLSG"),c("ALST","TXCB"),c("TXCB","ALST"),c("TXCB","TXCC"),
	c("TXCC","TXCB"),c("TXCC","TXSP"),c("TXSP","TXCC"))
pairwise.m<-do.call("rbind",apply(pairwise.pops,1,function(x){
	m[m$start==x[1]&m$end==x[2],] }))

add.arrows<-function(row,shift.x,shift.y){
	arrows(x0=mar.coor[mar.coor$site==pairwise.m$start[row],"lon"]+shift.x,
		y0=mar.coor[mar.coor$site==pairwise.m$start[row],"lat"]+shift.y,
		x1=mar.coor[mar.coor$site==pairwise.m$end[row],"lon"]+shift.x,
		y1=mar.coor[mar.coor$site==pairwise.m$end[row],"lat"]+shift.y,
		lwd=pairwise.m$Mean[row]/50,length=0.1)
}


png("migrate_map.png",height=7,width=14,units="in",res=300)
par(oma=c(0,0,0,0),mar=c(0,0,0,0),pin=c(7,7))
map("worldHires", "usa",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray90", mar=c(0,0,0,0),fill=TRUE, res=300,myborder=0)
map("worldHires", "mexico",xlim=c(-100,-76), ylim=c(24,32), 
	col="gray95", fill=TRUE, add=TRUE)
points(mar.coor$lon, mar.coor$lat,  col="black", cex=theta.cex*.5, pch=19)
abline(h=c(25,30),lty=3)
abline(v=c(-80,-85,-90,-95,-100),lty=3)
text(x=c(-99.5,-99.5,-99.5),y=c(25,30,35),c("25N","30N"))
text(x=c(-80,-85,-90,-95),y=rep(31.7,4),c("80W","85W","90W","95W"))
#text(y=26,x=-90,"Gulf of Mexico")
#text(x=-88,y=32,"USA")
#text(x=-78,y=29.5,"Atlantic Ocean")
text(x=-98,y=26,"TXSP",font=2)
text(x=-98.2,y=27,"TXCC",font=2)
text(x=-96,y=29,"TXCB",font=2)
text(x=-88,y=30.75,"ALST",font=2)
text(x=-85,y=30,"FLSG",font=2)
text(x=-82.75,y=29.75,"FLKB",font=2)
text(x=-82,y=28,"FLFD",font=2)
text(x=-81.5,y=26.5,"FLSI",font=2)
text(x=-80,y=24.5,"FLAB",font=2)
text(x=-80.6,y=26.7,"FLPB",font=2)
text(x=-80.75,y=27.2,"FLHB",font=2)
text(x=-81.3,y=28.2,"FLCC",font=2)
#FL Atlantic
add.arrows(1,0.5,0)
add.arrows(2,1,0)
add.arrows(3,0.5,0)
add.arrows(4,1,0)
add.arrows(5,0.5,0)
add.arrows(6,1,0)
#FL Gulf
add.arrows(7,-0.5,0)
add.arrows(8,-1,0)
add.arrows(9,-0.5,0)
add.arrows(10,-1,0)
add.arrows(11,-0.5,0)
add.arrows(12,-1,0)
add.arrows(13,0,0.9)
add.arrows(14,0,0.4)
add.arrows(15,0,-0.5)
add.arrows(16,0,-1)
#TX Gulf
add.arrows(17,0,0.25)
add.arrows(18,0,0.75)
add.arrows(19,0.5,0)
add.arrows(20,1,0)
add.arrows(21,0.5,0)
add.arrows(22,1,0)
dev.off()

pop.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB","FLFD","FLSI",
	"FLAB","FLPB","FLHB","FLCC")
migrate.matrix.dat<-m[,c("Mean","start","end")]
migrate.matrix<-as.data.frame(setNames(
	replicate(length(pop.list),numeric(0), simplify = F), pop.list))
for(i in 1:(length(pop.list)-1)){
	for(j in (i+1):length(pop.list)){
		mi<-migrate.matrix.dat[migrate.matrix.dat$start==pop.list[i] &
			migrate.matrix.dat$end==pop.list[j],"Mean"]
		mj<-migrate.matrix.dat[migrate.matrix.dat$start==pop.list[j] &
			migrate.matrix.dat$end==pop.list[i],"Mean"]
		migrate.matrix[pop.list[i],pop.list[j]]<-mi
		migrate.matrix[pop.list[j],pop.list[i]]<-mj
	}
}
migrate.matrix<-as.matrix(migrate.matrix)
diag(migrate.matrix)<-thetas$Mean
write.table(migrate.matrix,"migrate_matrix.txt",col.names=T,row.names=T,quote=F)


