#Author: Sarah P. Flanagan
#Date: 27 July 2015
#Purpose: Filter the environmental data from World Ocean Database
#The data have already been consolidated using parse_wod_data.cpp
#one file per collection region.

setwd("B:/ubuntushare/popgen/sw_results/environmental_assoc/environ_files")
wod.files<-list.files(pattern="_out.txt")

wod.dat<-data.frame()

for(i in 1:length(wod.files)){
	wod.dat<-rbind(wod.dat,read.delim(wod.files[i], sep="\t"))
}
#keep only those from 2004-2014 (2014 is newest)
wod.filt<-wod.dat[wod.dat$Year >= 2004,]


mar.coor<-read.csv("marine_coordinates_revised.csv")
adj.coor<-as.data.frame(cbind(site=as.character(mar.coor$site), 
	lat.l=as.numeric(mar.coor$lat-0.5), 
	lat.r=as.numeric(mar.coor$lat),
	lat.h=as.numeric(mar.coor$lat+0.5), 
	lon.l=as.numeric(mar.coor$lon-0.5),
	lon.r=as.numeric(mar.coor$lon),
	lon.h=as.numeric(mar.coor$lon+0.5)), stringsAsFactors = FALSE)
adj.coor[,2:7]<-as.data.frame(sapply(adj.coor[,2:7],as.numeric))
adj.coor[adj.coor$site=="FLCC",5]<-adj.coor[adj.coor$site=="FLCC",6]-0.5
adj.coor[adj.coor$site=="FLCC",7]<-adj.coor[adj.coor$site=="FLCC",6]+0.5
#restrict to the actual coordinates
wod.rest<-list()
for(i in 1:nrow(adj.coor)){
	wod.temp<-as.data.frame(subset(wod.filt,
		round(Long,1) <= adj.coor$lon.h[i] & round(Long,1) >= adj.coor$lon.l[i]))
	wod.rest[[i]]<-as.data.frame(subset(wod.temp,  
		Lat <= adj.coor$lat.h[i] & Lat >= adj.coor$lat.l[i]))
}
names(wod.rest)<-paste(mar.coor$site)
averages<-lapply(wod.rest, function(x){
	avgs<-apply(x,2,mean)
	n<-nrow(x)
	return(list("avgs"=avgs, "n"=n))
})

environ.file<-NULL
for(i in 1:length(averages)){
	temp.avg<-append(averages[[i]]$avgs,averages[[i]]$n)
	environ.file<-cbind(environ.file,temp.avg)
}
colnames(environ.file)<-mar.coor$site
rownames(environ.file)<-c(names(averages[[1]]$avgs), "n")
environ.file<-environ.file[-9,]#remove pH
environ.file<-environ.file[-9,]#remove Oxygen

#write environmental data to file
write.table(environ.file, "wod_data_bayenv.txt",sep='\t',eol='\n', quote=F)

##include temperature variance
temp.var<-lapply(wod.rest,function(x){
	tempvar<-var(x$Temp)
	return(tempvar)
})
temp.var<-t(do.call("rbind",temp.var))
rownames(temp.var)<-c("tempvar")
std.tempvar<-(temp.var-mean(temp.var))/sd(temp.var)
write.table(temp.var,"../new_bayenv/wod_tempvar_data_std.txt",sep='\t',
	eol="\t\n",quote=F)


