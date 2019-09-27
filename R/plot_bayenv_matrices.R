# script to generate a figure comparing Bayenv matrices

args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied filename", call.=FALSE)
} else if (length(args)==1) {
  # current directory is 
  directory <- getwd()
}
if(length(grep("png",args[1])>0)){
  name<-args[1]
}else{
  name<-paste(args[1],"png",sep=".")
}


files<-list.files(pattern="out.last",path = directory)
png(name,height=8,width=6,units="in",res=300,pointsize = 1)
par(mfrow=c(2,length(files)/2))
dat<-lapply(files,function(f){
  mat<-read.delim(f,skip = 1)
  image(mat)
})
dev.off()
