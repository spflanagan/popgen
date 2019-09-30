# script to generate a figure comparing Bayenv matrices

args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied filename", call.=FALSE)
} else if (length(args)==1) {
  # current directory is 
  directory <- getwd()
} else if (length(args)>1){
  directory<-args[2]
}
if(length(grep("png",args[1])>0)){
  name<-args[1]
}else{
  name<-paste(args[1],"png",sep=".")
}


files<-list.files(pattern="out.last",path = directory,full.names = TRUE)
png(name,height=6,width=6,units="in",res=300,pointsize = 1)
par(mfrow=c(2,length(files)/2))
dat<-lapply(files,function(f){
  mat<-as.matrix(read.delim(f,skip = 1,header = FALSE))
  mat<-mat[,colSums(is.na(mat))<nrow(mat)] # remove columns with NAs
  image(mat)
  return(mat)
})
dev.off()

matnum<-sample(x = 1:length(dat),size = 1)
write.table(dat[[matnum]],paste(directory,"matrix",sep="/"),col.names = FALSE,row.names = FALSE,quote=FALSE, sep="\t")
