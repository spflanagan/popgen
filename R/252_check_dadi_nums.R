pops=c("FLLG", "FLCC", "ALFW", "ALST", "LAFW", "TXFW", "TXCC")   
models=c("SI", "IM", "AM", "SC", "SI2N", "SIG", "SI2NG", "IMG", "IM2N", "IM2m", "IM2NG", "AM2N", "AMG", "AM2m", "AM2NG", "AM2N2m", "AM2mG", "AM2N2mG", "SCG", "SC2N", "SC2m", "SC2NG", "SC2N2m", "SC2mG", "SC2N2mG")

############ CREATE POP COMBOS TO RUN ############
tasks=NULL
for(i in 1:(length(pops)-1)){
  for(j in (i+1):length(pops)){
    for(mod in 1:length(models)){
      tasks<-c(tasks,paste(pops[i],pops[j],models[mod],sep="_"))
    }
  }
}

      

outs<-list.files(pattern="optimized.txt",path="dadi_results",recursive=TRUE)
outputs<-data.frame(dir=gsub("(\\w{4}_\\w{4})\\/.*","\\1",outs),
                    file=gsub("(\\w{4}_\\w{4})\\/(.*)","\\2",outs),
                    stringsAsFactors = FALSE)
outputs$pattern<-gsub("(\\w{4}_\\w{4})_\\d+.*\\.(\\w+.*)\\.optimized.txt","\\1_\\2",outputs$file)

counts<-as.data.frame(table(outputs$pattern))
mods<-as.data.frame(tasks)
all_counts<-merge(x=mods,y=counts,by.x="tasks",by.y="Var1",all.x=TRUE)
all_counts$Freq[is.na(all_counts$Freq)]<-0


sysinfo<-Sys.info()

outname<-paste0("model_counts",Sys.Date(),"_",sysinfo["nodename"],".csv")
write.csv(all_counts,outname, row.names = FALSE,quote=FALSE)
