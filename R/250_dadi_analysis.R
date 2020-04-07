parse_dadi_opt<-function(opt.file){
  model.opt<-read.delim(opt.file,sep='\t',header=TRUE,stringsAsFactors = FALSE)
  model.opt$params<-gsub("\\.\\.",",",gsub("optimized_params.(.*)","\\1",colnames(model.opt)[7]))
  colnames(model.opt)[7]<-"optimized_params"
  if(length(grep("Model",model.opt$Model))>0){  #THERE SHOULDN'T BE REPEATED HEADER ROWS
    model.opt<-model.opt[-which(model.opt$Model=="Model"),]
    warning('The model has been run multiple times and output has been concatenated') 
  }
  
  return(model.opt)
}
dadi.optimal<-function(opt.file){
  model.opt<-parse_dadi_opt(opt.file)
  opt<-model.opt[which.max(as.numeric(as.character(model.opt$log.likelihood))),]
  return(opt)
}

dadi.modelcomp<-function(path,pattern,id){
  opt.files<-list.files(path = path,pattern = pattern,full.names = TRUE)
  #find the best parameter set for each model using maximum log likelihood
  opts<-do.call(rbind,lapply(opt.files,dadi.optimal))
  #rank them by AIC
  opts$rank<-rank(opts$AIC)
  opts$id<-id
  opts<-opts[order(opts$rank),]
  return(opts)
}


dadi.nan<-function(opt.file){
  model.opt<-read.delim(opt.file,sep='\t',header=TRUE)
  if(length(grep("Model",model.opt$Model))>0){  #THERE SHOULDN'T BE REPEATED HEADER ROWS
    model.opt<-model.opt[-which(model.opt$Model=="Model"),]
    warning('The model has been run multiple times and output has been concatenated') 
  }
  nn<-nrow(model.opt[is.na(model.opt$AIC),])
  nr<-nrow(model.opt)
  return(data.frame(Model=unique(model.opt$Model),NumReps=nr,NumNan=nn))
}


dadi_warnings<-function(logname, mods=NA){
  log_txt<-scan(logname,what = "character",sep = '\n')
  # browser()
  models<-grep("Model",log_txt[grep("====",log_txt)+1],value = TRUE)
  if(length(mods)>1 & !is.na(mods[1])){ models<-models[models %in% mods] }
  model_locs<-unlist(lapply(models,grep,x=log_txt))
  names(model_locs)<-unlist(lapply(models,grep,x=log_txt,value=TRUE))
  model_locs<-model_locs[order(model_locs)]
  model_locs<-model_locs[!duplicated(model_locs)]
  model_locs<-c(model_locs,end=length(log_txt))
  
  prob_reps<-data.frame(Model=NULL,Replicate=NULL)
  for(i in 1:(length(model_locs)-1)){
    prob_rounds<-grep("Round",log_txt[model_locs[i]:model_locs[i+1]][grep("WARNING",log_txt[model_locs[i]:model_locs[i+1]])-1],value=TRUE)
    prob_rep<-gsub(".*(Round) (\\d+) (Replicate) (\\d+) of \\d+:","\\1_\\2_\\3_\\4",prob_rounds)
    if(length(prob_rep)>0){ prob_reps<-rbind(prob_reps,cbind(Model=gsub("Model (.*)","\\1",names(model_locs[i])),Replicate=prob_rep)) }
  }
  return(prob_reps)
}


plot_params<-function(dat,params,pt.col="dark grey"){ 
  plot(0:(length(params)*2),xlim=c(0,2*(length(params))),
       ylim=c(min(dat$log.likelihood),max(dat[,params])),
       bty="L",axes=FALSE,ann=FALSE,type='n')
  gwsca.vioplot(dat$log.likelihood,col="grey",add=TRUE,xpd=TRUE,at=0)
  points(jitter(rep(0,nrow(dat)),2.5),dat$log.likelihood,
         col=alpha(pt.col,0.5),pch=19)
  axis(2,las=1,pos=-0.5)
  axis(1,lwd=0,label="log likelihood",at=0)
  #mapply(function(param,counts) browser(),colnames(params),2:(ncol(params)+1))
  p<-mapply(function(param,counts){
    par(new=TRUE)
    plot(0:(length(params)*2),xlim=c(0,2*(length(params))),
         ylim=c(min(dat[,param]),max(dat[,param])),bty="L",axes=FALSE,ann=FALSE,type='n')
    gwsca.vioplot(dat[,param],col="grey",add=TRUE,xpd=TRUE,ylim=c(min(dat[,param]),max(dat[,param])),at=counts)
    points(jitter(rep(counts,nrow(dat)),2.5),dat[,param],col=alpha(pt.col,0.5),pch=19)
    axis(2,las=1,pos=counts-0.75)
    axis(1,lwd=0,label=param,at=counts)
  },params,seq(2,(length(params)*2),by=2))
  mtext(unique(dat$Model),3)
}

extract_dadi_params<-function(opt_dat){
  #need to do this separately for each model
  mod_params<-lapply(unique(opt_dat$Model), function(mod){ #browser()
    params<-do.call(rbind,mapply(function(param_est,params){ #browser()
      p<-data.frame(rbind(as.numeric(as.character(rbind(unlist(strsplit(param_est,",")))))))
      colnames(p)<-unlist(strsplit(unlist(strsplit(params,",")),"\\."))
      p<-p[order(colnames(p))]
      p$Model<-mod
      return(p)
    },param_est=opt_dat$optimized_params[opt_dat$Model %in% mod],params=opt_dat$params[opt_dat$Model %in% mod],SIMPLIFY = FALSE))
    this_dat<-cbind(opt_dat[opt_dat$Model %in% mod,],params)
    return(this_dat[-ncol(this_dat)])
  })
  return(mod_params)
}


# this function is meant to be used on a set of models each represented by its best replicate 
calc_dAIC<-function(mod){
  minAIC<-min(as.numeric(as.character(mod$AIC)))
  return(as.numeric(as.character(mod$AIC))-minAIC)
  
}

