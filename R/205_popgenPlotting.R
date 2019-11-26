

plot_multiple_LGs<-function(list_fsts,fst_name,chr_name,bp_name,lgs,plot_labs,pt_cols=NULL,plot_scaffs=NULL,addSmooth=TRUE,smoothFst="Smoothed.Fst",smoothcol="cornflowerblue",ncol=2,...){
  nrow<-length(list_fsts)/ncol
  
  # check the variables
  if(length(list_fsts)>1){
    if(length(fst_name)==1){
      fst_names<-rep(list(fst_name),length(list_fsts))
    }else if(length(list_fsts)==length(fst_name)){
      fst_names<-fst_name
    }else{
      print("ERROR: invalid fst_name")
      return(NULL)
    }
    if(length(chr_name)==1){
      chr_names<-rep(list(chr_name),length(list_fsts))
    }else if(length(list_fsts)==length(chr_name)){
      chr_names<-chr_name
    }else{
      print("ERROR: invalid chr_name")
      return(NULL)
    }
    if(length(bp_name)==1){
      bp_names<-rep(list(bp_name),length(list_fsts))
    }else if(length(list_fsts)==length(bp_name)){
      bp_names<-bp_name
    }else{
      print("ERROR: invalid bp_name")
      return(NULL)
    }
    if(length(smoothFst)==1){
      smoothFsts<-rep(list(smoothFst),length(list_fsts))
    }else if(length(list_fsts)==length(smoothFst)){
      smoothFsts<-smoothFst
    }else{
      print("ERROR: invalid smoothFst")
      return(NULL)
    }
    if(!is.null(pt_cols)){ #if it's not null, then need to check it's a list
      if(length(pt_cols)==1){
        pch_cols<-rep(list(pt_cols),length(list_fsts))
      }else if(length(list_fsts)==length(pt_cols)){
        pch_cols<-pt_cols
      }else{
        print("WARNING: invalid pt_cols, using defaults")
        pch_cols<-c("darkgrey","lightgrey")
      }
    }
  }
  if(length(plot_labs) != length(list_fsts) | is.null(plot_labs)){
    print("WARNING: invalid plot labels (plot_labs). Omitting plot labels")
    plot_labs<-rep(list(""),length(list_fsts))
  }
  
  # aggregate data
  all_chr<-data.frame(Chr=unlist(lapply(list_fsts,function(x){ as.character(x[,chr_name])})),
                      BP=unlist(lapply(list_fsts,function(x){ as.character(x[,bp_name])})),stringsAsFactors = F)
  bounds<-tapply(as.numeric(as.character(all_chr$BP)), all_chr$Chr,max)
  bounds<-data.frame(Chrom=dimnames(bounds),End=bounds)
  colnames(bounds)<-c("Chrom","End")
  if(is.null(plot_scaffs)){
    plot_scaffs<-levels(bounds$Chr)
    plot_scaffs[1:22]<-lgs
  }
  bounds<-bounds[match(plot_scaffs,bounds$Chrom),]
  
  #Plot
  if(nrow*ncol < length(list_fsts)) nrow<-nrow+1
  par(mfrow=c(nrow,ncol),mar=c(3,3,2,2),oma=c(2,2,2,2))
  fsts<-mapply(function(f, fst,bp,chr,cols, plot_lab,smF,plot_scaffs,bounds,smoothcol,...){
    
    fst<-fst.plot(f,fst.name = fst, bp.name = bp,chrom.name = chr, 
                  scaffs.to.plot=plot_scaffs, scaffold.widths = bounds,
                  pt.cols = cols,...) #pch=19,y.lim = c(0,1),pt.cex=1,axis.size = 1
    
    if(addSmooth==TRUE) points(fst$plot.pos,fst[,smF],col=smoothcol,type="l") 
    
    clip(0,max(fst$plot.pos),0,1)
    
    mtext(plot_lab,2,cex=0.75)#,line=-1)
    labs<-tapply(fst$plot.pos,fst[,chr],median)
    text(x=labs[lgs],y=-0.1,labels=lgn,xpd=TRUE)
    
    return(fst)
  },f=list_fsts,fst=fst_names,bp=bp_names,chr=chr_names,cols=pch_cols,plot_lab=plot_labs,smF=smoothFsts,MoreArgs = list(plot_scaffs=plot_scaffs,bounds=bounds, smoothcol=smoothcol,...))
  return(fsts)
}



perlg.add.lines<-function(fwsw.plot,lgs,width=NULL,lwds=4,color="cornflowerblue"){
  
  for(i in 1:length(lgs)){
    this.df<-fwsw.plot[fwsw.plot$Chr %in% lgs[i],]
    if(is.null(width)){
      width<-(nrow(this.df)*0.15)
    }
    this.smooth<-do.call("rbind",lapply(seq(1,nrow(this.df),width/5),sliding.avg,
                                        dat=data.frame(Pos=this.df$plot.pos,
                                                       Fst=this.df$Corrected.AMOVA.Fst),
                                        width=width))
    points(this.smooth,col=color,type="l",lwd=lwds)
  }
}



assign.plotpos<-function(df, plot.scaffs, bounds, df.chrom="Chrom", df.bp="BP"){
  colnames(bounds)<-c("Chrom","End")
  new.dat<-data.frame(stringsAsFactors = F)
  last.max<-0
  for(i in 1:length(plot.scaffs)){
    #pull out the data for this scaffold
    if(nrow(bounds[bounds$Chrom %in% plot.scaffs[i],])>0){ #sanity check
      chrom.dat<-df[df[[df.chrom]] %in% plot.scaffs[i],]
      if(nrow(chrom.dat)>0){
        chrom.dat$plot.pos<-as.numeric(as.character(chrom.dat[[df.bp]]))+last.max
        new.dat<-rbind(new.dat,chrom.dat)
        #last.max<-max(chrom.dat$plot.pos)+
        #               as.numeric(scaffold.widths[scaffold.widths[,1] %in% scaffs.to.plot[i],2])
      }
      last.max<-last.max+
        as.numeric(bounds[bounds$Chrom %in% plot.scaffs[i],2])
    }
  }
  #make sure everything is the correct class
  new.dat$plot.pos<-as.numeric(as.character(new.dat$plot.pos))
  return(new.dat)
}