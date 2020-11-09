# run from fwsw_results
source("../../gwscaR/R/gwscaR.R")
source("../../gwscaR/R/gwscaR_plot.R")
source("../../gwscaR/R/gwscaR_utility.R")
source("../../gwscaR/R/gwscaR_fsts.R")
source("../../gwscaR/R/gwscaR_popgen.R")
source("../../gwscaR/R/vcf2dadi.R")
permute.gwsca<-function(vcf,map1,nperms,z=1.96, maf.cutoff = 0.05,cov.thresh=0.2){
  # calculate the actuals
  actual_fsts<-gwsca(vcf,colnames(vcf)[1:9],
                     map1[map1[,2] %in% unique(map1[,2])[1],1],
                     map1[map1[,2] %in% unique(map1[,2])[2],1],
                     maf.cutoff=maf.cutoff,prop.ind.thresh=cov.thresh)
  # do the permutations
  perm_fsts<-lapply(1:nperms,function(i,vcf,map1){
    perm_map<-map1
    perm_map[,2]<-perm_map[,2][permute::shuffle(perm_map[,2])]
    perm_dat<-gwsca(vcf,colnames(vcf)[1:9],
                    perm_map[perm_map[,2] %in% unique(perm_map[,2])[1],1],
                    perm_map[perm_map[,2] %in% unique(perm_map[,2])[2],1],
                    maf.cutoff,cov.thresh)
    
    return(perm_dat)
  },vcf=vcf,map1=map1)
  
  # calculate stats
  fsts<-t(do.call(rbind,lapply(perm_fsts,'[[',"Fst"))) #extract permuted fsts
  perm_fst_mu<-rowMeans(fsts)
  perm_fst_sem<-apply(fsts,1, sem)
  perm_fst_t<-(perm_fst_mu-actual_fsts$Fst)/perm_fst_sem
  perm_fst_p<-pt(perm_fst_t,df=nperms-1,lower.tail = FALSE) # one-tailed
  perm_fst_padj<-p.adjust(perm_fst_p,method = "bonferroni")
  perm_fst_in<-NULL
  for(i in 1:nrow(actual_fsts)){
    pmax<-max(fsts[i,] )
    pmin<-min(fsts[i,] )
    if(actual_fsts[i,"Fst"] > pmax | actual_fsts[i,"Fst"] < pmin ){
      perm_fst_in[i]<-1
    }else{
      perm_fst_in[i]<-0
    }
  }
  
  fst_dat<-data.frame(cbind(actual_fsts,
                            n_perms=nperms,
                            mean_perm=perm_fst_mu,
                            sem_perm=perm_fst_sem,
                            tstat=perm_fst_t,
                            pval=perm_fst_p,
                            adjP=perm_fst_padj,
                            act_in_perm=perm_fst_in))
  return(fst_dat)
}

vcf<-parse.vcf("converted_subset.vcf")
popmap<-data.frame(inds=colnames(vcf)[10:ncol(vcf)],
                   pops=gsub("sample_(\\w{4}).*","\\1",colnames(vcf)[10:ncol(vcf)]),
                   stringsAsFactors = FALSE)
pwise_maps<-list(popmap[popmap$pops %in% c("TXFW","TXCC"),],
                 popmap[popmap$pops %in% c("FLLG","FLCC"),],
                 popmap[popmap$pops %in% c("ALFW","ALST"),],
                 popmap[popmap$pops %in% c("LAFW","ALST"),])

permuted_fsts<-lapply(pwise_maps,permute.gwsca,vcf=vcf,nperms=1000, maf.cutoff=0)
saveRDS(permuted_fsts,"permuted_fsts_22062020.RDS")
