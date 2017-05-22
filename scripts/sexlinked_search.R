source("../../gwscaR/R/gwscaR.R")

pop.list<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLLG")
pop.labs<-c("TXSP","TXCC","TXFW","TXCB","LAFW","ALST","ALFW","FLSG","FLKB",
            "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC","FLFW")
fw.list<-c("TXFW","LAFW","ALFW","FLLG")
sw.list<-c("TXSP","TXCC","TXCB","ALST","FLSG","FLKB",
           "FLFD","FLSI","FLAB","FLPB","FLHB","FLCC")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)
all.colors<-c(rep("black",2),"#2166ac","black","#2166ac","black","#2166ac",
              rep("black",8),"#2166ac")
grp.colors<-c('#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837')

vcf<-parse.vcf("stacks/fw-sw_populations/batch_2.vcf")
vcf$SNP<-paste(vcf$`#CHROM`,vcf$POS,sep=".")
scaffs<-levels(as.factor(vcf[,1]))
scaffs[1:22]<-lgs
scaff.starts<-tapply(vcf$POS,vcf$`#CHROM`,max)
scaff.starts<-data.frame(rbind(cbind(names(scaff.starts),scaff.starts)),stringsAsFactors = F)
locus.info<-colnames(vcf)[1:9]

sexlinked<-NULL
for(i in 1:length(sw.list)){
  pop.vcf<-cbind(vcf[,1:9],vcf[,grep(sw.list[i],colnames(vcf))])
  mal<-cbind(vcf[,1:9],
             pop.vcf[,grep(paste("sample_",sw.list[i],"[PN]\\w+",sep=""),
                           colnames(pop.vcf))])
  fem<-cbind(vcf[,1:9],
             pop.vcf[,grep(paste("sample_",sw.list[i],"[FD]\\w+",sep=""),
                           colnames(pop.vcf))])
  outliers<-gwsca(vcf = pop.vcf,locus.info = locus.info,
                  group1 = colnames(mal),
                  group2 = colnames(fem))
  print(max(outliers$NumAlleles))
  #sexlinked<-c(sexlinked,outliers[,"POS"])
}