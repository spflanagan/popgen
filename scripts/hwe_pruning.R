#Author: Sarah P. Flanagan
#Date: 4 March 2016
#Purpose: To test for hardy weinberg in loci from a plink file.
#NOTE:Did not actually use this in my analysis.

hwe.test<-function(df){ #df with two columns.
	#remove missing data
	loc.gen<-df[df[,1] != "0",]
	#calc allele freqs
	af<-table(c(as.character(loc.gen[,1]),as.character(loc.gen[,2])))/
		sum(table(c(as.character(loc.gen[,1]),as.character(loc.gen[,2]))))
	#make them genotypes
	loc.gen$AB<-paste(loc.gen[,1],loc.gen[,2],sep="/")
	gf<-table(loc.gen$AB)
	#calculate expected values
	result<-data.frame(genotype=character(),exp=numeric(),obs=numeric(),
		stringsAsFactors=F)
	for(i in 1:length(af)){ #calculate expected values
		for(j in 1:length(af)){
			names.sort<-sort(names(af)[c(i,j)])
			pq.name<-paste(names.sort[1],names.sort[2],sep="/")
			if(i == j){ 
				pq.exp<-af[i]*af[j]
			} else {	
				pq.exp<-2*af[i]*af[j]
			}
			if(length(gf[names(gf) %in% pq.name])>0){
				pq.obs<-gf[names(gf) %in% pq.name] 
				result[nrow(result)+1,]<-
					c(genotype=as.character(pq.name),exp=pq.exp,obs=pq.obs)
			} else {
				result[nrow(result)+1,]<-
					c(genotype=as.character(pq.name),exp=pq.exp,obs=rep(0,1))
			}
		}#end j
	}#end i
	#result has some duplicates because the order of alleles is not guaranteed.
	#test to see if expected = observed
	result<-result[!duplicated(result$genotype),]
	result$exp<-round(as.numeric(result$exp)*sum(as.numeric(result$obs)))
		result<-result[result$exp!="0",]
	result$chi<-((as.numeric(result$obs)-as.numeric(result$exp))^2)/
		as.numeric(result$exp)
	chi.result<-1-pchisq(sum(result$chi),length(af)-1)
	
	return(chi.result)
}

setwd("E:/ubuntushare/popgen/sw_results/stacks/populations")

map<-read.table("pruned.map",header=F)
ped<-read.table("pruned.ped",header=F)
hwe<-NULL
for(i in seq(7,length(ped),2)){
	hwe<-c(hwe,hwe.test(ped[,c(i,(i+1))]))
}
hwe<-data.frame(locus=paste(map$V1,map$V4,sep="."),hwe.result=hwe)
keep<-hwe[hwe$hwe.result > 0.001,]

