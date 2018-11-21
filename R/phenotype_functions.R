#Author: Sarah P. Flanagan
#Last updated: 9 June 2016
#Date: 9 June 2016
#Purpose: Functions for analyzing Pst and P-matrices

#********************************PST-FST*********************************#
pairwise.pst<-function(dat, pop.order){
	#first column must be pop id/grouping factor
	library(nlme)
	dat.split<-split(dat, factor(dat[,1]))
	dat.var<-as.data.frame(setNames(
		replicate(length(pop.order),numeric(0), simplify = F), pop.order))
	for(i in 1:(length(pop.order)-1)){
	  for(j in (i+1):length(pop.order)){
		temp.data<-rbind(as.data.frame(dat.split[[pop.order[i]]]),
			as.data.frame(dat.split[[pop.order[j]]]))
		colnames(temp.data)<-c("PopID","Var")
		temp.data$PopID<-factor(temp.data$PopID)
		anv <- lme(fixed=Var ~ 1, random=~1|PopID,data=temp.data)
		varcomp <- VarCorr(anv)
		v.btwn<- as.numeric(varcomp[1])
		v.wthn <- as.numeric(varcomp[2])
		pst <- v.btwn/(v.btwn+2*v.wthn)
		dat.var[pop.order[i],pop.order[j]]<-pst
	  }
	}
	dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
	rownames(dat.var)<-colnames(dat.var)
	return(dat.var)
}

all.traits.pst.mantel<-function(trait.df,comp.df,id.index){
	results.mantel<-data.frame()
	for(i in 3:ncol(trait.df)){
		res<-mantel.rtest(
			as.dist(t(pairwise.pst(trait.df[,c(id.index,i)],pop.order))),
			as.dist(t(comp.df)), nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	rownames(results.mantel)<-colnames(trait.df)[3:ncol(trait.df)]
	colnames(results.mantel)<-c("Obs","P")
	return(results.mantel)
}

fst.pst.byloc<-function(ped.file,trait.df,pop.order,trait.ind){
	results.list<-list()
	for(j in 3:ncol(trait.df)){
	results.mantel<-data.frame()
	for(i in seq(7,ncol(ped.file),2)){
		res<-mantel.rtest(
			as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
			as.dist(t(pairwise.pst(trait.df[,c(trait.ind,j)],pop.order))),
			 nrepet=9999)
		results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
	}
	results.mantel<-as.data.frame(results.mantel)
	colnames(results.mantel)<-c("Obs","P")
	results.list<-append(results.list,data.frame(results.mantel))
	}
	#names(results.list)<-colnames(trait.df)[3:ncol(trait.df)]
	return(results.list)
}

create.extract.sh<-function(df, output.dir){
	#writes a script to extract a genome region from a chrom,start and end bp
	#the df needs chrom, start, end as columns
	commands<-paste(
		"../SCA/programs/extract_sequence_part/extract_sequence_part",
		" -f ",output.dir,df[,1],".fasta -s ",
		df[,2], " -e ", df[,3],sep="")
	return(commands)
}

#********************************P-MATRIX*********************************#

calc.pmat<-function(phen.dat.list, dim1, dim2){
	#BETWEEN POPS
	pmat<-lapply(phen.dat.list, function(x){
		cov(x[,dim1:dim2])
	})
	return(pmat)
}

###############################BLOWS METHOD#################################
#A=first three eigenvectors of P for first species
#B=first three eigenvectors of P for second species
#S=t(A)Bt(B)A
#sim=sum(eigen(S))/3

#####FUNCTIONS#####
calc.sim<-function(P.one, P.two){
	A<-eigen(P.one)$vectors[,1:3]
	B<-eigen(P.two)$vectors[,1:3]
	S<-t(A)%*%B%*%t(B)%*%A
	sim<-sum(eigen(S)$values)/3
	return(sim)
}

generate.sim.mat<-function(list.pmatrices){
	sim.mat<-matrix(nrow=length(list.pmatrices), ncol=length(list.pmatrices))
	for(i in 1:length(list.pmatrices)){
		for(j in 1:length(list.pmatrices)){
			sim.mat[i,j]<-
				calc.sim(list.pmatrices[[i]], list.pmatrices[[j]])
		}
	}
	colnames(sim.mat)<-names(list.pmatrices)
	rownames(sim.mat)<-names(list.pmatrices)
	return(sim.mat)
}

#vector correlation of pmax
find.pmax<-function(pmatrix){
	pmax<-eigen(pmatrix)$vectors[,which.max(eigen(pmatrix)$values)]
	pmax<-pmax/sqrt(sum(pmax*pmax))
	return(as.vector(pmax))
}
vector.correlations<-function(list.pmatrices){
	pmax.mat<-matrix(nrow=length(list.pmatrices), ncol=length(list.pmatrices))
	for(i in 1:(length(list.pmatrices)-1)){
		for(j in (i+1):length(list.pmatrices)){
			p1<-find.pmax(list.pmatrices[[i]])
			p2<-find.pmax(list.pmatrices[[j]])
			pmax.mat[i,j]<-p1%*%p2
		}
	}
	colnames(pmax.mat)<-names(list.pmatrices)
	rownames(pmax.mat)<-names(list.pmatrices)
	return(pmax.mat)

}

calc.h<-function(P.list){
	H<-0
	for(i in 1:length(P.list)){
		A<-eigen(P.list[[i]])$vectors[,1:3]
		A.out<-A%*%t(A)
		H<-H+A.out
	}
	colnames(H)<-colnames(P.list[[1]])
	rownames(H)<-colnames(P.list[[1]])
	return(H)
}

pop.h.angle<-function(Pmatrix,H,h.eig){
	b<-eigen(H)$vectors[,h.eig]
	A<-eigen(Pmatrix)$vectors[,1:3]
	trans<-sqrt(t(b)%*%A%*%t(A)%*%b)
	delta<-acos(trans^0.5)
}



covtensor<-function(matrix.list){
#Adapted from Aguirre et al. 2013 supplemental material
	n.traits<-dim(matrix.list[[1]])[1]
	n.pops<-length(matrix.list)
	trait.names<-colnames(matrix.list[[1]])
	pop.names<-names(matrix.list)
	n.eigten<-n.traits*(n.traits+1)/2
	mat.var<-matrix(nrow=n.pops,ncol=n.traits)
	mat.cov<-matrix(nrow=n.pops,ncol=length(lowerTriangle(matrix.list[[1]])))
	for(i in 1:n.pops){
		mat.var[i,]<-diag(matrix.list[[i]])
		mat.cov[i,]<-lowerTriangle(matrix.list[[i]])
	}
	#create the S matrix
	s.mat<-matrix(nrow=n.eigten,ncol=n.eigten)
	colnames(s.mat)<-paste("E", 1:n.eigten, sep="")
	rownames(s.mat)<-paste("E", 1:n.eigten, sep="")
	#fill in upper left quadrant
	s.mat[1:n.traits,1:n.traits]<-cov(mat.var,mat.var)
	#fill in lower right quadrant
	s.mat[(n.traits+1):n.eigten,(n.traits+1):n.eigten]<-2*cov(mat.cov,mat.cov)
	#fill in upper right quadrant
	s.mat[1:n.traits,(n.traits+1):n.eigten]<-sqrt(2)*cov(mat.var, mat.cov)
	#fill in lower left quadrant
	s.mat[(n.traits+1):n.eigten,1:n.traits]<-sqrt(2)*cov(mat.cov, mat.var)
	
	#eigenvalues and vectors of S
	s.val<-eigen(s.mat)$values
	s.vec<-eigen(s.mat)$vectors

	#make the eigentensors matrix
	et.mat<-array(, c(n.traits, n.traits, n.eigten))
	dimnames(et.mat) <- list(trait.names, trait.names, 
		paste("E", 1:n.eigten, sep=""))  
	for(i in 1:n.eigten){
		e.mat <- matrix(0, n.traits, n.traits) 
		lowerTriangle(e.mat) <- 1/sqrt(2)*s.vec[(n.traits+1):n.eigten,i]
		e.mat <- e.mat + t(e.mat)
		diag(e.mat) <- s.vec[1:n.traits,i]
		et.mat[,,i] <- e.mat 
	}
	
	#eigenvectors of eigentensors
	et.eigen<-array(,c((n.traits+1),n.traits, n.eigten))
	for(i in 1:n.eigten){
		#eigenvalues of ith eigentensor
		et.eigen[1,,i]<-t(eigen(et.mat[,,i])$values)
		#eigenvectors of the ith eigentensor
		et.eigen[2:(n.traits+1),,i]<-eigen(et.mat[,,i])$vectors
		et.eigen [,,i]<-et.eigen[,order(abs(et.eigen[1,,i]),
			decreasing = T), i]
	}

	#project eigenvectors of s onto s to determine alpha
	s.alpha<-matrix(ncol=n.eigten)
	for(i in 1:n.eigten){
		s.alpha[,i]<-t(s.vec[,i]) %*% s.mat %*% s.vec[,i]
	}
	#distribution of phenotypic variance for eigenvectors of S
	p.coord<-matrix(nrow=n.pops, ncol=n.eigten)
	for(i in 1:n.eigten){
		p.coord[,i]<-unlist(lapply(matrix.list, frobenius.prod, y = et.mat[,,i]))
	}
	rownames(p.coord)<-pop.names
	colnames(p.coord)<-paste("E", 1:n.eigten, sep="")
	tensor.summary <- data.frame(rep(s.val,each=n.traits), 
		t(data.frame(et.eigen)))
	colnames(tensor.summary) <- c("S.eigval", "eT.val", trait.names)
	rownames(tensor.summary)<- paste(
		paste("e", rep(1:n.eigten, each=n.traits), sep=""), 
		rep(1:n.traits,n.eigten), sep=".")
	#maximum number of nonzero eigentensors
	nonzero<-min(n.traits*(n.traits+1)/2,n.pops-1)
	

	return(list(tensor.summary = tensor.summary, s.mat = s.mat, 
		et.mat = et.mat, p.coord = p.coord, s.val = s.val, 
		s.alpha=s.alpha, nonzero=nonzero))
	#tensor.summary: sumamry of the covariance tensor for S.
		#contains eigenvalues of the tensor and 
		#the eigenvalues of the eigentensors. And the eigenvectors of 
		#eigentensors with trait loadings identified by the column names.
	#s.mat: the S matrix
	#et.mat: the eigentensors of S. Rows and columns identify elements 
		#of an eigentensor. 3rd dimension identifies the eigentensor
	#p.coord: Coordinates of the P matrix in the space of the eigentensors 
		#of S
	#s.val: eigenvalues of S
	#s.alpha: projection of the eigenvectors of S on S.
	#nonzero: The maximum number of nonzero eigentensors.
}

max.eig.val<-function(eig.vec){
	max.perc<-max(abs(eig.vec))/sum(abs(eig.vec))
	max.loc<-which.max(abs(eig.vec))
	return(list(max.loc=max.loc, max.perc=max.perc))
}

#Function to do projection
proj<- function(G, b) t(b) %*% G %*% (b)

