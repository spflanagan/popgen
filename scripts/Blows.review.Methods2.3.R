MCMCsamp <- 1000 
#number of MCMC samples
n <- 8 
#number of traits 
m <- 6 
#number of matrices to compare
r <- 3
#number of random effects specified in the model. In our analyses these were animal, vial and residual effects.
traitnames <- c("lc2","lc3","lc4","lc5","lc6","lc7","lc8","lc9") 
#trait names
Gnames <- c("b1","b2","c1","c2","m1","m2")
#matrix labels 


MCMCarray <- array(,c(MCMCsamp,(n^2)*r,m)) 
#empty array 
MCMCarray[,,1] <- as.matrix(read.csv(file="C:/Users/sflanagan/Downloads/Aguirre_et_al_ data_ped_output/VC_b1_3MCMC.csv",header =T))
#G1 stored as the 1st element of dim[3] 
MCMCarray[,,2] <- as.matrix(read.csv(file="C:/Users/sflanagan/Downloads/Aguirre_et_al_ data_ped_output/VC_b2_3MCMC.csv",header =T))
#G2 stored as the 2nd element of dim[3]
MCMCarray[,,3] <- as.matrix(read.csv(file="C:/Users/sflanagan/Downloads/Aguirre_et_al_ data_ped_output/VC_c1_3MCMC.csv",header =T))
#G3 stored as the 3rd element of dim[3]
MCMCarray[,,4] <- as.matrix(read.csv(file="C:/Users/sflanagan/Downloads/Aguirre_et_al_ data_ped_output/VC_c2_3MCMC.csv",header =T)) 
#G4 stored as the 4th element of dim[3]
MCMCarray[,,5] <- as.matrix(read.csv(file="C:/Users/sflanagan/Downloads/Aguirre_et_al_ data_ped_output/VC_m1_3MCMC.csv",header =T)) 
#G5 stored as the 5th element of dim[3]
MCMCarray[,,6] <- as.matrix(read.csv(file="C:/Users/sflanagan/Downloads/Aguirre_et_al_ data_ped_output/VC_m2_3MCMC.csv",header =T)) 
#G6 stored as the 6th element of dim[3]

Garray <- array(,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)
Parray <- array(,c(n,n,m,MCMCsamp))
dimnames(Parray) <- list(traitnames,traitnames,Gnames)
    for (i in 1:m){
      for (j in 1:MCMCsamp){
        G <- matrix(MCMCarray[j,1:(n^2),i],ncol= n)
        CE <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
        R <- matrix(MCMCarray[j,(((n^2)*2)+1):((n^2)*3),i],ncol= n)
        Garray[,,i,j] <- G
        Parray[,,i,j] <- G + CE + R
      }
    }

#Method 2: Krzanowski's subspace
#START
kr.subspace <- function(Gs, vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]] 
  if(length(vec) != m){stop("vec must have length = m")}
  h <- function (g, v){
    AA <- array(, c(n, n, m))  
    for (k in 1:m){
      g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
      AA[,,k] <- g.vec %*% t(g.vec)
    }
    H <- apply(AA, 1:2, sum)
    list(H = H, AA = AA)
  }
  #internal function to calculate AA and H
  MCMC.H <- array(, c(n, n, MCMCsamp))
  dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
  MCMC.AA <- array(, c(n, n, m, MCMCsamp))
  dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
   for (i in 1:MCMCsamp){
     kr <- h(Gs[,,,i], v = vec)
     MCMC.H[,,i] <- kr$H
     MCMC.AA[,,,i] <- kr$AA
   }	
   #calculate AA and H for the ith MCMC sample of the G array		
   avH <- apply(MCMC.H, 1:2, mean)
   rownames(avH) <- dimnames(Gs)[[1]]
   colnames(avH) <- dimnames(Gs)[[1]]
   #calculate the posterior mean H
   avAA <- apply(MCMC.AA, 1:3, mean)
   dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
   #calculate the posterior mean AA
   avH.vec <- eigen(avH)$vectors
   #eigenanalysis of posterior mean H	
   proj<- function(a, b) t(b) %*% a %*% b
   #internal function to do projection
   avH.theta <- matrix(, n, m)
    for (i in 1:n){
     for (i in 1:n){
      avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
     }
    }
    #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
    MCMC.H.val <- matrix(, MCMCsamp, n)
    colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
     for (i in 1:n){
      MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
     }
     #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
     MCMC.H.theta <- array(, c(n, m, MCMCsamp))
     rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
     colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
      for(i in 1:n){
       for(j in 1:MCMCsamp){
         MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
      }
     }
     #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
  list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
}
#END


#Method 3: tensor
#required packages
library(gdata);library(matrixcalc);library(MCMCglmm)

#START
covtensor <- function(Gs){
    if (dim(Gs)[[1]] != dim(Gs)[[2]]){
      stop("G array must be of order n x n x m x MCMCsamp")
    }
    if (is.na(dim(Gs)[4])) {
      stop("There are no MCMCsamples")
    }
    neigten <- n*(n+1)/2 
    #Number of eigentensors
    MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
    dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
      for (k in 1:MCMCsamp){
        MCMCG <- Gs[,,,k] 
          MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
          #find the variances of the kth G and store them 
          MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
          #find the covariances of the kth G and store them
          MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
          #fill the upper left quadrant of the kth S
          MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
          #fill the lower right quadrant of the kth S
          MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
          #fill the upper right quadrant of the kth S
          MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
          #fill the lower left quadrant of the kthS
        }  
        av.S <- apply(MCMC.S, 1:2, mean)
        #posterior mean S
        av.S.val <- eigen(av.S)$values
        #eigenvalues of posterior mean S 
        av.S.vec <- eigen(av.S)$vectors
        #eigenvalues of posterior mean S
        eTmat <- array(, c(n, n, neigten))
        dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
          for (i in 1:neigten){
            emat <- matrix(0, n, n) 
            lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
            emat <- emat + t(emat)
            diag(emat) <- av.S.vec[1:n,i]
            eTmat[,,i] <- emat 
          }
          #construct the second-order eigentensors of posterior mean S
          eT.eigen <- array(, c(n+1, n, neigten))
            for (i in 1:neigten){
              eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
              #Eigenvalues of the ith eigentensor
              eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
              #Eigenvectors of the ith eigentensor
              eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
            }
            MCMC.S.val <- matrix(, MCMCsamp, neigten)
            colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
            for (i in 1:MCMCsamp){
              for(j in 1:neigten){
                MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
              }
            }
            #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
            av.G.coord <- array(, c(m, neigten, 1))
            dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
              for (i in 1:neigten){
                av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
              }
              #Coordinates of the jth avG for the eigentensors of posterior mean S
              MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
              dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
                for (i in 1:neigten){
                  MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
                }
                #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
        tensorsummary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
        colnames(tensorsummary) <- c("S.eigval", "eT.val", traitnames)
        rownames(tensorsummary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensorsummary = tensorsummary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}
#END


covtensor(Garray)

