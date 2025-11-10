# GMM with single cell data
# bhs: beta hat matrix of p by 2, pop1, pop2
# ld: ld matrix of p by 3, pop1, pop2, cross
# Omega: per-SNP heritability matrix 2 by 2
gmm0 <- function(bhs,ses,ld,Omega){
  
  p <- nrow(bhs)
  var_bluex <- beta_bluex <- rep(0,p)
  for(j in 1:p){
    ld1 <- ld[j,1]
    ld2 <- ld[j,2]
    ldx <- ld[j,3]
    Omegaj <- Omega*matrix(c(ld1,ldx,ldx,ld2),2,2)
    A <- diag(2)
    
    lambda <- Omegaj[,1]/Omegaj[1,1]
    invLambda <- Omegaj-outer(Omegaj[1,],Omegaj[1,])/Omegaj[1,1] + diag(ses[j,]^2)
    Lambda <- solve(invLambda)
    
    var_bluex[j] <- drop(1/(t(lambda) %*% Lambda %*% lambda))
    beta_bluex[j] <- var_bluex[j] * t(lambda) %*% Lambda %*% bhs[j,]
  }
  
  cbind(bh=beta_bluex,se=sqrt(var_bluex))
}


# Original version of GMM with out SigmaO
# bhs: beta hat matrix of p by 3, pop1, pop2, tissue
# bhses: se matrix of p by 3, pop1, pop2, tissue
# props: pic, pio
# ld: ld matrix of p by 3, pop1, pop2, cross
# Omega: per-SNP heritability matrix 3 by 3
gmm12 <- function(bhs,ses,props,ld,Omega){

  p <- nrow(bhs)
  var_bluex <- beta_bluex <- rep(0,p)
  for(j in 1:p){
    ld1 <- ld[j,1]
    ld2 <- ld[j,2]
    ldx <- ld[j,3]
    Omegaj <- Omega[1:2,1:2]*matrix(c(ld1,ldx,ldx,ld2),2,2)
    A <- rbind(diag(2),c(0,props[1]))
    
    lambda <- A %*% Omegaj[,1]/Omegaj[1,1]
    invLambda <- A %*% (Omegaj-outer(Omegaj[1,],Omegaj[1,])/Omegaj[1,1]) %*% t(A) + diag(ses[j,]^2)
    # invLambda[3,3] <- invLambda[3,3] + props[2]^2 * ld2 * Omega[3,3]
    Lambda <- solve(invLambda)
    
    var_bluex[j] <- drop(1/(t(lambda) %*% Lambda %*% lambda))
    beta_bluex[j] <- var_bluex[j] * t(lambda) %*% Lambda %*% bhs[j,]
  }
  
  cbind(bh=beta_bluex,se=sqrt(var_bluex))
}

# GMM with latent cell types
# bhs: beta hat matrix of p by 3, pop1, pop2, tissue
# bhses: se matrix of p by 3, pop1, pop2, tissue
# props: pic, pio
# ld: ld matrix of p by 3, pop1, pop2, cross
# Omega: per-SNP heritability matrix 3 by 3
gmm12o <- function(bhs,ses,props,ld,Omega){
  
  p <- nrow(bhs)
  var_bluex <- beta_bluex <- rep(0,p)
  for(j in 1:p){
    ld1 <- ld[j,1]
    ld2 <- ld[j,2]
    ldx <- ld[j,3]
    Omegaj <- Omega
    Omegaj[1:2,1:2] <- Omegaj[1:2,1:2]*matrix(c(ld1,ldx,ldx,ld2),2,2)
    Omegaj[3,3] <- Omegaj[3,3]*ld2
    Omegaj[2,3] <- Omegaj[3,2] <- Omegaj[3,2]*ld2
    Omegaj[1,3] <- Omegaj[3,1] <- Omegaj[3,1]*ldx
    A <- cbind(rbind(diag(2),c(0,props[1])),c(0,0,props[2]))
    
    lambda <- A %*% Omegaj[,1]/Omegaj[1,1]
    invLambda <- A %*% (Omegaj-outer(Omegaj[1,],Omegaj[1,])/Omegaj[1,1]) %*% t(A) + diag(ses[j,]^2)
    Lambda <- solve(invLambda)
    
    var_bluex[j] <- drop(1/(t(lambda) %*% Lambda %*% lambda))
    beta_bluex[j] <- var_bluex[j] * t(lambda) %*% Lambda %*% bhs[j,]
  }
  
  cbind(bh=beta_bluex,se=sqrt(var_bluex))
}