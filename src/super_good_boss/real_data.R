source("~/Documents/Research/Project/sc_eQTL/gmm_functions.R")
library(mvtnorm)
library(data.table)
library(susieR)
library(ggplot2)
library(XMAP)

ldsc_warp <- function(sumstat1,sumstat2,ldscore1,ldscore2,ldscorex,
                      reg_w1=1,reg_w2=1,reg_wx=1,constrain_intercept=F,int1=1,int2=1,intx=0,
                      nblocks=200,jknife=T){
  ld1_w <- 1 / sapply(ldscore1, function(x) max(x, 1))
  ld2_w <- 1 / sapply(ldscore2, function(x) max(x, 1))
  if(!constrain_intercept){
    fit_ldsc <- estimate_gc(sumstat1,sumstat2,
                            ldscore1, ldscore2, ldscorex,
                            reg_w1 = ld1_w, reg_w2 = ld2_w, reg_wx = sqrt(ld1_w * ld2_w),
                            constrain_intercept = F)
    int1 <- fit_ldsc$tau1$coefs[1]
    int2 <- fit_ldsc$tau2$coefs[1]
    intx <- fit_ldsc$theta$coefs[2]
  }
  
  
  fit_ldsc <- estimate_gc(sumstat1,sumstat2,
                          ldscore1, ldscore2, ldscorex,
                          reg_w1 = ld1_w, reg_w2 = ld2_w, reg_wx = sqrt(ld1_w * ld2_w),
                          constrain_intercept = T, int1, int2, intx)
  fit_ldsc$tau1$pval <- 2*pnorm(abs(fit_ldsc$tau1$coefs[2]/fit_ldsc$tau1$coefs_se[2]),lower.tail = F)
  fit_ldsc$tau2$pval <- 2*pnorm(abs(fit_ldsc$tau2$coefs[2]/fit_ldsc$tau2$coefs_se[2]),lower.tail = F)
  fit_ldsc$theta$pval <- 2*pnorm(abs(fit_ldsc$theta$coefs[2]/fit_ldsc$theta$coefs_se[2]),lower.tail = F)
  
  fit_ldsc
}

dat_tar <- fread("/Users/cmx/Downloads/TAR_chr2.csv")
dat_aux <- fread("/Users/cmx/Downloads/AUX_chr2.csv")
dat_tis <- fread("/Users/cmx/Downloads/Tissue_chr2.csv")

ld1 <- fread("/Users/cmx/Downloads/TAR_AUX_std_chr2_pop1.gz",data.table = F)
ld2 <- fread("/Users/cmx/Downloads/TAR_AUX_std_chr2_pop2.gz",data.table = F)
ldx <- fread("/Users/cmx/Downloads/TAR_AUX_std_chr2_te.gz",data.table = F)

pic <- 0.15
pio <- 1-pic

genes <- colnames(ld1)[-1:-4]
G <- length(genes)

g_sk <- 0
out <- data.frame()
for(g in 1: G){
  snps <- subset(dat_tar,GENE==genes[g])$RSID
  bhs <- cbind(subset(dat_tar,GENE==genes[g])$BETA,subset(dat_aux,GENE==genes[g])$BETA,subset(dat_tis,GENE==genes[g])$BETA)
  ses <- cbind(subset(dat_tar,GENE==genes[g])$SE,subset(dat_aux,GENE==genes[g])$SE,subset(dat_tis,GENE==genes[g])$SE)
  lds  <- cbind(ld1[match(snps,ld1$SNP),genes[g]],ld2[match(snps,ld2$SNP),genes[g]],ldx[match(snps,ldx$SNP),genes[g]])
  
  if(length(snps)<500){
    cat(genes[g],"skipped \n")
    g_sk <- g_sk + 1
    next
  }
  
  fit12 <- ldsc_warp(data.frame(Z = bhs[,1]/ses[,1], N = 105), data.frame(Z = bhs[,2]/ses[,2], N = 424),
                     lds[,1],lds[,2],lds[,3],constrain_intercept = T)
  fitt <- ldsc_warp(data.frame(Z = bhs[,3]/ses[,3], N = 670), data.frame(Z = bhs[,3]/ses[,3], N = 670),
                    lds[,2],lds[,2],lds[,2],constrain_intercept = T)
  fit2t <- ldsc_warp(data.frame(Z = bhs[,2]/ses[,2], N = 424), data.frame(Z = bhs[,3]/ses[,3], N = 670),
                     lds[,2],lds[,2],lds[,2],constrain_intercept = T)
  
  fit1t <- ldsc_warp(data.frame(Z = bhs[,1]/ses[,1], N = 105), data.frame(Z = bhs[,3]/ses[,3], N = 670),
                     lds[,1],lds[,2],lds[,3],constrain_intercept = T)
  
  o1 <- fit12$tau1$coefs[2]
  oc <- fit12$tau2$coefs[2]
  o1c <- fit12$theta$coefs[2]
  oco <- fit2t$theta$coefs[2]/pio
  oo <- (fitt$tau1$coefs[2] - pic^2*oc -2*pic*pio*oco)/pio^2
  oxo <- (fit1t$theta$coefs[2] - pic*o1c)/pio
  
  
  if(o1<0 | fit12$tau1$pval>0.05){
    cat(g,"-th gene",genes[g],"skipped \n")
    g_sk <- g_sk + 1
    next
  }
  
  if(oc<0 | fit12$tau2$pval>0.05){
    oc <- 1e-16
    o1c <- 1e-30*sqrt(oc*o1)
  }
  
  if(oo<0 | fitt$tau1$pval>0.05){
    oo <- 1e-16
    oco <- 1e-30*sqrt(oc*oo)
    oxo <- 1e-30*sqrt(o1*oo)
  }
  
  r1c <- o1c/sqrt(o1*oc)
  if(abs(r1c)>1){
    r1c <- 0.99 * sign(r1c)
    o1c <- r1c*sqrt(o1*oc)
  } 
  
  if(fit12$theta$pval>0.05){
    r1c <- 1e-30
    o1c <- r1c*sqrt(o1*oc)
  }
  
  rco <- oco/sqrt(oc*oo)
  if(abs(rco)>1){
    rco <- 0.99 * sign(rco)
    oco <- rco*sqrt(oc*oo)
  }
  
  if(fit2t$theta$pval>0.05){
    rco <- 1e-30
    oco <- rco*sqrt(oc*oo)
  }
  
  rxo <- oxo/sqrt(o1*oo)
  if(abs(rxo)>1){
    rxo <- 0.99 * sign(rxo)
    oxo <- rxo*sqrt(o1*oo)
  }
  
  if(fit1t$theta$pval>0.05){
    rxo <- 1e-30
    oxo <- rxo*sqrt(o1*oo)
  }
  
  Omega12o <- diag(c(o1,oc,oo))
  Omega12o[1,2] <- Omega12o[2,1] <- o1c
  Omega12o[2,3] <- Omega12o[3,2] <- oco
  Omega12o[1,3] <- Omega12o[3,1] <- oxo

  # GMM no SigO (original version)
  out12_noSigO <- gmm12(bhs,ses,c(pic,pio),lds,Omega12o)
  
  # GMM with latent cell types
  out12o <- gmm12o(bhs,ses,c(pic,pio),lds,Omega12o)
  
  # GMM with latent cell types and zero correlation between latent cell type and others
  Omega0 <- Omega12o
  Omega0[2,3] <- Omega0[3,2] <- Omega0[1,3] <- Omega0[3,1] <- 0
  out12o_noCor <- gmm12o(bhs,ses,c(pic,pio),lds,Omega0)
  
  # GMM using only target and auxiliary sc data
  Omega12 <- Omega12o[1:2,1:2]
  out0 <- gmm0(bhs[,1:2],ses[,1:2],lds,Omega12)

  
  out <- rbind(out,data.frame(SNP=snps,
                              GENE=rep(genes[g],length(snps)),
                              TAR_BETA=bhs[,1],
                              TAR_SE=ses[,1],
                              TAR_P=2*pnorm(abs(bhs[,1]/ses[,1]),lower.tail = F),
                              AUX_BETA=bhs[,2],
                              AUX_SE=ses[,2],
                              AUX_P=2*pnorm(abs(bhs[,2]/ses[,2]),lower.tail = F),
                              TISSUE_BETA=bhs[,3],
                              TISSUE_SE=ses[,3],
                              TISSUE_P=2*pnorm(abs(bhs[,3]/ses[,3]),lower.tail = F),
                              GMM0_BETA=out0[,1],
                              GMM0_SE=out0[,2],
                              GMM0_P=2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F),
                              GMM12O_BETA=out12o[,1],
                              GMM12O_SE=out12o[,2],
                              GMM12O_P=2*pnorm(abs(out12o[,1]/out12o[,2]),lower.tail = F),
                              GMM12O_noCor_BETA=out12o_noCor[,1],
                              GMM12O_noCor_SE=out12o_noCor[,2],
                              GMM12O_noCor_P=2*pnorm(abs(out12o_noCor[,1]/out12o_noCor[,2]),lower.tail = F),
                              GMM_noSigO_BETA=out12_noSigO[,1],
                              GMM_noSigO_SE=out12_noSigO[,2],
                              GMM_noSigO_P=2*pnorm(abs(out12_noSigO[,1]/out12_noSigO[,2]),lower.tail = F)
  ))
  
  cat(g,"/",G," genes finished.\n")
  # par(mfrow=c(2,3))
  # plot(-log10(2*pnorm(abs(bhs[,1]/ses[,1]),lower.tail = F)),main="target original")
  # plot(-log10(2*pnorm(abs(bhs[,2]/ses[,2]),lower.tail = F)),main="auxiliary original")
  # plot(-log10(2*pnorm(abs(bhs[,3]/ses[,3]),lower.tail = F)),main="tissue original")
  # plot(-log10(2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F)),main="target+aux sc")
  # plot(-log10(2*pnorm(abs(out12o[,1]/out12o[,2]),lower.tail = F)),main="target+aux+tissue")
  # 
  # 
  # par(mfrow=c(2,4))
  # qq(2*pnorm(abs(bhs[,1]/ses[,1]),lower.tail = F),main="original")
  # qq(2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F),main="gmm0")
  # qq(2*pnorm(abs(out12o[,1]/out12o[,2]),lower.tail = F),main="gmm12o")
  # qq(2*pnorm(abs(out12[,1]/out12[,2]),lower.tail = F),main="gmm12")
  # plot(-log10(2*pnorm(abs(bhs[,1]/ses[,1]),lower.tail = F)),-log10(2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F)))
  # abline(a=0,b=1)
  # plot(-log10(2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F)),-log10(2*pnorm(abs(out12o[,1]/out12o[,2]),lower.tail = F)))
  # abline(a=0,b=1)
  # plot(-log10(2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F)),-log10(2*pnorm(abs(out120[,1]/out120[,2]),lower.tail = F)))
  # abline(a=0,b=1)
  # plot(-log10(2*pnorm(abs(out0[,1]/out0[,2]),lower.tail = F)),-log10(2*pnorm(abs(out12[,1]/out12[,2]),lower.tail = F)))
  # abline(a=0,b=1)
  # par(mfrow=c(1,1))
}

pval_gene <- sapply(out[,c("TAR_P","AUX_P","TISSUE_P","GMM0_P","GMM12O_P","GMM12O_noCor_P","GMM_noSigO_P")],function(x) tapply(x,out$GENE,min))
colSums(pval_gene<1e-5,na.rm = T)
