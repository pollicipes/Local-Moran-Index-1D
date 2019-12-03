#!/usr/bin/Rscript
####!/apps/R/2.14.2/bin/Rscript
args<-commandArgs(TRUE)
res <- as.character(args[1])
snp <- as.character(args[2]);
setwd(res)
loc <- read.table(paste(snp,".1KG.sumstats.ld",sep=''))
#loc=read.table("/home/jrodriguez/scratch/new_pancreas_HiC/snp_files/rs9908467.1KG.sumstats.ld")
colnames(loc) <- c("chr","position","snp","beta","p","maf","ld")

### Select Beta for our seed snp;
Bi=loc[which(loc$snp == snp),]$beta
### REMOVE ROWS THAT CONTAIN INFINITE VALUES
loci <- loc[!is.infinite(loc$beta),]
### REMOVE OUR SEED SNP FROM THE DATA FRAME;
a=loci[-which(loc$snp == snp),];
# Sum of all the r2 with our seed SNP
Er.ij <- sum(a$ld)
# Vector of R2
r2.ij <- a$ld
# Weights' vector for all the other SNPs (weighted by the sum of all R2)
W <- r2.ij/Er.ij
Bj <- a$beta
####
# Each of the weights times the beta of the j-SNP 
W.betaj <- W*Bj
EW.betaj <- sum(W.betaj)
# Get local Moran I
localMI <- EW.betaj*Bi
cat(localMI);
################################################

### Resampling-pvalue betas

# 
# moran <- function (W,Bi){
#     Bji <- sample(a$beta,length(a$beta));
#     W.betaj <- W*Bji;
#     EW.betaj <- sum(W.betaj);
#     lmsi <- EW.betaj*Bi;
# #    cat(localMI);
#     return(lmsi);
# }
# 
# ppt <- function(localMI,N=100000,...) {
#   stat <- localMI
#   sim <- replicate(N, moran(W,Bi))  
#   p.value <- mean((all <- c(stat, sim)) >= stat)
# #  hist(sim, sub=paste("p =", round(p.value, 4)), xlim=range(all), main="Null Distribution of Moran's I", xlab="I")
# #  abline(v = stat, col="#903030", lty=3, lwd=2)
#   return(p.value)  
# }
# 
# pv <- ppt(localMI)
# cat(localMI,"\t",pv);
