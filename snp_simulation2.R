## /mnt/software/R-3.2.4/bin/R
## /mnt/software/R-3.2.4/bin/Rscript

## options(download.file.method = "wget")
## options(download.file.method = "curl")

## install.packages("plyr")
## install.packages("echoseq-master",repos=NULL, type="source") 
## install.packages("R.matlab", repos="http://cran.cnr.berkeley.edu")
## library("devtools")
## install_github("hruffieux/echoseq")

## install.packages("devtools", repos="http://cran.cnr.berkeley.edu") 
library(echoseq)
library("R.matlab")


## user_seed <- 123; set.seed(user_seed)

Args<-commandArgs()
n <- as.numeric(Args[6]); p <- as.numeric(Args[7])  # arguments start from the 6 and 7

cor_type <- "autocorrelated"; 
set.seed(13)
vec_rho <- runif(1000, min = 0.25, max = 0.95)
set.seed(38)
vec_maf <- runif(p, min = 0.05, max = 0.5)

##vec_rho <- as.numeric(Args[8])
##vec_maf <- as.numeric(Args[9])
list_snps <- generate_snps(n, p, cor_type, vec_rho, vec_maf, n_cpus = 4)

## writeMat("C:/Users/Dell/Desktop/aa.mat", A=list_snps$snps)
writeMat("/DATA/249/xli/snp2.mat", snp=list_snps$snps)

