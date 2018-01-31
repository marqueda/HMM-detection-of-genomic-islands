#! /usr/bin/Rscript

#########################################
## Detection of genomic islands        ##
## from outlier probabilities [0..1]   ##
## using a Hidden Markov model with    ##
## 2 normally distributed states       ##
#########################################

# based on Hofer et al. 2012 BMC Genomics
# adapted by David Alexander Marques
# please cite both Hofer et al. 2012 BMC Genomics
# and Marques et al. 2016 PLoS Genetics

# Changelog
# 2018-01-08: output filenames adapted

# USAGE
# HMM_outlierprob_2norm.R filename.suffix ncores [nstart]
# filename.suffix is a file with (outlier) probabilities between 0 and 1 (quantile/q-values from outlier analysis) and a header 'x'
# ncores is the number of cores allowed for parallel computation
# nstart (optional) is the number of random starting parameters in the parameter estimation [default=1000]

# REQUIRED
# R-libraries HiddenMarkov, foreach, doParallel, PLIS

# OUTPUT
# - filename_inputdata_histogram.pdf : Histogram of the z-transformed q-values (input data), should look approximately normally distributed and continuous, otherwise consider data filtering (minor allele frequency filter) or using an other appraoch.
# - filename_2state_HMMparest.txt : Text file with Baum-Welch algorithm optimized parameters from nstart [default=1000] random starting values.
# - filename_2state_parest_plot.pdf : Scatterplot of the maximized parameter estimates against their likelihood
# - filename_zscore_vs_2statedistributions.pdf : Histogram of the z-transformed q-values with the distribution and parameters of the three states.
# - filename_2state_bestparameters.txt : Text file with the best HMM parameters used for the state reconstruction
# - filename_2state_HMMstates.txt : Text file with the HMM states from the state reconstruction using the best parameters
# - filename_2state_PLISsignificantsites.txt : Text file with multiple testing corrected significant sites coded as 0 and 1 (significant SNPs after multiple testing)

##################
# Load libraries #
##################

library("HiddenMarkov")
library("foreach")
library("doParallel")
library("PLIS")

##################################
# Read input parameters and data #
##################################

# Reads input parameters
args<-commandArgs(trailingOnly=T)

# Checks if enough input parameters were given
if(length(args)<2){
  stop("Usage: HMM_outlierprob_2norm.R filename.suffix ncores [nstart]\n\t\t- filename.suffix is a file with (outlier) probabilities between 0 and 1 (quantile/q-values from outlier analysis) and a header 'x'\n\t\t- ncores is the number of cores allowed for parallel computation\n\t\t- nstart (optional) is the number of random starting parameters in the parameter estimation [default=1000]") 
}

# Defines prefix for output files
elements<-strsplit(args[1],split="\\.")
base<-paste(elements[[1]][-length(elements[[1]])],collapse=".")
rm(elements)

# Loads input data
data<-read.table(file=args[1], header=T, quote="\"")

########################
# Z-transform q-values #
########################

qval<-data$x
qval[qval==1]<-0.999999
qval[qval==0]<-0.000001
zscore<-qnorm(qval,lower.tail=F)

# Prints distribution of z-transformed q-values
pdf(paste(base,"_inputdata_histogram.pdf",sep=""))
  hist(zscore,breaks=100,col=1,xlab="zscore (z-transformed q-value)",main="HMM input data distribution",freq=F)
dev.off()

#############################
## PARAMETER ESTIMATION BY ##
## BAUM-WELCH ALGORITHM    ##
#############################

# Parallelizes parameter estimation from args[3] random starting parameters

# Defines the number of cores to be used
registerDoParallel(cores=args[2])

# Default number of random starting parameters
if(is.na(args[3])){args[3]=1000}

# Runs parameter estimation in parallel on args[2] cores
parout<-foreach(i=1:as.integer(args[3]),.packages="HiddenMarkov",.combine='c') %dopar% {
  print(i);set.seed(i)
  # Samples random initial parameters for the transition matrix (trans),
  #   marginal/initial probabilities (init), for the state means (means)
  #   and state standard deviations (sdevs)
  prob<-runif(3,min=0,max=1)
  mean2<-runif(1,min=quantile(zscore,c(0.5)),max=quantile(zscore,c(0.8)))
  mean1<-runif(1,min=quantile(zscore,c(0.2)),max=quantile(zscore,c(0.5)))
  sdevs<-runif(2,min=0.5,max=2)
  trans<-matrix(c(prob[1],1-prob[1],
                  1-prob[2],prob[2]),byrow=T,nrow=2)
  init<-c(prob[3],1-prob[3])
  means<-c(mean1,mean2)

  # Builds Hidden Markov Model with random initial parameters
  myhmm<-dthmm(zscore,trans,init,"norm",list(mean=means,sd=sdevs),discrete=F)

  # Optimizes parameters with Baum-Welch algorithm, with 3 additional runs to find maximal estimates
  # Baum-Welch configuration
  a<-bwcontrol(maxiter=1000,tol=1e-07,prt=F,posdiff=F)
  bwhmm2<-try(BaumWelch(myhmm,control=a),silent=T)
  bwhmm1<-try(BaumWelch(bwhmm2,control=a),silent=T)
  bwhmm<-try(BaumWelch(bwhmm1,control=a),silent=T)
  
  # Outputs parameters
  if(length(bwhmm)>1){
    c(bwhmm$Pi,bwhmm$delta,bwhmm$pm$mean,bwhmm$pm$sd,bwhmm$LL,bwhmm$iter)
  }else{
    rep(NA,20)
  }
}

hmmpar<-as.data.frame(matrix(parout,ncol=12,byrow=T))
names(hmmpar)<-c("a11","a21","a12","a22","i1","i2","m1","m2","sd1","sd2","LL","iter")

# Writes all parameter estimates to a single file
write.table(hmmpar,file=paste(base,"_2state_HMMparest.txt",sep=""))
hmmpar<-read.table(paste(base,"_2state_HMMparest.txt",sep=""),stringsAsFactors=F)

########################################
## EVALUATION OF PARAMETER ESTIMATION ##
########################################

# Plots parameter estimates against their likelihood
pdf(file=paste(base,"_2state_parest_plot.pdf",sep=""),height=16,width=10)
par(mfrow=c(5,2))
attach(hmmpar)
for(k in colnames(hmmpar)[-(length(colnames(hmmpar))-c(0,1))]){
  plot(hmmpar$LL,get(k),xlab="Likelihood",ylab=k,pch=20,col="#88888888",xlim=c(min(hmmpar$LL,na.rm=T),max(hmmpar$LL,na.rm=T)))
  abline(h=with(hmmpar,get(k)[LL==max(LL,na.rm=T)]),col="red")
}
dev.off()

# Stores the best initial parameters (the one with the highest likelihood)
bestpar<-head(hmmpar[with(hmmpar,order(-hmmpar$LL)),],1)

# Runs another parameter optimization round to ensure maximum likelihood was reached
{
  # Builds Hidden Markov Model with best parameters
  myhmm<-dthmm(zscore,matrix(c(bestpar$a11,bestpar$a12,
                            bestpar$a21,bestpar$a22),byrow=T,nrow=2),
                c(bestpar$i1,bestpar$i2),"norm",
                list(mean=c(bestpar$m1,bestpar$m2),sd=c(bestpar$sd1,bestpar$sd2)),discrete=F)
  
  # Runs Baum-Welch algorithm 3 times in a row to maximize parameter estimates
  # Baum-Welch configuration
  a<-bwcontrol(maxiter=1000,tol=1e-07,prt=F,posdiff=F)
  bwhmm2<-try(BaumWelch(myhmm,control=a),silent=T)
  bwhmm1<-try(BaumWelch(bwhmm2,control=a),silent=T)
  bwhmm<-try(BaumWelch(bwhmm1,control=a),silent=T)  
}

# Saves new best parameters to file
bestpar<-c(bwhmm$Pi,bwhmm$delta,bwhmm$pm$mean,bwhmm$pm$sd,bwhmm$LL,bwhmm$iter)
names(bestpar)<-c("a11","a21","a12","a22","i1","i2","m1","m2","sd1","sd2","LL","iter")
write.table(bestpar,file=paste(base,"_2state_bestparameters.txt",sep=""))

# Plots the data (z-scores) against the inferred state distributions
pdf(file=paste(base,"_zscore_vs_2statedistributions.pdf",sep=""))
par(mfrow=c(1,1))
hist(zscore,breaks=50,col="grey",main="Data vs. Emission & Transition Probabilities",border=F,xlab="z-transformed outlier probability",freq = F)
for(i in 1:length(bwhmm$pm$mean)){
  points(seq(par("usr")[1],par("usr")[2],length.out=1000),dnorm(seq(par("usr")[1],par("usr")[2],length.out=1000),mean=bwhmm$pm$mean[i],sd=bwhmm$pm$sd[i]),type="l",col=c(1,2)[i],lwd=2)
}
legend("topright",legend=c(paste("state1"," mean=",round(bwhmm$pm$mean[1],2)," sd=",round(bwhmm$pm$sd[1],2)),
                           paste("state2"," mean=",round(bwhmm$pm$mean[2],2)," sd=",round(bwhmm$pm$sd[2],2)),
                           paste("1>1: ",round(bwhmm$Pi[1,1],2),"; 1>2: ",round(bwhmm$Pi[1,2],2)),
                           paste("2>1: ",round(bwhmm$Pi[2,1],2),"; 2>2: ",round(bwhmm$Pi[2,2],2))),
       fill=c(1,2,NA,NA),border=c(1,1,0,0),bty="n")
dev.off()

##############################
## RECONSTRUCTION OF STATES ##
##############################

# Reconstruction of states with the Viterbi algorithm
states<-Viterbi(bwhmm)

# Output state reconstruction results
write.table(states,file=paste(base,"_2state_HMMstates.txt",sep=""))

#####################################
## CORRECTION FOR MULTIPLE TESTING ##
#####################################

# Correct for multiple testing using the PLIS method as described in Hofer et al. 2012 BMC Genomics
eml <- em.hmm(zscore, L=1)
plis_res <- plis(eml$LIS,fdr=0.001,adjust=F)

# Output information on significant sites after multiple testing coded as 0 / 1
write.table(plis_res$States,file=paste(base,"_2state_PLISsignificantsites.txt",sep=""))
