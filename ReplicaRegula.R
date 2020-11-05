##################################################
### 
### Replication-based regularization Bayesian approaches 
### to diagnose Reinke's edema by using recorded speech 
###  
### Lizbeth Naranjo (1), 
### Carlos J. Perez (2), 
### Yolanda Campos-Roca (3), 
### Mario Madruga (2).
###
### Submitted
### 
### (1) Departamento de Matematicas, Facultad de Ciencias, 
###     Universidad Nacional Autonoma de Mexico, Ciudad de 
###     Mexico, Mexico.
### (2) Departamento de Matemáticas, Facultad de Veterinaria, 
###     Universidad de Extremadura,  Cáceres, Spain.
### (3) Departamento de Tecnologias de los Computadores y 
###     las Comunicaciones, Escuela Politecnica, Universidad 
###     de Extremadura, Caceres, Spain
##################################################

##################################################
### R packages
###
### Instructions: 
### Load the R libraries
##################################################
library(rjags)   # JAGS
library(MCMCpack)   # MCMC methods
library(mnormt)   # Multivariate normal
library(coda)   # CODA
library(dclone)   # To run MCMC in
library(snow)   # several cores
library(mnormt) 
library(corrplot) 
library(pROC) 
library(xlsx) 
##################################################

##################################################
### Instructions: 
### Change the address where the data and codes are located. 
### setwd("HERE")
##################################################
setwd("~/FileDirectory/")
setwd("/Users/EstadisticaCiencias/Documents/Articulos/ReplicacionesORL/Reinke/Codes/")
getwd()
##################################################

##################################################
### Simulate mimic STdata
##################################################

##################################################
### Data, simulate covariates and response

N = 60   # subjects
J = 4   # replicates
K = 4   # covariates with replications
H = 1   # covariates 

### Parameters
BetaX = as.vector(c(1,-2, 1.6,-2.6))   
GammaZ = c(1)
Alpha = 1

### Covariates 
Z = array(NA,dim=c(N,H))
Z = rbinom(N,size=1,prob=0.5)

### Covariates and Replications
meanW = rep(2,K)
corW = matrix(0.1,K,K)    # correlations 
diag(corW) = 1   
varW = diag(0.5,K)   # variances 
covW = matrix(NA,K,K)
for(k1 in 1:K){
	for(k2 in 1:K){
		covW[k1,k2] = corW[k1,k2]*sqrt(varW[k1,k1])*sqrt(varW[k2,k2])
	}
} 


Xrep = array(NA,dim=c(N,K,J))
Wrep = array(NA,dim=c(N,K))


Wrep[,1:K] = rmnorm(N,mean=meanW,varcov=covW)   # generate W
for(i in 1:N){   # generate X
	for(j in 1:J){
		Xrep[i,,j] = Wrep[i,] + rmnorm(1,mean=rep(0,K),varcov=diag(0.1,K))
	}
}
etarep = Alpha + Wrep%*%BetaX + Z*GammaZ   	
theta = plogis(etarep)   # logistic

### Response
VV = runif(N)
Y = as.vector(ifelse(theta>VV,1,0))
plot(etarep,theta,ylim=c(0,1))
points(etarep,Y)
table(Y)
##################################################

##################################################
### Data

data1 <- list(
  Y = Y ,
  X = Xrep , 
  Z = Z ,
  N = N , 
  J = J ,
  K = K ,
  X.new = Xrep , 
  Z.new = Z ,
  N.new = N 
)
N.new = N
Y.new = Y

##################################################
### Parameters, Initials

paramN <- c(
  "Beta" ,
  "Gamma" , 
  "Alpha" ,
  "P.new"
)

initsN <- function(){	list(
  "Beta" = rnorm(K,0,0.01) , 
  "Gamma" = rnorm(1,0,0.01) ,
  "Alpha" = rnorm(1,0,0.01) , 
  "deltaX" = rep(1,K) , 
  "muW" = rnorm(K,0,0.01) , 
  "tauW" = rep(1,K) ,  
  "W" = matrix(rnorm(N*K),N,K) ,
  "W.new" = matrix(rnorm(N.new*K),N.new,K)  
)}


paramL <- c(
  "Beta" ,
  "Gamma" , 
  "Alpha" ,
  "lambda1" ,
  "P.new"
)

initsL <- function(){	list(
  "Beta" = rnorm(K,0,0.01) , 
  "Gamma" = rnorm(1,0,0.01) ,
  "Alpha" = rnorm(1,0,0.01) , 
  "deltaX" = rep(1,K) , 
  "muW" = rnorm(K,0,0.01) , 
  "tauW" = rep(1,K) ,  
  "W" = matrix(rnorm(N*K),N,K) ,
  "W.new" = matrix(rnorm(N.new*K),N.new,K)  ,
  "lambda1square" = 1 
)}


paramR <- c(
  "Beta" ,
  "Gamma" , 
  "Alpha" ,
  "lambda2" , 
  "P.new"
)

initsR <- function(){	list(
  "Beta" = rnorm(K,0,0.01) , 
  "Gamma" = rnorm(1,0,0.01) ,
  "Alpha" = rnorm(1,0,0.01) , 
  "deltaX" = rep(1,K) , 
  "muW" = rnorm(K,0,0.01) , 
  "tauW" = rep(1,K) , 
  "W" = matrix(rnorm(N*K),N,K) ,
  "W.new" = matrix(rnorm(N.new*K),N.new,K) ,  
  "lambda2" = 1 
)}


paramE <- c(
  "Beta" ,
  "Gamma" , 
  "Alpha" , 
  "lambda1" ,
  "lambda2", 
  "P.new"
)

initsE <- function(){	list(
  "Beta" = rnorm(K,0,0.01) , 
  "Gamma" = rnorm(1,0,0.01) ,
  "Alpha" = rnorm(1,0,0.01) , 
  "deltaX" = rep(1,K) , 
  "muW" = rnorm(K,0,0.01) , 
  "tauW" = rep(1,K) , 
  "W" = matrix(rnorm(N*K),N,K) ,
  "W.new" = matrix(rnorm(N.new*K),N.new,K)  , 
  "lambda1square" = 1 ,
  "lambda2" = 1
)}

##################################################
### Normal model 

fitN <- jags.model("ReplicaRegulaNorm.bug", data1, initsN, n.chains=3)

update(fitN,1000)
sampleN <- coda.samples(fitN, paramN, n.iter=1000, thin=5)

plot(sampleN)
summary(sampleN)
sumN <- summary(sampleN) 
tableN <- cbind(sumN$statistics[,1], sumN$quantiles[,"50%"],
                sumN$statistics[,2], sumN$quantiles[,c("2.5%","97.5%")])
colnames(tableN)[1:3] <- c("Mean", "Median", "SD")
write.table(tableN, file=paste0("tableNorm.txt"), quote=FALSE, sep=" ", row.names=TRUE)
round(tableN[,c(1,4,5)],5)

postN <- summary(sampleN)
PpredN <- postN$statistics[(K+H+2):(K+H+2+N.new-1),1]
YpredN <- ifelse(PpredN>0.5,1,0)
(tabN <- table(Y.new,YpredN))

TP = tabN[2,2]
TN = tabN[1,1]
FP = tabN[1,2]
FN = tabN[2,1]
TP/(TP+FN) # Sensitivity
TN/(TN+FP) # Specificity
TP/(TP+FP) # Precision
(TP+TN)/N.new # Accuracy rate
pROC:: auc(Y.new,PpredN) ### (response, predictor)

##################################################
### Ridge Model

fitR <- jags.model("ReplicaRegulaRidge.bug", data1, initsR, n.chains=3)

update(fitR,2000)
sampleR <- coda.samples(fitR, paramR, n.iter=2000, thin=10)

plot(sampleR)
summary(sampleR)
sumR <- summary(sampleR) 
tableR <- cbind(sumR$statistics[,1], sumR$quantiles[,"50%"],
                sumR$statistics[,2], sumR$quantiles[,c("2.5%","97.5%")])
colnames(tableR)[1:3] <- c("Mean", "Median", "SD")
write.table(tableR, file=paste0("tableRidge.txt"), quote=FALSE, sep=" ", row.names=TRUE)
round(tableR[,c(1,4,5)],5)

postR <- summary(sampleR)
PpredR <- postR$statistics[(K+H+2):(K+H+2+N.new-1),1]
YpredR <- ifelse(PpredR>0.5,1,0)
(tabR <- table(Y.new,YpredR))

TP = tabR[2,2]
TN = tabR[1,1]
FP = tabR[1,2]
FN = tabR[2,1]
TP/(TP+FN) # Sensitivity
TN/(TN+FP) # Specificity
TP/(TP+FP) # Precision
(TP+TN)/N.new # Accuracy rate
pROC:: auc(Y.new,PpredR) ### (response, predictor)


##################################################
### Lasso Model

fitL <- jags.model("ReplicaRegulaLasso.bug", data1, initsL, n.chains=3)

update(fitL,2000)
sampleL <- coda.samples(fitL, paramL, n.iter=2000, thin=10)

plot(sampleL)
summary(sampleL)
sumL <- summary(sampleL) 
tableL <- cbind(sumL$statistics[,1], sumL$quantiles[,"50%"],
                sumL$statistics[,2], sumL$quantiles[,c("2.5%","97.5%")])
colnames(tableL)[1:3] <- c("Mean", "Median", "SD")
write.table(tableL, file=paste0("tableLasso.txt"), quote=FALSE, sep=" ", row.names=TRUE)
round(tableL[,c(1,4,5)],5)

postL <- summary(sampleL)
PpredL <- postL$statistics[(K+H+2):((K+H+2)+N.new-1),1]
YpredL <- ifelse(PpredL>0.5,1,0)
(tabL <- table(Y.new,YpredL))

TP = tabL[2,2]
TN = tabL[1,1]
FP = tabL[1,2]
FN = tabL[2,1]
TP/(TP+FN) # Sensitivity
TN/(TN+FP) # Specificity
TP/(TP+FP) # Precision
(TP+TN)/N.new # Accuracy rate
pROC:: auc(Y.new,PpredL) ### (response, predictor)

##################################################
### Elastic Net Model

fitE <- jags.model("ReplicaRegulaElasticNet.bug", data1, initsE, n.chains=3)

update(fitE,2000)
sampleE <- coda.samples(fitE, paramE, n.iter=2000, thin=10)

plot(sampleE)
summary(sampleE)
sumE <- summary(sampleE) 
tableE <- cbind(sumE$statistics[,1], sumE$quantiles[,"50%"],
                sumE$statistics[,2], sumE$quantiles[,c("2.5%","97.5%")])
colnames(tableE)[1:3] <- c("Mean", "Median", "SD")
write.table(tableE, file=paste0("tableElasticNet.txt"), quote=FALSE, sep=" ", row.names=TRUE)
round(tableE[,c(1,4,5)],5)

postE <- summary(sampleE)
PpredE <- postE$statistics[(K+H+2):(K+H+2+N.new-1),1]
YpredE <- ifelse(PpredE>0.5,1,0)
(tabE <- table(Y.new,YpredE))

TP = tabE[2,2]
TN = tabE[1,1]
FP = tabE[1,2]
FN = tabE[2,1]
TP/(TP+FN) # Sensitivity
TN/(TN+FP) # Specificity
TP/(TP+FP) # Precision
(TP+TN)/N.new # Accuracy rate
pROC:: auc(Y.new,PpredE) ### (response, predictor)

##################################################

##################################################
