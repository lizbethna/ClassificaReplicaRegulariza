##################################################

Paper: 
Replication-based regularization Bayesian approaches to diagnose Reinke’s edema by using recorded speech 

Authors:
Lizbeth Naranjo (1), Carlos J. Perez (2), Yolanda Campos-Roca (3), Mario Madruga (2)

Journal:


(1) Departamento de Matematicas, Facultad de Ciencias, Universidad Nacional Autonoma de Mexico, Ciudad de Mexico, Mexico.
(2) Departamento de Matematicas, Facultad de Veterinaria, Universidad de Extremadura,  Caceres, Spain.
(3) Departamento de Tecnologias de los Computadores y las Comunicaciones, Escuela Politecnica, Universidad de Extremadura, Caceres, Spain

##################################################

Instructions to run the codes in R and JAGS are provided. 
The codes are applied to obtain a similar analysis as in Section 4 ‘Results’, but with simulated data and without cross-validation. 

##################################################
FILES 

The file ‘ReplicaRegulaReinke.R’ contains the R code. The JAGS code is run from this R file.

The files ‘ReplicaRegulaElasticNet.bug’, ‘ReplicaRegulaLasso.bug' and ‘ReplicaRegulaRidge.bug' contains the JAGS models. 

##################################################

To run the files, do the following.
 
1.- Download JAGS from www.mcmc-jags.sourceforge.net/

2.- Install the packages necessary to run the R file. 
These are indicated in the R file. 

3.- Change the address indicated in ‘setwd()’. 
setwd("HERE"). This is the address where the files ‘ReplicaRegulaxxxxx.bug’ are in.

4.- Run the R file, to simulate the data and to estimate the parameters of the  model.

##################################################
