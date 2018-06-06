# use simulated data to demonstrate the codes

# set working directory the path to the directory of 'demo_files'
setwd('demo_files')

# load required library for data generation
library(plyr)   # split 3D array to list of matrices
library(expm)   # calculate square root of a matrix

# load all R function codes
files = list.files('Code_Biostat/')
sapply(files,function(f) source(file.path('Code_Biostat',f)))
   

#############################
### Generat Complete Data ###
#############################

load('true_para.rda')

S=42; J=2; p=5; T=7; n=50; # no spatial covariate in simulation

true$nu = 20;

X = matrix(rnorm(p*n),p,n)
true$xbeta = t(true$beta)%*%X

true$gamma = rgamma(n,true$nu/2,rate=true$nu/2)
true$z = abs(sapply(true$gamma, function(x) matrix(rnorm(S*J)/sqrt(x),S,J), simplify='array'))
true$dzv = sapply(alply(true$z,3),function(x) diag(true$eta)%*%x%*% sqrtm(true$V), simplify='array')
true$U = replicate(n,true$u) + aperm(replicate(S,true$xbeta),c(3,1,2)) 
true$y = mapply(function(x1,x2) t(chol(solve(true$Omega)/x1))%*%x2%*%chol(true$V), true$gamma, 
                alply(array(rnorm(S*J*n),c(S,J,n)),3),SIMPLIFY='array')+true$U


########################
# Generate missingness #
########################

meanop =  kronecker(diag(1,T),rep(1,S/T)/(S/T)) 
true$miss_quant = sapply(alply(true$y,3),function(x) t(meanop)%*%x%*%true$misspara[2:3]+true$misspara[1]) 
true$miss_prob = apply(true$miss_quant,1:2,pnorm);

delta = 1*( matrix(runif(T*n),T,n) < true$miss_prob )
true$yMISS = aperm(replicate(2,apply(delta,2,rep,each=S/T)),c(1,3,2))
y = true$y; y[true$yMISS==1] = NA;


####################
### Run Analysis ###
####################

nu = true$nu;

N=2000; burn=1000; thin=1; update=100; ncore=8;

postsamp = mvst.info.miss(y,delta,X,missing.type='probit',nu=nu,runs=N,burn=burn,thin=thin,update=update,ncore=ncore)




