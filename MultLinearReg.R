###############################################################################
##In this code I develop a Markov chain Monte Carlo based normal linear
##regression model using non informative priors, at the end there is an
##applicationg with real data relating the birth rate of females between
##15 and 17 years old to a poverty index for all 51 states of the USA.
###############################################################################
###Simulating Simple Linear Regression
n <- 20
epsilon <- rnorm(n,0,1)
x <- cbind(rep(1,n),rnorm(n,0,4))
beta <- c(13,-4)
y <- (x%*%beta) + epsilon

hist(y)

###############################################################################
library(MASS) #This package needs to be load, so we can sample from a multivatiate
              #normal distribution

updateBeta <- function(Sig2Be,X,Y,Sig2){ 
  #This function provides a Gibbs sampler for the regression parameters, 
  # and it can adapt to simple and multiple linear regression, all that is
  #need is to adjust the dimension of the design matrix
  n <- dim(X)[2]
  B0 <- diag(Sig2Be,n,n)
  Sigma <- Sig2 * diag(1,ncol=n,nrow=n)
  sig2.post <- solve((t(X)%*%X)%*%Sigma + B0)
  m.post <- (t((t(X)%*%Y))%*%Sigma) %*% sig2.post
  
  return(mvrnorm(1,m.post,sig2.post))
}

updateSig2 <- function(a,b,N,Y,X,Beta){
  #This function provides a Gibbs sampler for the contant variance in a normal
  #linear regression model
  a.post <- N/2 + a
  b.post <- b + (t(Y - X%*%Beta)%*%(Y - X%*%Beta))*0.5
  
  return(1/rgamma(1,a.post,b.post))
}

############################################################################
Niter <- 10000
Beta.out <- array(NA, dim = c(Niter,dim(x)[2]))
Sig2.out <- array(NA,dim=Niter)

Beta.out[1,] <- beta
Sig2.out[1] <- 1
X <- x
Y <- y
N <- length(Y)

for(i in 2:Niter){
  Beta.out[i,] <- updateBeta(0.001,X,Y,Sig2.out[i-1])
  Sig2.out[i] <- updateSig2(0.01,0.01,N,Y,X,Beta.out[i,])
  print(i)
}

plot(Beta.out[,1],type='l')
mean(Beta.out[,1])
hist(Beta.out[,1])

plot(Beta.out[,2],type='l')
mean(Beta.out[,2])
hist(Beta.out[,2])

plot(Sig2.out,type='l')
mean(Sig2.out)
hist(Sig2.out)

#############################################################################
#I got the data from https://newonlinecourses.science.psu.edu/stat462/node/101/
#The response is birth rate per 1000 females 15 to 17 years old, for all
#51 states of the United States and the regressor is the poverty rate,
#which is the percent of the state's population living in households with
# incomes below the federally defined poverty level.
#############################################################################

#####Loading and checking the data
poverty <- as.matrix(read.csv("Poverty.csv", sep = ';')[-1])
head(poverty)
plot(poverty[,1],poverty[,2])
abline(lm(poverty[,2] ~ poverty[,1]),col='red')

hist(poverty[,1])
shapiro.test(poverty[,1])
############################################################################
####Variables
X <- cbind(rep(1,length(poverty[,1])), poverty[,1])
Y <- poverty[,2]
N <- length(Y)

##MCMC
Niter <- 100000
Beta.out <- array(NA, dim = c(Niter,dim(X)[2]))
Sig2.out <- array(NA,dim=Niter)

Beta.out[1,] <- rep(1,1)
Sig2.out[1] <- 1

for(i in 2:Niter){
  Beta.out[i,] <- updateBeta(0.001,X,Y,Sig2.out[i-1])
  Sig2.out[i] <- updateSig2(0.01,0.01,N,Y,X,Beta.out[i,])
  print(i)
}

#####Briefly checking convergence of parameters 
plot(Beta.out[10001:100000,1],type='l')
mean(Beta.out[10001:100000,1])
hist(Beta.out[10001:100000,1])

plot(Beta.out[10001:100000,2],type='l')
mean(Beta.out[10001:100000,2])
hist(Beta.out[10001:100000,2])

plot(Sig2.out[10001:100000],type='l')
mean(Sig2.out[10001:100000])
hist(Sig2.out[10001:100000])

##################################################################################
#####Calculating the Fitted Values

Nburn <- 10000 #Burn-in phase 
Thinning <- 10 #Thinning to reduce autocorrelation
Beta <- Beta.out[seq(Nburn+1,Niter,by=Thinning),]
Sig2 <- Sig2.out[seq(Nburn+1,Niter,by=Thinning)]

####Fitting the values
Fitted <- array(NA, dim = c(nrow(Beta),N))

for(i in 1:nrow(Fitted)){
    Fitted[i,] <- X %*% Beta[i,] + sqrt(Sig2[i])*rnorm(N)
}
medians <- apply(Fitted,2,median)

####Calculating the residuals 
Residuals <- Y - medians
plot(Residuals)
####Testing for normality
shapiro.test(Residuals) #Apparently the residuals follow a normal distribution