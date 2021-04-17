#Clear
#...variables
rm(list=ls())
#...console
cat("\014\n")
#...graphs
tryCatch(dev.off(), error = function(e) {NULL})
dev.new()

library(tidyverse)
library(MAST)

#Trying the example of running glm with only matrices
#...not SingleCellAssay object

#sample size of 500
#x is continuous 
#z is bernoulli
data<- data.frame(x=rnorm(500), z=rbinom(500, 1, .3))

#Probability distribution based on covariates
logit.y <- with(data, x*2 + z*2)
mu.y <- with(data, 10+10*x+10*z + rnorm(500))

#Output based on the covariates
y <- (runif(500)<exp(logit.y)/(1+exp(logit.y)))*1
y[y>0] <- mu.y[y>0]
data$y <- y

fit <- zlm(y ~ x+z, data)

#Get B_x = 2, B_z = 2, B_cons = 0
summary.glm(fit$disc)
#Get B_x = 10, B_z = 10, B_cons = 10
summary.glm(fit$cont)
