#=========================================================================
# Econ 613 A2
# Yuqi Zhou
#=========================================================================
install.packages("fastDummies")
install.packages("tidyverse")
library(tidyverse)
library(readr)
library(ggplot2)
library(fastDummies)
setwd("/Users/jennyzhou/Desktop/ECON 613/a1/Data")
datind = list.files(pattern="datind")
# Read all the datind files by their names.
for (i in 1:16) {
  assign(datind[i], read.csv(datind[i]))
}

# ================= Exercise 1 OLS estimate ================
#(a)
# Get rid of all NA wage and NA age 
datind2009 <- subset(datind2009.csv, select = c("empstat", "age", "wage"))
datind2009 <- na.omit(datind2009) %>% filter(age > 0, wage != 0)
Y <- datind2009$wage
X <- datind2009$age
corr <- cor(Y,X, method = "pearson", use = "complete.obs")
# The correlation between Y and X is 0.143492.

#(b)
X <- cbind(1, datind2009$age)
beta_hat <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)
# beta_hat = 230.9923.

#(c)
# Using standard formulas of OLS: e = Y - beta*X
residuals <- Y - X %*% beta_hat
residuals_square <- t(residuals) %*% residuals/(length(Y) - 2)
beta_var <- residuals_square[1] * diag(2) %*% solve(t(X) %*% X)
std <- sqrt(diag(beta_var))
# The standard error is 14.8774.
reg = lm(wage ~ age, data = datind2009)
summary(reg)


#Using bootstrap with 49 replications
num_ind = nrow(datind2009) # number of individuals is 11469
num_var = length(reg$coefficients)  # number of variables is 2
set.seed(123)
bootf <- function(R,X,Y) {
  output = mat.or.vec(R, num_var)
  betas = c()
  for (i in 1:R) {
    samp = sample(1:num_ind, num_ind, replace = TRUE)
    sample_data = datind2009[samp, ]
    reg = lm(wage ~ age, data = datind2009)
    output[i,] = reg$coefficients
  }
 # mean_est = apply(output,2,mean)
  #sd_est = apply(output,2,sd)
  est = cbind(summary(reg)$coefficients[,1], summary(reg)$coefficients[,2])
  return(est)
}
output1 = bootf(49,X,Y)

# Using bootstrap with 499 replications

# ================= Exercise 2 Detrend Data =================
datind_05_to_18 <- rbind(datind2005.csv, datind2006.csv, datind2007.csv,
                         datind2008.csv, datind2009.csv, datind2010.csv, 
                         datind2011.csv, datind2012.csv, datind2013.csv, 
                         datind2014.csv, datind2015.csv, datind2016.csv, 
                         datind2017.csv, datind2018.csv)
#(a)
datind_05_to_18 <- subset(datind_05_to_18, select = c("year","empstat", "age", "wage"))
datind_05_to_18 <- na.omit(datind_05_to_18) %>% filter(age > 0, wage != 0)
# Get rid of people who are under 18 or wage is equal to 0.
ag <- data.frame(datind_05_to_18, bin=cut(datind_05_to_18$age, 
                                          c(18,25,30,35,40,45,50,55,60,120),
                                          include.lowest = TRUE))

#(b)
ag_plot <- ag %>% group_by(year) %>% summarise(mean_wage = mean(wage, na.rm = TRUE))
ggplot(data = ag_plot, mapping = aes(x = year, y = mean_wage)) + geom_point()

#(c)
reg3 = lm(wage ~ age + year, data = ag)
reg3_sum = summary(reg3)


# ============ Exercise 3 Numerical Optimization ============== Finished 
#(a)
# Import the data of 2007 and get rid of all Inactive/Retired labors.
datind2007 <- datind2007.csv %>% filter(age > 0, wage != 0, empstat != "Inactive", empstat != "Retired") 

#(b)
datind2007$empstat[which(datind2007$empstat == "Employed")] = 1 
datind2007$empstat[which(datind2007$empstat == "Unemployed")] = 0 
age <- datind2007$age
empstat <- as.numeric(datind2007$empstat)

flike <- function(par, age, empstat) {
  xbeta = par[1] + par[2]*age
  pr = pnorm(xbeta) 
  pr[pr>0.999999] = 0.999999 
  pr[pr<0.000001] = 0.000001
  like = empstat * log(pr) + (1-empstat) * log(1-pr)
  return(-sum(like)) 
}

#(c)
output3 <- mat.or.vec(100, 3) 
for (i in 1:100) {
  searchv = runif(2, -5, 5) 
  result  = optim(searchv, fn = flike, method = "BFGS", 
                  control = list(trace = 6, maxit = 3000),
                  age = age, empstat = empstat)
  output3[i,] = c(result$par, result$value) 
}
output3 <- as.data.frame(output3)
output3[which(output3$V3 == min(output3$V3)),] 

#(d)
wage <- datind2007$wage
flike2 <- function(par, age, wage, empstat) {
  xbeta = par[1] + par[2]*age + par[3]*wage
  pr = pnorm(xbeta) 
  pr[pr>0.999999] = 0.999999 
  pr[pr<0.000001] = 0.000001
  like = empstat * log(pr) + (1-empstat) * log(1-pr)
  return(-sum(like)) 
}
 
output4 <- mat.or.vec(100, 4) 
for (i in 1:100) {
  searchv = runif(3, -5, 5) 
  result  = optim(searchv, fn = flike2, method = "BFGS", 
                  control = list(trace = 6, maxit = 3000),
                  age = age, wage = wage, empstat = empstat)
  output4[i,] = c(result$par, result$value) 
}

output4 <- as.data.frame(output4)
output4[which(output4$V4 == min(output4$V4)), ] 


# ============== Exercise 4 Discrete choice ================
datind_05_to_15 <- rbind(datind2005.csv, datind2006.csv, datind2007.csv,
                         datind2008.csv, datind2009.csv, datind2010.csv, 
                         datind2011.csv, datind2012.csv, datind2013.csv, 
                         datind2014.csv, datind2015.csv)
#(a)
datind_05_to_15 <- na.omit(datind_05_to_15) %>% filter(age > 0, wage != 0, empstat != "Inactive", empstat != "Retired")

#(b)
# Make dummy variables for employment status
datind_05_to_15$empstat[which(datind_05_to_15$empstat == "Employed")] = 1
datind_05_to_15$empstat[which(datind_05_to_15$empstat == "Unemployed")] = 0
age2 <- datind_05_to_15$age
empstat2 <- as.numeric(datind_05_to_15$empstat)
datind_05_to_15 <- dummy_cols(datind_05_to_15, select_columns ='year',
                              remove_first_dummy = TRUE)
 
# =============== Probit Model ===================

flike2 <- function(par, age, year, empstat) {
  xbeta = par[1] + par[2]*age + par[3]*year5 + par[4]*year6 + par[5]*year7 + par[6]*year8 +
          par[7]*year9 + par[8]*year10 + par[9]*year11 + par[10]*year12 + par[11]*year13 + 
          par[12]*year14 + par[13]*year15 
  pr = pnorm(xbeta)
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like = empstat*log(pr) + (1-empstat)*log(1-pr)
  return(-sum(like))
}


# =============== Logit Model ================
logitlike <- function(par, age, year5, year6, year7, year8, year9, year10, year11,
                   year12, year13, year14, year15, empstat) {
  xbeta = par[1] + par[2]*age + par[3]*year5 + par[4]*year6 + par[5]*year7 + par[6]*year8 +
          par[7]*year9 + par[8]*year10 + par[9]*year11 + par[10]*year12 + par[11]*year13 + 
          par[12]*year14 + par[13]*year15 
  pr = exp(xbeta)/(1+exp(xbeta))
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like = empstat*log(pr) + (1-empstat)*log(1-pr)
  return(-sum(like))
}



# =============== Linearity Model ================
linearlike <- function(par, age, year5, year6, year7, year8, year9, year10, year11,
                       year12, year13, year14, year15, empstat) {
  xbeta = par[1] + par[2]*age + par[3]*year5 + par[4]*year6 + par[5]*year7 + par[6]*year8 +
          par[7]*year9 + par[8]*year10 + par[9]*year11 + par[10]*year12 + par[11]*year13 + 
          par[12]*year14 + par[13]*year15 
  pr = xbeta
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like = empstat*log(pr) + (1-empstat)*log(1-pr)
  return(-sum(like))
}


# =============== Exercise 5 Marginal Effects =================
#(a)



#  ========== Probit Model =============
probit_model <- function(beta,x1,x2){
  xbeta <- beta[1] + beta[2]*x1 + beta[3]*x2 
  pbeta <- matrix(dnorm(xbeta)) %*% t(beta)
  return(pbeta)
}



#  ========== Logit Model=============

logit_model <- function(beta,x1,x2){
  xbeta = beta[1] + beta[2]*x1 + beta[3]*x2 
  ebeta <- (exp(xbeta)/(1+exp(xbeta))^2) %*% t(beta)
  return(ebeta)
}

