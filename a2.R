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

# ================= Exercise 1 OLS estimate ================ Finished
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
beta_var <- residuals_square[1,1] * diag(2) %*% solve(t(X) %*% X)
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
    reg = lm(wage ~ age, data = sample_data)
    output[i,] = reg$coefficients
  }
  mean_est = apply(output,2,mean)
  sd_est   = apply(output,2,sd)
  
  est = cbind(summary(reg)$coefficients[,1], summary(reg)$coefficients[,2],
              mean_est, sd_est)
  colnames(est) = c("CF: est","CF: sd","BT: est","BT: sd")
  return(est)
}
output1 = bootf(49,X,Y)

# Using bootstrap with 499 replications
output2 = bootf(499,X,Y)

# ================= Exercise 2 Detrend Data ================= Finished
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
ag_data <- ag %>% group_by(year, bin) %>% summarise(average_wage = mean(wage, na.rm = TRUE))
ggplot(data = ag_data, mapping = aes(x = year, y = average_wage, color = bin)) + geom_line()

#(c)
reg3 = lm(wage ~ age + year, data = ag)
reg3_sum = summary(reg3)


# ============ Exercise 3 Numerical Optimization ============== Finished 
#(a)
# Import the data of 2007 and get rid of all Inactive/Retired labors.
datind2007 <- subset(datind2007.csv, select = c("empstat", "age", "wage"))
datind2007 <- datind2007 %>% filter(age > 0, wage != 0, empstat != "Inactive", empstat != "Retired") 

#(b)
datind2007$empstat[which(datind2007$empstat == "Employed")] = 1 
datind2007$empstat[which(datind2007$empstat == "Unemployed")] = 0 
age <- as.numeric(datind2007$age)
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
output3[which(output3$V2 == min(output3$V2)),] 

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
output4[which(output4$V3 == min(output4$V3)), ] 


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
year <- as.numeric(datind_05_to_15$empstat)

 
# =============== Probit Model ===================

fprobit <- function(par, age2, year, empstat2) {
  xbeta = par[1] + par[2] * age + par[3] * year
  pr = pnorm(xbeta)
  pr[pr > 0.999999] = 0.999999
  pr[pr < 0.000001] = 0.000001
  like = empstat * log(pr) + (1 - empstat) * log(1 - pr)
  return(-sum(like))
}

output5 <- mat.or.vec(100, 5)
 for (i in 1:100) {
  start_point = runif(4, -10, 10)
  result = optim(start_point, fn = fprobit, method = "BFGS", 
                 control = list(trace = 6, maxit = 3000), 
                 age = datind_05_to_15$age, year = datind_05_to_15$year, 
                 empstat = as.numeric(datind_05_to_15$empstat), hessian = TRUE)
  output5[i, ] = c(result$par, result$value)

}

output5 <- as.data.frame(output5)
output5[which.max(output5$V5 == min(output5$V5)), ]
fisher_info = solve(result$hessian)      
prop_sigma  = sqrt(diag(fisher_info))
est = cbind(par,summary(reg1)$coefficients[, 1],summary(reg1)$coefficients[, 2],res$par,prop_sigma)
colnames(est) = c("True parameter","R: GLM : est","R: GLM :se","R: own : est","R: own :se")
est

# =============== Logit Model ================

flogit <- function(par, age2, year, empstat2) {
  xbeta = par[1] + par[2] * age + par[3] * year
  pr = 1/(1+exp(-xbeta))
  pr[pr > 0.999999] = 0.999999
  pr[pr < 0.000001] = 0.000001
  like = empstat * log(pr) + (1 - empstat) * log(1 - pr)
  return(-sum(like))
}

output6 <- mat.or.vec(100, 5)
hessian_list <- list()
for (i in 1:100) {
  start_point = runif(4, -10, 10)
  result = optim(start_point, fn = flogit, method = "BFGS", 
                 control = list(trace = 6, maxit = 3000), 
                 age = datind_05_to_15$age, year = datind_05_to_15$year, 
                 empstat = as.numeric(datind_05_to_15$empstat), hessian = TRUE)
  output6[i, ] = c(result$par, result$value)
}

output6 <- as.data.frame(output6)
output6[which.max(output6$V5 == min(output6$V5)), ]


# =============== Linearity Model ================
flinear <- function(par, age2, year, empstat2) {
  xbeta = par[1] + par[2] * age + par[3] * year
  pr = xbeta
  pr[pr > 0.999999] = 0.999999
  pr[pr < 0.000001] = 0.000001
  like = empstat * log(pr) + (1 - empstat) * log(1 - pr)
  return(-sum(like))
}

output7 <- mat.or.vec(100, 5)
for (i in 1:100) {
  start_point = runif(4, -10, 10)
  result = optim(start_point, fn = flinear, method = "BFGS", 
                 control = list(trace = 6, maxit = 3000), 
                 age = datind_05_to_15$age, year = datind_05_to_15$year, 
                 empstat = as.numeric(datind_05_to_15$empstat), hessian = TRUE)
  output7[i, ] = c(result$par, result$value)
}

output7 <- as.data.frame(output7)
output7[which.max(output7$V5 == min(output7$V5)), ]



# =============== Exercise 5 Marginal Effects =================

#  ========== Probit Model =============
set.seed(123)

m_probit <- function(fun, data, reps = 100, digits = 3) {
  x <- glm(fun, data, family = binomial(link = "probit"))
  pdf_probit <- mean(dnorm(predict(x, type = "link")))
  m_probit <- pdf_probit * coef(x)
  outputs6 <- matrix(rep(NA, reps * length(coef(x))), nrow = reps)
  for(i in 1:reps) {
    samp <- sample(1:dim(data)[1], dim(data)[1], rep = TRUE)
    data_probit <- data[samp, ]
    reg_probit <- glm(fun, data_probit, family = binomial(link = "probit"))
    pr <- mean(dnorm(predict(reg_probit, type = "link")))
    outputs6[i, ] <- pr * coef(reg_probit)
  }
  probit <- cbind(m_probit, apply(outputs6, 2, sd))
  colnames(probit) <- c("Marginal Effect", "SE")  
  return(probit)
}

m_probit(fun <- empstat2 ~ age2 + year, data <- datind_05_to_15)

#  ========== Logit Model =============
m_logit <- function(fun, data, reps = 100, digits = 3) {
  x <- glm(fun, data, family = binomial(link = "logit"))
  pdf_logit <- mean(dlogis(predict(x, type = "link")))
  m_logit <- pdf_logit * coef(x)
  output7 <- matrix(rep(NA, reps * length(coef(x))), nrow = reps)
  for(i in 1:reps) {
    samp = sample(1:dim(data)[1], dim(data)[1], rep = TRUE)
    data_logit = data[samp, ]
    reg_logit <- glm(fun, data_logit, family = binomial(link = "logit"))
    pr <- mean(dlogis(predict(reg_logit, type = "link")))
    output7[i, ] <- pr * coef(reg_logit)
  }
  logit() <- cbind(m_logit, apply(output7, 2, sd))
  colnames(logit) <- c("Marginal Effect", "SE")  
  return(logit)
}

m_logit(fun = empstat2 ~ age2 + year, data = datind_05_to_15)

