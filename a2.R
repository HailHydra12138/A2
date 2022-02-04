#=========================================================================
# Econ 613 A1
# Yuqi Zhou
#=========================================================================
install.packages("tidyverse")
library(tidyverse)
library(readr)
library(ggplot2)
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
residuals <- Y - beta_hat * X
residuals_square <- t(residuals) %*% residuals/(nrow(X) - ncol(X))
beta_var <- residuals_square %*% solve(t(X) %*% X)
std <- sqrt(diag(beta_var))
# The standard deviation is 509.2344
reg = lm(wage ~ age, data = datind2009)
summary(reg)


#Using bootstrap with 49 replications
R = 49 # number of bootstraps 
num_ind = nrow(datind2009) # number of individuals
num_var = length(reg$coefficients)  # number of variables
output = mat.or.vec(R, num_var)
set.seed(123)

for (i in 1:R)
{
  samp = sample(1:num_ind, num_ind, rep = TRUE)
  sample_data = datind2009[samp, ]
  reg = lm(wage ~ age, data = datind2009)
  output[i,] = reg$coefficients
}

mean1 = apply(output, 2, mean)
sd1 = apply(output, 2, sd)
est1 = cbind(summary(reg)$coefficients[ , 1], summary(reg)$coefficients[ , 2], 
             mean1, sd1)
colnames(est1) = c("CF: est","CF: std","BT (49): estimate","BT (49): std dev")
est1

# Using bootstrap with 499 replications
R2 = 499 # number of bootstraps 
num_ind = nrow(datind2009) # number of individuals
num_var = length(reg$coefficients)  # number of variables
output2 = mat.or.vec(R2, num_var)
set.seed(123)

for (i in 1:R2) {
  samp = sample(1:num_ind, num_ind, rep = TRUE)
  sample_data = datind2009[samp, ]
  reg = lm(wage ~ age, data = datind2009)
  output2[i,] = reg$coefficients
}

mean2 = apply(output2, 2, mean)
sd2 = apply(output2, 2, sd)
est2 = cbind(summary(reg)$coefficients[ , 1], summary(reg)$coefficients[ , 2], 
             mean2, sd2)
colnames(est1) = c("CF: est","CF: std","BT (499): estimate","BT (499): std dev")
est2

# ================= Exercise 2 Detrend Data =================
datind_05_to_18 <- rbind(datind2005.csv, datind2006.csv, datind2007.csv,
                         datind2008.csv, datind2009.csv, datind2010.csv, 
                         datind2011.csv, datind2012.csv, datind2013.csv, 
                         datind2014.csv, datind2015.csv, datind2016.csv, 
                         datind2017.csv, datind2018.csv)
#(a)
datind_05_to_18 <- subset(datind_05_to_18, select = c("year","empstat", "age", "wage"))
datind_05_to_18 <- na.omit(datind_05_to_18) %>% filter(age > 0, wage != 0)
ag <- data.frame(datind_05_to_18, bin=cut(datind_05_to_18$age, 
                                          c(18,25,30,35,40,45,50,55,60,120),
                                          include.lowest = TRUE))

#(b)
ag_plot <- ag %>% group_by(year) %>% summarise(mean_wage = mean(wage, na.rm = TRUE))
ggplot(data = ag_plot, mapping = aes(x = year, y = mean_wage)) + geom_point()

#(c)
reg3 = lm(wage ~ age + year, data = ag)
summary(reg3)

# ============ Exercise 3 Numerical Optimization ==============
#(a)
# Import the data of 2007 and get rid of all Inactive labors.
datind2007 <- subset(datind2007.csv, select = c("empstat", "age", "wage"))
datind2007 <- na.omit(datind2007) %>% filter(age > 0, wage != 0, empstat != "Inactive", empstat != "Retired")

#(b)
datind2007$empstat[which(datind2007$empstat == "Employed")] = 1
datind2007$empstat[which(datind2007$empstat == "Unemployed")] = 0
age <- datind2007$age
empstat <- as.numeric(datind2007$empstat)
flike <-function(par, age, empstat) {
  xbeta = par[1] + par[2]*age
  pr = pnorm(xbeta)
  # pr = exp(beta)/(1+exp(beta)) logit
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like = empstat*log(pr) + (1-empstat)*log(1-pr)
  return(-sum(like))
}

#(c)
output3 <- mat.or.vec(100, 3) 
for (i in 1:100) {
  start = runif(2, -5, 5)
  result  = optim(start, fn = flike, method = "BFGS", control = list(trace = 6, maxit = 3000),
                         age = age, empstat = empstat)
  
  output3[i,] = c(result$par, result$value)
}


#(d)
wage <- datind2007$wage
flike2 = function(par, age, wage, empstat) {
  xbeta = par[1] + par[2]*age +par[3]*wage
  pr = pnorm(xbeta)
  # pr = exp(beta)/(1+exp(beta)) logit
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like = empstat*log(pr) + (1-empstat)*log(1-pr)
  return(-sum(like))
}

output4 <- mat.or.vec(100, 4) 
for (i in 1:100) {
  start = runif(2, -5, 5)
  result  = optim(start, fn = flike, method = "BFGS", control = list(trace = 6, maxit = 3000),
                  age = age, empstat = empstat)
  
  output3[i,] = c(result$par, result$value)
}


# ============== Exercise 4 Discrete choice ================
datind_05_to_15 <- rbind(datind2005.csv, datind2006.csv, datind2007.csv,
                         datind2008.csv, datind2009.csv, datind2010.csv, 
                         datind2011.csv, datind2012.csv, datind2013.csv, 
                         datind2014.csv, datind2015.csv)
#(a)
datind_05_to_15 <- subset(datind_05_to_15, select = c("year","empstat", "age", "wage"))
datind_05_to_15 <- na.omit(datind_05_to_18) %>% filter(age > 0, wage != 0, empstat != "Inactive", empstat != "Retired")

#(b)
# Make dummy variables for employment status
datind_05_to_15$empstat[which(datind_05_to_15$empstat == "Employed")] = 1
datind_05_to_15$empstat[which(datind_05_to_15$empstat == "Unemployed")] = 0
age2 <- datind_05_to_15$age
empstat2 <- as.numeric(datind_05_to_15$empstat)



# =============== Probit Model ================
flike2 <- function(par, age, year5, year6, year7, year8, year9, year10, year11,
                   year12, year13, year14, year15, empstat) {
  xbeta = par[1] + par[2]*age + par[3]*year5 + par[4]*year6 + par[5]*year7 + par[6]*year8 +
          par[7]*year9 + par[8]*year10 + par[9]*year11 + par[10]*year12 + par[11]*year13 + 
          par[12]*year14 + par[13]*year15 
  pr = pnorm(xbeta)
  # pr = exp(beta)/(1+exp(beta)) logit
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
  # pr = exp(beta)/(1+exp(beta)) logit
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
  # pr = exp(beta)/(1+exp(beta)) logit
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like = empstat*log(pr) + (1-empstat)*log(1-pr)
  return(-sum(like))
}



# =============== Exercise 5 Marginal Effects =================



