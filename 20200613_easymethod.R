rm(list = ls())
library("xlsx")
require("xlsx")

library("stats")

library("grDevices")
library("actuar")

library("copula")
library("psych")
library("dplyr")

library("FAdist")
require("grid")

library("vcd")

library("MASS")

library("fitdistrplus")
require("survival")
require("npsurv")
require("lsei")
library("EnvStats")
library("ismev")
library("mgcv")
library("nlme")
library("reliaR")
library("lestat") #inverse CDF

library("goftest")
library("goft")
library("gumbel")


# Read data from excel flie
data <- read.xlsx(file="E:\\R_coding\\test_2.xlsx",
                  sheetIndex="DATA",startRow = 1,endRow = 960,
                  header = T,colIndex =2:3,
                  encoding = "UTF-8")

attach(data)
#
# ----------------------------- Remove missing value ------------------------------
#
summary(data)
sat <- describe(data)

# ?ˬd???S??missing value
sum(is.na(data$Q)) #data have missing value
sum(is.na(data$S)) #data have missing value

# ????missing value
complete.cases(data)
rm.data <- data[complete.cases(data), ]

# ?w?q?s??Q?MS
Q <- rm.data$Q
S <- rm.data$S

# weibull distribution
Q.par.weibull   <- fitdist(Q, distr = "weibull", method = "mle") # Q parameter
Q.shape.weibull <- Q.par.weibull[["estimate"]][["shape"]]
Q.scale.weibull <- Q.par.weibull[["estimate"]][["scale"]]
S.par.weibull   <- fitdist(S, distr = "weibull", method = "mle",lower=c(0, 0))# s parameter
S.shape.weibull <- S.par.weibull[["estimate"]][["shape"]]
S.scale.weibull <- S.par.weibull[["estimate"]][["scale"]]
## gamma distribution
Q.par.gamma   <- fitdist(Q, distr = "gamma", method = "mle",lower=c(0, 0)) # Q parameter
Q.shape.gamma <- Q.par.gamma[["estimate"]][["shape"]]
Q.rate.gamma  <- Q.par.gamma[["estimate"]][["rate"]]
S.par.gamma   <- fitdist(S, distr = "gamma", method = "mle",lower=c(0, 0))# s parameter
S.shape.gamma <- S.par.gamma[["estimate"]][["shape"]]
S.rate.gamma  <- S.par.gamma[["estimate"]][["rate"]]
## log-normal distribution
Q.par.lnorm       <- fitdist(Q,"lnorm",method = "mle") # Q parameter
Q.meanlog.lnorm   <- Q.par.lnorm[["estimate"]][["meanlog"]]
Q.sdlog.lnorm     <- Q.par.lnorm[["estimate"]][["sdlog"]]
S.par.lnorm       <- fitdist(S,"lnorm",method = "mle") # S parameter
S.meanlog.lnorm   <- S.par.lnorm[["estimate"]][["meanlog"]]
S.sdlog.lnorm     <- S.par.lnorm[["estimate"]][["sdlog"]]
# pearson tpye three distribution
# Q.par.lgamma      <- fitdist(Q,"lgamma",method = "mle") # Q parameter
# Q.shapelog.lgamma <- Q.par.lgamma[["estimate"]][["shapelog"]]
# Q.ratelog.lgamma  <- Q.par.lgamma[["estimate"]][["ratelog"]]
# S.par.lgamma      <- fitdist(S,"lgamma",method = "mle") # S parameter
# S.shapelog.lgamma <- S.par.lgamma[["estimate"]][["shapelog"]]
# S.ratelog.lgamma  <- S.par.lgamma[["estimate"]][["ratelog"]]
# gumbel distribution
Q.par.gumbel      <- eevd(Q, method = "mle") #EnvStats # Q parameter
Q.location.gumbel <- Q.par.gumbel[["parameters"]][["location"]]
Q.scale.gumbel    <- Q.par.gumbel[["parameters"]][["scale"]]
S.par.gumbel      <- eevd(S, method = "mle") #EnvStats # S parameter
S.location.gumbel <- S.par.gumbel[["parameters"]][["location"]]
S.scale.gumbel    <- S.par.gumbel[["parameters"]][["scale"]]

plot(Q.par.weibull)
plot(Q.par.gamma)
plot(Q.par.lnorm)
#plot(Q.par.lgamma)


plot(S.par.weibull)
plot(S.par.gamma)
plot(S.par.lnorm)
#plot(S.par.lgamma)