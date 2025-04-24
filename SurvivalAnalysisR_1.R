library(survival)
library(parmsurvfit)

########EX B###############

#### insert data ############
cbdat <- read.table("C:\\Users\\Admin\\Desktop\\8ο ΕΞΑΜΗΜΟ\\cbb.txt", header = TRUE)
attach(cbdat) 
cbdat


##Kaplan Meier Estimates

outp<-survfit(Surv(t, censor)~1) 
summary(outp) 
plot(outp, main=expression(paste("Kaplan-Meier-estimate ", hat(S)(t))), 
conf.int=FALSE,xlab="t", ylab="Survival")
________________________________________________________________________


##METHOD KAPLAN MEIER

##EXPONENTIAL



outp$surv
SKM<-outp$surv [outp$n.event ==1] 
SKM 
Utime<-outp$time[outp$n.event ==1] 
Utime

#----------------------------------------------------------------------------
##Examine whether data follows Exponential distribution
 
plot(Utime,-log(SKM), main=expression(paste("Exponential Plot")), xlab="t", 
ylab="-ln S(t)", pch=19) 
abline(lm(-log(SKM)~Utime))


#----------------------------------------------------------------------------
##Examine whether data follows Weibull distribution

plot(log(Utime),log(-log(SKM)), main=expression(paste("Weibull Plot")), xlab="ln + t", ylab="ln-(ln S(t))", pch=19)
abline(lm(log(-log(SKM))~log(Utime)))
#----------------------------------------------------------------------------
##Examine whether data follows LogNormal distribution 

norm_quant<-qnorm(1-SKM,0,1) 
norm_quant 

plot(log(Utime), norm_quant, pch=19, main=expression(paste("LogNormal Plot"))) 
abline(lm(norm_quant ~ log(Utime)))


#----------------------------------------------------------------------------
##Examine whether data follows LogLogistic distribution 

plot(log(Utime),log((1-SKM)/SKM), pch=19, main=expression(paste("LogLogistic Plot"))) 
abline(lm(log((1-SKM)/SKM) ~ log(Utime)))


#________________________________________________________________________________
#iii)
##METHOD OF MAXIMUM LIKELIHOOD

#EXPONENTIAL
mod1<-survreg(Surv(t, censor)~1, data=cbdat, dist="exponential")
mod1
summary(mod1)

----------------------------------------------------

#install.packages("parmsurvfit")
#library(parmsurvfit)


#WEIBULL
modweib<-fit_data(data = cbdat, dist = "weibull",    time = "t",   censor = "censor")
summary(modweib)

#1.
outp<-survfit(Surv(t,censor)~1)
SKM<-outp$surv [outp$n.event ==1]
SKM
Utime<-outp$time[outp$n.event ==1]
Utime
#2.
Weibml<- pweibull(Utime, shape= modweib$estimate[1], scale = modweib$estimate[2],
lower.tail = FALSE, log.p = FALSE)
Weibml
#3.
plot(Weibml~SKM, pch=19, main=expression(paste("S-ML Weibull vs S-KM")))
abline(0,1)
plot(Weibml~Utime,type="l",col="red", main=expression(paste("S-ML Weibull and S-KM vs Time")))
lines(SKM~Utime,col="black")

#----------------------------------------------------
#LONORMAL
modlnorm<-fit_data(data = cbdat, dist = "lnorm",    time = "t",   censor = "censor")
summary(modlnorm)

lnorm.ml<-plnorm(Utime, meanlog=modlnorm$estimate[1] , sdlog = modlnorm$estimate[2],
lower.tail = FALSE, log.p = FALSE)
lnorm.ml

plot(lnorm.ml~Utime,type="l",col="red", main=expression(paste("S-ML Weibull,S-ML LogNormal and S-KM vs Time")))
lines(SKM~Utime,col="black")
lines(Weibml~Utime,col="green")
#-----------------------------------------------------
#LOGLOGISTIC
install.packages("actuar")
library(actuar)

mod.log.logis<-fit_data(data = cbdat, dist = "llogis",    time = "t",   censor = "censor")
summary(mod.log.logis)

log.logis.ml <-pllogis(Utime, shape = mod.log.logis$estimate[1], scale =
mod.log.logis$estimate[2], lower.tail = FALSE, log.p = FALSE)
log.logis.ml

plot(lnorm.ml~Utime,type="l", col="red", main=expression(paste("S-ML Weibull,S-ML LogNormal,S-ML LogLogistic  and S-KM vs Time")))
lines(SKM~Utime,col="black")
lines(Weibml~Utime,col="green")
lines(log.logis.ml~Utime,col="blue")

#______________________________________________________________
#iv )
compute_AD (data = cbdat, dist = "exp",    time = "t",   censor = "censor")
compute_AD (data = cbdat, dist = "weibull",    time = "t",   censor = "censor")
compute_AD (data = cbdat, dist = "lnorm",    time = "t",   censor = "censor")
compute_AD (data = cbdat, dist = "llogis",    time = "t",   censor = "censor")

#______________________________________________________________
#v )

summary(mod1)

mod2<- survreg(Surv(t, censor)~1, data=cbdat, dist="weibull")
mod3<- survreg(Surv(t, censor)~1, data=cbdat, dist="lognormal")
mod4<- survreg(Surv(t, censor)~1, data=cbdat, dist="loglogistic")
AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)
AIC(mod1)




#___________________________________________________________________________


#########################################################
################ ΕΧ C ###################################
#########################################################

library(splines)
library(survival)

cc <- read.table("C:\\Users\\Admin\\Desktop\\8ο ΕΞΑΜΗΜΟ\\lung-patients.txt", header = TRUE)
attach(cc) 
cc
#____________________________________________________________________________
#i)
#KAPLAN-MEIER
outp<-survfit(Surv(t,c)~Group)
summary(outp)

plot(outp, lty = 1:2, main=expression(paste("Kaplan-Meier-estimate ", hat(S)(t))),
xlab="t", ylab="Survival", lwd=2)
legend("topright", c("Group 1", "Group 2"), lty = 1:2)

#____________________________________________________________________________
#ii)

#Log-rank test  (Mantel-Haenszel test) with R
out1<-survdiff(Surv(t, c) ~ Group)
out1

#_____________________________________________________________________________

#iii) FIRST APPROACH TO SOLVING THE EXERCISE

##########Group0###############################

group0<- Surv(t[Group=="0"],c[Group=="0"])
outp0<-survfit(group0~1, type="kaplan-meier",data=cc)
summary(outp0)



#Extract values when Event=1

Utime0<- outp0$time[outp0$n.event==1]
Utime0
SKM0<-outp0$surv[outp0$n.event==1]  
SKM0

#Weibull
plot(log(Utime0),log(-log(SKM0)) , main=expression(paste("Weibull Plot – Group0")),
xlab="ln t", ylab="ln-(ln S(t))", pch=19)
abline(lm(log(-log(SKM0))~ log(Utime0)))

#Log-Logistic
plot(log(Utime0),log((1-SKM0)/SKM0), main=expression(paste("Log-Logistic Plot–
Group0")), xlab="ln t", ylab=expression(ln (F(t)/S(t))), pch=19)
abline(lm(log((1-SKM0)/SKM0) ~ log(Utime0)))

##########Group1###############################

group1<- Surv(t[Group=="1"],c[Group=="1"])
outp1<-survfit(group1~1, type="kaplan-meier",data=cc)
summary(outp1)


#Extract values when Event=1

Utime1<- outp1$time[outp1$n.event==1]
Utime1
SKM1<-outp1$surv[outp1$n.event==1] 
SKM1 

#Weibull
plot(log(Utime1),log(-log(SKM1)) , main=expression(paste("Weibull Plot – Group1")),
xlab="ln t", ylab="ln-(ln S(t))", pch=19)
abline(lm(log(-log(SKM1))~ log(Utime1)))


#Log-Logistic
plot(log(Utime1),log((1-SKM1)/SKM1), main=expression(paste("Log-Logistic Plot–
Group1")), xlab="ln t", ylab=expression(ln (F(t)/S(t))), pch=19)
abline(lm(log((1-SKM1)/SKM1) ~ log(Utime1)))


##########################################################################
##########################################################################
##########################################################################

#iii) SECOND APPROACH TO SOLVING THE EXERCISE - BETTER VERSSION

library(survival)
library(parmsurvfit)

cc <- read.table("C:\\Users\\Admin\\Desktop\\8ο ΕΞΑΜΗΜΟ\\lung-patients.txt", header = TRUE)
attach(cc) 
cc

group0<- Surv(t[Group=="0"],c[Group=="0"])
outp0<-survfit(group0~1, type="kaplan-meier",data=cc)
summary(outp0)

Utime0<- outp0$time[outp0$n.event==1]
SKM0<-outp0$surv[outp0$n.event==1]

##########################################################################

###Split data into two groups

ccdat0 <-subset(cc, Group ==0)
ccdat1 <-subset(cc, Group ==1)

#-------------------------------------------------------------------------
###Parametric Results for 1st group

###Fitting Weibull distribution

modweib0<-fit_data(data = ccdat0, dist = "weibull", time = "t", censor = "c")
summary(modweib0)

SWeib.ml0<- pweibull(Utime0, shape= modweib0$estimate[1], scale =
modweib0$estimate[2], lower.tail = FALSE, log.p = FALSE)
SWeib.ml0

plot_ppsurv(ccdat0, "weibull", time = "t",  censor = "c")

###Fitting Log-logistic distribution

library(actuar)

modllog0<-fit_data(data = ccdat0, dist = "llogis", time = "t", censor = "c")
summary(modllog0)

Sllogis.ml0<-pllogis(Utime0, shape = modllog0$estimate[1], scale =
modllog0$estimate[2], lower.tail = FALSE, log.p = FALSE)
Sllogis.ml0

plot_ppsurv(ccdat0, "llogis", time = "t",  censor = "c")


plot(SKM0~ Utime0,type="l", col="black", main=expression(paste("Group 0  - 
Estimated Survival Functions vs Time")), xlab="Time", ylab="S(t)")
lines(SWeib.ml0~ Utime0,col="green", lwd=2)
lines(Sllogis.ml0~ Utime0,col="blue", lwd=2)  
legend("topright", c("S-KM", "Sweib", "Sllogis"), col=c("black", "green","blue"), lty=c(1,2,3) )


#-------------------------------------------------------------------------
###Parametric Results for 2nd group

###Fitting Weibull distribution

group1<- Surv(t[Group=="1"],c[Group=="1"])
outp1<-survfit(group1~1, type="kaplan-meier",data=cc)
summary(outp1)

Utime1<- outp1$time[outp1$n.event==1]
SKM1<-outp1$surv[outp1$n.event==1] 

modweib1<-fit_data(data = ccdat1, dist = "weibull", time = "t", censor ="c")
summary(modweib1)

SWeib.ml1<- pweibull(Utime1, shape= modweib1$estimate[1], scale =
modweib1$estimate[2], lower.tail = FALSE, log.p = FALSE)
SWeib.ml1

plot_ppsurv(ccdat1, "weibull", time = "t",  censor = "c")


###Fitting Log-logistic distribution

library(actuar)

modllog1<-fit_data(data = ccdat1, dist = "llogis", time = "t", censor = "c")
summary(modllog1)

Sllogis.ml1<-pllogis(Utime1, shape = modllog1$estimate[1], scale =
modllog1$estimate[2], lower.tail = FALSE, log.p = FALSE)
Sllogis.ml1

plot_ppsurv(ccdat1, "llogis", time = "t",  censor = "c")

plot(SKM1~ Utime1,type="l", col="black", main=expression(paste("Group 1 (NHL) - 
Estimated Survival Functions vs Time")), xlab="Time", ylab="S(t)")
lines(SWeib.ml1~ Utime1,col="green", lwd=2)
lines(Sllogis.ml1~ Utime1,col="blue", lwd=2)  
legend("topright", c("S-KM", "Sweib", "Sllogis"), col=c("black", "green","blue"), lty=c(1,2,3) )
############################################

plot_surv(cc, "weibull",  time = "t", censor = "c", by = "Group")

#_____________________________________________________________________________

#iV) 
modllog0<-fit_data(data = ccdat0, dist = "llogis", time = "t", censor = 
"c") 
summary(modllog0) 

##########################################################################
##########################################################################
##########################################################################






