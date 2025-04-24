
##SEIRA 3 ###

#EXERCISE A

#--------------------------------------PART 1----------------------------------------
library(survival) 
cdat<-read.table("C:\\Users\\Admin\\Desktop\\8ο ΕΞΑΜΗΜΟ\\ΜΟΝΤΕΛΑ ΑΞΙΟΠΙΣΤΙΑΣ ΚΑΙ ΕΠΙΒΙΩΣΣ\\σειρα 3\\lung-cancer.txt",header=TRUE) 
attach(cdat) 
cdat 

#Kaplan Meier

outp<-survfit(Surv(time,status)~cell.type) 
summary(outp)

#Kaplan-Meier GRAPΗ

# Sample colors to use for the lines
colors <- c("green", "blue","red", "purple")

# Plot the Kaplan-Meier estimates with different line types and colors
plot(outp, lty = 1:4, col = colors, main = expression(paste("Kaplan-Meier-estimate ", hat(S)(t))),
     xlab = "t", ylab = "Survival", lwd = 2)

# Add a legend with different line types and colors
legend("topright", c("Cell Type 1", "Cell Type 2", "Cell Type 3", "Cell Type 4"), lty = 1:4, col = colors)

#ελέγχους log-rank και Wilcoxon

out1<-survdiff(Surv(time,status) ~ cell.type)
out1

out2<-survdiff(Surv(time, status) ~ cell.type, rho=1)
out2
#----------------- END  of PART 1 --------------------------------------

#--------------------------------------PART 2----------------------------------------

###Kaplan-Meier for 1st group
group1<- Surv(time[cell.type=="1"],status[cell.type=="1"]) 
outp1<-survfit(group1~1, type="kaplan-meier", data=cdat)
#Extract values when status=1
Utime1<- outp1$time[outp1$n.event==1] 
SKM1<-outp1$surv[outp1$n.event==1] 

###Kaplan-Meier for 2nd group
group2<- Surv(time[cell.type=="2"],status[cell.type=="2"]) 
outp2<-survfit(group2~1, type="kaplan-meier", data=cdat)
#Extract values when status=2
Utime2<- outp2$time[outp2$n.event==1] 
SKM2<-outp2$surv[outp2$n.event==1] 

###Kaplan-Meier for 3rd group
group3<- Surv(time[cell.type=="3"],status[cell.type=="3"]) 
outp3<-survfit(group3~1, type="kaplan-meier", data=cdat)

#Extract values when status=3
Utime3<- outp3$time[outp3$n.event==1] 
SKM3<-outp3$surv[outp3$n.event==1] 

###Kaplan-Meier for 4th group
group4<- Surv(time[cell.type=="4"],status[cell.type=="4"]) 
outp4<-survfit(group4~1, type="kaplan-meier", data=cdat)
#Extract values when status=4
Utime4<- outp4$time[outp4$n.event==1] 
SKM4<-outp4$surv[outp4$n.event==1] 


Utime1
SKM1
Utime2
SKM2
Utime3
SKM3
Utime4
SKM4

# Remove the last element from Utime1 and SKM1
Utime1 <- head(Utime1, -1)
SKM1 <- head(SKM1, -1)

# Remove the last element from Utime2 and SKM2
Utime2 <- head(Utime2, -1)
SKM2 <- head(SKM2, -1)

# Remove the last element from Utime3 and SKM3
Utime3 <- head(Utime3, -1)
SKM3 <- head(SKM3, -1)

# Remove the last element from Utime4 and SKM4
Utime4 <- head(Utime4, -1)
SKM4 <- head(SKM4, -1)

#PH assumption plot
plot(log(Utime1),log(-log(SKM1)) , main=expression(paste("PH assumption")), xlab="ln t", 
ylab="ln-(ln S(t))", col="blue",pch=19) 
abline(lm(log(-log(SKM1))~ log(Utime1)), col="blue")
 
points(log(Utime2),log(-log(SKM2)),col="red",pch=19) 
abline(lm(log(-log(SKM2))~ log(Utime2)), col="red")

points(log(Utime3),log(-log(SKM3)),col="purple",pch=19) 
abline(lm(log(-log(SKM3))~ log(Utime3)), col="purple")

points(log(Utime4),log(-log(SKM4)),col="green",pch=19) 
abline(lm(log(-log(SKM4))~ log(Utime4)), col="green")
 
legend("bottomright", c("Cell Type 1", "Cell Type 2","Cell Type 3", "Cell Type 4"), col=c("blue","red","purple","green"), 
lty=1:4) 

#------------------------------------------------------------------------------------
#Examine the AL assumption 
plot(log(Utime1), SKM1 , main=expression(paste("AL assumption")), xlab="ln t", 
ylab="S(t)", col="blue",pch=19,lty = 1:4) 
points(log(Utime2), SKM2,col="red",pch=19)
points(log(Utime3), SKM3,col="purple",pch=19) 
points(log(Utime4), SKM4,col="green",pch=19) 
 
lines(log(Utime1), SKM1, col="blue") 
lines(log(Utime2), SKM2, col="red")
lines(log(Utime3), SKM3, col="purple") 
lines(log(Utime4), SKM4, col="green")  
legend("bottomleft", c("Cell Type 1", "Cell Type 2","Cell Type 3", "Cell Type 4"), col=c("blue","red","purple","green"), 
lty=1:4)


#----------------- END  of PART 2 --------------------------------------

##########################################################################

#ASKISI B
#--------------------------B1--------------------------
#Να προσαρμοστεί ένα μοντέλο παλινδρόμησης της κατανομής Weibull

mod <- survreg(Surv(time,status)~ treatment + factor(cell.type)+ karno.score + months + age + prior.therapy + id, data=cdat, dist="weibull") 
mod
summary(mod) 

#

step(mod, direction="backward", test="Chisq") 


#FINAL

modfinal <- survreg(Surv(time,status)~ factor(cell.type)+ karno.score, data=cdat, dist="weibull") 
modfinal 
summary(modfinal) 


#--------------------------B1 ENDS---------------------

#--------------------------B2--------------------------

#ΔΕ
confint.default(modfinal) 
exp(confint.default(modfinal)) 

#--------------------------B2 ENDS---------------------
#--------------------------B3--------------------------
#COX-SNELL

standres<-(log(time)-modfinal$linear.predictors)/modfinal$scale 
csres<-exp(standres) 
csres 

Cox.Snell<- csres [status == 1] 

library(EnvStats) 

qqPlot(Cox.Snell, distribution="exp", param.list = list(rate=1), add.line=TRUE, pch=19) 

#--------------------------B3 ENDS---------------------
#############################################################################

#EXERCISE C

#--------------------------C1--------------------------
mod <- survreg(Surv(time,status)~ treatment + factor(cell.type)+ karno.score + months + age + prior.therapy + id, data=cdat, dist="weibull") 
mod
mod2<- coxph(Surv(time,status)~treatment +factor(cell.type)+karno.score + months + age + prior.therapy + id,ties=c("breslow")) 
mod2
summary(mod2)

#Backward τεχνική με βήματα
mod3<-step(mod2, direction="backward", test="Chisq")
#--------------------------C1 ENDS---------------------

#--------------------------C2--------------------------

#Τελικό μοντέλο 
 
mod4<- coxph(Surv(time,status)~factor(cell.type)+karno.score , ties=c("breslow")) 
mod4 
summary(mod4)

#--------------------------C2 ENDS---------------------

#--------------------------C3--------------------------

#Κλιμακοποιημένα υπόλοιπα   Schoenfeld  για το τελικό μοντέλο γαι τις τρεις μεταβλητές  
#cell.type1, cell.type2, cell.type3 και cell.type4 , και karno.score    
 
sresid<-resid(mod4,type="scaledsch") 
sresid 

#Η  παρακάτω εντολή εξετάζει το typef σαν μια ενιαία μεταβλητή 
 
vv<-cox.zph(mod4,transform='identity', terms=FALSE) 
vv 
par(mfrow=c(2,2)) 
plot(vv) 
#--------------------------C3 ENDS---------------------

#--------------------------C4--------------------------
mod4$linear.predictor

library(risksetROC) 
eta<-mod4$linear.predictor 

ROC10=risksetROC(Stime= time, status=status, marker=eta, predict.time=60, method="Cox", main="ROC Curve", lty=2, col="red", ylab="True Positive", xlab="False Positive")

ROC50=risksetROC(Stime= time, status=status,marker=eta, predict.time=80, method="Cox", plot=FALSE) 

lines(ROC50$FP,ROC50$TP, lty=3,col="darkblue") 
legend(.6,.25,lty=c(2,3),col=c("red","darkblue"), 
legend=c("t=60","t=80"), bty="n") 
#--------------------------C4 ENDS----------------------
