
# last edited 12-05-17

# load packages
library(nlme) # 3.1-131
library(visreg) # 2.2-2
library(MuMIn) # 1.15.6

# make sure folder data is working directory
rh <- read.csv("traindisplaydataRH2.csv")
head(rh)
dim(rh);summary(rh)
names(rh)[10:11] <- c('dewpoint','relhumid')
head(rh)
rh$sqrtYmod <- sqrt(3.671-0.022*rh$relhumid) # sqare root of Y modulus
plot(jitter(sqrtYmod)~relhumid, rh, pch=16, cex=0.5, las=1, bty='l')
hist(rh$relhumid[!duplicated(rh$boutID)], breaks=30) 
hist(rh$temperatureC[!duplicated(rh$boutID)], breaks=30) 

# median conditions during peafowl vibration displays recorded in the field:
median(rh$temperatureC[!duplicated(rh$boutID)]) # air temperature
median(rh$relhumid[!duplicated(rh$boutID)]) # rel. huidity

# avg. lab conditions = 21.1C temperature and 74.8% rel. humidity
labT <- 21.1
labH <- 74.8

x <- seq(0,10,by=0.01)
for(i in 1:length(x)){
  mypercent <- sum(rh$temperatureC[!duplicated(rh$boutID)] < (labT+x[i]) & rh$temperatureC[!duplicated(rh$boutID)] > (labT-x[i]))/length(rh$temperatureC[!duplicated(rh$boutID)])
  if(mypercent < 0.5){
    next
  } else {
    print(c(mypercent*100, x[i]))
    break
  }
} # therefore, over 50% within 2.2ºC of lab

x <- seq(0,20,by=0.01)
for(i in 1:length(x)){
  mypercent <- sum(rh$relhumid[!duplicated(rh$boutID)] < (labH+x[i]) & rh$relhumid[!duplicated(rh$boutID)] > (labH-x[i]))/length(rh$relhumid[!duplicated(rh$boutID)])
  if(mypercent < 0.5){
    next
  } else {
    print(c(mypercent*100, x[i]))
    break
  }
} # over 50% within 13.8% of lab

plot(freq ~ relhumid, rh, pch=16, cex=0.5, las=1, bty='l', col=type)
legend('bottomleft', col=1:3, legend=levels(factor(rh$type)), pch=16, bty='n', cex=0.5)
rh$sample <- factor(rh$sample, levels=c('pre','peak','post'))
rh$displaydateDOY <- as.numeric(substr(rh$date, 3,4)) + 59

rh.mod <- lme(freq ~ sample + displaydateDOY + timeH + relhumid + trainL + traindate, random=~1|id, data=subset(rh, type='male'), na.action=na.omit, method='REML')
summary(rh.mod) # ns. relationship with relative humidity, accounting for date, time, and morphology
plot(rh.mod)
dev.off(); hist(residuals(rh.mod))
par(mfrow=c(3,3), mar=c(4,4,0.1,0.1), mgp=c(2.5,1,0)); visreg(rh.mod)
r.squaredGLMM(rh.mod)
rh.mod2 <- lme(freq ~ sample + displaydateDOY + timeH + trainL + traindate, random=~1|id, data=subset(rh, type='male'), na.action=na.omit, method='REML') # remove humidity term
r.squaredGLMM(rh.mod2) # virtually no change in variance explained when it is removed


