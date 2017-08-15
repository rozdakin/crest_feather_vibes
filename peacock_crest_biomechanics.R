


(vib <- read.csv('crest_vibration_data2.csv'))
(morph <- read.csv('crest_samples.csv'))
morph$crest_number <- ifelse(nchar(as.character(morph$crest_number))<2, paste('0',as.character(morph$crest_number),sep=''), as.character(morph$crest_number))
morph$crest_number <- paste('crest', morph$crest_number, sep='_')

vib$feather_ID <- gsub('Crest', 'crest', as.character(vib$feather_ID))
vib$crest_number <- unlist(strsplit(as.character(vib$feather_ID), split='_'))[seq(2,480,by=4)]
vib$crest_number <- ifelse(nchar(as.character(vib$crest_number))<2, paste('0',as.character(vib$crest_number),sep=''), as.character(vib$crest_number))
vib$crest_number <- paste('crest', vib$crest_number, sep='_')
vib$run <- unlist(strsplit(as.character(vib$feather_ID), split='_'))[seq(4,480,by=4)]
vib$exp <- unlist(strsplit(as.character(vib$feather_ID), split='_'))[seq(3,480,by=4)]

library(nlme)
library(visreg)
library(MuMIn)
library(lme4)

head(vib); dim(vib)
head(morph); dim(morph)
vib <- merge(vib, morph, by='crest_number'); dim(vib)
vib$run <- as.numeric(vib$run)

head(vib) # omit the 0-15 sweeps
vib <- subset(vib, exp!='0-15')

# unique identifier for single feathers
vib$single_id <- ifelse(vib$whole_single=='singl', paste(vib$crest_number, vib$exp, sep='.'), NA)

dim(vib) # 112 rows (not all complete)
vib$sweep <- factor(ifelse(vib$whole_single=='whole', vib$exp, NA))
summary(vib$sweep)

summary(vib)
vib$orientation <- factor(vib$orient, levels=c('standard','sideways'))
vib$sex <- factor(vib$sex, levels=c('male','female'))
vib$sex_col <- ifelse(vib$sex=='female', 'green', 'blue')

# add quality factor, Q
vib$q <- vib$f_res/vib$del_f

# whole crest

# quality factor.
mod.q <- lme(q ~ orientation + sweep + sex + width + height + top_area + nfeathers + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
# first we check run and sweep. both ns therefore eliminate.
mod.q <- lme(q ~ orientation + sex + width + height + top_area + nfeathers, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.q)
dev.new(); plot(mod.q) # fine, no outliers
dev.new(); hist(residuals(mod.q)) # good
dev.new(width=4,height=6); par(mfrow=c(3,2), bty='l'); visreg(mod.q, cond=list(orientation='standard'), ylab='Q', ylim=c(3,9))

r.squaredGLMM(mod.q) # 46% exp lained, very high
mod.q.rep <- lme(q ~ 1, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
r.squaredGLMM(mod.q.rep) # 50% repeatability

# resonant freq.
mod.f <- lme(f_res ~ orientation + sweep + sex + width + height + top_area + nfeathers + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
# first we check run and sweep. both ns therefore eliminate.
mod.f <- lme(f_res ~ orientation + sex + width + height + top_area + nfeathers, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit, weights=varIdent(form=~1|crest_number))
summary(mod.f)
dev.new(); plot(mod.f) # some outliers, but better in unequal variance model
dev.new(); hist(residuals(mod.f)) # ok
dev.new(width=4,height=6); par(mfrow=c(3,2), bty='l'); visreg(mod.f, cond=list(orientation='standard'), ylab='f_res', ylim=c(19,33))
# the more top area, the lower the f_res. b/c more massive, perhaps? confounded with sex
# ns. positive effect of n feathers, width. ns. negative effect of height
# higher in standard orientation

r.squaredGLMM(mod.f) # 34% exp lained, very high
mod.f.rep <- lme(f_res ~ 1, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
r.squaredGLMM(mod.f.rep) # 82% repeatability


# plot data and model fits for fr and Q

dev.new(width=7, height=2)
par(mfrow=c(2, 6), bty='l', las=1, mar=c(3,3,0.25,0.25), mgp=c(1.5,0.5,0), cex.axis=0.5, cex.lab=0.5, tck=-0.05)

plot(f_res ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('standard','sideways'))
axis(2, at=c(20,25,30))
legend('topright', col=c('green','blue'), pch=16, legend=c('female','male'), bty='n', cex=0.5)
segments(x0=0.75,x1=1.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
polygon(x=c(0,5,5,0), y=c(23.9,23.9,27.1,27.1), col=rgb(0,0,1,0.2), border=NA)
polygon(x=c(0,5,5,0), y=c(25.2,25.2,27.1,27.1), col=rgb(1,0,0,0.2), border=NA)
segments(x0=c(0), x1=c(5), y0=c(25.6), y1=c(25.6), lty=3, col=rgb(0,0,1))
segments(x0=c(0), x1=c(5), y0=c(26.1), y1=c(26.1), lty=3, col=rgb(1,0,0))

plot(f_res ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('male','female'))
axis(2, at=c(20,25,30))
segments(x0=0.75,x1=1.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(f_res ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='width (cm)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4,5.5,7))

plot(f_res ~ jitter(height,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='height (cm)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4.5,5.5,6.5))

plot(f_res ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='# feathers', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(20,25,30))

plot(f_res ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='top area (cm2)', ylab='f (Hz)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4,6.5,9))
segments(x0=c(4), x1=c(9), y0=predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=4, nfeathers=median(morph$nfeathers, na.rm=T)), level=0), y1=predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=9, nfeathers=median(morph$nfeathers, na.rm=T)), level=0), lty=3)
# sex effect is partially due to top area? one or the other is significant

plot(q ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(2,9), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('standard','sideways'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=((predict(mod.q, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.q, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.q, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.q, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(q ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(2,9), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('male','female'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=((predict(mod.q, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.q, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.q, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.q, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(q ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='width (cm)', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4,5.5,7))
plot(q ~ jitter(height,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='height (cm)', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4.5,5.5,6.5))
plot(q ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='# feathers', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(20,25,30))
plot(q ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='top area (cm2)', ylab='Q', ylim=c(2,9), cex=0.5, yaxt='n',xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4,6.5,9))

# male mean 25.6, range of bird means 23.9-27.1
# female mean 26.1, range of bird means 25.2-27.1

# plot (old)

dev.new(width=4,height=6)
par(mfrow=c(3,2), bty='l', las=1)
plot(f_res ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5)
axis(1, at=c(1,2), labels=c('standard','sideways'))
legend('topright', col=c('green','blue'), pch=16, legend=c('female','male'), bty='n', cex=0.5)
segments(x0=0.75,x1=1.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(f_res ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='top area (cm2)', ylab='f (Hz)', ylim=c(19,33), cex=0.5)
segments(x0=c(4), x1=c(9), y0=predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=4, nfeathers=median(morph$nfeathers, na.rm=T)), level=0), y1=predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=9, nfeathers=median(morph$nfeathers, na.rm=T)), level=0))
# sex effect is partially due to top area? one or the other is significant

plot(f_res ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5)
axis(1, at=c(1,2), labels=c('male','female'))
segments(x0=0.75,x1=1.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(f_res ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='width (cm)', ylim=c(19,33), cex=0.5)

plot(f_res ~ jitter(height,2), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='height (cm)', ylim=c(19,33), cex=0.5)

plot(f_res ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='# feathers', ylim=c(19,33), cex=0.5)



# FWHM
mod.FWHM <- lme(del_f ~ orientation + sweep + sex + width + height + top_area + nfeathers + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
# first we check run and sweep. both ns therefore eliminate.
mod.FWHM <- lme(del_f ~ orientation + sex + width + height + top_area + nfeathers, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.FWHM)
plot(mod.FWHM) 
hist(residuals(mod.FWHM))
shapiro.test(residuals(mod.FWHM)) # not fixed by varIdent
dev.new(width=3,height=6); par(mfrow=c(4,2), bty='l'); visreg(mod.FWHM, ylab='FWHM', ylim=c(1,9))
# male crests have lower FWHM; also lower in standard orientation
# do we care about FWHM for sideways orient? do we need to account for different variances?

r.squaredGLMM(mod.FWHM) # model explains 31% of the variance
mod.FWHM.rep <- lme(del_f ~ 1, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
r.squaredGLMM(mod.FWHM.rep) # rep is 62%

# plot and predictions

dev.new(width=4,height=6)
par(mfrow=c(3,2), bty='l', las=1)
plot(del_f ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='FWHM', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(1,9), cex=0.5)
axis(1, at=c(1,2), labels=c('standard','sideways'))
legend('bottomright', col=c('green','blue'), pch=16, legend=c('female','male'), bty='n', cex=0.5)
segments(x0=0.75,x1=1.25, y0=((predict(mod.FWHM, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.FWHM, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.FWHM, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.FWHM, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(del_f ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylab='FWHM', ylim=c(1,9), cex=0.5)
axis(1, at=c(1,2), labels=c('male','female'))
segments(x0=0.75,x1=1.25, y0=((predict(mod.FWHM, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.FWHM, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.FWHM, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.FWHM, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
plot(del_f ~ width, subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='width (cm)', ylab='FWHM', ylim=c(1,9), cex=0.5)
plot(del_f ~ height, subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='height (cm)', ylab='FWHM', ylim=c(1,9), cex=0.5)
plot(del_f ~ top_area, subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='top area (cm2)', ylab='FWHM', ylim=c(1,9), cex=0.5)
plot(del_f ~ nfeathers, subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='# feathers', ylab='FWHM', ylim=c(1,9), cex=0.5)




# single feather - not done

# freq.
mod.f <- lme(f_res ~ sex + rachis_l + run, random=~1|single_id, data=subset(vib, whole_single=='singl'))
summary(mod.f)
plot(mod.f) # note 0-15 data have outliers; suspect not supposed to include
par(mfrow=c(3,2)); visreg(mod.f) # the more area, the lower the f_res. does that make sense? more massive?
hist(residuals(mod.f))
shapiro.test(residuals(mod.f)) # outlier(s)?
r.squaredGLMM(mod.f) # >99% repeatability






