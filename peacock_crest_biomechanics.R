
# last edited 15-08-17

# load packages
library(nlme) # 1.1-12
library(visreg) # 2.2-2
library(MuMIn) # 1.15.6
library(lme4) # 3.1-131

(vib <- read.csv('./data/crest_vibration_data2.csv'))
(morph <- read.csv('./data/crest_samples.csv'))
morph$crest_number <- ifelse(nchar(as.character(morph$crest_number))<2, paste('0',as.character(morph$crest_number),sep=''), as.character(morph$crest_number)) # ID with consistent no. chars
morph$crest_number <- paste('crest', morph$crest_number, sep='_')

vib$feather_ID <- gsub('Crest', 'crest', as.character(vib$feather_ID)) # ID for feathers
vib$crest_number <- unlist(strsplit(as.character(vib$feather_ID), split='_'))[seq(2,480,by=4)]
vib$crest_number <- ifelse(nchar(as.character(vib$crest_number))<2, paste('0',as.character(vib$crest_number),sep=''), as.character(vib$crest_number))
vib$crest_number <- paste('crest', vib$crest_number, sep='_')
vib$run <- unlist(strsplit(as.character(vib$feather_ID), split='_'))[seq(4,480,by=4)]
vib$exp <- unlist(strsplit(as.character(vib$feather_ID), split='_'))[seq(3,480,by=4)] # Dan's experiment details. run is the sequence of repeated measures. exp gives the frequency sweep type for the whole crests

head(vib); dim(vib)
head(morph); dim(morph)
vib <- merge(vib, morph, by='crest_number'); dim(vib)
vib$run <- as.numeric(vib$run)

head(vib) # omit the 0-15 sweeps because why?
vib <- subset(vib, exp!='0-15')

# add unique identifier for single feathers
vib$single_id <- ifelse(vib$whole_single=='singl', paste(vib$crest_number, vib$exp, sep='.'), NA)

dim(vib) # 112 rows (not all complete)
vib$sweep <- factor(ifelse(vib$whole_single=='whole', vib$exp, NA))
summary(vib$sweep)

summary(vib)
vib$orientation <- factor(vib$orient, levels=c('standard','sideways')) # ask Dan to verify anatomical reference for the two orientations
vib$sex <- factor(vib$sex, levels=c('female','male'))
vib$sex_col <- ifelse(vib$sex=='female', 'green', 'blue')

# add quality factor, Q
vib$q <- vib$f_res/vib$del_f

# whole crest

# quality factor.
mod.q <- lme(q ~ orientation + sweep + sex + width + height + top_area + nfeathers + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.q)
# first we check run and sweep. both ns therefore remove those parameters
mod.q <- lme(q ~ orientation + sex + width + height + top_area + nfeathers, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.q)
dev.new(); plot(mod.q) # fine, no outliers
dev.new(); hist(residuals(mod.q)) # good

png(file='./figures/Qmod_resid_plot.png', width=4, height=5, res=300, units='in', bg='white')
par(mfrow=c(3,2), bty='l', mar=c(4,4,0.25,0.25), mgp=c(1.5,0.5,0)); visreg(mod.q, cond=list(orientation='standard'), ylab='Q', ylim=c(3,9))
dev.off()

r.squaredGLMM(mod.q) # 46% var explained in Q by fixed effects, very high
mod.q.rep <- lme(q ~ 1, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
r.squaredGLMM(mod.q.rep) # 50% measurement repeatability

# resonant freq.
mod.f <- lme(f_res ~ orientation + sweep + sex + width + height + top_area + nfeathers + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.f) # first we check run and sweep. both ns therefore eliminate. residuals indicate that variance differs among the crests, so we furthermore use a model with different variances:
mod.f <- lme(f_res ~ orientation + sex + width + height + top_area + nfeathers, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit, weights=varIdent(form=~1|crest_number))
summary(mod.f)
dev.new(); plot(mod.f) # good
dev.new(); hist(residuals(mod.f)) # ok

png(file='./figures/fmod_resid_plot.png', width=4, height=5, res=300, units='in', bg='white')
par(mfrow=c(3,2), bty='l', mar=c(4,4,0.25,0.25), mgp=c(1.5,0.5,0)); visreg(mod.f, cond=list(orientation='standard'), ylab='f_res', ylim=c(19,33))
dev.off()
# the more top area, the lower the f_res. b/c more massive, perhaps? confounded with sex
# ns. positive effect of n feathers, width. ns. negative effect of height
# higher in standard orientation

r.squaredGLMM(mod.f) # 34% ver explained in f res. by fixed effects, high
mod.f.rep <- lme(f_res ~ 1, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
r.squaredGLMM(mod.f.rep) # 82% measurement repeatability, very high


# plot data and some model fits for both fr and Q, for Figure 2 of manuscript
# add:
# male mean display freq = 25.6, range of bird means 23.9-27.1; range 22.2 to 28.2
# female mean display freq = 26.1, range of bird means 25.2-27.1; range 25.2 to 28.5

png(file='./figures/fig2_vib_results.png', width=7, height=2, res=300, units='in', bg='white')

par(mfrow=c(2, 6), bty='l', las=1, mar=c(3,3,0.25,0.25), mgp=c(1.5,0.5,0), cex.axis=0.5, cex.lab=0.5, tck=-0.05)

plot(f_res ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('standard','sideways'))
axis(2, at=c(20,25,30))
legend('topright', col=c('green','blue'), pch=16, legend=c('female','male'), bty='n', cex=0.5)
segments(x0=0.75,x1=1.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

plot(f_res ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('female','male'))
axis(2, at=c(20,25,30))
segments(x0=0.75,x1=1.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.f, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.f, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

polygon(x=c(0,5,5,0), y=c(22.2,22.2,28.5,28.5), col=rgb(0,0,0,0.2), border=NA) # range of frequencies recorded over all shaking displays 
segments(x0=c(0), x1=c(5), y0=c(25.6), y1=c(25.6), lty=3, col=rgb(0,0,1)) # peacock shaking display shaking avg
segments(x0=c(0), x1=c(5), y0=c(26.1), y1=c(26.1), lty=3, col=rgb(1,0,0)) # peahen shaking display avg

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
axis(1, at=c(1,2), labels=c('female','male'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=((predict(mod.q, newdata=data.frame(orientation=c('standard'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.q, newdata=data.frame(orientation=c('sideways'), sex='female',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)
segments(x0=1.75,x1=2.25, y0=((predict(mod.q, newdata=data.frame(orientation=c('standard'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)) + (predict(mod.q, newdata=data.frame(orientation=c('sideways'), sex='male',width=median(morph$width), height=median(morph$height), top_area=median(morph$top_area), nfeathers=median(morph$nfeathers, na.rm=T)), level=0)))/2)

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

dev.off()


