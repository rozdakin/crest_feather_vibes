
# last edited 09-13-17

# load packages
library(nlme) # 1.1-12
library(visreg) # 2.2-2
library(MuMIn) # 1.15.6
library(lme4) # 3.1-131
library(dplyr) # 0.5.0

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
vib$orientation <- factor(vib$orient, levels=c('standard','sideways')) # change to anatomical reference for the two orientations
vib$orientation <- factor(ifelse(vib$orient=='standard','A_out-plane','B_in-plane'))
vib$sex <- factor(vib$sex, levels=c('female','male'))
vib$sex_col <- ifelse(vib$sex=='female', 'green', 'blue')

# add quality factor, Q
vib$q <- vib$f_res/vib$del_f

# whole crest

# quality factor.
mod.q <- lme(q ~ orientation + sex + sweep + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.q)
# first we check run and sweep. both ns therefore remove those parameters
# next we assess morphology. with only 15 crests, we consider only one trait at a time
mod.q <- lme(q ~ orientation + sex, random=~1|crest_number, data=subset(vib, whole_single=='whole'&!is.na(nfeathers)), method='ML') # no trait model
mod.q1 <- update(mod.q, .~.+top_area)
mod.q2 <- update(mod.q, .~.+width)
mod.q3 <- update(mod.q, .~.+height)
mod.q4 <- update(mod.q, .~.+nfeathers)
mod.q5 <- update(mod.q, .~.+percent_unalign)
mod.q6 <- update(mod.q, .~.+percent_short)
AICc(mod.q1, mod.q2, mod.q3, mod.q4, mod.q5, mod.q6) # mod.q1 is lowest
mod.q1 <- update(mod.q1, method='REML')
summary(mod.q1)
dev.new(); plot(mod.q1) # fine, no outliers
dev.new(); hist(residuals(mod.q1)) # good

png(file='./figures/Qmod_resid_plot.png', width=4, height=5, res=300, units='in', bg='white')
par(mfrow=c(3,2), bty='l', mar=c(4,4,0.25,0.25), mgp=c(1.5,0.5,0)); visreg(mod.q1, cond=list(orientation='A_out-plane'), ylab='Q', ylim=c(3,9))
dev.off()

r.squaredGLMM(mod.q1) # 49% var explained in Q by fixed effects, fairly high
vcomp <- as.numeric(VarCorr(mod.q1)[,1])
vcomp[1]/(vcomp[1]+vcomp[2]) # 47% measurement repeatability

vib <- group_by(vib, crest_number, sex_col, orientation)
summarize(group_by(summarize(subset(vib, whole_single=='whole'), q=mean(q)), sex_col, orientation), mean(q))

# resonant freq.
mod.f <- lme(f_res ~ orientation + sex + sweep + run, random=~1|crest_number, data=subset(vib, whole_single=='whole'), na.action=na.omit)
summary(mod.f) # first we check run and sweep. both ns therefore eliminate. residuals indicate that variance differs among the crests, so we furthermore use a model with different variances:
mod.f <- lme(f_res ~ orientation + sex, random=~1|crest_number, data=subset(vib, whole_single=='whole'&!is.na(nfeathers)), weights=varIdent(form=~1|crest_number), method='ML') # no trait model
mod.f1 <- update(mod.f, .~.+top_area)
mod.f2 <- update(mod.f, .~.+width)
mod.f3 <- update(mod.f, .~.+height)
mod.f4 <- update(mod.f, .~.+nfeathers)
mod.f5 <- update(mod.f, .~.+percent_unalign)
mod.f6 <- update(mod.f, .~.+percent_short)
AICc(mod.f1, mod.f2, mod.f3, mod.f4, mod.f5, mod.f6) # mod.f1 is lowest
mod.f1 <- update(mod.f1, method='REML')
summary(mod.f1)
dev.new(); plot(mod.f1) # good
dev.new(); hist(residuals(mod.f1)) # ok

png(file='./figures/fmod_resid_plot.png', width=4, height=5, res=300, units='in', bg='white')
par(mfrow=c(3,2), bty='l', mar=c(4,4,0.25,0.25), mgp=c(1.5,0.5,0)); visreg(mod.f1, cond=list(orientation='A_out-plane'), ylab='f_res', ylim=c(19,33))
dev.off()
# higher in standard orientation, larger top area, tend to have lower res. freq

r.squaredGLMM(mod.f1) # 28% var explained in f res. by fixed effects
# repeatability
vcomp <- as.numeric(VarCorr(mod.f1)[,1])
vcomp[1]/(vcomp[1]+vcomp[2]) # 94% repeatability, very high

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
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='A_out-plane')$f_res))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='B_in-plane')$f_res))

plot(f_res ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('female','male'))
axis(2, at=c(20,25,30))
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&sex=='female')$f_res))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&sex=='male')$f_res))

polygon(x=c(0,5,5,0), y=c(22.2,22.2,28.5,28.5), col=rgb(0,0,0,0.2), border=NA) # range of freq for all shaking displays 
segments(x0=c(0), x1=c(5), y0=c(25.6), y1=c(25.6), lty=3, col=rgb(0,0,1)) # peacock shaking display shaking avg
segments(x0=c(0), x1=c(5), y0=c(26.1), y1=c(26.1), lty=3, col=rgb(1,0,0)) # peahen shaking display avg

plot(f_res ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='top area (cm2)', ylab='f (Hz)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4,6.5,9))
# sex effect is partially due to top area?

plot(f_res ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='width (cm)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4,5.5,7))

plot(f_res ~ jitter(height,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='height (cm)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4.5,5.5,6.5))

plot(f_res ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='f (Hz)', xlab='# feathers', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(20,25,30))

plot(q ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(2,9), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('standard','sideways'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='A_out-plane')$q))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='B_in-plane')$q))

plot(q ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(2,9), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('female','male'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&sex=='female')$q))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&sex=='male')$q))

plot(q ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', xlab='top area (cm2)', ylab='Q', ylim=c(2,9), cex=0.5, yaxt='n',xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4,6.5,9))
plot(q ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='width (cm)', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4,5.5,7))
plot(q ~ jitter(height,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='height (cm)', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4.5,5.5,6.5))
plot(q ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=16, col=sex_col, bty='l', ylab='Q', xlab='# feathers', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(20,25,30))

dev.off()

# plot morphology of dried crests, as compared to measurements on live birds from Dakin (2011)
# Dakin data:
# lengths: male 5.68 [5.52, 5.81] cm; female 5.69 [5.63, 5.75] cm.
# widths: Male 7.77 [7.43, 8.11] cm; female 6.62 [6.36, 6.98] cm  
morph$sex <- factor(ifelse(morph$color=='blue','b_male','a_female'))

png(file='./figures/fig1_morph_compare.png', width=4, height=2, res=300, units='in', bg='white')
par(mfrow=c(1,2), bty='l', las=1, mar=c(3,3,0.25,0.25), mgp=c(1.5,0.5,0), cex.axis=0.5, cex.lab=0.5, tck=-0.05, cex=0.5)
plot(c(5.69,5.68) ~ c(1,2), pch=16, ylim=c(3.25,8.25), xlim=c(0.5,2.5), ylab='length (cm)', xaxt='n', xlab='', col=c('green','blue'))
points(height ~ c(as.numeric(sex)-0.25), pch=1, col=as.numeric(sex)+2, data=morph)
axis(1, at=c(1,2), labels=c('female','male'))
segments(x0=1:2, x1=c(1:2), y0=c(5.63,5.52), y1=c(5.75,5.81))
segments(x0=c(0.5,1.5), x1=c(1,2), y0=c(mean(subset(morph,sex=='a_female')$height), mean(subset(morph,sex=='b_male')$height)))
plot(c(6.62,7.77) ~ c(1,2), pch=16, ylim=c(3.25,8.25), xlim=c(0.5,2.5), ylab='width (cm)', xaxt='n', xlab='', col=c('green','blue'))
points(width ~ c(as.numeric(sex)-0.25), pch=1, col=as.numeric(sex)+2, data=morph)
axis(1, at=c(1,2), labels=c('female','male'))
segments(x0=1:2, x1=c(1:2), y0=c(6.36,7.43), y1=c(6.98,8.11))
segments(x0=c(0.5,1.5), x1=c(1,2), y0=c(mean(subset(morph,sex=='a_female')$width), mean(subset(morph,sex=='b_male')$width)))
dev.off()

# plot example vibrational spectrum and fit
spec <- read.csv('data/Crest_transfer_function_plot_-_Sheet1.csv')[,1:3] # crest 1
names(spec) <- c('freq', 'transferfcn', 'fit')

png(file='./figures/fig2_sample_spectrum.png', width=4, height=4, res=300, units='in', bg='white')
par(bty='l', las=1)
plot(transferfcn ~ freq, spec, type='l', xlab='Frequency (Hz)', ylab='Transfer function (arb. units)', lty=3, ylim=c(0,8))
points(fit ~ freq, spec, type='l')
legend('topright', lty=c(3,1), legend=c('data','L. fit'), bty='n', cex=0.75)
text(x=spec[spec$fit==max(spec$fit),]$freq, y=spec[spec$fit==max(spec$fit),]$fit, pos=4, 'fr = 25.5', cex=0.75)
dev.off()

# plot example audio results
audio <- read.csv('data/Audio_file_sample_data_-_Sheet1.csv')[,1:2] # crest 13, 25Hz, 30 cm distance
names(audio) <- c('freq','power')
audiofit <- read.csv('data/Audio_file_sample_data_-_Sheet1.csv')[,4:5] # fit
names(audiofit) <- c('freq','fit')

png(file='./figures/fig3_audio_result.png', width=4, height=4, res=300, units='in', bg='white')
par(bty='l', las=1, mgp=c(3,0.5,0))
plot(power ~ freq, data=audio, type='l', ylab='Power as MSA (mm2)', xlab='Frequency (Hz)', lty=3, xlim=c(0,110))
points(fit ~ freq, data=audiofit, type='l')
legend('topright', lty=c(3,1), legend=c('data','L. fit'), bty='n', cex=0.75)
dev.off()


head.csv('data/Vortex_response_figure_-_Sheet1.csv')




