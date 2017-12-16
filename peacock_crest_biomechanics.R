
# last edited 12-05-17

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

head(vib) # omit the 0-15 sweeps
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

par(mfrow=c(3,2), bty='l', mar=c(4,4,0.25,0.25), mgp=c(1.5,0.5,0)); visreg(mod.q1, cond=list(orientation='A_out-plane'), ylab='Q', ylim=c(3,9))
dev.off()

r.squaredGLMM(mod.q1) # 49% var explained in Q by fixed effects, fairly high
vcomp <- as.numeric(VarCorr(mod.q1)[,1])
vcomp[1]/(vcomp[1]+vcomp[2]) # 47% measurement repeatability

vib <- group_by(vib, crest_number, sex_col, orientation)
summarize(group_by(summarize(subset(vib, whole_single=='whole'), q=mean(q)), sex_col, orientation), mean(q))
range(summarize(group_by(subset(vib, whole_single=='singl'), feather_ID), f_res=mean(f_res))$f_res)

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

# higher in standard orientation, larger top area, tend to have lower res. freq

r.squaredGLMM(mod.f1) # 28% var explained in f res. by fixed effects

# repeatability
vcomp <- as.numeric(VarCorr(mod.f1)[,1])
vcomp[1]/(vcomp[1]+vcomp[2]) # 94% repeatability, very high

# plot data and some model fits for both fr and Q, for Figure 2 of manuscript
# add:
# male mean display freq = 25.6, range of bird means 23.9-27.1; range 22.2 to 28.2
# female mean display freq = 26.1, range of bird means 25.2-27.1; range 25.2 to 28.5

mean(vib$f_upper - vib$f_lower, na.rm=T) # the average span of error for f is 0.07 Hz, too small to plot
vib$q_lower <- vib$f_lower/vib$del_upper
vib$q_upper <- vib$f_upper/vib$del_lower
mean(vib$q_upper - vib$q_lower, na.rm=T) # the average span of error for Q is 0.23

# add symbols to denote unique crests. there are 7 male and 8 female crests
# add legend as well. use pch=1:8
lookup <- data.frame(table(vib$crest_number, vib$sex))
lookup <- lookup[lookup$Freq>0,-3]
lookup$pch <- c(1:8,1:7)
vib$pch <- lookup$pch[match(vib$crest_number, lookup$Var1)]

png(file='./figures/fig2_vib_results.png', width=7, height=2, res=300, units='in', bg='white')

par(mfrow=c(2, 6), bty='l', las=1, mar=c(3,3,0.25,0.25), mgp=c(1.5,0.5,0), cex.axis=0.5, cex.lab=0.5, tck=-0.05)

plot(f_res ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='f (Hz)', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('standard','sideways'))
axis(2, at=c(20,25,30))
legend('topright', col=c('green','blue'), pch=16, legend=c('female','male'), bty='n', cex=0.5)
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='A_out-plane')$f_res))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='B_in-plane')$f_res))
points(f_res ~ jitter(rep(1.5,33), 1.5), subset(vib, whole_single=='singl')[,c('f_res','orientation','sex_col','pch')], cex=0.5, pch=pch, col=sex_col)

plot(f_res ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='f (Hz)', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(19,33), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('female','male'))
axis(2, at=c(20,25,30))
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&sex=='female')$f_res))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&sex=='male')$f_res))

polygon(x=c(0,5,5,0), y=c(22.2,22.2,28.5,28.5), col=rgb(0,0,0,0.2), border=NA) # range of freq for all shaking displays 
segments(x0=c(0), x1=c(5), y0=c(25.6), y1=c(25.6), lty=3, col=rgb(0,0,1)) # peacock shaking display shaking avg
segments(x0=c(0), x1=c(5), y0=c(26.1), y1=c(26.1), lty=3, col=rgb(1,0,0)) # peahen shaking display avg

plot(f_res ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', xlab='top area (cm2)', ylab='f (Hz)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4,6.5,9))
# sex effect is partially due to top area?

plot(f_res ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='f (Hz)', xlab='width (cm)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4,5.5,7))

plot(f_res ~ jitter(height,1), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='f (Hz)', xlab='height (cm)', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(4.5,5.5,6.5))

plot(f_res ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='f (Hz)', xlab='# feathers', ylim=c(19,33), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(20,25,30))
axis(1, at=c(20,25,30))

plot(q ~ jitter(as.numeric(orientation), 0.5), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='Q', xlab='orientation', xaxt='n', xlim=c(0.5,2.5), ylim=c(2,9), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('standard','sideways'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='A_out-plane')$q))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&orientation=='B_in-plane')$q))

plot(q ~ jitter(as.numeric(sex), 0.5), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='Q', xlab='sex', xaxt='n', xlim=c(0.5,2.5), ylim=c(2,9), cex=0.5, yaxt='n')
axis(1, at=c(1,2), labels=c('female','male'))
axis(2, at=c(3,6,9))
segments(x0=0.75,x1=1.25, y0=mean(subset(vib, whole_single=='whole'&sex=='female')$q))
segments(x0=1.75,x1=2.25, y0=mean(subset(vib, whole_single=='whole'&sex=='male')$q))

plot(q ~ jitter(top_area, 3), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', xlab='top area (cm2)', ylab='Q', ylim=c(2,9), cex=0.5, yaxt='n',xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4,6.5,9))
plot(q ~ jitter(width,2), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='Q', xlab='width (cm)', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4,5.5,7))
plot(q ~ jitter(height,1), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='Q', xlab='height (cm)', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(4.5,5.5,6.5))
plot(q ~ jitter(nfeathers,1), subset(vib, whole_single=='whole'), pch=pch, col=sex_col, bty='l', ylab='Q', xlab='# feathers', ylim=c(2,9), cex=0.5, yaxt='n', xaxt='n')
axis(2, at=c(3,6,9))
axis(1, at=c(20,25,30))

dev.off()

png(file='./figures/fig2_legend.png', width=7, height=2, res=300, units='in', bg='white') # legend for crest IDs
par(mfrow=c(2, 6), bty='l', las=1, par(mar=c(0,0,0,0)), mgp=c(1.5,0.5,0))
plot.new()
legend('topleft', legend=lookup$Var1, pch=lookup$pch, bty='n', col=c('green','blue')[c(rep(1,8),rep(2,7))], cex=0.4)
dev.off()

# plot morphology of dried crests, as compared to measurements on live birds from Dakin (2011)
# Dakin data:
# lengths: male 5.68 [5.52, 5.81] cm; female 5.69 [5.63, 5.75] cm.
# widths: Male 7.77 [7.43, 8.11] cm; female 6.62 [6.36, 6.98] cm  
morph$sex <- factor(ifelse(morph$color=='blue','b_male','a_female'))
morph <- group_by(morph, sex)
morphsumm <- summarize(morph, L=mean(height), L_SE=sd(height)/sqrt(n()), L_lower=L-1.96*L_SE, L_upper=L+1.96*L_SE, W=mean(width), W_SE=sd(width)/sqrt(n()), W_lower=W-1.96*W_SE, W_upper=W+1.96*W_SE)

morph$pch <- lookup$pch[match(morph$crest_number, lookup$Var1)]

png(file='./figures/fig1_morph_compare.png', width=4, height=2, res=300, units='in', bg='white')
par(mfrow=c(1,2), bty='l', las=1, mar=c(3,3,0.25,0.25), mgp=c(1.5,0.5,0), cex.axis=0.5, cex.lab=0.5, tck=-0.05, cex=0.5)
plot(c(5.69,5.68) ~ c(1,2), pch=16, ylim=c(3.25,8.25), xlim=c(0.5,2.5), ylab='length (cm)', xaxt='n', xlab='', col=c('green','blue'))
points(height ~ c(as.numeric(sex)-0.25), pch=pch, col=as.numeric(sex)+2, data=morph)
axis(1, at=c(1,2), labels=c('female','male'))
segments(x0=1:2, x1=c(1:2), y0=c(5.63,5.52), y1=c(5.75,5.81))
segments(x0=c(0.5,1.5), x1=c(1,2), y0=c(mean(subset(morph,sex=='a_female')$height), mean(subset(morph,sex=='b_male')$height)))
segments(x0=c(0.75,1.75), x1=c(0.75,1.75), y0=c(morphsumm$L_lower), y1=c(morphsumm$L_upper))
segments(x0=c(0.95,1.95), x1=c(1.05,2.05), y0=c(5.69,5.68), y1=c(5.69,5.68))
segments(x0=0.75, x1=0.75, y0=max(morph$height+0.1), y1=max(morph$height-0.1))
plot(c(6.62,7.77) ~ c(1,2), pch=16, ylim=c(3.25,8.25), xlim=c(0.5,2.5), ylab='width (cm)', xaxt='n', xlab='', col=c('green','blue'))
points(width ~ c(as.numeric(sex)-0.25), pch=pch, col=as.numeric(sex)+2, data=morph)
axis(1, at=c(1,2), labels=c('female','male'))
segments(x0=1:2, x1=c(1:2), y0=c(6.36,7.43), y1=c(6.98,8.11))
segments(x0=c(0.5,1.5), x1=c(1,2), y0=c(mean(subset(morph,sex=='a_female')$width), mean(subset(morph,sex=='b_male')$width)))
segments(x0=c(0.75,1.75), x1=c(0.75,1.75), y0=c(morphsumm$W_lower), y1=c(morphsumm$W_upper))
segments(x0=c(0.95,1.95), x1=c(1.05,2.05), y0=c(6.62,7.77), y1=c(6.62,7.77))
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

# plot wingflap experiment resuls
flap <- read.csv('data/wingflap_power_spectrum.csv')
head(flap)

png(file='./figures/fig4_wingflap_result.png', width=4, height=4, res=300, units='in', bg='white')
par(mfrow=c(1,1), bty='l', las=1, mgp=c(3,0.5,0))
plot(FFTpower ~ frequency, type='l', flap, xlim=c(0,30), yaxt='n', ylim=c(0,0.000002), xlab='Frequency (Hz)', ylab='FFT spectral power (arb. units)')
axis(2, at=c(0,0.000001,0.000002))
dev.off()

# plot vortex gun results
vortex <- read.csv('data/Vortex_response_figure_-_Sheet1.csv')[,1:2]
names(vortex) <- c('t','disp')
vortexfit <- read.csv('data/Vortex_response_figure_-_Sheet1.csv')[,4:5]
names(vortexfit) <- c('t','fit')
vortexall <- read.csv('data/vortex_results_natural_frequency.csv')

png(file='./figures/fig5_vortex.png', width=8, height=4, res=300, units='in', bg='white')
par(mfrow=c(1,2), mar=c(4,4,0.25,0.25), bty='l', las=1, mgp=c(2.5,0.5,0))
plot(disp ~ t, data=vortex, type='l', ylab='Displacement (mm)', xlab='Time (s)', lty=3)
points(fit/5.9 ~ t, data=vortexfit, type='l')
legend('topright', lty=c(3,1), legend=c('data','fit'), bty='n', cex=0.75)
plot(vortexall$f_res, pch=16, ylim=c(20,30), ylab='Frequency (Hz)', xlab='Crest', xaxt='n', xlim=c(0.9,3.5))
segments(x0=1:3, y0=vortexall$res_lower, y1=vortexall$res_upper)
points(vortexall$f_vortex_response ~c(1:3+0.25), pch=1)
segments(x0=1:3+0.25, y0=vortexall$vortex_lower, y1=vortexall$vortex_upper)
axis(1, at=1:3, labels=c('A','B','C'))
legend('bottomleft', pch=c(16,1), legend=c('fr','vortex response'), bty='n', cex=0.75)
dev.off()

# plot mechanical properties results
mech <- read.csv('./data/mechanical_crest_properties.csv') # note the # of feathers here is sometimes less than in the morph dataset because of removal
mech$color <- ifelse(mech$sex=='F','green','blue')
mech$pch <- c(15,16,15,17,16,17)[factor(mech$crestID)]
set.seed(1); mech$jitt <- as.numeric(mech$sex) + rnorm(18, 0, 0.2)
mech <- group_by(mech, sex, crestID)
mechsumm <- summarize(group_by(summarize(mech, k=mean(k_Nmm)), sex), k=mean(k))
range(mech$k_Nmm)

mechmod <- lme(k_Nmm ~ sex, random=~1|crestID, data=mech)
summary(mechmod)
as.numeric(VarCorr(mechmod)[,1])[1]/(as.numeric(VarCorr(mechmod)[,1])[1] + as.numeric(VarCorr(mechmod)[,1])[2])

mech$crestID <- ifelse(nchar(mech$crestID)<2, paste('crest_0', mech$crestID, sep=''), paste('crest_',mech$crestID, sep=''))
mech$pch <- vib$pch[match(mech$crestID, vib$crest_number)]
mech$crestID <- factor(mech$crestID)

mech2 <- read.csv('./data/Crest_static_force_vs_displacement.csv')
head(mech2)
range(mech2$displacement_mm)
range(mech2$force_N)
mech2$series <- paste('Crest', ifelse(nchar(mech2$crestID)<2, paste('0', mech2$crestID, sep=''), mech2$crestID), mech2$trial_number, sep='_')

png(file='./figures/fig6_mech.png', width=6, height=3, res=300, units='in', bg='white')
par(mar=c(4,4,0.25,0.25), mfrow=c(1,2), bty='l', las=1, mgp=c(2.5,0.5,0))
i=9
temp <- subset(mech2, crestID==i)
plot(force_N~displacement_mm, subset(temp, trial_number==1), type='n', ylim=c(0,0.063), xlim=c(0,13.5), main=paste('Crest_',i))
points(predict(lm(force_N~0+displacement_mm, subset(temp, trial_number==1))) ~ subset(temp, trial_number==1)$displacement_mm, type='l')
segments(x0=subset(temp, trial_number==1)$displacement_mm, y0=subset(temp, trial_number==1)$force_N-0.001, y1=subset(temp, trial_number==1)$force_N+0.001)
points(force_N~displacement_mm, subset(temp, trial_number==1), cex=1, type='p', col=as.character(color), pch=7)
plot(k_Nmm ~ as.numeric(crestID), data=mech, col=color, pch=pch, xlim=c(0.5,6.5), xaxt='n', xlab='Crest ID', cex=1, ylim=c(0.0015, 0.007))
# segments(x0=c(0.75,1.75), x1=c(1.25,2.25), y0=mechsumm$k, lwd=2)
# axis(1, at=1:2, labels=c('Female','Male'))
segments(x0=as.numeric(mech$crestID), y0=c(mech$k_Nmm+mech$SE), y1=c(mech$k_Nmm-mech$SE), col=mech$color)
dev.off()





