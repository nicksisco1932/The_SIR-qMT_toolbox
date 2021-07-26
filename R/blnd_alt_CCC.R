# Title     : A Script to Analyze SIR MC
# Objective : Analysis
# Created by: Nick Sisco
# Created on: 4/23/2021

library(blandr)
library(ggplot2)
library(oro.nifti)
library(reshape2)
library(BlandAltmanLeh)


## PSR
# a<-bland.altman.plot(PSR_mdata[,-1],PSRd_mdata[,-1],graph.sys = 'ggplot2',conf.int=0.95)
# print(a+xlab('Mean (PSR and PSR Default)') +
#         ggtitle("BA: PSR Forward Diff vs Default") +ylab('PSR-PSR(default)')+theme_classic())

known <- readNIfTI('./sim_knowns_MC.nii.gz',reorient = FALSE)
fit <- readNIfTI('./PSR_julia_MC.nii.gz',reorient = FALSE)

# known <- readNIfTI('/Users/nicks/Documents/MRI_data/DTI_VSI/PT16/reg2sage/SAGE_R2s.nii.gz',reorient = FALSE)
# fit <- readNIfTI('/Users/nicks/Documents/MRI_data/DTI_VSI/PT16/reg2sage/SAGE_R2s_fit.nii.gz',reorient = FALSE)

library(epiR)

tmp <- data.frame(m1, m2)
tmp.ccc <- epi.ccc(method1, method2, ci = "z-transform", conf.level = 0.95,
   rep.measure = FALSE)

tmp.lab <- data.frame(lab = paste("CCC: ",
   round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ",
   round(tmp.ccc$rho.c[,2], digits = 2), " - ",
   round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))

z <- lm(method2 ~ method1)
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

## Concordance correlation plot:

library(ggplot2)
# ggplot(tmp, aes(x = method1, y = method2)) +
ggplot(tmp, aes(x = m1, y = m2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta),
      linetype = "dashed") +
  scale_x_continuous(limits = c(0,30), name = "PSR Fit") +
  scale_y_continuous(limits = c(0,30), name = "PSR Known") +
  geom_text(data = tmp.lab, x = 10, y = 25, label = tmp.lab$lab) +
  coord_fixed(ratio = 1 / 1)+
  theme_classic()




#----------------------------------------------------------------------------
## R1f

known <- readNIfTI('./sim_knowns_MC.nii.gz',reorient = FALSE)
fit <- readNIfTI('./R1f_julia_MC.nii.gz',reorient = FALSE)


#
# To do some data reduction, can use this rolling window mean
library(zoo)

m1<-rollapply(method1[method1>0],width = 5,by=5,FUN=mean,align = 'left')
m2<-rollapply(method2[method2>0],width = 5,by=5,FUN=mean,align = 'left')

library(epiR)
# tmp <- data.frame(method1, method2)
tmp <- data.frame(m1, m2)
tmp.ccc <- epi.ccc(method1, method2, ci = "z-transform", conf.level = 0.95,
                   rep.measure = FALSE)

tmp.lab <- data.frame(lab = paste("CCC: ",
                                  round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ",
                                  round(tmp.ccc$rho.c[,2], digits = 2), " - ",
                                  round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))

z <- lm(method2 ~ method1)
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

## Concordance correlation plot:

library(ggplot2)
ggplot(tmp, aes(x = m1, y = m2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta),
              linetype = "dashed") +
  scale_x_continuous(limits = c(0,2), name = "R1f Fit (1/s)") +
  scale_y_continuous(limits = c(0,2), name = "R1f Known (1/s)") +
  geom_text(data = tmp.lab, x = 0.5, y = 1.5, label = tmp.lab$lab) +
  coord_fixed(ratio = 1 / 1)+
  theme_classic()


##-------------------------------------

known <- readNIfTI('./sim_knowns_MC.nii.gz',reorient = FALSE)
fit <- readNIfTI('./PSR_julia_MC.nii.gz',reorient = FALSE)
# fit <- readNIfTI('./R1f_julia_MC.nii.gz',reorient = FALSE)

tmp1<-known[,,1]*100 # for PSR
# tmp1<-known[,,2] # for R1f
tmp2<-fit

ind<-tmp2>0
a<-tmp1[ind]
b<-tmp2[ind]
b[b==0]<-NaN
a[a==0]<-NaN

df <- data.frame(A=a,B=b)
df$diff <- df$A-df$B
df$avg <- rowMeans(df)
mean_diff <- mean(df$diff)
lower <- mean_diff - 1.96*sd(df$diff)
upper <- mean_diff + 1.96*sd(df$diff)

p1<-ggplot(df, aes(x = avg, y = diff)) +
  geom_point(size=2,alpha=0.1) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ggtitle("Bland-Altman Plot") +
  ylab("Difference Between Fit and Known") +
  xlab("Average R_2 (1/s)")+
  xlim(c(0,1.25))+
  theme_bw()
ggMarginal(p1, type="histogram", bins = 50)



library(epiR)

tmp <- data.frame(a, b)
tmp.ccc <- epi.ccc(method1, method2, ci = "z-transform", conf.level = 0.95,
                   rep.measure = FALSE)

tmp.lab <- data.frame(lab = paste("CCC: ",
                                  round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ",
                                  round(tmp.ccc$rho.c[,2], digits = 2), " - ",
                                  round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))

z <- lm(method2 ~ method1)
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

## Concordance correlation plot:

library(ggplot2)
ggplot(tmp, aes(x = a, y = b)) +
  geom_point(alpha=0.1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta),
              linetype = "dashed") +
  # scale_x_continuous(limits = c(0,2), name = "R1f Fit (1/s)") +
  # scale_y_continuous(limits = c(0,2), name = "R1f Known (1/s)") +
  # geom_text(data = tmp.lab, x = 0.5, y = 1.5, label = tmp.lab$lab) +
  scale_x_continuous(limits = c(0,25), name = "PSR Fit") +
  scale_y_continuous(limits = c(0,25), name = "PSR Known") +
  geom_text(data = tmp.lab, x = 10, y = 25, label = tmp.lab$lab) +
  coord_fixed(ratio = 1 / 1)+
  theme_classic()

# dfa<-data.frame(x=a)
# dfb<-data.frame(x=b)
# 
# dfa$e <- 'Known'
# dfb$e <- 'Fit'
# 
# df<-rbind(dfa,dfb)
# 
# ggplot(df,aes(x,fill=e))+
#   # geom_histogram(aes(y=..density..),bins=100)+
#   geom_density(alpha = 0.5)+
#   # xlim(c(0,max(x)))+
#   theme_classic()+
#   xlab("Relaxation Rate (1/s)")
# 
# df3=data.frame(x=colMeans(rbind(a,b)) )
# ggplot(df3,aes(x))+
#   geom_histogram(aes(y=..density..),bins=100)+
#   geom_density(alpha = 0.5)+
#   theme_classic()+
#   xlab("Relaxation Rate (1/s)")
# 
# df2=data.frame(x=a-b)
# ggplot(df2,aes(x))+
#   geom_histogram(aes(y=..density..),bins=100)+
#   geom_density()+
#   theme_classic()+
#   xlab("Relaxation Rate (1/s)")
