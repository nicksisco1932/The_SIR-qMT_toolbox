# Title     : A Script to Analyze SIR MC
# Objective : Analysis
# Created by: Nick Sisco
# Created on: 8/17/2021

library(blandr)
library(ggplot2)
library(oro.nifti)
library(reshape2)
library(BlandAltmanLeh)
library(ggExtra)
library(epiR)
library(ggpubr)

LCCC <- function(a,b,l,u,x_lab,y_lab,legend_x,legend_y) {
  tmp <- data.frame(a, b)
  tmp.ccc <- epi.ccc(a, b, ci = "z-transform", conf.level = 0.95,
                     rep.measure = FALSE)
  
  tmp.lab <- data.frame(lab = paste("CCC: ",
                                    round(tmp.ccc$rho.c[,1], digits = 3), " (95% CI ",
                                    round(tmp.ccc$rho.c[,2], digits = 3), " - ",
                                    round(tmp.ccc$rho.c[,3], digits = 3), ")", sep = ""))
  z <- lm(b ~ a)
  alpha <- summary(z)$coefficients[1,1]
  beta <-  summary(z)$coefficients[2,1]
  tmp.lm <- data.frame(alpha, beta)
  p<-ggplot(tmp,aes(x=a,y=b)) +
    geom_point(size=2,alpha=0.1) +
    geom_abline(intercept = 0, slope = 1) +
    geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta),
                linetype = "dashed") +
    scale_x_continuous(limits = c(l,u), name = x_lab) +
    scale_y_continuous(limits = c(l,u), name = y_lab) +
    geom_text(data = tmp.lab, x = legend_x, y = legend_y, label = tmp.lab$lab) +
    coord_fixed(ratio = 1 / 1)+
    theme_classic()
  return(list(p,tmp,tmp.lm,tmp.ccc))
}

#----------------------------------------------------------------------------
# Edit, used percent difference and calculated in Julia. Loading in saved Niftis

# PSR_perc_diff <- readNIfTI('./Res_PSR_julia_MC.nii.gz')
# R1f_perc_diff <- readNIfTI('./Res_R1f_julia_MC.nii.gz')
PSR_perc_diff <- readNIfTI('/Users/nicks/Documents/Github/SIR_qMT_python/Res_PSR_julia_MC.nii.gz')
R1f_perc_diff <- readNIfTI('/Users/nicks/Documents/Github/SIR_qMT_python/Res_R1f_julia_MC.nii.gz')

ind<-PSR_perc_diff!=Inf
a<-PSR_perc_diff[TRUE]
b<-R1f_perc_diff[TRUE]

df<-data.frame(x=a,y=b)

p1<-ggplot(df,aes(x=x))+
    geom_histogram(aes(y=..ndensity..),bins=100)+
    # geom_density(aes(y=..ndensity..))+
    theme_classic()+
    xlab("PSR Percent Difference (%)")
p2<-ggplot(df,aes(x=y))+
  geom_histogram(aes(y=..ndensity..),bins=100)+
  # geom_density(aes(y=..ndensity..))+
  theme_classic()+
  xlab("R1f Percent Difference (%)")

ggarrange(p1,p2,
          labels = c("A","B"),
          align = 'hv')

# ## R1f
# 
known <- readNIfTI('/Users/nicks/Documents/Github/SIR_qMT_python/sim_knowns_MC.nii.gz',reorient = FALSE)
fit <- readNIfTI('/Users/nicks/Documents/Github/SIR_qMT_python/R1f_julia_MC.nii.gz',reorient = FALSE)

tmp2<-known[,,2] # for R1f
tmp1<-fit

ind<-tmp2!=Inf
a<-tmp1[TRUE]
b<-tmp2[TRUE]
# b[b==0]<-NaN # only for brain masks
# a[a==0]<-NaN

df<-data.frame(x=a-b)
# p1<-ggplot(df,aes(x))+
#   geom_histogram(aes(y=..ndensity..),bins=100)+
#   # geom_density(aes(y=..ndensity..))+
#   theme_classic()+
#   xlab("R1f Difference, Fit-Known (%)")

A<-LCCC(b,a,0,2,"R1f Known (s-1)","R1f Fit (s-1)",0.5,2)
A[[1]]

##-------------------------------------
# PSR

fit <- readNIfTI('/Users/nicks/Documents/Github/SIR_qMT_python/PSR_julia_MC.nii.gz')

tmp2<-known[,,1]*100 # for PSR
tmp1<-fit

a<-tmp1[TRUE]
b<-tmp2[TRUE]
# b[b==0]<-NaN
# a[a==0]<-NaN


df<-data.frame(x=a-b)
# p2<-ggplot(df,aes(x))+
#   geom_histogram(aes(y=..ndensity..),bins=100)+
#   # geom_density(aes(y=..ndensity..))+
#   theme_classic()+
#   xlab("PSR Difference, Fit-Known (%)")

# B<-LCCC(a,b,0,25,"PSR Fit (%)","PSR Known (%)",8,23)
B<-LCCC(b,a,0,25,"PSR Known (%)","PSR Fit (%)",8,23)
B[[1]]
# 
ggarrange(p2,p1,A[[1]],B[[1]],
          labels = c("A","B","C","D"),align = 'hv')
