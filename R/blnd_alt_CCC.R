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



m_known<-melt(known[,,1]*100)
m_fit<-melt(fit)
b<-bland.altman.plot(m_known[,-1][,-1],m_fit[,-1][,-1],graph.sys = 'ggplot2',conf.int=0.95)
print(b+xlab('Mean (Matlab vs Julia Default)') +
        ggtitle("BA: PSR ") +ylab('Matlab-Julia')+theme_classic())




method1=m_fit[,-1][,-1]
method2=m_known[,-1][,-1]

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
## Not run:
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


diff<-(m1-m2)
diffp<-(diff)/mean(c(m1,m2))*100
sd.diff <- sd(diff)
sd.diffp <- sd(diffp)
my.data<-data.frame(m1,m2,diff,diffp)
my.formula<-y~x
library(ggExtra)
library(ggpmisc)
library(devtools)
diffplot <- ggplot(my.data, aes(m1, diff)) + 
  geom_point(size=2, colour = rgb(0,0,0, alpha = 0.5)) + 
  # Line Fit
  # geom_smooth(method = "lm", se=FALSE, color="red", formula = my.formula) +
  # stat_poly_eq(formula = my.formula,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE)+
  # Line Fit END  
  theme_bw() + 
  #when the +/- 2SD lines will fall outside the default plot limits 
  #they need to be pre-stated explicitly to make the histogram line up properly. 
  #Thanks to commenter for noticing this.
  ylim(mean(my.data$diff) - 3*sd.diff, mean(my.data$diff) + 3*sd.diff) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_hline(yintercept = mean(my.data$diff)) +
  geom_hline(yintercept = mean(my.data$diff) + 2*sd.diff, linetype = 2) +
  geom_hline(yintercept = mean(my.data$diff) - 2*sd.diff, linetype = 2) +
  ylab("Difference PSR from Known and Fit") +
  xlab("Average PSR Values")+
  xlim(c(0,25))

ggMarginal(diffplot, type="histogram", bins = 50)


#----------------------------------------------------------------------------
## R1f
# a<-bland.altman.plot(PSR_mdata[,-1],PSRd_mdata[,-1],graph.sys = 'ggplot2',conf.int=0.95)
# print(a+xlab('Mean (PSR and PSR Default)') +
#         ggtitle("BA: PSR Forward Diff vs Default") +ylab('PSR-PSR(default)')+theme_classic())

known <- readNIfTI('./sim_knowns_MC.nii.gz',reorient = FALSE)
fit <- readNIfTI('./R1f_julia_MC.nii.gz',reorient = FALSE)



m_known<-melt(known[,,2])
m_fit<-melt(fit)
b<-bland.altman.plot(m_known[,-1][,-1],m_fit[,-1][,-1],graph.sys = 'ggplot2',conf.int=0.95)
print(b+xlab('Mean (Known vs Fit (R1f)') +
        ggtitle("BA: R1f ") +ylab('Known-Fit')+theme_classic())




method1=m_fit[,-1][,-1]
method2=m_known[,-1][,-1]

#
# To do some data reduction, can use this rolling window mean
library(zoo)
# TS<-zoo(method1)
# TS_y<-zoo(method2)
# method1<-rollapply(TS[TS>20],width = 300,by=200,FUN=mean,align = 'left')
# method2<-rollapply(TS_y[TS_y>20],width = 300,by=200,FUN=mean,align = 'left')

# m1<-rollapply(TS[TS>0],width = 5,by=5,FUN=mean,align = 'left')
# m2<-rollapply(TS_y[TS_y>0],width = 5,by=5,FUN=mean,align = 'left')

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
## Not run:
library(ggplot2)
# ggplot(tmp, aes(x = method1, y = method2)) +
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


diff<-(m1-m2)
diffp<-(diff)/mean(c(m1,m2))*100
sd.diff <- sd(diff)
sd.diffp <- sd(diffp)
my.data<-data.frame(m1,m2,diff,diffp)
my.formula<-y~x
library(ggExtra)
library(ggpmisc)
library(devtools)
diffplot <- ggplot(my.data, aes(m1, diff)) + 
  geom_point(size=2, colour = rgb(0,0,0, alpha = 0.5)) + 
  # Line Fit
  # geom_smooth(method = "lm", se=FALSE, color="red", formula = my.formula) +
  # stat_poly_eq(formula = my.formula,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE)+
  # Line Fit END  
  theme_bw() + 
  #when the +/- 2SD lines will fall outside the default plot limits 
  #they need to be pre-stated explicitly to make the histogram line up properly. 
  #Thanks to commenter for noticing this.
  ylim(mean(my.data$diff) - 3*sd.diff, mean(my.data$diff) + 3*sd.diff) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_hline(yintercept = mean(my.data$diff)) +
  geom_hline(yintercept = mean(my.data$diff) + 2*sd.diff, linetype = 2) +
  geom_hline(yintercept = mean(my.data$diff) - 2*sd.diff, linetype = 2) +
  ylab("Difference R1f from Known and Fit") +
  xlab("Average R1f Values (1/s)")+
  xlim(c(0,2))

ggMarginal(diffplot, type="histogram", bins = 50)

