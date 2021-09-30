library(blandr)
library(ggplot2)
library(oro.nifti)
library(reshape2)
library(BlandAltmanLeh)
library(ggExtra)
library(epiR)
library(ggpubr)


BSA_percent<-c(1.25,2.5,5,10,20)
BSA_frac<-BSA_percent/100

R1f<-c(0.41,0.46,0.51,0.63,0.83) # measured from phantoms
R1f_stdev<-c(0.01,0.03,0.02,0.02,0.05)

PSR<-c(0.9,1.9,3.9,6.5,13.2) # measured from phantoms
PSR_frac<-PSR/100

PSR_stdev<-c(0.7,1.5,1.1,1.1,2.9)

PSR_stdev_frac<-PSR_stdev/100

df<-data.frame(x=BSA_frac,y=PSR_frac,err=PSR_stdev_frac,y2=R1f,err2=R1f_stdev)

z <- lm( PSR_frac ~ BSA_frac)
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lab <- data.frame(lab = paste("y =  ",
                                  round(beta, digits = 3), "x + ",
                                  round(alpha, digits = 3)))


p1<-ggplot(df,aes(x=x,y=y))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=-err+y,ymax=err+y), width=0.004,
                position=position_dodge(.9)
                )+
  # geom_abline(data = df, aes(intercept = alpha, slope = beta),
  #             linetype = "dashed")+
  geom_smooth(method = 'lm',se=FALSE,color='red')+
  theme_bw()+
  xlab('BSA Fraction (w/v)')+
  ylab('f (PSR/(1+PSR)')+
  geom_text(data = tmp.lab, x = 0.05, y = 0.1, label = tmp.lab$lab) 

z <- lm( R1f ~ BSA_frac)
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lab <- data.frame(lab = paste("y =  ",
                                  round(beta, digits = 3), "x + ",
                                  round(alpha, digits = 3)))


p2<-ggplot(df,aes(x=x,y=y2))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=-err2+y2,ymax=err2+y2), width=0.004,
                position=position_dodge(.9)
  )+
  # geom_abline(data = df, aes(intercept = alpha, slope = beta),
  #             linetype = "dashed")+
  geom_smooth(method = 'lm',se=FALSE,color='red')+
  theme_bw()+
  xlab('BSA Fraction (w/v)')+
  ylab('R1f (s-1)')+
  geom_text(data = tmp.lab, x = 0.05, y = 0.7, label = tmp.lab$lab) 
p2

ggarrange(p1,p2)

