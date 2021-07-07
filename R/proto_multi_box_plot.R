library(oro.nifti)
library(ggplot2)
library(ggpubr)
library(ppcor)

seg_fname<-'/Users/nicks/Documents/Github/SAGE_Structural_proc/tempDir/PT16_aseg_2_MNI.nii.gz'
st1t2_fname<-'/Users/nicks/Documents/Github/SAGE_Structural_proc/tempDir/sT1T2.nii.gz'
FA_fname<-'/Users/nicks/Documents/Github/SAGE_Structural_proc/PT16_FA_MNI.nii.gz'
CBF_fname<-'/Users/nicks/Documents/Github/SAGE_Structural_proc/PT16_CBF_MNI.nii.gz'

seg_mask<-readNIfTI(seg_fname)
st1t2_img<-readNIfTI(st1t2_fname)
fa_img<-readNIfTI(FA_fname)
cbf_img<-readNIfTI(CBF_fname)

cutoff<-500
pallidum_left<-st1t2_img[seg_mask==13]
pallidum_right<-st1t2_img[seg_mask==52]
pallidum_left[pallidum_left>cutoff]=0
pallidum_right[pallidum_right>cutoff]=0

data<-data.frame(pallidum_left=pallidum_left,group="1")

library(dplyr)
data<-dplyr::bind_rows(data.frame(pallidum_left=pallidum_left,group="Left"),
  data.frame(pallidum_right=pallidum_right,group="Right"))

p1<-ggplot(data,aes(x=group,y=pallidum_left,fill=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=pallidum_right),outlier.shape = NA)+
  ylab("GMT1T2 Ratio")+theme_bw()+
  ylim(c(0,cutoff))

left<-st1t2_img[seg_mask==2]
right<-st1t2_img[seg_mask==41]
left[left>cutoff]=0
right[right>cutoff]=0

library(dplyr)
data<-dplyr::bind_rows(data.frame(left=left,group="Left"),
                       data.frame(right=right,group="Right"))

p2<-ggplot(data,aes(x=group,y=left,fill=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=right),outlier.shape = NA)+
  ylab("GMT1T2 Ratio")+theme_bw()+
  ylim(c(0,cutoff))

hyop_WM<-st1t2_img[seg_mask==77]
hyop_WM[hyop_WM>cutoff]=0
# right[right>cutoff]=0

library(dplyr)
# data<-dplyr::bind_rows(data.frame(left=left,group="Left"),
#                        data.frame(right=right,group="Right"))
data<-data.frame(hyop_WM=hyop_WM,group="WM Hypo")

p3<-ggplot(data,aes(x=group,y=hyop_WM,fill=group))+
  geom_boxplot(outlier.shape = NA)+
  # geom_boxplot(aes(x=group,y=right),outlier.shape = NA)+
  ylab("GMT1T2 Ratio")+theme_bw()+
  ylim(c(0,cutoff))

ggarrange(p1,p2,p3)

cutoff<-500
pallidum_left<-st1t2_img[seg_mask==13]
pallidum_right<-st1t2_img[seg_mask==52]
pallidum_left[pallidum_left>cutoff]=0
pallidum_right[pallidum_right>cutoff]=0

putamen_left<-st1t2_img[seg_mask==12]
putamen_right<-st1t2_img[seg_mask==51]
putamen_left[putamen_left>cutoff]=0
putamen_right[putamen_right>cutoff]=0

WM_left<-st1t2_img[seg_mask==2]
WM_right<-st1t2_img[seg_mask==41]
WM_left[WM_left>cutoff]=0
WM_right[WM_right>cutoff]=0
hypo_WM<-st1t2_img[seg_mask==77]
hypo_WM[hypo_WM>cutoff]=0
post_CC<-st1t2_img[seg_mask==251]
post_CC[post_CC>cutoff]=0

data<-dplyr::bind_rows(
  data.frame(pallidum_left=pallidum_left,group="Pallidum Left"),
  data.frame(pallidum_right=pallidum_right,group="Pallidum Right"),
  data.frame(putamen_left=putamen_left,group="Putamen Left"),
  data.frame(putamen_right=putamen_right,group="Putamen Right"),
  data.frame(WM_left=WM_left,group="WM Left"),
  data.frame(WM_right=WM_right,group="WM Right"),
  data.frame(hypo_WM=hypo_WM,group="WM Hypo"),
  data.frame(post_CC=post_CC,group="Posterior CC")
)

ggplot(data,aes(x=group,y=pallidum_left,fill=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=pallidum_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=putamen_left),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=putamen_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=WM_left),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=WM_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=hypo_WM),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=post_CC),outlier.shape = NA)+
  ylab("GMT1T2 Ratio")+theme_bw()+
  ylim(c(0,cutoff))


cutoff<-1
pallidum_left<-fa_img[seg_mask==13]
pallidum_right<-fa_img[seg_mask==52]
pallidum_left[pallidum_left>cutoff]=0
pallidum_right[pallidum_right>cutoff]=0

putamen_left<-fa_img[seg_mask==12]
putamen_right<-fa_img[seg_mask==51]
putamen_left[putamen_left>cutoff]=0
putamen_right[putamen_right>cutoff]=0

WM_left<-fa_img[seg_mask==2]
WM_right<-fa_img[seg_mask==41]
WM_left[WM_left>cutoff]=0
WM_right[WM_right>cutoff]=0
hypo_WM<-fa_img[seg_mask==77]
hypo_WM[hypo_WM>cutoff]=0
post_CC<-fa_img[seg_mask==251]
post_CC[post_CC>cutoff]=0

data<-dplyr::bind_rows(
  data.frame(pallidum_left=pallidum_left,group="Pallidum Left"),
  data.frame(pallidum_right=pallidum_right,group="Pallidum Right"),
  data.frame(putamen_left=putamen_left,group="Putamen Left"),
  data.frame(putamen_right=putamen_right,group="Putamen Right"),
  data.frame(WM_left=WM_left,group="WM Left"),
  data.frame(WM_right=WM_right,group="WM Right"),
  data.frame(hypo_WM=hypo_WM,group="WM Hypo"),
  data.frame(post_CC=post_CC,group="Posterior CC")
)

ggplot(data,aes(x=group,y=pallidum_left,fill=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=pallidum_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=putamen_left),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=putamen_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=WM_left),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=WM_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=hypo_WM),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=post_CC),outlier.shape = NA)+
  ylab("FA")+theme_bw()+
  ylim(c(0,cutoff))


cutoff<-250
pallidum_left<-cbf_img[seg_mask==13]
pallidum_right<-cbf_img[seg_mask==52]
pallidum_left[pallidum_left>cutoff]=0
pallidum_right[pallidum_right>cutoff]=0

putamen_left<-cbf_img[seg_mask==12]
putamen_right<-cbf_img[seg_mask==51]
putamen_left[putamen_left>cutoff]=0
putamen_right[putamen_right>cutoff]=0

WM_left<-cbf_img[seg_mask==2]
WM_right<-cbf_img[seg_mask==41]
WM_left[WM_left>cutoff]=0
WM_right[WM_right>cutoff]=0
hypo_WM<-cbf_img[seg_mask==77]
hypo_WM[hypo_WM>cutoff]=0
post_CC<-cbf_img[seg_mask==251]
post_CC[post_CC>cutoff]=0

data<-dplyr::bind_rows(
  data.frame(pallidum_left=pallidum_left,group="Pallidum Left"),
  data.frame(pallidum_right=pallidum_right,group="Pallidum Right"),
  data.frame(putamen_left=putamen_left,group="Putamen Left"),
  data.frame(putamen_right=putamen_right,group="Putamen Right"),
  data.frame(WM_left=WM_left,group="WM Left"),
  data.frame(WM_right=WM_right,group="WM Right"),
  data.frame(hypo_WM=hypo_WM,group="WM Hypo"),
  data.frame(post_CC=post_CC,group="Posterior CC")
)

ggplot(data,aes(x=group,y=pallidum_left,fill=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=pallidum_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=putamen_left),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=putamen_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=WM_left),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=WM_right),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=hypo_WM),outlier.shape = NA)+
  geom_boxplot(aes(x=group,y=post_CC),outlier.shape = NA)+
  ylab("CBF")+theme_bw()+
  ylim(c(0,cutoff))

