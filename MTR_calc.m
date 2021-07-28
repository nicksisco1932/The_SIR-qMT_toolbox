
close all;clc;clear all

path = "/Users/nicks/Documents/MRI_data/Testing/20210720_DTI_SAGE_MTR_CTRL001/DICOM/";

S0 = sprintf("%s/MTR_0000.nii.gz",path);

Smt = sprintf("%s/RegWarped.nii.gz",path);
MASK = sprintf("%s/MTR_arm_mask.nii.gz",path);


addpath(genpath('C:\Users\nicks\Documents\Github\The_MRI_toolbox\Matlab\SAGE_DSC\mfiles'))


S0_img = loader(S0);
Smt_img = loader(Smt);
MASK_img = loader(MASK);

MTR = MASK_img.*( (S0_img - Smt_img)./S0_img);
MTR(isnan(MTR))=0;
MTR(isinf(MTR))=0;
MTR(MTR<0)=0;
%%

img1 = permute(MTR(:,:,1),[2 1 ]);
img2 = permute(MTR(:,:,36),[2 1 ]);
img3 = permute(MTR(:,:,36*2),[2 1 ]);
img4 = permute(MTR(:,:,36*3),[2 1 ]);
multi = cat(3,img1,img2,img3,img4);
montage(multi);
%%
info = niftiinfo(Smt);

niftiwrite( MTR,sprintf("%s/MTR.nii.gz",path),info)
