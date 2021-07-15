# The MRI toolbox

<p align="center">
  <img src="https://github.com/nicksisco1932/The_MRI_toolbox/blob/master/Images/MR_logo_big.png" alt="drawing" width="400"/>
</p>

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nicksisco1932.github.io/The_MRI_toolbox/docs)

Welcome to the magnetic resonance toolbox. This is intended to become a repository for magnetic resonance imaging processing pipelines that do not already have dedicated and extensive software written. The main pipelines shown here are for processing quantitative magnetization transfer imaging using selective inversion recovery (SIR-qMT), spin- and gradient-echo non-linear fitting (SAGE), and multi-flip angle (MFA) T1 maps. All the fitting is implemented using Julia, which is computationally fast but still readable.

Updates to this will include a tutorial on how to use the code and more detailed documentation. 


**Command Line**

Make sure you change the <PATH> to the absolute path of your SIR data and brain mask.
```Bash
julia ./SIR_fit.jl --TI 15 15 278 1007 --TD 684 4121 2730 10 --SIR_Data <PATH>/SIR_DATA.nii.gz --SIR_brainMask <PATH>/brain_mask.nii.gz --kmf 14.5 --Sm 0.83
```


**Results of a Successful SIR-qMT fit**
<p align="center">
  <img src="https://github.com/nicksisco1932/The_MRI_toolbox/blob/master/Images/Brain_Figure.png" alt="drawing" width="400"/>
</p>
Representative SIR-qMT on a healthy volunteer. A represents the first data point corresponding to t<sub>I</sub>,t<sub>D</sub> = 278,2730 ms. B, C, and D are maps from the fit parameters pool size ratio, R<sub>1f</sub>, and S<sub>f</sub> (B<sub>1</sub> inhomogeneity), respectively. These images are consistent with published parameters, white matter have the highest relative PSR and R<sub>1f</sub>, while S<sub>f</sub> remains relatively flat at 3T with slight increases near the posterior of this map.
  
  
# Using Matlab or Python to call Julia examples.
  

# SAGE Fitting call
```bash
julia ./SAGE_biexp_fit.jl --TE_nii_names <PATH>/PT1319001_TE1_img_w_Skull.nii.gz <PATH>/PT1319001_TE2_img_w_Skull.nii.gz <PATH>/PT1319001_TE3_img_w_Skull.nii.gz <PATH>/PT1319001_TE4_img_w_Skull.nii.gz <PATH>/PT1319001_TE5_img_w_Skull.nii.gz --SAGE_nii_brainMask <PATH>/bPT1319001_preb_mask.nii.gz --echos 7.82 28.8 60.7 81.6 102.6  
```
  
**Python**
```Python
# A function wrapper to call Julia to fit spin- and gradient-echo signal to a piecewise function using Julia
import nibabel as nib
import os

def loader(path):
  headerinfo = nib.load(path)
  vol = headerinfo.get_fdata()
  return vol
end

def sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname):
  julia_cmd = "<PATH TO JULIA>/julia";
  SAGE_Fit_julia = "<PATH TO MRI TOOLBOX>/The_MRI_toolbox/SAGE_biexp_fit.jl"
  cmd = ‘{} {} {} {} {} {} {} {}’.format(julia_cmd,SAGE_Fit_julia,te1,te2,te3,te4,te5,b_fname)  
  os.system(cmd)
  r2s = os.path.join(‘{}/SAGE_R2s_Brain_julia.nii.gz’.format(base_path))
  r2 = os.path.join (‘{}/SAGE_R2_Brain_julia.nii.gz’.format(base_path))  
  r2simg = loader(r2s)
  r2img = loader(r2)
  return r2simg,r2img

r2simg,r2img = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname) # this should return two arrays

```

**MATLAB**

```MATLAB
% A function wrapper to call Julia to fit spin- and gradient-echo signal to a piecewise function using Julia
[r2simg,r2img] = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname);

function [vol, headerinfo] = loader(path)
  % requires MATLAB: image processing toolbox
  headerinfo = niftiinfo(path);
  vol = niftiread(info);
end

function [r2simg,r2img] = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname)
  julia_cmd = "<PATH TO JULIA>/julia";
  SAGE_Fit_julia = "<PATH TO MRI TOOLBOX>/The_MRI_toolbox/SAGE_biexp_fit.jl";
  unix(sprintf("%s %s %s %s %s %s %s %s",julia_cmd,SAGE_Fit_julia,te1,te2,te3,te4,te5,b_fname))  
  r2s = sprintf("%s/SAGE_R2s_Brain_julia.nii.gz",base_path);
  r2 = sprintf("%s/SAGE_R2_Brain_julia.nii.gz",base_path);  
  r2simg = loader(r2s);  
  r2img = loader(r2);
end

```


