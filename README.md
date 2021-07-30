# The SIR-qMT toolbox

<p align="center">
  <img src="https://github.com/nicksisco1932/The_MRI_toolbox/blob/master/Images/MR_logo_big.png" alt="drawing" width="400"/>
</p>

Welcome to the magnetic resonance toolbox. The main pipeline shown here are for processing quantitative magnetization transfer imaging using selective inversion recovery (SIR-qMT). All the fitting is implemented using Julia, which is computationally fast but still readable.

The first thing you will have to do is download and install Julia! 
https://julialang.org/downloads/

# Jupyter Notebook with Julia kernel
If you have not cloned this repository, you can download the zip in the top right of this page or if you are familiar with git, use this. 

```bash
git clone https://github.com/nicksisco1932/The_SIR-qMT_toolbox
```

After that:
  Do these steps
   1) Download and install Julia 1.6
   2) Then open the Julia app.
   3) Using the Julia REPL command line
```Julia
      using Pkg 
      Pkg.add("IJulia") 
      Pkg.precompile() 
      notebook()
```
Then open and run the SIR_qMT_test.ipynb from your local device. 

**Compare your results to our results**
<p align="center">
  <img src="https://github.com/nicksisco1932/The_SIR-qMT_toolbox/blob/master/Images/Figure%201.png" alt="drawing" width="600"/>
</p>

# Command Line

Make sure you change the <PATH> to the absolute path of your SIR data and brain mask.
  
_The dimensions of your SIR_Data file must be __[nx,ny,nz,nd]__ where nd is the dynamics and nx,ny,nz are 3D matrix dims_

  Set your file path first
```Bash
path=<PATH TO FILES>
```
  Then copy and paste this.
```Bash
brain_mask=$path/brain_mask.nii.gz
SIR_4D_DATA=$path/SIR_DATA.nii.gz
julia ./SIR_fit.jl --TI 15 15 278 1007 --TD 684 4121 2730 10 --SIR_Data $SIR_4D_DATA --SIR_brainMask $brain_mask --kmf 14.5 --Sm 0.83
```


**Results of a Successful SIR-qMT fit**
<p align="center">
  <img src="https://github.com/nicksisco1932/The_MRI_toolbox/blob/master/Images/Brain_Figure.png" alt="drawing" width="400"/>
</p>
Representative SIR-qMT on a healthy volunteer. A represents the first data point corresponding to t<sub>I</sub>,t<sub>D</sub> = 278,2730 ms. B, C, and D are maps from the fit parameters pool size ratio, R<sub>1f</sub>, and S<sub>f</sub> (B<sub>1</sub> inhomogeneity), respectively. These images are consistent with published parameters, white matter have the highest relative PSR and R<sub>1f</sub>, while S<sub>f</sub> remains relatively flat at 3T with slight increases near the posterior of this map.
  
  
# Using Matlab or Python to call Julia examples. 

**Shell Script**
```tcsh
#!/bin/tcsh
path=<FULL PATH TO FILES>           # Full path of image directory
brain_mask=$path/brain_mask.nii.gz  # Mask name
SIR_4D_DATA=$path/SIR_DATA.nii.gz   # 4D dataset
julia ./SIR_fit.jl --TI 15 15 278 1007 --TD 684 4121 2730 10 --SIR_Data $SIR_4D_DATA --SIR_brainMask $brain_mask --kmf 12.5 --Sm 0.83 

```
  
**Python**
```Python
#!/usr/local/bin/Python3.8              # or this can be your virtual environment
import os                               # A core library, no virtual environment needed
def pj(in1, in2):
    return os.path.join(in1, in2)

def julia_call(ti,td,kmf,sm,data,mask): # A function to call julia
    cmd = 'julia ./SIR_fit.jl --TI {} --TD {} --kmf {} --Sm {} --SIR_Data {} --SIR_brainMask {}'.format(ti, td, kmf, sm, data, mask)
    print(cmd)
    os.system(cmd)
  
path=<FULL PATH TO FILES>               # Full path of image directory
brain_mask=pj(path, brain_mask.nii.gz)  # Mask name
SIR_4D_DATA=pj(path, SIR_DATA.nii.gz)   # 4D dataset

ti_values = [15 15 278 1007]
ti_values = [684 4121 2730 10]
kmf = 12.5
sm = 0.83 
julia_call(ti_values, td_values, kmf, sm, SIR_4D_DATA, brain_mask)
```

**MATLAB**

```MATLAB
path=<FULL PATH TO FILES>                     % Full path of image directory
brain_mask=fullfile(path, brain_mask.nii.gz)  % Mask name
SIR_4D_DATA= fullfile(path, SIR_DATA.nii.gz)  % 4D dataset

ti_values = [15,15,278,1007];
ti_values = [684,4121,2730,10];
kmf = 12.5;
sm = 0.83;
julia_call(ti_values, td_values, kmf, sm, SIR_4D_DATA, brain_mask)
  
function julia_call(ti,td,kmf,sm,data,mask)
  cmd = sprintf(“julia ./SIR_fit.jl --TI %s --TD %s --kmf %s --Sm %s --SIR_Data %s --SIR_brainMask %s”, ti, td, kmf, sm, data, mask)
  disp(cmd)
  system(cmd)                                   % alternatively, use unix(cmd)
end
  
```
