# The MRI Toolbox

## Table of contents

## Introduction

This is toolbox for qunatitative MRI processing. It is intended as a one stop shop for DSC fitting (non-linear fitting and log ratio calculations), SIR-qMT, and other useful tools for analyzing MRI data. 

As of 20210701, it is just getting off the ground. Improvements will include documented pages for each method used here. There will be MATLAB, Julia, Python, and R code that will be used for fitting and analysis.

## Work in progress

Temp math
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\space\frac{M_f(t)}{M_{f\infty}}=b_f^+exp(-R_1^+t)+b_f^-exp(-R_1^-t)+1" title="\Equation 1" />
</p>

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\space\2R_1^\pm=R_{1f}+R_{1m}+k_{fm}+k_{mf}\pm\sqrt{(R_{1f}-R_{1m}+k_{fm}-k_{mf})^2+4k_{fm}k_{mf}}" title="\Equation 2" />
</p>

and

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\space\b_f^\pm=\pm\frac{[\frac{s_fM_f(0^-)}{M_{f\infty}}-1](R_{1f}-R_1^\pm)+[\frac{S_fM_f(0^-)}{M_{f\infty}}-\frac{S_mM_m(0^-)}{M_{m\infty}}]k_{fm}}{R_1^+-R_1^-}" title="\Equation 3" />
</p>