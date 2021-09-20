# The MRI Toolbox

## Table of contents
- [Introduction](https://github.com/nicksisco1932/The_SIR-qMT_toolbox/blob/master/docs/index.md#introduction)
- [Selective Inversion Recovery](https://github.com/nicksisco1932/The_SIR-qMT_toolbox/blob/master/docs/index.md#selective-inversion-recovery)
- [Matrix formulation](https://github.com/nicksisco1932/The_SIR-qMT_toolbox/blob/master/docs/index.md#matrix-formulation)
- [General formulation](https://github.com/nicksisco1932/The_SIR-qMT_toolbox/blob/master/docs/index.md#general-formulation)

## Introduction

This is toolbox for qunatitative MRI processing. It is intended as a one stop shop for DSC fitting (non-linear fitting and log ratio calculations), SIR-qMT, and other useful tools for analyzing MRI data. 

As of 20210701, it is just getting off the ground. Improvements will include documented pages for each method used here. There will be MATLAB, Julia, Python, and R code that will be used for fitting and analysis.

### Selective Inversion Recovery

The selective inversion recovery experiment relies on the following scheme for inversion of the spin polarization:

<img src="https://latex.codecogs.com/svg.image?pi&space;\rightarrow&space;(\pi/2&space;\rightarrow&space;\pi)&space;\rightarrow&space;ECHO" title="pi \rightarrow (\pi/2 \rightarrow \pi) \rightarrow ECHO" />

where the first $\pi$ proceeds the Hahn echo, or spin echo, sequence by a time delay called TI (inversion time). The effect of this is to polarize the spin vectors into a high energy state, which is against the static magnetic field. The polarized spin vectors are then allowed to recover for a specific amount of time, which is the inversion time. At the end of this time, the spins are irradiated with the Hahn echo sequence to polarize the spins into the transverse plane and then to refocus the dephasing spins. The relaxation from this final refocusing pulse is detected and thus the spins magnetization are separated according to their inversion recovery from the initial $\pi$ pulse. The overall effect is that the spins that have high degrees of freedom to relax tend to have short relaxation times. 


### Matrix formulation
The time evolution of the magnetization vector in _x_, _y_, and _z_ planes can be expressed as a matrix: 

<img src="https://latex.codecogs.com/svg.image?A&space;=&space;\begin{bmatrix}&space;&space;&space;-R_{2f}&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;&space;\Delta\omega&space;&space;&space;&space;&space;&space;&space;&space;&&space;-\omega_1\sin\phi&space;&&space;0\\&space;&space;&space;-\Delta\omega&space;&space;&space;&space;&space;&space;&space;&&space;&space;-R_{2f}&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;\omega_1\cos\phi&space;&&space;&space;0&space;\\&space;&space;&space;\omega_1\sin\phi&space;&space;&&space;&space;-\omega_1\cos\phi&space;&&space;-(R_{1f}&space;&plus;&space;k_{mf})&space;&space;&&space;k_{mf}&space;\\&space;&space;&space;0&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;0&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;&space;k_{mf}&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;-(R_{1m}&plus;R_{RF}&plus;k_{mf})&space;\\&space;\end{bmatrix}" title="A = \begin{bmatrix} -R_{2f} & \Delta\omega & -\omega_1\sin\phi & 0\\ -\Delta\omega & -R_{2f} & \omega_1\cos\phi & 0 \\ \omega_1\sin\phi & -\omega_1\cos\phi & -(R_{1f} + k_{mf}) & k_{mf} \\ 0 & 0 & k_{mf} & -(R_{1m}+R_{RF}+k_{mf}) \\ \end{bmatrix}" />

<img src="https://latex.codecogs.com/svg.image?B&space;=&space;\begin{bmatrix}&space;0&space;\\&space;0&space;\\&space;R_{1f}M_{0f}&space;\\&space;R_{1m}M_{0m}&space;\end{bmatrix}" title="B = \begin{bmatrix} 0 \\ 0 \\ R_{1f}M_{0f} \\ R_{1m}M_{0m} \end{bmatrix}" />

$\Delta\omega$ is the frequency offset from resonance from the radio frequency pulse, $\omega_1$ is the frequency of precession about the RF pulse, and $\phi$ is the phase of the RF pulse in the transverse plane. The standard Block equations are not valid for the macromlecular pool due to the fact that the signal for this pool does not reflect a  Lorentzian lineshape, rather a Super-Lorentzian(cite). Given this, equation 2 has been augmented to include by a single longitudinal component whose saturation is governed by the rate $R_{RF} = \pi\omega^2_1g_m(\Delta\omega)$, where $g_m$ is the lineshape function of the macromolecular pool. 

The general solution to the system of equations, is given by:

<img src="https://latex.codecogs.com/svg.image?M(t)&space;=&space;exp(At)M(0)&space;&plus;[exp(At)-I]A^{-1}B" title="M(t) = exp(At)M(0) +[exp(At)-I]A^{-1}B" />

where $M(0)$ is the initial condition of the system and I is the identity matrix. This expression can be used to describe a system when it is in free precession; e.g., when $\omega_1 = 0$. In that case, equation 4 can be reduced by noting that the z-component is decoupled from the $x$- and $y$-components, which results in the following
equation. 

<img src="https://latex.codecogs.com/svg.image?M_z&space;=&space;\begin{bmatrix}&space;M_{zf}&space;&&space;&space;M_{zm}&space;\\&space;\end{bmatrix}^T&space;" title="M_z = \begin{bmatrix} M_{zf} & M_{zm} \\ \end{bmatrix}^T " />

therefore:

<img src="https://latex.codecogs.com/svg.image?M_z(t)&space;=&space;exp(A_zt)M(0)&space;&plus;[I-exp(A_zt)]M(0)" title="M_z(t) = exp(A_zt)M(0) +[I-exp(A_zt)]M(0)" />

and:

<img src="https://latex.codecogs.com/svg.image?M_0&space;=&space;\begin{bmatrix}&space;M_{0f}&space;&&space;&space;M_{0m}&space;\\&space;\end{bmatrix}^T" title="M_0 = \begin{bmatrix} M_{0f} & M_{0m} \\ \end{bmatrix}^T" />

Given that the z-component of A comprise the following with $R_{RF} = 0$:

A_z = 
\begin{bmatrix}
    -(R_{1f} + k_{mf})  & k_{mf} \\
   k_{mf}                    & -(R_{1m}+R_{RF}+k_{mf}) \\
 \end{bmatrix}
 
 Expanding exponentials, the biexponential becomes:

<img src="https://latex.codecogs.com/svg.image?A_z&space;=&space;\begin{bmatrix}&space;&space;&space;exp(-\lambda^&plus;t)&space;&space;&&space;0&space;\\&space;&space;0&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;exp(-\lambda^-t)&space;\\&space;\end{bmatrix}&space;U^{-1}M(0)&plus;&space;\begin{pmatrix}I-U\begin{bmatrix}&space;&space;&space;exp(-\lambda^&plus;t)&space;&space;&&space;0&space;\\&space;&space;0&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;exp(-\lambda^-t)&space;\\&space;\end{bmatrix}&space;&space;\end{pmatrix}&space;M_0" title="A_z = \begin{bmatrix} exp(-\lambda^+t) & 0 \\ 0 & exp(-\lambda^-t) \\ \end{bmatrix} U^{-1}M(0)+ \begin{pmatrix}I-U\begin{bmatrix} exp(-\lambda^+t) & 0 \\ 0 & exp(-\lambda^-t) \\ \end{bmatrix} \end{pmatrix} M_0" />

where $\lambda^{+/-}$ are the negative eigenvalues of $A_z$ and $U$ is a matrix whose columns are the corresponding eigenvectors. Equation \ref{equ:final} tells us that $M_{zf}$ recovers as a bi-exponential function that is determined by the fast and slow rate constants $\lambda^{+}$ and $\lambda^{-}$, respectively. Using this formulations, we can determine the qMT parameters for PSR and $k_{mf}$.

The final signal solution is:

<img src="https://latex.codecogs.com/svg.image?M_z(t_I,t_D)&space;=&space;[exp(At_I)&space;\times&space;S&space;(I&space;-&space;exp(At_D))&plus;(I&space;-&space;exp(At_I)]M_0" title="M_z(t_I,t_D) = [exp(At_I) \times S (I - exp(At_D))+(I - exp(At_I)]M_0" />

where $M_0$ is the longitudinal magnetization vector, $I$ is the identity matrix, $M_0$ is the vector of equilibrium magnetizations, and S accounts for the effect of the inversion pulse on each pool. 

<img src="https://latex.codecogs.com/svg.image?S&space;=&space;diag(S_f,S_m)" title="S = diag(S_f,S_m)" />

When $S_f =1$ there is complete inversion of $M_{fm}$ and when $S_m = -1$ there is no saturation of the magnetization $M_{zm}$

### General formulation
_This form was used in the Julia code_


<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\space\frac{M_f(t)}{M_{f\infty}}=b_f^+exp(-R_1^+t)+b_f^-exp(-R_1^-t)+1" title="\Equation 1" />
</p>

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\space\2R_1^\pm=R_{1f}+R_{1m}+k_{fm}+k_{mf}\pm\sqrt{(R_{1f}-R_{1m}+k_{fm}-k_{mf})^2+4k_{fm}k_{mf}}" title="\Equation 2" />
</p>

and

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=b_f^\pm=\pm\frac{\begin{bmatrix}\frac{M_f(0)}{M_{0f}}-1\end{bmatrix}(R_1^--R_1^\pm)&plus;\begin{bmatrix}\frac{M_f(0)-M_m(0)}{M_{0f}-M_{0f}}\end{bmatrix}k_{fm}}{R_1^&plus;-R_1^-}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?b_f^\pm=\pm\frac{\begin{bmatrix}\frac{M_f(0)}{M_{0f}}-1\end{bmatrix}(R_1^--R_1^\pm)&plus;\begin{bmatrix}\frac{M_f(0)-M_m(0)}{M_{0f}-M_{0f}}\end{bmatrix}k_{fm}}{R_1^&plus;-R_1^-}" title="\Equation 3" /></a>
</p>

_Notes on markdown_
To embed math, easiest to go to this website. https://www.codecogs.com/latex/eqneditor.php
