# DEcomposition and Component Analysis of Exponential Signals (DECAES)

## Table of contents

## Introduction

DECAES provides tools for decomposing multi-exponential signals which arise from multi spin-echo magnetic resonance imaging (MRI) scans into exponential components.
The main decomposition method used is an inverse Laplace transform-based technique which involves solving the regularized nonnegative least squares (NNLS) inverse problem

```math
X = \mathrm{argmin}_{x \ge 0} ||Cx - d||_2^2 + \mu^2 ||x||_2^2
```
