#!/bin/bash

path="Data_used_for_publication/SIR_brains/"
SIR_4D_DATA=$path/SIR_DATA.nii.gz
brain_mask=$path/brain_mask.nii.gz
julia ./SIR_fit.jl --TI 15 15 278 1007 --TD 684 4121 2730 10 --SIR_Data $SIR_4D_DATA --SIR_brainMask $brain_mask --kmf 12.5 --Sm 0.83
