#!/bin/bash

path=/mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618/
SIR_4D_DATA=/mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618//SIR_DATA.nii.gz
brain_mask=/mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618//brain_mask.nii.gz
julia ./SIR_fit.jl --TI 15 15 278 1007 --TD 684 4121 2730 10 --SIR_Data /mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618//SIR_DATA.nii.gz --SIR_brainMask /mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618//brain_mask.nii.gz --kmf 12.5 --Sm 0.83
