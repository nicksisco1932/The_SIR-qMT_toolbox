function [r2simg,r2img] = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname,DSC)
    TE1 = DSC.Parms.TE1;
    TE2 = DSC.Parms.TE2;
    TE3 = DSC.Parms.TE3;
    TE4 = DSC.Parms.TE4;
    TE5 = DSC.Parms.TE5;
    julia_cmd = '"/Users/nicks/AppData/Local/Programs/Julia 1.5.3/bin/julia.exe"';
    SAGE_Fit_julia = "C:\Users\nicks\Documents\Github/The_MRI_toolbox/SAGE_biexp_fit.jl";
    cmd = sprintf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s",julia_cmd,SAGE_Fit_julia,te1,te2,te3,te4,te5,b_fname,TE1,TE2,TE3,TE4,TE5);
    system(cmd)
    r2s = sprintf("%s/SAGE_R2s_Brain_julia.nii.gz",base_path);
    r2 = sprintf("%s/SAGE_R2_Brain_julia.nii.gz",base_path);
    r2simg = loader(r2s);
    r2img = loader(r2);
end