#=
simulating_R2R2s
=#

using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Optim

include("./utils.jl")


te1=0.007653;
te2=0.027969;
te3=0.059074;
te4=0.07939;
te5=0.099706;
echos=[te1,te2,te3,te4,te5]


using Plots
X0 = [1,150,80,500];Yy = SAGE_biexp3p_d(echos,X0);scatter(echos,Yy);plot!(echos,Yy)

a = LinRange(50,200,128);
b = LinRange(40,150,128);
mat_img = zeros(128,128,5);
R2_map = zeros(128,128);
R2s_map = zeros(128,128);
SNR=150;
for ii in 1:128
    for jj in 1:128
        X0 = [1,a[ii],b[jj],500];
        Yy = SAGE_biexp3p_d(echos,X0);
        R2s_map[ii,jj] = a[ii]
        R2_map[ii,jj] = b[ii]
        mat_img[ii,jj,:] = Yy+rand(5)/SNR
    end
end

niwrite("./R2_R2s_sim.nii.gz",NIVolume(mat_img))
niwrite("./R2s_sim.nii.gz",NIVolume(R2s_map))
niwrite("./R2_sim.nii.gz",NIVolume(R2_map))

function optim_fitty(f::Function, xdata::Vector{Float64},data::Vector{Float64},x0::Vector{Float64})
    # fitdata = Optim.minimizer( optimize(b -> sqerrorSAGE(b, xdata,  data), x0))
    fitdata = curve_fit(f,xdata,data,x0;autodiff=:finiteforward) # this is faster and has good results
    return fitdata
end

x0 = [10.0,150,80,10];
R2sfitdata = zeros(128,128);
R2fitdata = zeros(128,128);
for ii in 1:128
    for jj in 1:128
        fitY = optim_fitty(SAGE_biexp3p_d,echos,mat_img[ii,jj,:],x0)
        parm = fitY.param;
        R2sfitdata[ii,jj] = parm[2]
        R2fitdata[ii,jj] = parm[3]
    end
end
niwrite("./R2s_fit.nii.gz",NIVolume(R2sfitdata))
niwrite("./R2_fit.nii.gz",NIVolume(R2fitdata))
niwrite("./R2s_diff.nii.gz",NIVolume(R2s_map - R2sfitdata))
niwrite("./R2_diff.nii.gz",NIVolume(R2_map - R2fitdata))

