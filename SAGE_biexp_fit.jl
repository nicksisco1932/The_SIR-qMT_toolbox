#=
    Author: Nicholas J. Sisco, Ph.D. 
    Email: X@barrowneuro.org where X = nicholas.sisco
    Affiliation: Barrow Neurological Institue in Phoenix, AZ

    Title: SAGE R2star and R2 fitting using least squares and Levenberg-Marquardt 
    - This is a basic script written in Julia for fitting signal acquired from SAGE to a piecewise exponential 
    equation described in (1). 

    Requirements:
    - Julia 1.5

    Input and Output Parameters for Code

    utils.jl needs to be in the same directory as this script

    User supplied data:
        - Command line usage
        $ julia -t <threads #> SAGE_biexp_fit.jl <path to Python processed data>

        - The echo times are still hard coded

    Output:
        - The output from this script are two nifti files for the fitting parameters R2star and R2.


    Future Updates TO Do:
        1) user defined echo times
        2) user defined preprocessing niftis, i.e. flexible input names
        3) main function within utils
        4) module creation
        5) unit testing
        6) depolyable docker

    References:
    1. H. Schmiedeskamp, M. Straka, R. D. Newbould, G. Zaharchuk, J. B. Andre, J. M. Olivot, M. E. Moseley, G. W. Albers, R. Bammer, Combined spin- and gradient-echo perfusion-weighted imaging. Magn. Reson. Med. 68, 30–40 (2012).

    Change log:
    "history (of nifti library changes):\n"
    "\n",
    "0.0  June 20, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - initial version \n"

=#
using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Optim

include("./utils.jl")

parsed_args = parse_commandline();
for (out, val) in parsed_args
    println(" $out => $val")
end

base = parsed_args["input"]

paths = [ @sprintf("%sPT1319001_TE1_img_w_Skull.nii.gz",base),
            @sprintf("%sPT1319001_TE2_img_w_Skull.nii.gz",base),
            @sprintf("%sPT1319001_TE3_img_w_Skull.nii.gz",base), 
            @sprintf("%sPT1319001_TE4_img_w_Skull.nii.gz",base),
            @sprintf("%sPT1319001_TE5_img_w_Skull.nii.gz",base)
            ]

# paths = [ @sprintf("%sMFA_10_angles.nii.gz",base)]
b_fname = @sprintf("%sbPT1319001_preb_mask.nii.gz",base)

DATA,MASK,nx,ny,nz,ne,nt=load_sage_data(paths,b_fname);
tot_voxels = nx*ny*nz
vec_mask = reshape(MASK[:,:,:,1,1],tot_voxels,1);

vec_data = reshape(DATA,tot_voxels,ne,nt)


thr = Threads.nthreads()
if thr==16
    Threads.nthreads() = 8
    return thr = 8
end
println("Fitting with $thr threads")

te1=0.007653;
te2=0.027969;
te3=0.059074;
te4=0.07939;
te5=0.099706;
echos=[te1,te2,te3,te4,te5]
TR = 1800/1000;

X0=[1,100.0,50,1]
IND=findall(x->x.>0,vec_mask[:,1]);

# function fitty(f::Function,xData::Vector{Float64},Ydata::Array{Float64},x0::Vector{Float64})
#     # fitY = Optim.minimizer( optimize(b -> f(b, xData,  Ydata), x0))
    
#     ub=[Inf,500.0,490,Inf];
#     lb=[0.0,0,0,0];
#     fitparm = curve_fit(f,xData,Ydata,x0;autodiff=:finiteforward)
#     fitY = fitparm.param
#     return fitY #array of 4
# end

# function fitSAGE(xData::Vector{Float64},mask::Array{Bool},Ydata::Array{Float64},n::Int64,nx::Int64,ny::Int64,nz::Int64,ne::Int64,nt::Int64,x0::Vector{Float64})
# function fitSAGE(xData::Vector{Float64},ind::Vector{Int64},Ydata::Array{Float64},n::Int64,nx::Int64,ny::Int64,nz::Int64,ne::Int64,NT::Int64,x0::Vector{Float64})
#     out = zeros(n,4,NT);
    
#     @time for jj in 1:NT
#         # jj=1
#         for ii in ind
#                 # fitY = fitty(sqerrorSAGE,xData,Ydata[ii,:,jj],x0)
#                 fitY = fitty(SAGE_biexp3p_d,xData,Ydata[ii,:,jj],x0)
#                 out[ii,:,jj] = fitY;
#         end
#     end    
#     # @time for jj in 1:nt
#     #     for ii = 1:n
#     #         if mask[ii]
#     #             # fitY = fitty(sqerrorSAGE,xData,Ydata[ii,:,jj],x0)
#     #             fitY = fitty(SAGE_biexp3p_d,xData,Ydata[ii,:,jj],x0)
#     #             out[ii,:,jj] = fitY;
#     #         end
#     #     end
#     # end
#     return out
# end


function optim_fitty(f::Function, xdata::Vector{Float64},data::Vector{Float64},x0::Vector{Float64})
    # fitdata = Optim.minimizer( optimize(b -> sqerrorSAGE(b, xdata,  data), x0))
    fitdata = curve_fit(f,xdata,data,x0;autodiff=:finiteforward) # this is faster and has good results
    return fitdata
end

fitdata = zeros(tot_voxels,4,nt);
@time for JJ=1:150;
    for ii in IND
        temp = vec_data[ii,:,JJ];
        fitY = optim_fitty(SAGE_biexp3p_d, echos, temp, X0)
        # fitY = optim_fitty(echos, temp, X0)
        fitdata[ii,:,JJ] = fitY.param
    end
end
#-------------------------------------  
Xv = Array{Float64}(fitdata);
Xv = reshape(Xv,nx,ny,nz,4,nt);
R2s = zeros(nx,ny,nz,nt);
R2 = zeros(nx,ny,nz,nt);
R2s = Xv[:,:,:,2,:];
R2 = Xv[:,:,:,3,:];
R2s[findall(x->x.>1000,R2s)].=0
R2s[findall(x->x.<0,R2s)].=0
R2[findall(x->x.>1000,R2)].=0
R2[findall(x->x.<0,R2)].=0
#-------------------------------------

d=niread(paths[1]);

tmp1 = voxel_size(d.header)[1]
tmp2 = voxel_size(d.header)[2]
tmp3 = voxel_size(d.header)[3]

R2s_fname = @sprintf("%sSAGE_R2s_Brain_julia.nii.gz",base)
R2_fname = @sprintf("%sSAGE_R2_Brain_julia.nii.gz",base)
niwrite(R2s_fname,NIVolume(R2s;voxel_size=(tmp1,tmp2,tmp3)))
niwrite(R2_fname,NIVolume(R2;voxel_size=(tmp1,tmp2,tmp3)))

println("The R2s is saved at $R2s_fname")
println("The R2 is saved at $R2_fname")
