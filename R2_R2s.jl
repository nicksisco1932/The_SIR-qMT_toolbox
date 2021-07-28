using Pkg	
try
    println("If is first time you ran the code. It will take a minute to precompile.")
    @eval using NIfTI; 
    @eval using LsqFit;
    @eval using Printf
    @eval using ArgParse;
    @eval using Statistics;
    @eval using Optim;
    @eval using LinearAlgebra;
    Pkg.precompile()
catch e
    # not found; install and try loading again
    Pkg.add("NIfTI")
    Pkg.add("LsqFit")
    Pkg.add("Printf")
    Pkg.add("ArgParse")
    Pkg.add("Optim")
    Pkg.add("LinearAlgebra")
    @eval using NIfTI; 
    @eval using LsqFit;
    @eval using Printf
    @eval using ArgParse;
    @eval using Statistics;
    @eval using Optim;
    @eval using LinearAlgebra;
end
using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Statistics;
using Optim;
using LinearAlgebra;

SAGE_base="/mnt/box/MSPerfusion_Data/SAGE_dcm_converted/SAGE_niftis/"
SAGE_TE1_fname = @sprintf("%sPT1319016_TE1_img_w_Skull.nii.gz",SAGE_base)
SAGE_TE1_img = niread(SAGE_TE1_fname);
SAGE_TE2_fname = @sprintf("%sPT1319016_TE2_img_w_Skull.nii.gz",SAGE_base)
SAGE_TE2_img = niread(SAGE_TE2_fname);
SAGE_TE3_fname = @sprintf("%sPT1319016_TE3_img_w_Skull.nii.gz",SAGE_base)
SAGE_TE3_img = niread(SAGE_TE3_fname);
SAGE_TE4_fname = @sprintf("%sPT1319016_TE4_img_w_Skull.nii.gz",SAGE_base)
SAGE_TE4_img = niread(SAGE_TE4_fname);
SAGE_TE5_fname = @sprintf("%sPT1319016_TE5_img_w_Skull.nii.gz",SAGE_base)
SAGE_TE5_img = niread(SAGE_TE5_fname);
MASK_fname = @sprintf("%s/bPT1319016_preb_mask.nii.gz",SAGE_base)
MASK_img = niread(MASK_fname);
MASK = Array{Bool}(MASK_img.raw);
ADC_path = "/mnt/c/Users/nicks/Documents/MRI_data/DTI_VSI/PT16/reg2sage/"
ADC_fname = @sprintf("%s/Reg2Warped.nii.gz",ADC_path)
ADC_img = niread(ADC_fname);
adc=ADC_img.raw;

echos = [8.8;26.0;49.0;66.0;88.0]*1E-3

function r2_r2s(echos::Array{Float64},data::Array{Float64},mask::Array{Bool})
    te1=echos[1];
    te2=echos[2];
    te3=echos[3];
    te4=echos[4];
    te5=echos[5];
    
    Y = zeros(4,4);
    Y[1,:] = [1  0  -te1 0]
    Y[2,:] = [1  0  -te2 0]
    Y[3,:] = [1 -1  -te5+te3 te5-2*te3]
    Y[4,:] = [1 -1  0 -te5]

    Y = zeros(5,4);
    Y[1,:] = [1  0  -te1 0 ]
    Y[2,:] = [1  0  -te2 0 ]
    Y[3,:] = [1 -1  -te5+te3 te5-2*te3 ]
    Y[4,:] = [1 -1  -te5+te4 te5-2*te4 ]
    Y[5,:] = [1 -1  0 -te5 ]

    
    nx,ny,nz,nt,ne = size(data)
    S = zeros(5,nx*ny*nz,nt);
    vec_mask = reshape(mask,nx*ny*nz)

    vec_te1 = reshape(data[:,:,:,:,1],nx*ny*nz,nt);
    vec_te2 = reshape(data[:,:,:,:,2],nx*ny*nz,nt);
    vec_te3 = reshape(data[:,:,:,:,3],nx*ny*nz,nt);
    vec_te4 = reshape(data[:,:,:,:,4],nx*ny*nz,nt);
    vec_te5 = reshape(data[:,:,:,:,5],nx*ny*nz,nt);
    S[1,:,:]=log.(vec_te1);
    S[2,:,:]=log.(vec_te2);
    S[3,:,:]=log.(vec_te3);
    S[4,:,:]=log.(vec_te4);
    S[5,:,:]=log.(vec_te5);

    
    
    A = zeros(nx*ny*nz,4,nt)

    for ii in nx*ny*nz
        # if vec_mask[ii]
            for jj = 1:nt
                A[ii,:,jj] = Y \ S[:,ii,jj]
            end
        # end
    end

    

    A = zeros(4,nx*ny*nz,nt);
    for ii in 1:nt
        A[:,:,ii] = Y \ S[:,:,ii]
    end

    Am = reshape(A,4,nx,ny,nz,150)
    S0I = Am[1,:,:,:,:];
    δ = Am[2,:,:,:,:];
    R₂star = Am[3,:,:,:,:];
    R₂= Am[4,:,:,:,:];

    return S0I,δ,R₂star,R₂
end

function SAGE_biexp4p_NS(te,p)
    τ=te;
    TEₛₑ = τ[end];
    SI_sage = similar(τ);
    for (n,ii) in enumerate(τ)
        if ii<TEₛₑ/2
            SI_sage[n] = p[1]*exp(-ii*p[3]);
        else
            SI_sage[n] = (p[2])*exp(-TEₛₑ*(p[3]-p[4]))*exp(-ii*(2*p[4]-p[3]));
        end
    end
    return SI_sage
end
function nlsfit(f::Function, xvalues::Array{Float64},yvalues::Array{Float64},guesses::Array{Float64})
    # fit = curve_fit(f,xvalues,yvalues,guesses,lower=[-Inf, -Inf, 0.0, 0],upper=[Inf, Inf, 200, 200];autodiff=:finiteforward)
    fit = curve_fit(f,xvalues,yvalues,guesses;autodiff=:finiteforward)
    return fit.param
end



ste1 = Array{Float64}(SAGE_TE1_img.raw);
ste2 = Array{Float64}(SAGE_TE2_img.raw);
ste3 = Array{Float64}(SAGE_TE3_img.raw);
ste4 = Array{Float64}(SAGE_TE4_img.raw);
ste5 = Array{Float64}(SAGE_TE5_img.raw);
nx,ny,nz,nt=size(ste1)
ne=5;
data=zeros(nx,ny,nz,nt,ne);
data[:,:,:,:,1]=ste1;
data[:,:,:,:,2]=ste2;
data[:,:,:,:,3]=ste3;
data[:,:,:,:,4]=ste4;
data[:,:,:,:,5]=ste5;

S0I,δ,R₂star,R₂ = r2_r2s( echos,data,MASK);

# tot_voxels = nx*ny*nz*nt;
# vec_data = reshape(data,tot_voxels,ne);
# vec_a = reshape(S0I,tot_voxels);
# vec_b = reshape(δ,tot_voxels);
# vec_c = reshape(R₂star,tot_voxels);
# vec_d = reshape(R₂,tot_voxels);
# rmask = repeat(MASK,1,1,1,150);
# vec_mask = reshape(rmask,tot_voxels);
# vec_out = zeros(tot_voxels,4);
# IND=findall(x->x==true,vec_mask);
# l=length(IND)
# @time Threads.@threads for ii in IND
#         vec_out[ii,:] = nlsfit(SAGE_biexp4p_NS,echos,vec_data[ii,:],[vec_a[ii],vec_b[ii],vec_c[ii],vec_d[ii]])
#         # println("Done with $ii of $l")        
#         #  vec_data[ii,:];[vec_a[ii];vec_b[ii];vec_c[ii];vec_d[ii]]
# end

IND=findall(x->x==true,MASK);

p0 = [1000.0,0,25,10]
l=length(IND)
OUT=zeros(nx,ny,nz,nt,4);
residuals = similar(data);
@time for jj in 1:nt
    for ii in IND
        # OUT[ii,jj,:],residuals[ii,jj,:] = nlsfit(SAGE_biexp4p_NS,echos,data[ii,jj,:],[ S0I[ii,jj],δ[ii,jj], R₂star[ii,jj],R₂[ii,jj] ])
        OUT[ii,jj,:] = nlsfit(SAGE_biexp4p_NS,echos,data[ii,jj,:],[ S0I[ii,jj],δ[ii,jj], R₂star[ii,jj],R₂[ii,jj] ])
    end
    println("Done with $jj of $nt")
end

R₂star[findall(x->x.<0,R₂star)].=0;
R₂star[findall(x->isnan(x),R₂star)].=0;
R₂star[findall(x->isinf(x),R₂star)].=0;
R₂[findall(x->isnan(x),R₂)].=0;
R₂[findall(x->isinf(x),R₂)].=0;
R₂[findall(x->x.<0,R₂)].=0;

IND=findall(x->x==1,MASK);

for ii in IND
    for jj in 1:150
        if R₂star[ii,jj] > 1000 || R₂star[ii,jj] < 0
            R₂star[ii,jj]=0
        end
        if R₂[ii,jj] > 1000 ||  R₂[ii,jj] < 0
            R₂[ii,jj]=0
        end
    end
end

dynamic=27

maxR2 = R₂[:,:,:,dynamic];
maxR2s = R₂star[:,:,:,dynamic];

function vsi(r2s::Array{Float64},r2::Array{Float64},adc::Array{Float64},mask::Array{Bool})

    Δχ = 0.264e-6; # 0.264 ppm Spees et al. for RBC 
    γ = 267.522E6;  # rad/s/T
    B₀ = 3               # T

    nx,ny,nz = size(r2s)
    vsi_img = similar(r2s)
    for ii in 1:nx
        for jj in 1:ny
            for kk in 1:nz
                if mask[ii,jj,kk]
                    vsi_img[ii,jj,kk] = 0.425 * sqrt(adc[ii,jj,kk] / (γ*Δχ*B₀) ) * (r2s[ii,jj,kk] / r2[ii,jj,kk])^(3/2)
                end
            end
        end
    end
    return vsi_img
end

vsi_img = vsi(maxR2s,maxR2,adc,MASK)

tmp1 = voxel_size(SAGE_TE1_img.header)[1]
tmp2 = voxel_size(SAGE_TE1_img.header)[2]
tmp3 = voxel_size(SAGE_TE1_img.header)[3]

R2s_fname_log = @sprintf("%s/SAGE_R2s.nii.gz",ADC_path)
R2_fname_log = @sprintf("%s/SAGE_R2.nii.gz",ADC_path)
temp1 = NIVolume(MASK.* R₂star;voxel_size=(tmp1,tmp2,tmp3));
niwrite(R2s_fname_log,temp1)
temp2 = NIVolume(MASK.*R₂;voxel_size=(tmp1,tmp2,tmp3));
niwrite(R2_fname_log,temp2)

tmp1 = voxel_size(ADC_img.header)[1]
tmp2 = voxel_size(ADC_img.header)[2]
tmp3 = voxel_size(ADC_img.header)[3]


ADC_fname = @sprintf("%s/ADC_2_SAGE.nii.gz",ADC_path)
temp1 = NIVolume(adc;voxel_size=(tmp1,tmp2,tmp3));
niwrite(ADC_fname,temp1)

vsi_img[findall(x->x.>100,vsi_img)].=0;
vsi_img[findall(x->x.<0,vsi_img)].=0;

VSI_fname = @sprintf("%s/VSI_SAGE.nii.gz",ADC_path)
temp1 = NIVolume(MASK.* vsi_img;voxel_size=(tmp1,tmp2,tmp3));
niwrite(VSI_fname,temp1)


tmp1 = voxel_size(SAGE_TE1_img.header)[1];tmp2 = voxel_size(SAGE_TE1_img.header)[2];tmp3 = voxel_size(SAGE_TE1_img.header)[3]
r2star = OUT[:,:,:,:,3];
r2 = OUT[:,:,:,:,4];
res1 = residuals[:,:,:,:,1];
res2 = residuals[:,:,:,:,2];
res3 = residuals[:,:,:,:,3];
res4 = residuals[:,:,:,:,4];
res5 = residuals[:,:,:,:,5];



R2s_fname_log = @sprintf("%s/SAGE_R2s_fit.nii.gz",ADC_path)
R2_fname_log = @sprintf("%s/SAGE_R2_fit.nii.gz",ADC_path)
# res1_fname = @sprintf("%s/SAGE_res1.nii.gz",ADC_path)
# res2_fname = @sprintf("%s/SAGE_res2.nii.gz",ADC_path)
# res3_fname = @sprintf("%s/SAGE_res3.nii.gz",ADC_path)
# res4_fname = @sprintf("%s/SAGE_res4.nii.gz",ADC_path)
# res5_fname = @sprintf("%s/SAGE_res5.nii.gz",ADC_path)
# temp1 = NIVolume(res1;voxel_size=(tmp1,tmp2,tmp3));
# niwrite(res1_fname,temp1)
# temp2 = NIVolume(res2;voxel_size=(tmp1,tmp2,tmp3));
# niwrite(res2_fname,temp2)
# temp1 = NIVolume(res3;voxel_size=(tmp1,tmp2,tmp3));
# niwrite(res3_fname,temp1)
# temp2 = NIVolume(res4;voxel_size=(tmp1,tmp2,tmp3));
# niwrite(res4_fname,temp2)
# temp1 = NIVolume(res5;voxel_size=(tmp1,tmp2,tmp3));
# niwrite(res5_fname,temp1)


temp1 = NIVolume(r2star;voxel_size=(tmp1,tmp2,tmp3));
niwrite(R2s_fname_log,temp1)
temp2 = NIVolume(r2;voxel_size=(tmp1,tmp2,tmp3));
niwrite(R2_fname_log,temp2)



R2s_fname_log = @sprintf("%s/SAGE_R2s_fit_v_linearAlgebra.nii.gz",ADC_path)
R2_fname_log = @sprintf("%s/SAGE_R2_fit_v_linearAlgebra.nii.gz",ADC_path)
temp1 = NIVolume(MASK.* R₂star .- r2star ;voxel_size=(tmp1,tmp2,tmp3));
niwrite(R2s_fname_log,temp1)
temp2 = NIVolume(MASK.* R₂ .- r2 ;voxel_size=(tmp1,tmp2,tmp3));
niwrite(R2_fname_log,temp2)