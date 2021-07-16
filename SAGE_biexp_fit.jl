#=
    Author: Nicholas J. Sisco, Ph.D. 
    Email: X@barrowneuro.org where X = nicholas.sisco
    Affiliation: Barrow Neurological Institue in Phoenix, AZ

    Title: SAGE R2star and R2 fitting using least squares and Levenberg-Marquardt 
    - This is a basic script written in Julia for fitting signal acquired from SAGE to a piecewise exponential 
    equation described in (1). 

    Requirements:
    - Julia 1.5
    - utils.jl
    - five echo times with corresponding images
    - one brain image

    Input and Output Parameters for Code

    utils.jl needs to be in the same directory as this script

    User supplied data:
        - Command line usage
        $ julia -t <threads #> SAGE_biexp_fit.jl <path to Python processed data>

        - The echo times are still hard coded

    To Run this within MATLAB:
        % CHANGE PATHS AND NAMES ACCORDINGLY. NOT TESTED ON WINDOWS
        % To set: Julia path, This scrip path, your files pathes.
        base_path = "/Volumes/GraceHopper/HFTT/1420021/NS/";
        te1 = sprintf("%s1420021_LEFT_e1_tshift_bet.nii.gz",base_path);
        te2 = sprintf("%s1420021_LEFT_e2_tshift_bet.nii.gz",base_path);
        te3 = sprintf("%s1420021_LEFT_e3_tshift_bet.nii.gz",base_path);
        te4 = sprintf("%s1420021_LEFT_e4_tshift_bet.nii.gz",base_path);
        te5 = sprintf("%s1420021_LEFT_e5_tshift_bet.nii.gz",base_path);
        b_fname = sprintf("%s1420021_LEFT_e2_tshift_tmean_bet_mask.nii.gz",base_path);
        R2s_fname = sprintf("%sSAGE_R2s_Brain_julia.nii.gz",path);
        R2_fname = sprintf("%sSAGE_R2s_Brain_julia.nii.gz",path);
        %% Function definition from here.
        [r2simg,r2img] = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname);
            function [r2simg,r2img] = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname)
            julia_cmd = "/usr/local/bin/julia";
            SAGE_Fit_julia = "/Users/StokesLab/Documents/EK/Code/The_MRI_toolbox/SAGE_biexp_fit.jl";
            unix(sprintf("%s %s %s %s %s %s %s %s",julia_cmd,SAGE_Fit_julia,te1,te2,te3,te4,te5,b_fname))
            r2s = sprintf("%s/SAGE_R2s_Brain_julia.nii.gz",base_path);
            r2 = sprintf("%s/SAGE_R2_Brain_julia.nii.gz",base_path);
            info = niftiinfo(r2s);
            r2simg = niftiread(info);
            info = niftiinfo(r2);
            r2img = niftiread(info);
        end

    Output:
        - The output from this script are two nifti files for the fitting parameters R2star and R2.


    Future Updates TO Do:
        1) user defined echo times
        2) user defined preprocessing niftis, i.e. flexible input names
        3) module creation
        4) unit testing
        5) depolyable docker

    References:
    1. H. Schmiedeskamp, M. Straka, R. D. Newbould, G. Zaharchuk, J. B. Andre, J. M. Olivot, M. E. Moseley, G. W. Albers, R. Bammer, Combined spin- and gradient-echo perfusion-weighted imaging. Magn. Reson. Med. 68, 30–40 (2012).

    Change log:
    "history (of nifti library changes):\n"
    "\n",
    "0.0  June 20, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - initial version \n"

    "\n",
    "0.1  June 21, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - Updated I/O \n"
    "\n",
    "0.2  June 23, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - Added a progress meter and tested multi threading \n"
    "\n",
    "0.3  June 24, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - Added some MATLAB coding implimentation \n"

=#

using Pkg	
try
    println("If is first time you ran the code. It will take a minute to precompile.")
    @eval using NIfTI; 
    @eval using LsqFit;
    @eval using Printf
    @eval using ArgParse;
    @eval using Statistics;
    @eval using Optim;
    Pkg.precompile()
catch e
    # not found; install and try loading again
    Pkg.add("NIfTI")
    Pkg.add("LsqFit")
    Pkg.add("Printf")
    Pkg.add("ArgParse")
    @eval using NIfTI; 
    @eval using LsqFit;
    @eval using Printf
    @eval using ArgParse;
    @eval using Statistics;
    @eval using Optim;
end

using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Statistics;
using Optim;

include("./utils.jl")

function commandline()
        

    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "--TE_nii_names"
            nargs=5
            arg_type = String
            required = true
        # "--TE_nii_2"
        #     required = true
        # "--TE_nii_3"
        #     required = true
        # "--TE_nii_4"
        #     required = true
        # "--TE_nii_5"
        #     required = true
        "--SAGE_nii_brainMask"
            required = true
        "--echos"
            nargs = 5
            arg_type = Float64
            help = "Echo times, 5 required"
    end

    println(parse_args(settings))

    # parsed_args = TE_parse_commandline();
    for (out, val) in parse_args(settings)
        println(" $out => $val")
    end
    return parse_args(settings)
end

A = commandline()

function fguess(guess::Vector{Float64},y::Vector{Float64},vec_r2s_log::Float64,vec_r2_log::Float64)
    guess[1] = maximum(y)
    guess[2] = vec_r2s_log
    if vec_r2_log<0
        guess[3]=10
    else
        guess[3] = vec_r2_log # sometimes the log ratio makes this negative.
    end
    guess[4] = 0
    return guess
end

function sage_ns(echotimes::Vector{Float64},p::Vector{Float64})
    # this is the fastest way to get the signal 
    τ=echotimes;
    TEₛₑ=τ[end];
    S₀I = p[1]
    R₂star = p[2]
    R₂ = p[3]
    S₀II = p[4]

    M = similar(τ)
    for (n,k) in enumerate(τ)
        if k < TEₛₑ/2
            M[n] = S₀I * exp(-k * R₂star)
        elseif   k >= TEₛₑ/2
            M[n] = S₀II * exp(-TEₛₑ * (R₂star - R₂) * exp(-k * (2*R₂-R₂star) ) )
        end
    end
    return M
end

function sqerrorSAGE2(betas::Vector{Float64}, X::Vector{Float64}, Y::Vector{Float64})
    err = 0.0
    pred_i = sage_ns(X,betas)
    for (n,ii) in enumerate( Y)
        
        err += (ii -pred_i[n]).^2
    end
    return err
end

function linearize(nx::Int64,ny::Int64,nz::Int64,ne::Int64,nt::Int64,data::Array{Float64},mask::Array{Bool},echos::Vector{Float64})
    te1=echos[1];
    te2=echos[2];
    te3=echos[3];
    te4=echos[4];
    te5=echos[5];
    tot_voxels=nx*ny*nz
        
    newDATA=zeros(nx,ny,nz,nt,ne);

    for ii in 1:nx
        for jj in 1:ny
            for kk in 1:nz
                for tt in 1:nt
                    if mask[ii,jj,kk]
                        newDATA[ii,jj,kk,tt,:] = data[ii,jj,kk,:,tt]
                    end
                end
            end
        end
    end


    # STE1=zeros(nx,ny,nz);
    # STE2=similar(STE1);
    # STE5=similar(STE1);
    # begin
    #     for ii in 1:nx
    #         for jj in 1:ny
    #             for kk in 1:nz
    #                 STE1[ii,jj,kk] = mean(data[ii,jj,kk,1,1:15]);
    #                 STE2[ii,jj,kk] = mean(data[ii,jj,kk,2,1:15]);
    #                 STE5[ii,jj,kk] = mean(data[ii,jj,kk,5,1:15]);
    #             end
    #         end
    #     end
    # end

    # STE1_pre = repeat(STE1,1,1,1,150);
    # STE2_pre = repeat(STE1,1,1,1,150);
    # STE5_pre = repeat(STE1,1,1,1,150);

    ste1 = data[:,:,:,1,:];
    ste2 = data[:,:,:,2,:];
    ste5 = data[:,:,:,5,:];


    
    R₂star_log = similar(ste1);
    R₂_log = similar(ste1);

    IND=findall(x->x==true,mask);
    begin
        for ii in 1:nx
            for jj in 1:ny
                for kk in 1:nz
                    if mask[ii,jj,kk]
                        for bb = 1:nt
                            a = ste1[ii,jj,kk,bb];
                            # b = STE1_pre[ii,jj,kk,bb];
                            c = ste2[ii,jj,kk,bb];
                            # d = STE2_pre[ii,jj,kk,bb];
                            e = ste5[ii,jj,kk,bb];
                            # f = STE5_pre[ii,jj,kk,bb];
                            ste0 = a*(a/c)^(te1/(te2-te1));
                            R₂star_log[ii,jj,kk,bb] = log( (a) / (c)) /(te2-te1);
                            R₂_log[ii,jj,kk,bb] = log( (ste0 / e)) /te5
                        end
                    end
                end
            end
        end
    end


    vec_R2fit = zeros(tot_voxels,nt);
    vec_R2sfit = zeros(tot_voxels,nt);
    vec_R2_log = reshape(R₂_log,tot_voxels,nt);
    vec_R2s_log = reshape(R₂star_log,tot_voxels,nt);
    tmpOUT=Array{Float64}(zeros(nx*ny*nz,4,nt));
    # vec_R2fit = zeros(nx*ny*nz*nt,4);
    # vec_R2sfit = zeros(nx*ny*nz*nt,4);
    # vec_R2_log = reshape(R₂_log,nx*ny*nz*nt);
    # vec_R2s_log = reshape(R₂star_log,nx*ny*nz*nt);
    # tmpOUT=Array{Float64}(zeros(nx*ny*nz*nt,4));

    return vec_R2_log,vec_R2s_log,vec_R2sfit,vec_R2fit,tmpOUT,R₂star_log,R₂_log,newDATA
end

function nlsworker_sage(tot_voxels::Int64,vec_data::Array{Float64},vec_mask::Array{Bool},vec_R2s_log::Array{Float64},vec_R2_log::Array{Float64},x0::Vector{Float64},echos::Vector{Float64})
    tmpOUT=zeros(tot_voxels,4)
    @time for vox in 1:tot_voxels
        if vec_mask[vox]
            initial = fguess(x0,vec_data[vox,:],vec_R2s_log[vox],vec_R2_log[vox])
            tmpOUT[vox,:] = nlsfit(SAGE_biexp4p_d,echos,vec_data[vox,:],initial)
        end
    end
    return tmpOUT
end

function main(a)
    


    base = dirname(a["TE_nii_names"][1])
    # paths = [
    #     joinpath(base,basename(a["TE_nii_1"])),
    #     joinpath(base,basename(a["TE_nii_2"])),
    #     joinpath(base,basename(a["TE_nii_3"])),
    #     joinpath(base,basename(a["TE_nii_4"])),
    #     joinpath(base,basename(a["TE_nii_5"]))
    #     ]
    paths = [
        joinpath(base,basename(a["TE_nii_names"][1])),
        joinpath(base,basename(a["TE_nii_names"][2])),
        joinpath(base,basename(a["TE_nii_names"][3])),
        joinpath(base,basename(a["TE_nii_names"][4])),
        joinpath(base,basename(a["TE_nii_names"][5]))
        ]
    b_fname = joinpath(base,basename(a["SAGE_nii_brainMask"]))

    DATA,MASK,nx,ny,nz,ne,nt=load_sage_data(paths,b_fname);
    tot_voxels = nx*ny*nz
    vec_mask = reshape(MASK[:,:,:,1,1],tot_voxels,1);
    vec_data = reshape(DATA,tot_voxels,ne,nt);

    # MASK = Array{Bool}(niread(b_fname))
    # MASK = MASK.raw
    # vec_mask = reshape(MASK,tot_voxels,1)
    

#= Long comment
    te1=0.00782;
    te2=0.028769;
    te3=0.060674;
    te4=0.081622;
    te5=0.102571;
    echos=Vector{Float64}([te1,te2,te3,te4,te5])
 
    te1=parse(Float64,a["TE1"])
    te2=parse(Float64,a["TE2"])
    te3=parse(Float64,a["TE3"])
    te4=parse(Float64,a["TE4"])
    te5=parse(Float64,a["TE5"])
    echos=Array{Float64}([te1,te2,te3,te4,te5])

=# 
    echos = Array{Float64}(a["echos"])
    if echos[1]>1
        println("Converting to seconds")
        echos = echos.* 1E-3
    else
        echos = echos;
    end

    TR = 1800/1000;

    X0=[1000,100.0,50,100]
    IND=findall(x->x.>0,vec_mask[:,1]);

    vec_R2_log,vec_R2s_log,vec_R2sfit,vec_R2fit,tmpOUT,R₂star_log,R₂_log,newDATA = linearize(nx,ny,nz,ne,nt,DATA,MASK,echos)

    x0 = Vector{Float64}(zeros(4))
#= Trie multithreading and got weird maps
    # tot_voxels = nx*ny*nz*nt
    # vec_data = reshape(newDATA,tot_voxels,ne)
    # vec_mask = repeat(vec_mask,1,150);
    # @time begin
    #     Threads.@threads for vox in 1:tot_voxels
    #         if vec_mask[vox]
    #             initial = fguess(x0,vec_data[vox,:],vec_R2s_log[vox],vec_R2_log[vox])
    #             tmpOUT[vox,:] = nlsfit(SAGE_biexp4p_d,echos,vec_data[vox,:],initial)
    #         end
    #     end
    # end
=#


    # tmpOUT = Array{Float64}(zeros(tot_voxels,4,nt));

    begin
        @time @inbounds for dynamics in 1:nt
            tmpOUT[:,:,dynamics] = nlsworker_sage(tot_voxels,vec_data[:,:,dynamics],vec_mask,vec_R2s_log,vec_R2_log,x0,echos)
            println("Done with Fit $dynamics of $nt")
        end    
    end
    
    vec_R2sfit = tmpOUT[:,2,:]
    vec_R2fit = tmpOUT[:,3,:]

    # vec_R2sfit = tmpOUT[:,2]
    # vec_R2fit = tmpOUT[:,3]


    #=
    @time for JJ=1:nt;
        temp = vec_data[:,:,JJ]
        tempFit_data = fitdata[:,:,JJ]
        @time fitdata[:,:,JJ] = work(IND,SAGE_biexp3p_d,echos,temp,X0,tempFit_data)
        
        println("Done with Fit $JJ of $nt")
    end
    =#
    #-------------------------------------  
    # Xv = Array{Float64}(fitdata);
    # Xv = reshape(Xv,nx,ny,nz,4,nt);
    # R2s = zeros(nx,ny,nz,nt);
    # R2 = zeros(nx,ny,nz,nt);
    # R2s = Xv[:,:,:,2,:];
    # R2 = Xv[:,:,:,3,:];
    R2s = reshape(vec_R2sfit,nx,ny,nz,nt);
    R2 = reshape(vec_R2fit,nx,ny,nz,nt);
    R2s[findall(x->x.>1000,R2s)].=0;
    R2s[findall(x->x.<0,R2s)].=0;
    R2[findall(x->x.>1000,R2)].=0;
    R2[findall(x->x.<0,R2)].=0;
    #-------------------------------------

    d=niread(paths[1]);

    tmp1 = voxel_size(d.header)[1]
    tmp2 = voxel_size(d.header)[2]
    tmp3 = voxel_size(d.header)[3]

    R2s_fname = @sprintf("%s/SAGE_R2s_Brain_julia.nii.gz",base)
    R2_fname = @sprintf("%s/SAGE_R2_Brain_julia.nii.gz",base)

   
    temp1 = NIVolume(R2s;voxel_size=(tmp1,tmp2,tmp3))
    niwrite(R2s_fname,temp1)
    temp2 = NIVolume(R2;voxel_size=(tmp1,tmp2,tmp3))
    niwrite(R2_fname,temp2)

    R2s_fname_log = @sprintf("%s/SAGE_R2s_Brain_julia_log.nii.gz",base)
    R2_fname_log = @sprintf("%s/SAGE_R2_Brain_julia_log.nii.gz",base)
    temp1 = NIVolume(R₂star_log;voxel_size=(tmp1,tmp2,tmp3));
    niwrite(R2s_fname_log,temp1)
    temp2 = NIVolume(R₂_log;voxel_size=(tmp1,tmp2,tmp3));
    niwrite(R2_fname_log,temp2)

    R2s_fname_log = @sprintf("%s/SAGE_R2s_Brain_residual_fit_m_log.nii.gz",base)
    R2_fname_log = @sprintf("%s/SAGE_R2_Brain_residual_fit_m_log.nii.gz",base)
    temp1 = NIVolume(R2s .- R₂star_log;voxel_size=(tmp1,tmp2,tmp3));
    niwrite(R2s_fname_log,temp1)
    temp2 = NIVolume(R2 .- R₂_log;voxel_size=(tmp1,tmp2,tmp3));
    niwrite(R2_fname_log,temp2)


    println("The R2s is saved at $R2s_fname")
    println("The R2 is saved at $R2_fname")

    for ii in 1:ne
        temp = NIVolume(DATA[:,:,:,ii,:];voxel_size=(tmp1,tmp2,tmp3))
        fname = @sprintf("%s/Input_%s.nii.gz",base,ii)
        niwrite(fname, temp)
    end


end

main(A)