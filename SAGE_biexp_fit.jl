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

thr = Threads.nthreads()
if thr==16
    Threads.nthreads() = 8
    return thr = 8
end
function commandline()
        

    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "TE_nii_1"
        required = true
        "TE_nii_2"
        required = true
        "TE_nii_3"
        required = true
        "TE_nii_4"
        required = true
        "TE_nii_5"
        required = true
        "SAGE_nii_brainMask"
        required = true
        "TE1"
        required = true
        "TE2"
        required = true
        "TE3"
        required = true
        "TE4"
        required = true
        "TE5"
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

a = commandline()

function optim_fitty(f::Function, xdata::Vector{Float64},data::Vector{Float64},x0::Vector{Float64})
    # fitdata = Optim.minimizer( optimize(b -> sqerrorSAGE(b, xdata,  data), x0))
    fitdata = curve_fit(f,xdata,data,x0) # this is faster and has good results
    return fitdata
end

function work(IND::Vector{Int64},f::Function,echos::Vector{Float64},vec_data::Array{Float64},X0::Vector{Float64},fitdata::Array{Float64})
    # Threads.@threads for ii in IND # surprisingly, threads is not faster
    for ii in IND
        temp = vec_data[ii,:];
        fitY = optim_fitty(f, echos, temp, X0);fitdata[ii,:] = fitY.param
        # fitY = optim_fitty(echos, temp, X0);fitdata[ii,:] = fitY
    end
    
    return fitdata
end

function main()
    


    base = dirname(a["TE_nii_1"])
    paths = [
        joinpath(base,basename(a["TE_nii_1"])),
        joinpath(base,basename(a["TE_nii_2"])),
        joinpath(base,basename(a["TE_nii_3"])),
        joinpath(base,basename(a["TE_nii_4"])),
        joinpath(base,basename(a["TE_nii_5"]))
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
    echos = Array{Float64}(parse(Float64,a["echos"]))
    if echos[1]>1
        echos = echos.* 1E-3
    else
        echos = echos;
    end

    TR = 1800/1000;

    X0=[1000,100.0,50,100]
    IND=findall(x->x.>0,vec_mask[:,1]);

    STE1=zeros(nx,ny,nz)
    STE2=similar(STE1)
    STE5=similar(STE1)
    begin
        for ii in 1:nx
            for jj in 1:ny
                for kk in 1:nz
                    STE1[ii,jj,kk] = mean(DATA[ii,jj,kk,1,1:15]);
                    STE2[ii,jj,kk] = mean(DATA[ii,jj,kk,2,1:15]);
                    STE5[ii,jj,kk] = mean(DATA[ii,jj,kk,5,1:15]);
                end
            end
        end
    end

    STE1_pre = repeat(STE1,1,1,1,150);
    STE2_pre = repeat(STE1,1,1,1,150);
    STE5_pre = repeat(STE1,1,1,1,150);

    ste1 = DATA[:,:,:,1,:];
    ste2 = DATA[:,:,:,2,:];
    ste5 = DATA[:,:,:,5,:];


    
    R₂star_log = similar(STE1_pre);
    R₂_log = similar(STE1_pre);

    IND=findall(x->x==true,MASK);
    begin
        for ii in 1:nx
            for jj in 1:ny
                for kk in 1:nz
                    if MASK[ii,jj,kk]
                        for bb = 1:nt
                            a = ste1[ii,jj,kk,bb];
                            # b = STE1_pre[ii,jj,kk,bb];
                            c = ste2[ii,jj,kk,bb];
                            # d = STE2_pre[ii,jj,kk,bb];
                            e = ste5[ii,jj,kk,bb];
                            # f = STE5_pre[ii,jj,kk,bb];
                            ste0 = a*(a/c)^(te1/(te2-te1))
                            R₂star_log[ii,jj,kk,bb] = log( (a) / (c)) /(te2-te1)
                            R₂_log[ii,jj,kk,bb] = log( (ste0 / e)) /te5
                        end
                    end
                end
            end
        end
    end


    fitdata = zeros(tot_voxels,4,nt);
#= Testing
    # Threads.@threads for ii in IND
        
    #     fitY = optim_fitty(SAGE_biexp3p_d, echos, vec_data[ii,:,1], X0)
    # end
    # tempFit_data = fitdata[:,:,1]
    # @time work(IND,SAGE_biexp3p_d,echos,vec_data[:,:,1],X0,tempFit_data,1)
=#

    R2fit = zeros(nx,ny,nz,nt);
    R2sfit = zeros(nx,ny,nz,nt);
    
    newDATA = zeros(nx,ny,nz,nt,ne);
    
    function sage_ns(echotimes,p)
        x=echotimes;
        TE=x[end];
        R₂star = p[2]
        R₂ = p[3]
        S₀I = p[1]
        S₀II = p[4]

        M = similar(echotimes)
        for k in 1:length(echotimes)
            if x[k]<TE/2
                M[k] = S₀I * exp(-x[k] * R₂star)
            elseif x[k]>TE/2
                M[k] = S₀II * exp(-TE * (R₂star - R₂) * exp(-x[k] * (2*R₂-R₂star)))
            end
        end
        return M
    end

    function sqerrorSAGE(betas::Vector{Float64}, X::Vector{Float64}, Y::Vector{Float64}) 
        err = 0.0
        pred_i = sage_ns(X,betas)
        for ii in 1:length(Y)
            err += (abs.(pred_i[ii])-Y[ii]).^2
        end
        return err
    end
    
    #= Too slow
    x0 = zeros(4,1)
    @time for ii in 1:nx
        for jj in 1:ny
            for kk in 1:nz
                if MASK[ii,jj,kk]
                    Threads.@threads for bb = 1:nt
                        x0[1] = maximum(DATA[ii,jj,kk,:,bb])
                        x0[2] = R₂star_log[ii,jj,kk,bb]
                        x0[3] = R₂_log[ii,jj,kk,bb]
                        x0[4] = 1
                        model(echos,x0) = sage_ns(echos,x0)
                        # fit = nlsfit(model,echos,DATA[ii,jj,kk,:,bb],x0)
                        fit = Optim.minimizer(optimize(b -> sqerrorSAGE(b, echos,  DATA[IND[1],:,1]), X0))
                        R2sfit[ii,jj,kk,bb] = fit[2]
                        R2fit[ii,jj,kk,bb] = fit[3]
                    end
                end
            end
        end
    end
=#

    vec_R2fit = zeros(tot_voxels,nt);
    vec_R2sfit = zeros(tot_voxels,nt);
    vec_R2_log = reshape(R₂_log,tot_voxels,nt);
    vec_R2s_log = reshape(R₂star_log,tot_voxels,nt);

    x0 = zeros(4,1)
    for dynamics in 1:nt
        @time for vox in 1:tot_voxels
            if vec_mask[vox]
                x0[1] = maximum(vec_data[vox,:,dynamics])
                x0[2] = vec_R2s_log[vox,dynamics]
                x0[3] = vec_R2_log[vox,dynamics]
                x0[4] = 1
                println(x0)
                model(xs,p) = sage_ns(echos,x0)
                fit = nlsfit(model,echos,vec_data[vox,:,dynamics],x0)
                # fit = Optim.minimizer(optimize(b -> sqerrorSAGE(b, echos,  DATA[IND[1],:,1]), x0))
                # vec_R2sfit[vox,:,dynamics] = fit[2];
                # vec_R2fit[vox,:,dynamics] = fit[3];
            end
        end
        println("Done with Fit $dynamics of $nt")
    end

    @time for JJ=1:nt;
        temp = vec_data[:,:,JJ]
        tempFit_data = fitdata[:,:,JJ]
        @time fitdata[:,:,JJ] = work(IND,SAGE_biexp3p_d,echos,temp,X0,tempFit_data)
        
        println("Done with Fit $JJ of $nt")
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

    println("The R2s is saved at $R2s_fname")
    println("The R2 is saved at $R2_fname")

    for ii in 1:ne
        temp = NIVolume(DATA[:,:,:,ii,:];voxel_size=(tmp1,tmp2,tmp3))
        fname = @sprintf("%s/Input_%s.nii.gz",base,ii)
        niwrite(fname, temp)
    end


end

main()