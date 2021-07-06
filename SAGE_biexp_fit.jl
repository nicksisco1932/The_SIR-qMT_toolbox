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
using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Optim

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
    
    # MASK = Array{Bool}(niread(b_fname))
    # MASK = MASK.raw
    # vec_mask = reshape(MASK,tot_voxels,1)
    

    vec_data = reshape(DATA,tot_voxels,ne,nt)



    # te1=0.00782;
    # te2=0.028769;
    # te3=0.060674;
    # te4=0.081622;
    # te5=0.102571;
    # echos=Vector{Float64}([te1,te2,te3,te4,te5])
    
    te1=parse(Float64,a["TE1"])
    te2=parse(Float64,a["TE2"])
    te3=parse(Float64,a["TE3"])
    te4=parse(Float64,a["TE4"])
    te5=parse(Float64,a["TE5"])
    echos=Vector{Float64}([te1,te2,te3,te4,te5])

    # echos=Vector{Float64}()
    TR = 1800/1000;

    X0=[1000,100.0,50,100]
    IND=findall(x->x.>0,vec_mask[:,1]);


    fitdata = zeros(tot_voxels,4,nt);
#= Testing
    # Threads.@threads for ii in IND
        
    #     fitY = optim_fitty(SAGE_biexp3p_d, echos, vec_data[ii,:,1], X0)
    # end
    # tempFit_data = fitdata[:,:,1]
    # @time work(IND,SAGE_biexp3p_d,echos,vec_data[:,:,1],X0,tempFit_data,1)
=#
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

    println("The R2s is saved at $R2s_fname")
    println("The R2 is saved at $R2_fname")

    for ii in 1:ne
        temp = NIVolume(DATA[:,:,:,ii,:];voxel_size=(tmp1,tmp2,tmp3))
        fname = @sprintf("%s/Input_%s.nii.gz",base,ii)
        niwrite(fname, temp)
    end


end

main()