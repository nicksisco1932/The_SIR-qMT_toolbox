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
    fitdata = curve_fit(f,xdata,data,x0;autodiff=:finiteforward) # this is faster and has good results
    return fitdata
end

function work(IND::Vector{Int64},f::Function,echos::Vector{Float64},vec_data::Array{Float64},X0::Vector{Float64},fitdata::Array{Float64})
    # Threads.@threads for ii in IND # surprisingly, threads is not faster
    for ii in IND
        temp = vec_data[ii,:];
        fitY = optim_fitty(f, echos, temp, X0)
        # fitY = optim_fitty(echos, temp, X0)
        fitdata[ii,:] = fitY.param
    end
    
    return fitdata
end

function main()
    


    base = dirname(a["TE_nii_1"])
    paths = [joinpath(base,basename(a["TE_nii_1"])),
        joinpath(base,basename(a["TE_nii_2"])),
        joinpath(base,basename(a["TE_nii_3"])),
        joinpath(base,basename(a["TE_nii_4"])),
        joinpath(base,basename(a["TE_nii_5"]))
        ]
    b_fname = joinpath(base,basename(a["SAGE_nii_brainMask"]))

    DATA,MASK,nx,ny,nz,ne,nt=load_sage_data(paths,b_fname);
    tot_voxels = nx*ny*nz
    vec_mask = reshape(MASK[:,:,:,1,1],tot_voxels,1);

    vec_data = reshape(DATA,tot_voxels,ne,nt)

    te1=0.007653;
    te2=0.027969;
    te3=0.059074;
    te4=0.07939;
    te5=0.099706;
    echos=[te1,te2,te3,te4,te5]
    TR = 1800/1000;

    X0=[1,100.0,50,1]
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

    R2s_fname = @sprintf("%sSAGE_R2s_Brain_julia.nii.gz",base)
    R2_fname = @sprintf("%sSAGE_R2_Brain_julia.nii.gz",base)
    niwrite(R2s_fname,NIVolume(R2s;voxel_size=(tmp1,tmp2,tmp3)))
    niwrite(R2_fname,NIVolume(R2;voxel_size=(tmp1,tmp2,tmp3)))

    println("The R2s is saved at $R2s_fname")
    println("The R2 is saved at $R2_fname")


end

main()