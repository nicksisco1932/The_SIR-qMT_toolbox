#=
    Author: Nicholas J. Sisco, Ph.D. 
    Email: X@barrowneuro.org where X = nicholas.sisco
    Affiliation: Barrow Neurological Institue in Phoenix, AZ

    Title: Selective inversion recovery fitting using least squares and Levenberg-Marquardt 
    - This is a basic script written in Julia for fitting signal acquired from SIR-qMT to a double exponential 
    equation described in (1,2,3,4). 


    Requirements:
    - Julia 1.5

    Input and Output Parameters for Code

    utils.jl needs to be in the same directory as this script

    User supplied data:
        - Command line usage
        $ julia ./SIR_fit.jl --TI 15 15 278 1007 --TD 684 4121 2730 10 --SIR_Data <PATH>/SIR_DATA.nii.gz --SIR_brainMask <PATH>/brain_mask.nii.gz --kmf 14.5 --Sm 0.83

    Output:
        - The output from this script are three nifti files for the fitting parameters pool size ratio, R1f (free transverse relaxation rate), and Sf (field inhomogeneity)


    Future Updates TO Do:
        6) unit testing
        7) depolyable docker

    References:
    1. R. D. Dortch, J. Moore, K. Li, M. Jankiewicz, D. F. Gochberg, J. A. Hirtle, J. C. Gore, S. A. Smith, Quantitative magnetization transfer imaging of human brain at 7T. Neuroimage. 64, 640–649 (2013).
    2. F. Bagnato, G. Franco, F. Ye, R. Fan, P. Commiskey, S. A. Smith, J. Xu, R. Dortch, Selective inversion recovery quantitative magnetization transfer imaging: Toward a 3 T clinical application in multiple sclerosis. Mult. Scler. J. 26, 457–467 (2020).
    3. R. D. Dortch, F. Bagnato, D. F. Gochberg, J. C. Gore, S. A. Smith, Optimization of selective inversion recovery magnetization transfer imaging for macromolecular content mapping in the human brain. Magn. Reson. Med. 80, 1824–1835 (2018).
    4. R. D. Dortch, K. Li, D. F. Gochberg, E. B. Welch, A. N. Dula, A. A. Tamhane, J. C. Gore, S. A. Smith, Quantitative magnetization transfer imaging in human brain at 3 T via selective inversion recovery. Magn. Reson. Med. 66, 1346–1352 (2011).

    Change log:
    "history (of nifti library changes):\n"
    "\n",
    "0.0  June 18, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - initial version \n"
    "\n",
    "0.0  June 21, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - Updated I/O \n"
    "\n",
    "1.0  June 14, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - Major updated I/O 
        - Change fitting function and for loop 
        - Dramatic improvement in speed 
        \n"

=#
using NIfTI; 
using LsqFit;
using Printf
using ArgParse;

#=----------------------------------------------
    Commandline Arguments
----------------------------------------------
=#
include("./utils.jl")
function commandline()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "--SIR_Data"
            required = true
        "--SIR_brainMask"
            required = true
        "--TD"
            nargs = 4
            arg_type = Float64
            help = "tD values, 4 required"
        "--TI"
            nargs = 4
            arg_type = Float64
            help = "tI values, 4 required"
        "--kmf"
            nargs = 1
            arg_type = Float64
            help = "kmf value, for phantoms this is 35.0 and brains 14.5"
        "--Sm"
            nargs = 1
            arg_type = Float64
            help = "Sm value, unless >3T this is 0.83"
    end

    println(parse_args(settings))

    # parsed_args = SIR_parse_commandline();
    for (out, val) in parse_args(settings)
        println(" $out => $val")
    end
    return parse_args(settings)
end

#=----------------------------------------------
    Commandline Arguments END
----------------------------------------------
=#

#=
test data /mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618/SIR_mo_corr.nii.gz
test mask /mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618/brain_mask.nii.gz
temp = "/mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618/SIR_mo_corr.nii.gz"
base = dirname(temp)
paths = [joinpath(base,"SIR_mo_corr.nii.gz")]
b_fname = joinpath(base,"brain_mask.nii.gz")
=#


function main(a)
    
    base = dirname(a["SIR_Data"])
    paths = a["SIR_Data"]
    b_fname = joinpath(base,basename(a["SIR_brainMask"]))
 
    #=----------------------------------------------
    Variable Definitions
    ----------------------------------------------
    =#
    # If Error in TD unit handling. If the user puts in ms, convert to s
    if a["TI"][1]>1
        ti_times = a["TI"]* 1E-3
        td_times = a["TD"]* 1E-3
    else
        ti_times = a["TI"]
        td_times = a["TD"]
    end
    kmfmat = a["kmf"][1]
    Smmat=a["Sm"][1]
    X0 = Array{Float64}([0.07, 1, -1, 1.5]);

    #=----------------------------------------------
    Variable Definitions END
    ----------------------------------------------
    =#

    thr = Threads.nthreads()
    println("Fitting with $thr threads")

    
    #=----------------------------------------------
    Setting up data to be fit
    ----------------------------------------------
    =#
    nx,ny,nz,nt = size(niread(paths))
    temp=niread(paths)
    temp2=Array{Float64}(temp.raw)
    tempdata=Array{Float64}(zeros(nt,nx,ny,nz))
    for ii in 1:nt
        tempdata[ii,:,:,:] = temp2[:,:,:,ii];
    end
    DATA=tempdata;
    MASKtmp = niread(b_fname);
    MASK=Array{Bool}(MASKtmp.raw);

    tot_voxels,times,Yy = reshape_and_normalize( DATA,ti_times,td_times,nx,ny,nz,nt);

    tmpOUT = Array{Float64}(zeros(4,tot_voxels));
    vec_mask = reshape(MASK,(tot_voxels));


    #=----------------------------------------------
    Setting up data to be fit END
    ----------------------------------------------
    =#


    #=----------------------------------------------
    Fit END
    ----------------------------------------------
    =#
    
    mag=true;
    model(x,p) = SIR_Mz0(x,p,kmfmat,Sm=Smmat,R1m=NaN,mag=true) 

    function f() # anonymous function for fitting
        begin
            Threads.@threads for ii in 1:tot_voxels
                if vec_mask[ii]
                    tmpOUT[:,ii] = nlsfit(model,times,Yy[:,ii],X0)
                end #if
            end #for
        end #begin
    end #f()

    @time f() # Fitting call

    #=----------------------------------------------
    Fit END
    ----------------------------------------------
    =#

    #= 
    Finishing up 
    =#
    Xv = tmpOUT;
    Xv = reshape(Xv,4,nx,ny,nz);
    PSR = zeros(nx,ny,nz);
    R1f = zeros(nx,ny,nz);
    Sf = zeros(nx,ny,nz);
    PSR = Xv[1,:,:,:]*100;
    R1f = Xv[2,:,:,:];
    Sf = Xv[3,:,:,:];


    PSR[findall(x->x.>100,PSR)].=0
    PSR[findall(x->x.<0,PSR)].=0
    R1f[findall(x->x.>50,R1f)].=0
    R1f[findall(x->x.<0,R1f)].=0

    Sf[findall(x->x.>10,Sf)].=0
    Sf[findall(x->x.<-10,Sf)].=0

    d=niread(paths);

    tmp1 = voxel_size(d.header)[1]
    tmp2 = voxel_size(d.header)[2]
    tmp3 = voxel_size(d.header)[3]

    PSR_fname = @sprintf("%s/SIR_Brain_PSR_julia.nii.gz",base)
    R1f_fname = @sprintf("%s/SIR_Brain_R1f_julia.nii.gz",base)
    SF_fname = @sprintf("%s/SIR_Brain_Sf_julia.nii.gz",base)

    fname_1 = @sprintf("%s/SIR_DATA.nii.gz",base)
    niwrite(PSR_fname,NIVolume(PSR;voxel_size=(tmp1,tmp2,tmp3)))
    niwrite(R1f_fname,NIVolume(R1f;voxel_size=(tmp1,tmp2,tmp3)))
    niwrite(SF_fname,NIVolume(Sf;voxel_size=(tmp1,tmp2,tmp3)))

    niwrite(fname_1,NIVolume(temp2;voxel_size=(tmp1,tmp2,tmp3)))

    println("The PSR is saved at $PSR_fname")
    println("The R1f is saved at $R1f_fname")
    println("The Sf is saved at $SF_fname")
    println("The original data is saved at $fname_1")
end

A = commandline()
@time main(A)
