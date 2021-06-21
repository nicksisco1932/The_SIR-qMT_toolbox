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
        $ julia -t <threads #> ./SIR_fit.jl <path to Python processed data>

        - In the Python processed directory, there should be nt number of preprocessed nifti files, corresponding 
        to the number of dynamic time points in the data, e.g. nt = 4 means there are ti = [15,15,...,...], etc.

    Output:
        - The output from this script are three nifti files for the fitting parameters pool size ratio, R1f (free transverse relaxation rate), and Sf (field inhomogeneity)


    Future Updates TO Do:
        1) user defined ti and td
        2) user defined kmf
        3) user defined preprocessing niftis, i.e. flexible input names
            - I think we need== ("-i" action => :append_arg)
            - then called on command line "-i" <value> "-i" <value2>, etc. 
        4) main function within utils
        5) module creation
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

=#
using NIfTI; 
using LsqFit;
using Printf
using ArgParse;

include("./utils.jl")
function commandline()
        

    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "SIR_nii_1"
        required = true
        "SIR_nii_2"
        required = true
        "SIR_nii_3"
        required = true
        "SIR_nii_4"
        required = true
        "SIR_nii_brainMask"
        required = true
    end

    println(parse_args(settings))

    # parsed_args = SIR_parse_commandline();
    for (out, val) in parse_args(settings)
        println(" $out => $val")
    end
    return parse_args(settings)
end

a = commandline()

function main()
    

    base = dirname(a["SIR_nii_1"])
    paths = [joinpath(base,basename(a["SIR_nii_1"])),
        joinpath(base,basename(a["SIR_nii_2"])),
        joinpath(base,basename(a["SIR_nii_3"])),
        joinpath(base,basename(a["SIR_nii_4"]))
        ]
    b_fname = joinpath(base,basename(a["SIR_nii_brainMask"]))
    ti_times,td_times,X0 = params();
    DATA,MASK,nx,ny,nz,nt=load_data(paths,b_fname);

    tot_voxels,vec_mask,pre_zeros,times,Yy = reshape_and_normalize(DATA,MASK,ti_times,td_times,nx,ny,nz,nt);

    # It is faster for Julia as it is column major, but I suck at writing in column major. This means I transpose it and retranspose. :(
    thr = Threads.nthreads()
    println("Fitting with $thr threads")
    ROW_MAJOR = true
    if ROW_MAJOR
        fitOUT = timed(tot_voxels,vec_mask,pre_zeros,times,Yy,X0); # row major
    else
        colY = Array{Float64}(Yy');
        col_pre_zeros = Array{Float64}(pre_zeros');
        fitOUT = timed(tot_voxels,vec_mask,col_pre_zeros,times,colY,X0); # column major
    end

    # param, resids = nlsfit(SIR_signal_absolute,yy,vec_mask,times,X0)

    if ROW_MAJOR
        Xv = Array{Float64}(fitOUT);
    else
        Xv = Array{Float64}(fitOUT');
    end
    Xv = reshape(Xv,nx,ny,nz,4);
    PSR = zeros(nx,ny,nz);
    R1f = zeros(nx,ny,nz);
    Sf = zeros(nx,ny,nz);
    PSR = Xv[:,:,:,1]*100;
    R1f = Xv[:,:,:,2];
    Sf = Xv[:,:,:,3];


    PSR[findall(x->x.>100,PSR)].=0
    PSR[findall(x->x.<0,PSR)].=0
    R1f[findall(x->x.>50,R1f)].=0
    R1f[findall(x->x.<0,R1f)].=0

    Sf[findall(x->x.>10,Sf)].=0
    Sf[findall(x->x.<-10,Sf)].=0

    d=niread(b_fname);
    d=niread(paths[1]);

    tmp1 = voxel_size(d.header)[1]
    tmp2 = voxel_size(d.header)[2]
    tmp3 = voxel_size(d.header)[3]

    PSR_fname = @sprintf("%sSIR_Brain_PSR_julia.nii.gz",base)
    R1f_fname = @sprintf("%sSIR_Brain_R1f_julia.nii.gz",base)
    SF_fname = @sprintf("%sSIR_Brain_Sf_julia.nii.gz",base)
    niwrite(PSR_fname,NIVolume(PSR;voxel_size=(tmp1,tmp2,tmp3)))
    niwrite(R1f_fname,NIVolume(R1f;voxel_size=(tmp1,tmp2,tmp3)))
    niwrite(SF_fname,NIVolume(Sf;voxel_size=(tmp1,tmp2,tmp3)))

    println("The PSR is saved at $PSR_fname")
    println("The R1f is saved at $R1f_fname")
    println("The Sf is saved at $SF_fname")
end

main()