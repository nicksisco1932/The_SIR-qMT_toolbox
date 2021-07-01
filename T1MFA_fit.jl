#=
    Author: Nicholas J. Sisco, Ph.D. 
    Email: X@barrowneuro.org where X = nicholas.sisco
    Affiliation: Barrow Neurological Institue in Phoenix, AZ

    Title: Multi-flip angle T1 fitting using least squares and Levenberg-Marquardt 
    - This is a basic script written in Julia for fitting signal acquired multi-flip angle MRI. It is based on the Ernst angle.

    Requirements:
    - Julia 1.5

    Input and Output Parameters for Code

    utils.jl needs to be in the same directory as this script

    User supplied data:
        - julia needs to be in your $PATH for this command structure to work
        - Command line usage
        $ julia -t <threads #> T1MFA_fit.jl <path to Python processed data>

        - In the Python processed directory, there should be nt number of preprocessed nifti files, corresponding 
        to the number of dynamic time points in the data, e.g. nt = 4 means there are ti = [15,15,...,...], etc.

    Output:
        - The output from this script are three nifti files for the fitting parameters pool size ratio, R1f (free transverse relaxation rate), and Sf (field inhomogeneity)


    Future Updates TO Do:
        1) add in user defined B1 files
        2) other things

    References:
    

    Change log:
    "history (of nifti library changes):\n"
    "\n",
    "0.0  June 18, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - initial version \n"
    "\n",
    "0.1  June 22, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - Changed I/O \n"

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

function commandline()
        

    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "MFA_nii"
        required = true
        "nii_brainMask"
        required = false
    end

    println(parse_args(settings))

    # parsed_args = TE_parse_commandline();
    for (out, val) in parse_args(settings)
        println(" $out => $val")
    end
    return parse_args(settings)
end

a = commandline()



function main()
    
    base = dirname(a["MFA_nii"])
    try
        b_fname = joinpath(base,basename(a["nii_brainMask"]))
        brain_data = niread(b_fname);
        MASK=Array{Bool}(brain_data.raw);
        NO_MASK = false
    catch
        println("No mask, fitting all voxels")
        NO_MASK = true
    end

    b_fname = joinpath(base,basename(a["nii_brainMask"]))
    
    DATA,MASK,nx,ny,nz,nt=load_data(paths,b_fname);

    data1 = niread(paths[1]);
    
    nx,ny,nz=size(data1);
    nt = length(paths);
    DATA = Array{Float64}(zeros(nx,ny,nz,nt));
    try
        for (n,ii) in enumerate(paths)
            DATA[:,:,:,n] = loader(ii);
        end    
    catch
        DATA = loader(paths[1]);
        nx,ny,nz,nt = size(DATA)
    end


    tot_voxels = nx*ny*nz
    if NO_MASK
        vec_mask = ones(tot_voxels)
    else
        vec_mask = dropdims(reshape(MASK[:,:,:,1,1],tot_voxels,1);dims=2) #change this eventually
    end

    vec_data = reshape(DATA,tot_voxels,nt)

    thr = Threads.nthreads()
    println("Fitting with $thr threads")
    
    TR = 7.59/1000;
    alpha_default = 1.0;
    FA = Array{Float64}([20,18,16 ,14,12,10,8,6,4,2.0]);
    LSQNLIN = false
    if LSQNLIN
        fitdata = nlsFit_MFA(vec_mask,TR,FA,vec_data,alpha_default);
    else

        fitdata = alt_fit(vec_mask,TR,FA,vec_data,alpha_default)
    end
    # m0 = fitdata[:,1];
    # t1=fitdata[:,2];

    Xv = Array{Float64}(fitdata);
    Xv = reshape(Xv,nx,ny,nz,2);
    M0 = zeros(nx,ny,nz);
    T1 = zeros(nx,ny,nz);
    M0 = Xv[:,:,:,1];
    T1 = Xv[:,:,:,2];
    m=1E9;
    T1[findall(x->x.>10,T1)].=0
    T1[findall(x->x.<0,T1)].=0
    M0[findall(x->x.>m,M0)].=0
    M0[findall(x->x.<0,M0)].=0

    # Sf[findall(x->x.>10,Sf)].=0
    # Sf[findall(x->x.<-10,Sf)].=0

    d=niread(b_fname);
    d=niread(paths[1]);

    tmp1 = voxel_size(d.header)[1]
    tmp2 = voxel_size(d.header)[2]
    tmp3 = voxel_size(d.header)[3]

    T1_fname = @sprintf("%sMFA_T1_Brain_julia.nii.gz",base)
    M0_fname = @sprintf("%sMFA_M0_Brain_julia.nii.gz",base)
    niwrite(T1_fname,NIVolume(T1;voxel_size=(tmp1,tmp2,tmp3)))
    niwrite(M0_fname,NIVolume(M0;voxel_size=(tmp1,tmp2,tmp3)))

    println("The T1 is saved at $T1_fname")
    println("The M0 is saved at $M0_fname")
end

main()