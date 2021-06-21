#=
    Author: Nicholas J. Sisco, Ph.D. 
    Email: X@barrowneuro.org where X = nicholas.sisco
    Affiliation: Barrow Neurological Institue in Phoenix, AZ

    Title: Multi-flip angle T1 fitting using least squares and Levenberg-Marquardt 
    - This is a basic script written in Julia for fitting signal acquired from SIR-qMT to a double exponential 
    equation described in (1,2,3,4). 


    Requirements:
    - Julia 1.5

    Input and Output Parameters for Code

    utils.jl needs to be in the same directory as this script

    User supplied data:
        - Command line usage
        $ julia -t <threads #> SIR_PING_Brain_20210617_v2.jl <path to Python processed data>

        - In the Python processed directory, there should be nt number of preprocessed nifti files, corresponding 
        to the number of dynamic time points in the data, e.g. nt = 4 means there are ti = [15,15,...,...], etc.

    Output:
        - The output from this script are three nifti files for the fitting parameters pool size ratio, R1f (free transverse relaxation rate), and Sf (field inhomogeneity)


    Future Updates TO Do:
        1) user defined ti and td
        2) user defined kmf
        3) user defined preprocessing niftis, i.e. flexible input names
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

function main()
    

    paths = [ @sprintf("%sMFA_10_angles.nii.gz",base)]
    b_fname = @sprintf("%sMFA_brain_mask.nii.gz",base)
    
    DATA,MASK,nx,ny,nz,nt=load_data(paths,b_fname);
    tot_voxels = nx*ny*nz
    vec_mask = dropdims(reshape(MASK[:,:,:,1,1],tot_voxels,1);dims=2) #change this eventually

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