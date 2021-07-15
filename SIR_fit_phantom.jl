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

base="/mnt/c/Users/nicks/Documents/MRI_data/Phantoms/20210713_Phantom_Agarose/DICOM/SIR/"
paths = [@sprintf("%sDICOM_WIP_2D_SIR_TSE_Log_16_TI_2500_TD_20210713135226_901.nii.gz",base)]

# base = dirname(a["SIR_nii_1"])
# paths = [joinpath(base,basename(a["SIR_nii_1"])),
#     joinpath(base,basename(a["SIR_nii_2"])),
#     joinpath(base,basename(a["SIR_nii_3"])),
#     joinpath(base,basename(a["SIR_nii_4"]))
#     ]
# b_fname = joinpath(base,basename(a["SIR_nii_brainMask"]))
ti_times = Array{Float64}([10,15,22,33,48,72,107,158,235,348,516,766,1136,1685,2500,10000] * 1E-3);
td_times = Array{Float64}(repeat([2500], 16, 1)*1E-3);
X0 = Array{Float64}([0.07, 1.5, -1, 1.5]);

function load_data_phantom(paths)
    data1 = niread(paths[1]);
    
    nx,ny,nz,nt = size(data1);
    if nz ==1
        temp=Array{Float64}(data1.raw)
        DATA = dropdims(temp;dims=3);
    else
        temp=Array{Float64}(data1.raw)
        DATA=temp;
    end
    return DATA,nx,ny,nz,nt    
end

DATA,nx,ny,nz,nt=load_data_phantom(paths);

function reshape_and_normalize_phantom(data_4d::Array{Float64},TI::Array{Float64},TD::Array{Float64},NX::Int64,NY::Int64,NZ::Int64,NT)
    # [pmf R1f Sf M0f]

    tot_voxels = NX*NY*NZ
    X=hcat(TI,TD)
    tmp = reshape(data_4d,(NT,tot_voxels));
    Yy=Array{Float64}(zeros(size(tmp)))
    
    for ii in 1:tot_voxels
        
        Yy[:,ii] = tmp[:,ii] / (tmp[end, ii])
        
    end

    Yy[findall(x->isinf(x),Yy)].=0
    Yy[findall(x->isnan(x),Yy)].=0

    return tot_voxels,X,Yy
end

tot_voxels,times,Yy = reshape_and_normalize_phantom( DATA,ti_times,td_times,nx,ny,nz,nt);

thr = Threads.nthreads()
println("Fitting with $thr threads")

tmpOUT = Array{Float64}(zeros(4,tot_voxels));
# IND = findall(x->x>0,Yy[:,1])

function nlsfit(f::Function, xvalues::Array{Float64},yvalues::Array{Float64},guesses::Array{Float64})
    fit = curve_fit(f,xvalues,yvalues,guesses)
    return fit.param
end

kmfmat=35.0;
Smmat=0.83;
mag=true;
model(x,p) = SIR_Mz0(x,p,kmfmat,Sm=Smmat,R1m=NaN,mag=true) 

function f()
    begin
        Threads.@threads for ii in 1:tot_voxels
            tmpOUT[:,ii] = nlsfit(model,times,Yy[:,ii],X0)
        end 
    end
end

@time f()


if nz == 1
    Xv = tmpOUT;
    Xv = reshape(Xv,nx,ny,4);
    PSR = zeros(nx,ny);
    R1f = zeros(nx,ny);
    Sf = zeros(nx,ny);
    PSR = Xv[:,:,1]*100;
    R1f = Xv[:,:,2];
    Sf = Xv[:,:,3];
else
    Xv = tmpOUT;
    Xv = reshape(Xv,nx,ny,nz,4);
    PSR = zeros(nx,ny,nz);
    R1f = zeros(nx,ny,nz);
    Sf = zeros(nx,ny,nz);
    PSR = Xv[:,:,:,1]*100;
    R1f = Xv[:,:,:,2];
    Sf = Xv[:,:,:,3];
end


PSR[findall(x->x.>100,PSR)].=0;
PSR[findall(x->x.<0,PSR)].=0;
R1f[findall(x->x.>50,R1f)].=0;
R1f[findall(x->x.<0,R1f)].=0;

Sf[findall(x->x.>10,Sf)].=0;
Sf[findall(x->x.<-10,Sf)].=0;

d=niread(paths[1]);

tmp1 = voxel_size(d.header)[1];
tmp2 = voxel_size(d.header)[2];
tmp3 = voxel_size(d.header)[3];

PSR_fname = @sprintf("%s/SIR_Brain_PSR_julia.nii.gz",base)
R1f_fname = @sprintf("%s/SIR_Brain_R1f_julia.nii.gz",base)
SF_fname = @sprintf("%s/SIR_Brain_Sf_julia.nii.gz",base)

# fname_1 = @sprintf("%s/SIR_T1_1.nii.gz",base)
# fname_2 = @sprintf("%s/SIR_T1_2.nii.gz",base)
# fname_3 = @sprintf("%s/SIR_T1_3.nii.gz",base)
# fname_4 = @sprintf("%s/SIR_T1_4.nii.gz",base)
niwrite(PSR_fname,NIVolume(PSR;voxel_size=(tmp1,tmp2,tmp3)))
niwrite(R1f_fname,NIVolume(R1f;voxel_size=(tmp1,tmp2,tmp3)))
niwrite(SF_fname,NIVolume(Sf;voxel_size=(tmp1,tmp2,tmp3)))

niwrite(fname_1,NIVolume(DATA[:,:,:,1];voxel_size=(tmp1,tmp2,tmp3)))
niwrite(fname_2,NIVolume(DATA[:,:,:,2];voxel_size=(tmp1,tmp2,tmp3)))
niwrite(fname_3,NIVolume(DATA[:,:,:,3];voxel_size=(tmp1,tmp2,tmp3)))
niwrite(fname_4,NIVolume(DATA[:,:,:,4];voxel_size=(tmp1,tmp2,tmp3)))

println("The PSR is saved at $PSR_fname")
println("The R1f is saved at $R1f_fname")
println("The Sf is saved at $SF_fname")

