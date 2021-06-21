#=
    Selective inversion recovery fitting using least squares and Levenberg-Marquardt 
    - This is a collection of helper functions written in Julia for fitting signal acquired from SIR-qMT to a double exponential 
    equation described in (1,2,3,4). 


    Input and Output Parameters for Code

    implimented within a Julia script or in REPL
    > include('./utils.jl')

        Future Updates TO Do:
        1) user defined ti and td, misc variables, B1 map for T1 MFA, etc.
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

################################################################################
##########################  User Supplied Arguments   ##########################
################################################################################
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin

        "input"
        help = "Input path for data. Currently, files should be read in separately."
        required = true
        
    end
    return parse_args(s)
end

################################################################################
##########################  User Supplied Arguments END  #######################
################################################################################

function reshape_and_normalize(data_4d::Array{Float64},mask_4d::Array{Bool},TI::Array{Float64},TD::Array{Float64},NX,NY,NZ,NT)
    # [pmf R1f Sf M0f]
    X=hcat(TI,TD)
    tmp = reshape(data_4d,(NX*NY*NZ,NT));
    tmpOUT = Array{Float64}(zeros(NX*NY*NZ,NT));
    vec_mask=reshape(mask_4d,(NX*NY*NZ))
    Yy=Array{Float64}(zeros(size(tmp)))

    tot_voxels = NX*NY*NZ

    mvec = vec_mask
    for ii = 1:tot_voxels
        if mvec[ii]
            Yy[ii,:] = tmp[ii, :] / (tmp[ii, end])
        end
    end

    Yy[findall(x->isinf(x),Yy)].=0
    Yy[findall(x->isnan(x),Yy)].=0

    return tot_voxels,vec_mask,tmpOUT,X,Yy
end

function params()
    inversion_time = Array{Float64}([15, 15, 278, 1007] * 1E-3);
    delay_time = Array{Float64}([684, 4121, 2730, 10] * 1E-3);
    X0 = Array{Float64}([0.1, 1, -0.95, 1]);
    return inversion_time, delay_time, X0
end

function load_sage_data(paths,b_fname)
    data1 = niread(paths[1]);
    brain_data = niread(b_fname);
    ne = size(paths)[1];
    nx,ny,nz,nt=size(data1);
    
    DATA = Array{Float64}(zeros(nx,ny,nz,ne,nt));
    for (n,ii) in enumerate(paths)
            DATA[:,:,:,n,:] = loader(ii);
    end    
    
    MASK=Array{Bool}(brain_data.raw);
    return DATA,MASK,nx,ny,nz,ne,nt
end

function load_data(paths,b_fname)
    data1 = niread(paths[1]);
    brain_data = niread(b_fname);
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
    MASK=Array{Bool}(brain_data.raw);
    return DATA,MASK,nx,ny,nz,nt
end

function loader(path)
    temp = niread(path)
    return Array{Float64}(temp.raw)
end

# function timed(tot_voxels,mvec::Vector{Bool},tmpOUT::Matrix{Float64},x::Array{Float64},yy::Array{Float64},x0::Array{Float64})::Array{Float64}
    
#     ind = Vector{Int64}(findall(x->x==Bool(1),mvec)); 
#     @time Threads.@threads for ii in ind
#             fit = curve_fit(SIR_signal_absolute,x,yy[:,ii],x0;autodiff=:finiteforward)
#             tmpOUT[:,ii] = fit.param
#         end
#     return tmpOUT
# end

# function alt_fit(tot_voxels,mvec::Vector{Bool},tmpOUT::Matrix{Float64},x::Array{Float64},yy::Array{Float64},x0::Array{Float64})
#         ## Optim Tricks
#         function sqerror(betas, X, Y)
#             err = 0.0
#             pred_i = SIR_signal(X,betas)
#             for ii in 1:length(Y)
#                 err += (abs.(pred_i[ii])-Y[ii]).^2
#             end
#             return err
#         end
#         using Optim
#         @time Threads.@threads for ii in 1:tot_voxels
#         if mvec[ii]
#             tmpOUT[ii,:]=Optim.minimizer(optimize(b -> sqerror(b, X,  yData[ii,:]), X0))
#         end
#     end
# end

function timed(tot_voxels,mvec::Vector{Bool},tmpOUT::Matrix{Float64},x::Array{Float64},yy::Array{Float64},x0::Array{Float64})
    @time if Threads.nthreads() == 1
        for ii in 1:tot_voxels
            if mvec[ii]
                fit = curve_fit(SIR_signal_absolute,x,yy[ii,:],x0;autodiff=:finiteforward)
                tmpOUT[ii,:] = fit.param
            end
        end
    else
    @time Threads.@threads for ii in 1:tot_voxels
            if mvec[ii]
                fit = curve_fit(SIR_signal_absolute,x,yy[ii,:],x0;autodiff=:finiteforward)
                tmpOUT[ii,:] = fit.param
            end
        end
    end
    return tmpOUT
end

function R1_minus_plus(R1f::Float64,R1m::Float64,pmf::Float64,kmf::Float64)
    # Calculate R1+/- in Eq. 4
    R1diff = sqrt((R1f - R1m + (pmf - 1) * kmf)^2 + 4 * pmf * kmf^2)
    R1plus = (R1f + R1m + (1 + pmf) * kmf + R1diff) / 2
    R1minus = R1plus - R1diff
    return R1minus,R1plus,R1diff
end
function amplitude(R1f::Float64,R1m::Float64,R1minus::Float64,R1diff::Float64,R1plus::Float64)
    # Component amplitude terms for td terms (Eq. 5)
    bftdplus = -(R1f - R1minus) / R1diff
    bftdminus = (R1f - R1plus) / R1diff
    bmtdplus = -(R1m - R1minus) / R1diff
    bmtdminus = (R1m - R1plus) / R1diff
    return bftdplus,bftdminus,bmtdminus,bmtdplus
end
function signal_recovery(R1plus::Float64,R1minus::Float64,td::Array{Float64},bftdplus::Float64,bftdminus::Float64,bmtdminus::Float64,bmtdplus::Float64)
    # Signal recovery during td (Eq. 5)
    Mftd = bftdplus .* exp.(-R1plus * td) + bftdminus * exp.(-R1minus * td) .+ 1
    Mmtd = bmtdplus .* exp.(-R1plus * td) + bmtdminus * exp.(-R1minus * td) .+ 1
    return Mftd,Mmtd
end
function amplitude_ti(Sf::Float64,Sm::Float64,R1f::Float64,R1minus::Float64,R1plus::Float64,R1diff::Float64,Mftd::Array{Float64},Mmtd::Array{Float64},pmf::Float64,kmf::Float64)
    # Component amplitude terms for ti terms (Eq. 5)
    bfplus = (
            (Sf * Mftd .- 1) * (R1f .- R1minus) .+
            (Sf .* Mftd .- Sm * Mmtd) .* pmf .* kmf  ) ./ R1diff
    bfminus =
        -(
            (Sf * Mftd .- 1) * (R1f .- R1plus) .+
            (Sf .* Mftd .- Sm * Mmtd) .* pmf .* kmf
        ) ./ R1diff
    return bfminus,bfplus
end
function signal_equation(R1plus::Float64,R1minus::Float64,bfplus::Array{Float64},bfminus::Array{Float64},ti::Array{Float64},Mfinf::Float64)
     # Signal equation (Eq. 3)
     M =
     (
         bfplus .* exp.(-R1plus .* ti) .+ bfminus .* exp.(-R1minus .* ti) .+
         1
     ) .* Mfinf
     return M
end

function SIR_signal_absolute(x::Matrix{Float64}, p::Vector{Float64})
#=
    # Selective inversion recovery signal for fitting absolute value siganl

    - Requirements
    1) R1_minus_plus
    2) amplitude
    3) signal_recovery
    4) amplitude_ti
    5) signal_equation

    - Constants
        * currently there are only two constants that can be fixed in later versions.
            - kmf is held constant at 14.5 for human brains and 35 for phantoms
            - Sm is held constant at 0.83 
=#
    ti = x[:,1]
    td = x[:,2]
    kmf = 14.5; # not sure how to pass kwargs to this function using curve fit
    # kmf = 35.0
    M0 = [1, p[1]] * p[4]
    R1 = [p[2], p[2]]
    S = [p[3], 0.83]
    # Get pmf and Mfinf
    pmf = M0[2] / M0[1]
    Mfinf = M0[1]

    # Get R1f/R1m and Sf/Sm
    R1f = R1[1]
    R1m = R1[2]
    Sf = S[1]
    Sm = S[2]

    R1minus,R1plus,R1diff = R1_minus_plus(R1f,R1m,pmf,kmf)
    bftdplus,bftdminus,bmtdminus,bmtdplus = amplitude(R1f,R1m,R1minus,R1diff,R1plus)
    Mftd,Mmtd = signal_recovery(R1plus,R1minus,td,bftdplus,bftdminus,bmtdminus,bmtdplus)
    bfminus,bfplus = amplitude_ti(Sf,Sm,R1f,R1minus,R1plus,R1diff,Mftd,Mmtd,pmf,kmf)
    Mm = signal_equation(R1plus,R1minus,bfplus,bfminus,ti,Mfinf)
    M = abs.(Mm)
    return M 
end

function SIR_signal(x::Matrix{Float64},p::Vector{Float64})

    # Selective inversion recovery signal for fitting/simulating  siganl

    # THis is more or less for plotting
    ti = x[:,1]
    td = x[:,2]
    kmf = 14.5
    M0 = [1, p[1]] * p[4]
    R1 = [p[2], p[2]]
    S = [p[3], 0.83]
    # Get pmf and Mfinf
    pmf = M0[2] / M0[1]
    Mfinf = M0[1]

    # Get R1f/R1m and Sf/Sm
    R1f = R1[1]
    R1m = R1[2]
    Sf = S[1]
    Sm = S[2]

    # Calculate R1+/- in Eq. 4
    R1diff = sqrt((R1f - R1m + (pmf - 1) * kmf)^2 + 4 * pmf * kmf^2)
    R1plus = (R1f + R1m + (1 + pmf) * kmf + R1diff) / 2
    R1minus = R1plus - R1diff

    # Component amplitude terms for td terms (Eq. 5)
    bftdplus = -(R1f - R1minus) / R1diff
    bftdminus = (R1f - R1plus) / R1diff
    bmtdplus = -(R1m - R1minus) / R1diff
    bmtdminus = (R1m - R1plus) / R1diff

    # Signal recovery during td (Eq. 5)
    Mftd = bftdplus .* exp.(-R1plus * td) + bftdminus * exp.(-R1minus * td) .+ 1
    Mmtd = bmtdplus .* exp.(-R1plus * td) + bmtdminus * exp.(-R1minus * td) .+ 1

    # Component amplitude terms for ti terms (Eq. 5)
    bfplus =
        (
            (Sf * Mftd .- 1) * (R1f .- R1minus) .+
            (Sf .* Mftd .- Sm * Mmtd) .* pmf .* kmf
        ) ./ R1diff
    bfminus =
        -(
            (Sf * Mftd .- 1) * (R1f .- R1plus) .+
            (Sf .* Mftd .- Sm * Mmtd) .* pmf .* kmf
        ) ./ R1diff

    # Signal equation (Eq. 3)
    M =
        (
            bfplus .* exp.(-R1plus .* ti) .+ bfminus .* exp.(-R1minus .* ti) .+
            1
        ) .* Mfinf
    return M
end

function nlsFit_MFA(mask::Vector{Bool},tr::Float64,flip_angles::Vector{Float64},yData::Array{Float64},alpha::Float64)
    n,nt = size(yData)
    params = zeros(nt,2)
    params[:,1] = flip_angles
    params[1,2] = tr
    params[2,2] = alpha

    # tot_voxels = n;
    ind = findall(x->x==Bool(1),mask)
    x0 = [maximum(yData),1.5]
    parm_num = size(x0)[1]
    out = zeros(n,parm_num)
    
    if Threads.nthreads() == 1
        @time for ii in ind
            fit = curve_fit(MFA,params,yData[ii,:],x0)
            out[ii,:] = fit.param
        end
    else
        @time Threads.@threads for ii in ind
            fit = curve_fit(MFA,params,yData[ii,:],x0)
            out[ii,:] = fit.param
        end
    end

    return out
end


# function MFA(fa::Vector{Float64},p::Vector{Float64},TR::Float64,alpha::Float64)::Array{Float64}
function MFA(x::Matrix{Float64},p::Vector{Float64})
    #p[1] = M0
    #p[2] = T1
    
    # FA = [20 18 16 14 12 10 8 6 4 2]; %default vals
    # TR = 7.590/1000; %default val

    fa = x[:,1]
    TR = x[1,2]
    alpha = x[2,2]

    A = sind.(alpha.*fa)
    B = 1-exp.(-TR/p[2])
    C = cosd.(alpha.*fa)

    # S_mfa = p[1].*((A.*(B))./(B.*C))
    S_mfa = p[1] .* (sind.(fa.*alpha) .* (1-exp.(-TR/p[2]))) ./ (1 .-cosd.(fa.*alpha) .* exp.(-TR./p[2]))
    return S_mfa

end

## Optim Tricks
function sqerror(betas, X, Y)
    # fa = x[:,1]
    # TR = x[1,2]
    # alpha = x[2,2]
    err = 0.0
    pred_i = MFA(X,betas)
    for ii in 1:length(Y)
        err += (abs.(pred_i[ii])-Y[ii]).^2
    end
    return err
end

function alt_fit(mask::Vector{Bool},tr::Float64,flip_angles::Vector{Float64},yData::Array{Float64},alpha::Float64)
    

    n,nt = size(yData)
    params = zeros(nt,2)
    params[:,1] = flip_angles
    params[1,2] = tr
    params[2,2] = alpha

    # tot_voxels = n;
    ind = findall(x->x==Bool(1),mask)
    x0 = [maximum(yData),1.5]
    parm_num = size(x0)[1]
    out = zeros(n,parm_num)

    if Threads.nthreads()==1
        @time for ii in ind
        out[ii,:]=Optim.minimizer(optimize(b -> sqerror(b, params,  yData[ii,:]), x0))
        end
        else
    @time Threads.@threads for ii in ind
    
        out[ii,:]=Optim.minimizer(optimize(b -> sqerror(b, params,  yData[ii,:]), x0))
    
    end
end
    return out
end

function SAGE_biexp3p_d(te,x)
    dtemp=1;
    #x should be [SI_I R2s R2] where SI_I = SI_II; 
    tn=te;
    TE=te[end];
    SI_sage = zeros(size(tn));
    # %Create piece-wise function based on Schmiedeskamp H et al. MRM 2012
    # %67:378-388
    ind1 = findall(x->x<TE/2,tn)
    ind2 = findall(x->x>TE/2,tn)
    SI_sage[ind1] = x[1].*exp.(-tn[ind1].*x[2]);
    SI_sage[ind2] = (x[4]).*exp.(-TE*(x[2]-x[3])).*exp.(-tn[ind2]*(2*x[3]-x[2]));
    return SI_sage
end
## Optim Tricks
function sqerrorSAGE(betas::Vector{Float64}, X::Vector{Float64}, Y::Vector{Float64})
    err = 0.0
    pred_i = SAGE_biexp3p_d(X,betas)
    for ii in 1:length(Y)
        err += (abs.(pred_i[ii])-Y[ii]).^2
    end
    return err
end

