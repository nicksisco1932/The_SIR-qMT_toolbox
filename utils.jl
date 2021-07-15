#=
    
    - This is a collection of helper functions written in Julia for The MRI Toolbox

    Input and Output Parameters for Code

    implimented within a Julia script or in REPL
    > include('./utils.jl')

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

function reshape_and_normalize(data_4d::Array{Float64},TI::Array{Float64},TD::Array{Float64},NX::Int64,NY::Int64,NZ::Int64,NT)
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

function nlsfit(f::Function, xvalues::Array{Float64},yvalues::Array{Float64},guesses::Array{Float64})
    fit = curve_fit(f,xvalues,yvalues,guesses;autodiff=:finiteforward)
    return fit.param
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

#= TO BE DELETED
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

# function timed(tot_voxels,mvec::Vector{Bool},tmpOUT::Matrix{Float64},x::Array{Float64},yy::Array{Float64},x0::Array{Float64})
#     @time if Threads.nthreads() == 1
#         for ii in 1:tot_voxels
#             if mvec[ii]
#                 fit = curve_fit(SIR_signal_absolute,x,yy[ii,:],x0;autodiff=:finiteforward)
#                 tmpOUT[ii,:] = fit.param
#             end
#         end
#     else
#     @time Threads.@threads for ii in 1:tot_voxels
#             if mvec[ii]
#                 fit = curve_fit(SIR_signal_absolute,x,yy[ii,:],x0;autodiff=:finiteforward)
#                 tmpOUT[ii,:] = fit.param
#             end
#         end
#     end
#     return tmpOUT
# end

# function R1_minus_plus(R1f::Float64,R1m::Float64,pmf::Float64,kmf::Float64)
#     # Calculate R1+/- in Eq. 4
#     R1diff = sqrt((R1f - R1m + (pmf - 1) * kmf)^2 + 4 * pmf * kmf^2)
#     R1plus = (R1f + R1m + (1 + pmf) * kmf + R1diff) / 2
#     R1minus = R1plus - R1diff
#     return R1minus,R1plus,R1diff
# end
# function amplitude(R1f::Float64,R1m::Float64,R1minus::Float64,R1diff::Float64,R1plus::Float64)
#     # Component amplitude terms for td terms (Eq. 5)
#     bftdplus = -(R1f - R1minus) / R1diff
#     bftdminus = (R1f - R1plus) / R1diff
#     bmtdplus = -(R1m - R1minus) / R1diff
#     bmtdminus = (R1m - R1plus) / R1diff
#     return bftdplus,bftdminus,bmtdminus,bmtdplus
# end
# function signal_recovery(R1plus::Float64,R1minus::Float64,td::Array{Float64},bftdplus::Float64,bftdminus::Float64,bmtdminus::Float64,bmtdplus::Float64)
#     # Signal recovery during td (Eq. 5)
#     Mftd = bftdplus .* exp.(-R1plus * td) + bftdminus * exp.(-R1minus * td) .+ 1
#     Mmtd = bmtdplus .* exp.(-R1plus * td) + bmtdminus * exp.(-R1minus * td) .+ 1
#     return Mftd,Mmtd
# end
# function amplitude_ti(Sf::Float64,Sm::Float64,R1f::Float64,R1minus::Float64,R1plus::Float64,R1diff::Float64,Mftd::Array{Float64},Mmtd::Array{Float64},pmf::Float64,kmf::Float64)
#     # Component amplitude terms for ti terms (Eq. 5)
#     bfplus = (
#             (Sf * Mftd .- 1) * (R1f .- R1minus) .+
#             (Sf .* Mftd .- Sm * Mmtd) .* pmf .* kmf  ) ./ R1diff
#     bfminus =
#         -(
#             (Sf * Mftd .- 1) * (R1f .- R1plus) .+
#             (Sf .* Mftd .- Sm * Mmtd) .* pmf .* kmf
#         ) ./ R1diff
#     return bfminus,bfplus
# end
# function signal_equation(R1plus::Float64,R1minus::Float64,bfplus::Array{Float64},bfminus::Array{Float64},ti::Array{Float64},Mfinf::Float64)
#      # Signal equation (Eq. 3)
#      M =
#      (
#          bfplus .* exp.(-R1plus .* ti) .+ bfminus .* exp.(-R1minus .* ti) .+
#          1
#      ) .* Mfinf
#      return M
# end
=#

function SIR_Mz0(x::Matrix{Float64},p::Vector{Float64}, kmf::Float64;
    Sm::Float64=0.83, R1m::Float64=NaN, mag::Bool=true)

    # Extract ti and td values from x
    ti = x[:,1]
    td = x[:,2]

    # Define model parameters based on p
    pmf = p[1]
    R1f = p[2]
    Sf  = p[3]
    Mf∞ = p[4]

    # Define R1m based on user-defined value (=R1f when set to NaN)
    if isnan(R1m)
        R1m = R1f
    end

    # Define kfm based on kmf and pmf (assuming mass balance)
    kfm = kmf*pmf

    # Apparent rate constants
    ΔR1 = sqrt((R1f-R1m+kfm-kmf)^2.0 + 4.0*kfm*kmf)
    R1⁺ = (R1f + R1m + kfm + kmf + ΔR1) / 2.0
    R1⁻ = R1⁺ - ΔR1

    # Component amplitudes for td terms
    bf_td⁺ = -(R1f - R1⁻) / ΔR1
    bf_td⁻ =  (R1f - R1⁺) / ΔR1
    bm_td⁺ = -(R1m - R1⁻) / ΔR1
    bm_td⁻ =  (R1m - R1⁺) / ΔR1

    # Loop over ti/td values
    # make this a new function
    M = similar(ti)
    for k in 1:length(td)

        # Signal recovery during td
        E_td⁺ = exp(-R1⁺*td[k])
        E_td⁻ = exp(-R1⁻*td[k])
        Mf_td = bf_td⁺*E_td⁺ + bf_td⁻*E_td⁻ + 1.0
        Mm_td = bm_td⁺*E_td⁺ + bm_td⁻*E_td⁻ + 1.0

        # Component amplitude terms for ti terms
        a = Sf*Mf_td - 1.0
        b = (Sf*Mf_td - Sm*Mm_td) * kfm
        bf_ti⁺ =  (a*(R1f-R1⁻) + b) / ΔR1
        bf_ti⁻ = -(a*(R1f-R1⁺) + b) / ΔR1

        # Signal recovery during ti
        M[k] = (bf_ti⁺*exp(-R1⁺*ti[k]) + bf_ti⁻*exp(-R1⁻*ti[k]) + 1.0) * Mf∞

        # Take the magnitude of the signal
        if mag
            M[k] = abs(M[k])
        end
    end

    # Return signal
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

    R₂star = x[2]
    R₂ = x[3]
    S₀I=x[1]
    S₀II=x[4]
    
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

