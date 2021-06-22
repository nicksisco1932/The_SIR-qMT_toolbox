# % X: denoised matrix, s2: original noise variance, p: number of signal components, s2_after: noise variance after denoising
using LinearAlgebra
using Statistics
function  denoiseMatrix(X) 

    #= 
    From Mark Does Implimentation 
    =#

    M,N = size(X);
    minMN = minimum([M N]);
    Xm = mean(X,dims=2)
    X = X.-Xm

    if M<N
        Xin = Array{Float64}(X')
        # U,lambda = eigvals(X*Xin);
        # U = eigvecs(Xin*X);
        # lambda = eigvals(Xin*X);
        F = svd(Xin*X)
        U,S,V = F;
        lambda = S.^2/N;
    else
        Xin = Array{Float64}(X')
        
        # U = eigvecs(Xin*X);
        # lambda = eigvals(Xin*X);
        F = svd(X*Xin)
        U,S,V = F;
        lambda = S.^2/N;
    end
    
    p = 0;
    pTest = false;
    
    scaling = (M.-collect(0:minMN))/N;
    scaling[findall(x->x<1,scaling)] .= 1;
    # PSR[findall(x->x.>100,PSR)].=0

    while pTest == false
        sigma2 = (lambda[p+1] - lambda[minMN]) / (4*sqrt((M-p)/N))
        pTest = sum(lambda[p+1:minMN]/scaling[p+1]) >= (minMN-p)*sigma2
        if pTest == false
            p += 1
        end
    end
    sigma2 = sum(lambda[p+1:minMN]/ (minMN-p) /scaling[p+1])

    invV = Array{Float64}(V[:,1:p]')
    newX = U[:,1:p] * S[1:p,1:p] * invV + Xm

    return newX, sigma2, p

end

using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Optim

include("./utils.jl")


te1=0.007653;
te2=0.027969;
te3=0.059074;
te4=0.07939;
te5=0.099706;
echos=[te1,te2,te3,te4,te5]



a = LinRange(75,150,128);
b = LinRange(20,70,128);
# mat_img = zeros(128,128,15,5);
# R2_map = zeros(128,128);
# R2s_map = zeros(128,128);
# SNR=150;
# for ii in 1:128
#     for jj in 1:128
#         for kk in 1:15
#             X0 = [1,a[ii],b[jj],500];
#             Yy = SAGE_biexp3p_d(echos,X0);
#             mat_img[ii,jj,kk,:] = Yy+rand(5)/SNR
#         end
#     end
# end
X = rand(128,5)

denoiseMatrix(X)