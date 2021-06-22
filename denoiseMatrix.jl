# % X: denoised matrix, s2: original noise variance, p: number of signal components, s2_after: noise variance after denoising
using LinearAlgebra
using Statistics
function  denoiseMatrix(X) 

    #= 
    From Mark Does Implimentation 
    =#

    M,N = size(X);
    minMN = minimum([M N]);
    Xm = mean(X,dims=2)[:,1]
    X = X.-Xm

    F = svd(X)
    U,S,V = F;
    lambda = S.^2/N;

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

    newX = U[:,1:p] * diagm(S)[1:p,1:p] * transpose(V[:,1:p]) .+ Xm

    return newX, sigma2, p

end


X = rand(128,5)

dnX,sig,p = denoiseMatrix(X)