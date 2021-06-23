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

function denoiseVol(vol,w,mask=[])
    nx_orig,ny_orig,nz_orig,nt_orig = size(vol);
    nx,ny,nz,nt = size(vol);

    if mask == Any[]
        mask = ones(size(vol[:,:,:,1]));
    end

    @assert(length(w)>1 && length(w)<4)

    if length(w) == 2
        tmp= ones(size(x)[1]+1)
        tmp[1:2] = w
        w = tmp
    end

    M = nt;
    N = prod(w)
    denoised = zeros(size(vol));

    p = zeros(nx,ny,nz);
    S2 = zeros(nx,ny,nz);

    m = nx - w[1]+1;
    n = ny - w[2]+1;
    o = nz - w[3]+1;
    

    for (count,ii) in enumerate(1:m*n*o)
        
        k =  floor((1-1)/m/n)+1;
        j = floor((ii-1-(k-1)*m*n)/m )+1
        i = ii - (k-1)*m*n-(j-1)*m;

        rows = Array{Int64}(collect(i:i-1+w[1]))
        cols = Array{Int64}(collect(j:j-1+w[2]))
        slis = Array{Int64}(collect(k:k-1+w[3]))

        maskCheck = reshape(mask[rows,cols,slis],N,1)

        X = reshape(vol[rows,cols,slis,:],M,N)

    end

    
end

X = rand(128,5)

dnX,sig,p = denoiseMatrix(X)

w = 5

using RollingFunctions

X = rand(128,128,15,5);

rolling(denoiseMatrix,X,w)