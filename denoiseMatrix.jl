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
    
    newX = zeros(size(vol));
    sigs = zeros(nx,ny,nz);
    MP = zeros(nx,ny,nz);

    Threads.@threads for ii in 1:m*n*o
        
        k =  floor((ii-1)/m/n)+1;    
        j = floor((ii-1-(k-1)*m*n)/m )+1
        i = ii - (k-1)*m*n-(j-1)*m;

        rows = Array{Int64}(collect(i:i-1+w[1]))
        cols = Array{Int64}(collect(j:j-1+w[2]))
        slis = Array{Int64}(collect(k:k-1+w[3]))

        maskCheck = reshape(mask[rows,cols,slis],N,1)

        X = reshape(vol[rows,cols,slis,:],N,M)

        dnX,sig,p = denoiseMatrix(X)
        # @show size(sig)
        # @show sig
        rdnX = reshape(dnX,w[1],w[2],w[3],M)
        # rsig = reshape(sig,w[1],w[2],w[3])
        # rmp = reshape(p,w[1],w[2],w[3])
        
        newX[rows,cols,slis,:] = rdnX
        sigs[rows,cols,slis] .= sig
        MP[rows,cols,slis] .= p


    end

    return newX,sigs,MP
end

using NIfTI
using Printf
DATA = zeros(224,224,55,4);
for ii in 1:4
    fname = @sprintf("/mnt/c/Users/nicks/Documents/Github/SIR_qMT_python/inputs/SIR_T1_%s.nii.gz",ii)
    img = niread(fname);
    DATA[:,:,:,ii] = img.raw
end

mask=[]
dn_DATA,dnsigs,dnMP = denoiseVol(DATA,w,mask);

niwrite("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/scratch/in.nii.gz",NIVolume(DATA))
niwrite("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/scratch/out.nii.gz",NIVolume(dn_DATA))
niwrite("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/scratch/res.nii.gz",NIVolume(DATA-dn_DATA))
niwrite("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/scratch/out_sig.nii.gz",NIVolume(dnsigs))
niwrite("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/scratch/out_MP.nii.gz",NIVolume(dnMP))

X = rand(128,5)

dnX,sig,p = denoiseMatrix(X)

w = [5,5,5]

X = rand(128,128,15,5);

using IterTools

# for i in partition(1:128,3,1)
#     for j in partition(1:128,3,1)
#         @show X[i,j,:,1]
#     end
# end

for i in subsets(1:128,5)
    for j in subsets(1:128,5)
        for k in subsets(1:15,5)
            # s = X[i,j,k,:]
            @show i
            # denoiseMatrix(s)
        end
    end
end


n = 128
m = 15
k = 3
for (nn,ii) in enumerate(1:k:n-k+1)
    for (nm,jj) in enumerate(1:k:n-k+1)
        for (no,kk) in enumerate(1:k:m-k+1)
            @show size(X[nn:ii,nm:jj,no:kk,:])
        end
    end
end

function dnNS(X,n,m,k)
    nx,ny,nz,nt = size(X)
    L = size(X)[end]
    newX = zeros(size(X))
    sigs = zeros(nx,ny,nz)
    MP = zeros(nx,ny,nz)

    for ii in 1:k:n-k+1
        for kk in 1:k:n-k+1
            for jj in 1:k:m-k+1
                o1 = ii:ii+k-1
                o2 = kk:kk+k-1
                o3 = jj:jj+k-1
                vec_X = reshape(X[o1,o2,o3,:],k*k*k,L)
                dnX,sig,p = denoiseMatrix(vec_X)
                rdnX = reshape(dnX,k,k,k,L)
                # rdnSig = reshape(sig,k,k,k)
                # rdnP = reshape(p,k,k,k)
                newX[o1,o2,o3,:] = rdnX
                # sigs[o1,o2,o3] = rdnSig
                # MP[o1,o2,o3] = rdnP
            end
        end
    end
    return newX
end
x = rand(128,128,15,5);
dn_X = dnNS(x,n,m,k)