# % X: denoised matrix, s2: original noise variance, p: number of signal components, s2_after: noise variance after denoising
using LinearAlgebra

function  denoiseMatrix(X) 
    M = size(X,1);
    N = size(X,2);
    if M<N
        Xin = Array{Float64}(Yy')
        U,lambda = eigvals(X*Xin);
    else
        Xin = Array{Float64}(Yy')
        U,lambda = eigvals(Xin*X);
    end
    order_orig = sortperm(lambda);
    lambda = sort(lambda, rev=true);
    order = sortperm(lambda, rev=true);
    U = U[:,order];
    csum = cumsum(lambda[order_orig]);
    p = Array{Float64}([0:length(lambda)-1]');
    
    A = (lambda-lambda[end]).*(M-p).*(N-p)
    B = 4*csum*sqrt(M*N)
    p = -1 + Array{Int64}(findall(A < B));

    if p==0
        X = zeros(size(X));
    elseif M<N
        Binv = U[:,1:p]';
        X = U[:,1:p]*Binv*X;
    else
        Binv = U[:,1:p]';
        X = X*U[:,1:p]*Binv;
    end
    s2 = csum[p+1]/((M-p)*(N-p));
    s2_after = s2 - csum[p+1]/(M*N);
    return X,s2,p,s2_after 
end

denoiseMatrix