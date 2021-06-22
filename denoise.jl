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
    p = -1 + find((lambda-lambda[end]).*(M-p).*(N-p) < 4*csum*sqrt(M*N),1);
    if p==0
        X = zeros(size(X));
    elseif M<N
        X = U(:,1:p)*U(:,1:p)'*X;
    else
        X = X*U(:,1:p)*U(:,1:p)';
    end
    s2 = csum(p+1)/((M-p)*(N-p));
    s2_after = s2 - csum[p+1]/(M*N);
    return X,s2,p,s2_after 
end