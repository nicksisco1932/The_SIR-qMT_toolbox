function IMG = removeOutlierVolume(IMG,ns,K)

% Get size of image
[nr,nc,nz] = size(IMG);

% Find number of expected outliers from normal dist
nOutlierExp =  round((2*(1 - normcdf(ns)))*length(find(IMG)));

% Set filter kernel
H = fspecial3('gaussian',K);

% Loop until convergence
nOutlierCur  = nOutlierExp*2;
while (nOutlierCur > nOutlierExp) 
    
    % Median filter
    %IMGf = medfilt3(IMG,K);
    IMGf = imfilter(IMG,H,'replicate','conv');
    
    % Take difference, set cuoff as mean diff + 3SD
    D = abs(IMG - IMGf);
    cutoff = mean(D(find(D))) + ns*std(D(find(D))); 
    ind = find(D > cutoff);
    nOutlierCur = length(ind); 
    
    % Replace IMG w/ filtered value at these pixels
    ind = find(D > cutoff);
    IMG(ind) = IMGf(ind);
    % imagesc(myMontage(IMG),[0 .25]), pause
    
end
