% Parallel loop
X0 = [0.1  1   -0.95  1];  % [pmf R1f Sf M0f]
LB = [0    0.3 -1.05  0];
UB = [1    3    0    10];
% Xv = zeros(length(ind),4);
ti = x(:,1);
td = x(:,2);
tic
parfor ii = 1:length(yn)
%    disp(['Fitting SIR data: Voxel: ' num2str(ii) '/' num2str(length(ind))])
    fitSIR_fixedkmf(ti,td,yn(:,ii),Sm,X0,[],[],'y','n',kmf);
end
toc