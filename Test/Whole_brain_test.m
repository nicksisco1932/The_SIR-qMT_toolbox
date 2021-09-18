clear;clc;close all
addpath(genpath('/Users/nicks/Documents/Github/sisco_toolbox_goodness/Matlab/mfiles/'))

path='/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618/';

SIR_4D_DATA= fullfile(path,'SIR_DATA.nii.gz');
brain_mask= fullfile(path,'brain_mask.nii.gz'); %brain_mask=$path/brain_mask.nii.gz


data = niftiread(SIR_4D_DATA);
[nx,ny,nz,nt] = size(data);
brain = niftiread(brain_mask);


% Set data sampling parameters
ti = [10 10 278 1007 ]*1e-3;
td = [684 4171 2730 10]*1e-3;

% Set model parameters
kmf = 12.5;
R1 = [1 1];
PSR = 0.1;
M0 = [1 PSR];
S = [-0.95 0.83];

vec_data = reshape(data,[nx*ny*nz,nt]);
vec_mask = reshape(brain,[nx*ny*nz,1]);
ind = find(vec_mask);
%%
yn = zeros(size(vec_data));
for ii = 1:length(ind)
    yy= vec_data(ind(ii),:);
    yn(ind(ii),:) = yy ./ yy(end);
end
%%
% 596389 voxels in brain mask
% Initial guesses and bounds
p0 = [0.07    1.5  -1.0    1.5]
lb = [0.00   0.2   -1.05   0.0];
ub = [1.00   5.0    0.00   9.9];

Xv = zeros(size(vec_data));

run_single = 0; % run once Elapsed time is 1253.928551 seconds 476 voxels per second
if run_single
    tic
    for ii = 1:length(ind)
        yy = yn(ind(ii),:);
        Xv(ii,:) = fitSIR_fixedkmf(ti,td,yy,S(2),p0,[],[],'y','n',kmf);
    end
    toc
end

% Elapsed time is 224.189014 seconds. -> 2660 voxels per second
tic
parfor ii = 1:length(ind)
    yy = yn(ind(ii),:);
    Xv(ii,:) = fitSIR_fixedkmf(ti,td,yy,S(2),p0,[],[],'y','n',kmf);
end
toc