%% Sim signals, fit in MATLAB, and save for Julia
clear, close all, clc

% Create pool of workers
warning off
delete(gcp('nocreate'))
pool = parpool(4);
warning on

% Set data sampling parameters
ti = [10 10 278 1007 ]*1e-3
td = [684 4171 2730 10]*1e-3

% Set model parameters
kmf = 12.5;
R1 = [1 1];
PSR = 0.1;
M0 = [1 PSR];
S = [-0.95 0.83];

% Set SNR and number of trials
SNR = 2e2;
N = 5e5

% Calculate noise and add to signal
y = signalSIR(ti,td,kmf,R1,M0,S)'
n = randn(length(y),N)*(1/SNR);
yn = abs(repmat(y,[1 N]) + n);

% Initial guesses and bounds
p0 = [0.07    1.5  -1.0    1.5]
lb = [0.00   0.2   -1.05   0.0];
ub = [1.00   5.0    0.00   9.9];

% No threading
tic
pm1 = zeros(N,4);
for ii = 1:N
    pm1(ii,:) = fitSIR_fixedkmf(ti,td,yn(:,ii),S(2),p0,[],[],'y','n',kmf);
end
toc

% Parallel loop
tic
pm2 = zeros(N,4);
parfor ii = 1:N
    pm2(ii,:)= fitSIR_fixedkmf(ti,td,yn(:,ii),S(2),p0,[],[],'y','n',kmf);
end
toc

% Save data
Sm = S(2);
pact = [PSR R1(1) S(1) M0(1)];
x = [ti' td'];
save matData p0 pact pm1 pm2 x y yn kmf Sm

%% Run Julia
setenv('JULIADIR','/Applications/Julia-1.5.app/Contents/Resources/julia/bin');
setenv('PATH',[getenv('PATH') ':' getenv('JULIADIR')])
unix('julia -t 4 fitSIR_test.jl');

%% Compare MATLAB and Julia
clear, close all, clc

load matData
load juliaData
pj1 = pj1';
pj2 = pj2';

% Max diff - should be within convergence criteria
julMatMaxDiff = max(pm2(:)-pj2(:))

% Plot PSR for both to visualize any difference 
plot(pj2(1:100:end,1))
hold on
plot(pm2(1:100:end,1),'r--')

