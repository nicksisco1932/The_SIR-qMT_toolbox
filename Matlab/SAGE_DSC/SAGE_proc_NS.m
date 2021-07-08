% New proc script after misplacing data.

%===========================================================%
% Simple multi echo perfusion calculations using SAGE       %
% Much of this was adapted from Ashley M. Stokes, Ph.D.     %
%                                                           %
% Nicholas J. Sisco, Ph.D.                                  %
%===========================================================%

clc;clear;close all
script_path = "C:\Users\nicks\Documents\Github\The_MRI_toolbox\Matlab\SAGE_DSC\";
cd(script_path);
addpath(genpath('./mfiles/'))


% ptnums=1:100; %does not matter, keep the range
ptnums=1; %does not matter, keep the range

temp = uigetdir("C:\Users\nicks\Documents\GitHub\DSC_SAGE_python\");
out_path = "C:\Users\nicks\Documents\MRI_data\SAGE\DSC_standard_proc\";

f=dir(temp);
base_path=f(1).folder;

% base_path = uigetdir('C:/Users/nicks/Box/MSPerfusion Data/SAGE_dcm_converted/');
% brainMask_path = uigetdir('C:/Users/nicks/Box/MSPerfusion Data/SAGE_dcm_converted/');
% base_path = "/Users/nicks/Box/MSPerfusion Data/SAGE_dcm_converted/SAGE_niftis/";
% brainMask_path = "/Users/nicks/Box/MSPerfusion Data/SAGE_dcm_converted/SAGE_prebolus_TE5/";
if ~base_path
    return
end
%%
for index=1:ptnums(end)
    %%
    % if you want to use Volterra, put a 1 in the flag spot
    [DSC,CBF_map,CBFSE_map,CBV_all,CBV_SE,MTT,MTT_SE] = sage_proc_ns_func(base_path,index,1);
    
end

%%
MTT(MTT>100)=0;
MTT_SE(MTT_SE>100)=0;
close all
for ii = 1:DSC.Parms.nz
    subplot(2,1,1)
    imagesc(permute(MTT(:,:,ii),[2 1 3] ));colorbar
    subplot(2,1,2)
    imagesc(permute(MTT_SE(:,:,ii),[2 1 3] ));colorbar
    pause(0.1)
end

for ii = 1:DSC.Parms.nz
    subplot(2,1,1)
    imagesc(permute(CBV_all(:,:,ii),[2 1 3] ));colorbar
    subplot(2,1,2)
    imagesc(permute(CBV_SE(:,:,ii),[2 1 3] ));colorbar
    pause(0.1)
end

for ii = 1:DSC.Parms.nz
    subplot(2,1,1)
    imagesc(permute(CBF_map(:,:,ii),[2 1 3] ));colorbar
    subplot(2,1,2)
    imagesc(permute(CBFSE_map(:,:,ii),[2 1 3] ));colorbar
    pause(0.1)
end

subplot(2,1,1)
imagesc(permute(MTT(:,:,8),[2 1 3] ));colorbar
subplot(2,1,2)
imagesc(permute(MTT_SE(:,:,8),[2 1 3] ));colorbar