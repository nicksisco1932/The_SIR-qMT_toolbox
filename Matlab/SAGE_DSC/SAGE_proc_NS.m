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
for index=73
    %%
    % if you want to use Volterra, put a 1 in the flag spot
%     [DSC,CBF_map,CBFSE_map,CBV_all,CBV_SE,MTT,MTT_SE] = sage_proc_ns_func(base_path,index,0);
%     [DSC_volterra,CBF_map_volterra,CBFSE_map_volterra,CBV_all_volterra,CBV_SE_volterra,MTT_volterra,MTT_SE_volterra] = sage_proc_ns_func(base_path,index,1);
    [DSC_block,CBF_map_block,CBFSE_map_block,CBV_all_block,CBV_SE_block,MTT_block,MTT_SE_block,OEF_block,CMRO2_block,pO2_block] = sage_proc_ns_func(base_path,index,2);
    
end
%%
% CBV_all(CBV_all>100)=0;
CBV_all_volterra(CBV_all_volterra>100)=0;
CBV_all_block(CBV_all_block>100)=0;
% CBV_SE(CBV_SE>100)=0;
CBV_SE_volterra(CBV_SE_volterra>100)=0;
CBV_SE_block(CBV_SE_block>100)=0;

% CBF_map(CBF_map>1000)=0;
CBF_map_volterra(CBF_map_volterra>1000)=0;
CBF_map_block(CBF_map_block>1000)=0;
% CBFSE_map(CBFSE_map>1000)=0;
CBFSE_map_volterra(CBFSE_map_volterra>1000)=0;
CBFSE_map_block(CBFSE_map_block>1000)=0;

% MTT(MTT>100)=0;
MTT_volterra(MTT_volterra>100)=0;
MTT_block(MTT_block>100)=0;
% MTT_SE(MTT_SE>100)=0;
MTT_SE_volterra(MTT_SE_volterra>100)=0;
MTT_SE_block(MTT_SE_block>100)=0;

kernel=[5 5 5];
% CBF_map = removeOutlierVolume(CBF_map,2,kernel);
% CBFSE_map = removeOutlierVolume(CBFSE_map,2,kernel);
% CBV_all = removeOutlierVolume(CBV_all,2,kernel);
% CBV_SE = removeOutlierVolume(CBV_SE,2,kernel);
% MTT = removeOutlierVolume(MTT,2,kernel);
% MTT_SE = removeOutlierVolume(MTT_SE,2,kernel);

CBF_map_block = removeOutlierVolume(CBF_map_block,2,kernel);
CBFSE_map_block = removeOutlierVolume(CBFSE_map_block,2,kernel);
CBV_all_block = removeOutlierVolume(CBV_all_block,2,kernel);
CBV_SE_block = removeOutlierVolume(CBV_SE_block,2,kernel);
MTT_block = removeOutlierVolume(MTT_block,2,kernel);
MTT_SE_block = removeOutlierVolume(MTT_SE_block,2,kernel);

CBF_map_volterra = removeOutlierVolume(CBF_map_volterra,2,kernel);
CBFSE_map_volterra = removeOutlierVolume(CBFSE_map_volterra,2,kernel);
CBV_all_volterra = removeOutlierVolume(CBV_all_volterra,2,kernel);
CBV_SE_volterra = removeOutlierVolume(CBV_SE_volterra,2,kernel);
MTT_volterra = removeOutlierVolume(MTT_volterra,2,kernel);
MTT_SE_volterra = removeOutlierVolume(MTT_SE_volterra,2,kernel);

 
% CBF_outlier_volterra = removeOutlierVolume(CBF_map_volterra,2,kernel);
% CBF_se_outlier_volterra = removeOutlierVolume(CBFSE_map_volterra,2,kernel);
% CBV_outlier_volterra = removeOutlierVolume(CBV_all_volterra,2,kernel);
% CBV_se_outlier_volterra = removeOutlierVolume(CBV_SE_volterra,2,kernel);
% mtt_outlier_volterra = removeOutlierVolume(MTT_volterra,2,kernel);
% mtt_se_outlier_volterra = removeOutlierVolume(MTT_SE_volterra,2,kernel);


%%

close all
% for ii = 1:DSC.Parms.nz
%     subplot(2,1,1)
%     imagesc(permute(MTT(:,:,ii),[2 1 3] ));colorbar
%     subplot(2,1,2)
%     imagesc(permute(MTT_SE(:,:,ii),[2 1 3] ));colorbar
%     pause(0.1)
% end
% 
% for ii = 1:DSC.Parms.nz
%     subplot(2,1,1)
%     imagesc(permute(CBV_all(:,:,ii),[2 1 3] ));colorbar
%     subplot(2,1,2)
%     imagesc(permute(CBV_SE(:,:,ii),[2 1 3] ));colorbar
%     pause(0.1)
% end
% 
% for ii = 1:DSC.Parms.nz
%     subplot(2,1,1)
%     imagesc(permute(CBF_map(:,:,ii),[2 1 3] ));colorbar
%     subplot(2,1,2)
%     imagesc(permute(CBFSE_map(:,:,ii),[2 1 3] ));colorbar
%     pause(0.1)
% end

subplot(2,2,1)
imagesc(permute(CBF_map_volterra(:,:,8),[2 1 3] ));colorbar
subplot(2,2,2)
imagesc(permute(CBFSE_map_volterra(:,:,8),[2 1 3] ));colorbar
subplot(2,2,3)
imagesc(permute(CBF_map_block(:,:,8),[2 1 3] ));colorbar
subplot(2,2,4)
imagesc(permute(CBFSE_map_block(:,:,8),[2 1 3] ));colorbar

%%

temp1 = (MTT_volterra);
temp2 = (MTT_SE_volterra);
temp3 = (MTT_block);
temp4 = (MTT_SE_volterra);
base='C:\Users\nicks\Documents\MRI_data\SAGE\AIF_optimizations\';
save(sprintf('%sMTT_volterra.mat',base),'temp1','-v6')
save(sprintf('%sMTTSE_volterra.mat',base),'temp2','-v6')
save(sprintf('%sMTT_block.mat',base),'temp3','-v6')
save(sprintf('%sMTTSE_block.mat',base),'temp4','-v6')

temp1 = (CBF_map_volterra);
temp2 = (CBFSE_map_volterra);
temp3 = (CBF_map_block);
temp4 = (CBFSE_map_block);

save(sprintf('%sCBF_volterra.mat',base),'temp1','-v6')
save(sprintf('%sCBFSE_volterra.mat',base),'temp2','-v6')
save(sprintf('%sCBF_block.mat',base),'temp3','-v6')
save(sprintf('%sCBFSE_block.mat',base),'temp4','-v6')

temp1 = (CBV_all_volterra);
temp2 = (CBV_SE_volterra);
temp3 = (CBV_all_block);
temp4 = (CBV_SE_block);

save(sprintf('%sCBV_volterra.mat',base),'temp1','-v6')
save(sprintf('%sCBVSE_volterra.mat',base),'temp2','-v6')
save(sprintf('%sCBV_block.mat',base),'temp3','-v6')
save(sprintf('%sCBVSE_block.mat',base),'temp4','-v6')

%%
base='C:\Users\nicks\Documents\MRI_data\SAGE\AIF_optimizations\';
info = niftiinfo(sprintf('%s/bPT1319073_preb_mask.nii.gz',base_path));
info.Datatype='double';

% niftiwrite(CBF_map_volterra,sprintf('%s/PT73_CBF_volterra.nii',base),info,'Compressed',1)
niftiwrite(CBF_map_block,sprintf('%s/PT73_CBF_block.nii',base),info,'Compressed',1)


niftiwrite(OEF_block,sprintf('%s/PT73_OEF_block.nii',base),info,'Compressed',1)
niftiwrite(CMRO2_block,sprintf('%s/PT73_CMRO2_block.nii',base),info,'Compressed',1)
niftiwrite(pO2_block,sprintf('%s/PT73_pO2_block.nii',base),info,'Compressed',1)


%% OEF and pO2


