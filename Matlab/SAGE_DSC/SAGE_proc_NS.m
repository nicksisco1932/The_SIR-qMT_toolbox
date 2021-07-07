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
    [DSC,CBF_map,CBFSE_map,CBV_all,CBV_SE,MTT,MTT_SE] = sage_proc_ns_func(base_path,index);

end