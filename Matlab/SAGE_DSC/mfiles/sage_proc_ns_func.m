function [DSC,CBF_map,CBFSE_map,CBV_all,CBV_SE,MTT,MTT_SE,OEF,CMRO2,pO2] = sage_proc_ns_func(temp,index,FLAG)
%     Author: Nicholas J. Sisco, Ph.D. 
%     Email: X@barrowneuro.org where X = nicholas.sisco
%     Affiliation: Barrow Neurological Institue in Phoenix, AZ
%     Title: Processing script for Dynamic Susceptibility Contrast
%               -specifically for multiecho. Many other programs can do
%               single echo
%
%     Requirements:
%    AutoAIF_Brain.m       
%    CBF_calc_volterra.m   
%    adaptive_threshold.m  
%    circulant.m           
%    loader.m
%    CBF_calc.m
%    CBV_calc.m
%    aif_processing.m
%    linearizing_data.m
%    sage_proc_ns_func.m 
%     Output:
%         - The output from this script are CBF,CBV,MTT for SAGE. 

    
    fname = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(1)));
    info = niftiinfo(fname);
    tmp_data = niftiread(info);
    [nx,ny,nz,nt] = size(tmp_data);
    ne = 5;
    DATA = zeros(nx,ny,nz,ne,nt);
    for echos = 1:5
        fname = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(echos)));
        [tmp_data] = loader(fname);
        DATA(:,:,:,echos,:) = tmp_data;
    end
    DSC.Data = DATA;clear DATA;
    
    
    %% Brain Mask
    % If you need to make one, this will work.
    %     tmp = squeeze(DSC.Data(:,:,:,2,1));
    %     scale_fac = 1;
    %     tmp(tmp<threshold3(tmp)+scale_fac )=0;
    %     mask = tmp;
    %     mask(mask>0)=1;clear tmp;
    %     close all
    %     brain_img = erodeBWperim(mask,5); % Just to be conservative
    %     imagesc(brain_img(:,:,10)); % check to see if the mask is ok, otherwise change scale

    bfname = fullfile(temp,sprintf("bPT1319%03i_preb_mask.nii.gz",index));
    brain_img = double(niftiread(bfname));
    
    
    
    denoise=0;
    plots=0;
    
    if denoise~=0
        warning('Denoising is ON!')
        warning('Denoising has not yeilded good results with SAGE')
        DSC.DATAdn.DATA=zeros(nx,ny,nz,ne,nt);
        DSC.DATAdn.m=zeros(nx,ny,nz,nt);
        DSC.DATAdn.p=zeros(nx,ny,nz,nt);
        for ii = 1:nt
            [DSC.DATAdn.DATA(:,:,:,:,ii),DSC.DATAdn.m(:,:,:,ii),DSC.DATAdn.p(:,:,:,ii)] = denoiseCV(squeeze(DSC.Data(:,:,:,:,1)),[5 5 5],logical(brain_img));
            if plots==1
                subplot(2,3,1)
                imagesc(permute(DSC.DATAdn.DATA(:,:,9,5),[2 1 3]))
                title("Orig")
                subplot(2,3,2)
                imagesc(permute(DSC.Data(:,:,9,5),[2 1 3]))
                title("Denoised")
                subplot(2,3,3)
                imagesc(permute(DSC.DATAdn.DATA(:,:,9,5),[2 1 3])-permute(DSC.Data(:,:,9,5),[2 1 3]));colormap gray;colorbar
                title("Residual")
                subplot(2,3,4)
                imagesc(permute(DSC.DATAdn.m(:,:,9),[2 1 3]));colormap gray;colorbar
                title("M")
                subplot(2,3,5)
                imagesc(permute(DSC.DATAdn.p(:,:,9),[2 1 3]));colormap gray;colorbar
                title("P")
            end
        end
        DSC.Data = DSC.DATAdn.DATA;
    end
    
    %% SAGE DSC processing
    DSC.Parms.TR = 1800./1000;                          %TR in s
    DSC.Parms.time = DSC.Parms.TR:DSC.Parms.TR:DSC.Parms.TR*size(DSC.Data,5);
    DSC.Parms.TE1 = 8.8./1000;                          %TE in s
    DSC.Parms.TE2 = 26./1000;                           %TE in s
    DSC.Parms.TE3 = 49./1000;                           %TE in s
    DSC.Parms.TE4 = 66./1000;                           %TE in s
    DSC.Parms.TE5 = 88./1000;                           %TE in s
    DSC.Parms.flip= 90;                                 %FA in degrees
    [DSC.Parms.nx,DSC.Parms.ny,DSC.Parms.nz,DSC.Parms.ne,DSC.Parms.nt] = size(DSC.Data);

    %% process data - brain mask and dR2s
    filtered_image = abs(double(DSC.Data)); % I don't understand the point of this

    [DSC.Parms.ss_tp, DSC.Parms.gd_tp,DSC.Parms.pk_tp,~] = determineTimePoints_NS(filtered_image(:,:,:,1:2,:),brain_img,DSC.Parms.TR,0,0);
    prebolus_signal = squeeze(nanmean(filtered_image(:,:,:,:,DSC.Parms.ss_tp:DSC.Parms.gd_tp),5)); %average prebolus images together


    %% AIF and enhancing masks
    DSC.AIF = AutoAIF_Brain(filtered_image(:,:,:,1:2,:),[DSC.Parms.TR DSC.Parms.TE1 DSC.Parms.TE2],DSC.Parms.flip,0,1,double(brain_img),DSC.Parms.ss_tp,DSC.Parms.gd_tp,DSC.Parms.pk_tp);

    %%
    if FLAG==1
        AIF=1;
    elseif FLAG==2
        AIF=2;
    else
        AIF=0;
    end
    
    if AIF==1       
        [dR2star_all,dR2star_SE,CTC_all,CTC_SE] = linearizing_data(DSC,filtered_image,brain_img);
    elseif AIF==2
        [dR2star_all,dR2star_SE,CTC_all,CTC_SE] = linearizing_data(DSC,filtered_image,brain_img);
    else
        [dR2star_all,dR2star_SE,~,~] = linearizing_data(DSC,filtered_image,brain_img);
    end
    
    julia_fitting=1;
    if julia_fitting~=0
        out_path=sprintf("/Users/nicks/Documents/MRI_data/SAGE/nifti/PT%s/",num2str(index));
        if ~isfolder(out_path)
            mkdir(out_path)
        end
        te1 = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(1)));
        te2 = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(2)));
        te3 = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(3)));
        te4 = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(4)));
        te5 = fullfile(temp,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(5)));
        copyfile(te1,fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(1))))
        copyfile(te2,fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(2))))
        copyfile(te3,fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(3))))
        copyfile(te4,fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(4))))
        copyfile(te5,fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(5))))
        copyfile(bfname,fullfile(out_path,sprintf("PT1319%03i_brain_mask.nii.gz",index)));
        te1 = fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(1)));
        te2 = fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(2)));
        te3 = fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(3)));
        te4 = fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(4)));
        te5 = fullfile(out_path,sprintf("PT1319%03i_TE%s_img_w_Skull.nii.gz",index,num2str(5)));
        bfname=fullfile(out_path,sprintf("PT1319%03i_brain_mask.nii.gz",index));
        [dR2star_all,dR2star_SE] = sage_julia_fitty(out_path,te1,te2,te3,te4,te5,bfname,DSC);
    end
    
    %% CBV
    %--------------------------
    [CBV_all,CBV_SE,tidx] = CBV_calc(DSC,dR2star_all,dR2star_SE,brain_img);
    %--------------------------

    [threshold] = adaptive_threshold(DSC,filtered_image,prebolus_signal);

    %% CBF

    if AIF==1
        [CBF_map,CBFSE_map] = CBF_calc_volterra(DSC,brain_img,threshold,CTC_all,CTC_SE);
    elseif AIF==2
        [CBF_map,CBFSE_map] = CBF_calc_block_circulant(DSC,brain_img,threshold,CTC_all,CTC_SE);
    else
        
        [S_orig,maxS,U,S,V,AIFmatrixt,dtemp_all,dtemp_allSE] = aif_processing(DSC,dR2star_all,dR2star_SE,tidx);
        %--------------------------
        [CBF_map,CBFSE_map] = CBF_calc(dtemp_all,dtemp_allSE,DSC,U,S,V,maxS,S_orig,brain_img,threshold);
        %--------------------------
    end

    
    ind = find(brain_img);
    MTT = CBV_all./CBF_map;
    MTT(~ind)=0;
    MTT(isinf(MTT))=0;
    MTT_SE = CBV_SE./CBFSE_map;
    MTT_SE(~ind)=0;
    MTT_SE(isinf(MTT_SE))=0;
    
    CBV_all=0.733*100.*CBV_all;
    CBV_SE=0.733*100.*CBV_SE;
    CBF_map=6000.*CBF_map;
    CBFSE_map=6000.*CBFSE_map;
    
    %% OEF,CMRO2,pO2
    [OEF,CMRO2,pO2] = OEF_CMRO2_pO2(dR2star_all,dR2star_SE,CBF_map,CBV_all,DSC,brain_img);

    
   
end

function [r2simg,r2img] = sage_julia_fitty(base_path,te1,te2,te3,te4,te5,b_fname,DSC)
    TE1 = DSC.Parms.TE1;
    TE2 = DSC.Parms.TE2;
    TE3 = DSC.Parms.TE3;
    TE4 = DSC.Parms.TE4;
    TE5 = DSC.Parms.TE5;
    julia_cmd = '"/Users/nicks/AppData/Local/Programs/Julia 1.5.3/bin/julia.exe"';
    SAGE_Fit_julia = "C:\Users\nicks\Documents\Github/The_MRI_toolbox/SAGE_biexp_fit.jl";
    cmd = sprintf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s",julia_cmd,SAGE_Fit_julia,te1,te2,te3,te4,te5,b_fname,TE1,TE2,TE3,TE4,TE5);
    system(cmd)
    r2s = sprintf("%s/SAGE_R2s_Brain_julia.nii.gz",base_path);
    r2 = sprintf("%s/SAGE_R2_Brain_julia.nii.gz",base_path);
    r2simg = loader(r2s);
    r2img = loader(r2);
end