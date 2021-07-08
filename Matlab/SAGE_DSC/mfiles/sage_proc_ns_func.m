function [DSC,CBF_map,CBFSE_map,CBV_all,CBV_SE,MTT,MTT_SE] = sage_proc_ns_func(temp,index)

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
    %% SAGE DSC processing
    DSC.Parms.TR = 1800./1000;                          %TR in s
    DSC.Parms.time = DSC.Parms.TR:DSC.Parms.TR:DSC.Parms.TR*size(DSC.Data,5);
    DSC.Parms.TE1 = 8.8./1000;                          %TE in s
    DSC.Parms.TE2 = 26./1000;                           %TE in s
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
    volterra = 0
    if volterra
        [dR2star_all,dR2star_SE,CTC_all,CTC_SE] = linearizing_data(DSC,filtered_image,brain_img);
    else
        [dR2star_all,dR2star_SE,~,~] = linearizing_data(DSC,filtered_image,brain_img);
    end
    %% CBV
    %--------------------------
    [CBV_all,CBV_SE,tidx] = CBV_calc(DSC,dR2star_all,dR2star_SE,brain_img);
    %--------------------------

    [threshold] = adaptive_threshold(DSC,filtered_image,prebolus_signal);

    %% CBF

    if volterra
        [S_orig,maxS,U,S,V,AIFmatrixt,dtemp_all,dtemp_allSE] = aif_processing(DSC,dR2star_all,dR2star_SE,tidx);
    else
        [CBF_map,CBFSE_map] = CBF_calc_volterra(DSC,brain_img,threshold)
    end

    %--------------------------
    [CBF_map,CBFSE_map] = CBF_calc(dtemp_all,dtemp_allSE,DSC,U,S,V,maxS,S_orig,brain_img,threshold);
    %--------------------------
    MTT = CBV_all./CBF_map;
    MTT_SE = CBV_SE./CBFSE_map;
end