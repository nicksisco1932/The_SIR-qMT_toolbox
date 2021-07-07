function [DSC,CBF_map,CBFSE_map,CBV_all,CBV_SE,MTT,MTT_SE] = sage_proc_ns(temp,index)

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
    [dR2star_all,dR2star_SE] = linearizing_data(DSC,filtered_image,brain_img);
    %% CBV
    %--------------------------
    [CBV_all,CBV_SE,tidx] = CBV_calc(DSC,dR2star_all,dR2star_SE,brain_img);
    %--------------------------

    [threshold] = adaptive_threshold(DSC,filtered_image,prebolus_signal);

    %% CBF

    [S_orig,maxS,U,S,V,AIFmatrixt,dtemp_all,dtemp_allSE] = aif_processing(DSC,dR2star_all,dR2star_SE,tidx);

    %--------------------------
    [CBF_map,CBFSE_map] = CBF_calc(dtemp_all,dtemp_allSE,DSC,U,S,V,maxS,S_orig,brain_img,threshold);
    %--------------------------
    MTT = CBV_all./CBF_map;
    MTT_SE = CBV_SE./CBFSE_map;
end

function [dR2star_all,dR2star_SE] = linearizing_data(DSC,filtered_image,brain_img)

    ss_tp = DSC.Parms.ss_tp;
    gd_tp = DSC.Parms.gd_tp;
    pk_tp = DSC.Parms.pk_tp;
    nt=DSC.Parms.nt;
    STE1=squeeze(filtered_image(:,:,:,1,:));
    STE2=squeeze(filtered_image(:,:,:,2,:));
    STE5=squeeze(filtered_image(:,:,:,5,:));
    STE1_pre = repmat(nanmean(squeeze(STE1(:,:,:,ss_tp:gd_tp)),4),[1 1 1 nt]);
    STE2_pre = repmat(nanmean(squeeze(STE2(:,:,:,ss_tp:gd_tp)),4),[1 1 1 nt]);
    STE5_pre = repmat(nanmean(squeeze(STE5(:,:,:,ss_tp:gd_tp)),4),[1 1 1 nt]);

    TE1=DSC.Parms.TE1;
    TE2=DSC.Parms.TE2;
    TE5=DSC.Parms.TE5;
    STE0 = squeeze(STE1.*(STE1./STE2).^(TE1/(TE2-TE1))); %from Vonken paper
    DSC.dR2s.STE0 = squeeze(filtered_image(:,:,:,1,:).*(filtered_image(:,:,:,1,:)./filtered_image(:,:,:,2,:)).^(DSC.Parms.TE1/(DSC.Parms.TE2-DSC.Parms.TE1))); %from Vonken paper
    STE0_fraction = STE0./repmat(nanmean(STE0(:,:,:,ss_tp:gd_tp),4),[1 1 1 size(STE0,4)]);

    % Calculate delta_R2star
%     prebolus_signal5d = repmat(prebolus_signal,[1 1 1 1 DSC.Parms.nt]);
%     TE_fraction = filtered_image./prebolus_signal5d;

    rtissue = 87;
    rtissueSE = 20.4;
    dR2star_all = (log((STE1./STE1_pre)./(STE2./STE2_pre))./(TE2-TE1));
    dR2star_SE = (log(STE5_pre./STE5))./TE5;
    CTC_all = dR2star_all./rtissue;
    CTC_SE = dR2star_SE./rtissueSE;
    CTC_all(isinf(CTC_all)) = 0;
    CTC_all(isnan(CTC_all)) = 0;
    CTC_SE(isinf(CTC_SE)) = 0;
    CTC_SE(isnan(CTC_SE)) = 0;
    ind=find(brain_img);
    dR2star_all(isnan(dR2star_all))=0;
    dR2star_SE(isnan(dR2star_SE))=0;
    dR2star_all(isinf(dR2star_all))=0;
    dR2star_SE(isinf(dR2star_SE))=0;
    dR2star_all(ind)=0;
    dR2star_SE(ind)=0;

    clear STE1 STE2 STE0 STE0_fraction STE5 STE1_pre STE2_pre CTC_all CTC_SE
end

function [S_orig,maxS,U,S,V,AIFmatrixt,dtemp_all,dtemp_allSE] = aif_processing(DSC,dR2star_all,dR2star_SE,tidx)
    rtissue = 87;
    rtissueSE = 20.4;
    disp('Doing Ashley Original Code')
    dtemp_all = dR2star_all(:,:,:,tidx)./rtissue;
    dtemp_allSE = dR2star_SE(:,:,:,tidx)./rtissueSE;
    zpad = 2;
    DR2sAIF = DSC.AIF(2,tidx)';
    DR2sAIF((length(DR2sAIF)+1):zpad*length(DR2sAIF)) = 0;
    if size(DR2sAIF,1)<size(DR2sAIF,2)
        warning('Incorrect - AIF row matrix!');
        DR2sAIF = DR2sAIF';
    end
    dtemp_all(:,:,:,(size(dtemp_all,4)+1):zpad*size(dtemp_all,4))  = 0;
    dtemp_allSE(:,:,:,(size(dtemp_allSE,4)+1):zpad*size(dtemp_allSE,4))  = 0;
    AIFmatrixt = DSC.Parms.TR*circulant(DR2sAIF,1); %Create circular matrix
    %Performs standard SVD on the DeltaR2star values for the AIF
    [U,S,V]=svd(AIFmatrixt);
    maxS = max(diag(S));
    S_orig = S;
end

function [threshold] = adaptive_threshold(DSC,filtered_image,prebolus_signal)
    %% Adaptive Threshold
    Smin = squeeze(nanmin(filtered_image(:,:,:,:,DSC.Parms.ss_tp:DSC.Parms.gd_tp+40),[],5)); %find min signal (including during bolus passage)
    prebolus_std = squeeze(nanstd(filtered_image(:,:,:,:,DSC.Parms.ss_tp:DSC.Parms.gd_tp),[],5)); %find prebolus standard deviation
    % % SNRc(j,k)=(S0_TE2/std0_TE2)*(Smin/S0_TE2)*log(S0_TE2/Smin);
    SNRc2 = squeeze((prebolus_signal./prebolus_std).*(Smin./prebolus_signal).*log(prebolus_signal./Smin));
    SNRc2(isnan(SNRc2)) = 0;
    C1 = SNRc2;
    C2 = SNRc2;
    C3 = SNRc2;
    C1(C1>=4) = 0;
    C1(C1<2)=0;
    threshold_C1=(100./(1.7082*SNRc2-1.0988))./100;
    C2(C2<4)=0;
    threshold_C2 = (100./(-0.0061.*(SNRc2.^2)+0.7454*SNRc2+2.8697))./100; %AMS remove >37 criteria! C2(C2>37)=0;
    C3(C1~=0)=0;
    C3(C2~=0)=0;
    threshold_C3 = 0.15;
    threshold = (threshold_C1.*logical(C1))+(threshold_C2.*logical(C2))+(threshold_C3.*logical(C3));
end

function [CBV_all,CBV_SE,tidx] = CBV_calc(DSC,dR2star_all,dR2star_SE,brain_img)
    TR = DSC.Parms.TR;
    aif = DSC.AIF(2,:);

    % AIF = aif.fit.gv; % if run by DSC toolbox

    a1=0.493;
    a2=2.62;
    htcvar=1.1378*(0.4)/(1-0.4)^2; % 3T
    for ii = 1:size(aif,2)
        F=[a2*htcvar,a1,-aif(ii)];
        C=roots(F);
        if isreal(C(1))
            CTC(ii)=-C(1); % is this always a negative curve?
        end
    end
    AIF = CTC;
    
    gd_tp = DSC.Parms.gd_tp;
    
    ind=find(brain_img);
    tidx = max(1,round(gd_tp-60/TR)):round(gd_tp+60/TR);
    rtissue = 87;
    temp = dR2star_all(:,:,:,tidx)./rtissue;
    % temp = dR2star_all(:,:,:,:)./rtissue;
    temp(temp<0) = 0;
    CBV_all = squeeze(trapz(temp,4))./trapz(AIF(tidx));
    CBV_all(~isfinite(CBV_all)) = 0;
    CBV_all(~ind) = 0;

    rtissueSE = 20.4;
    temp = dR2star_SE(:,:,:,tidx)./rtissueSE;
    % temp = dR2star_SE(:,:,:,:)./rtissueSE;
    temp(temp<0) = 0;
    CBV_SE = squeeze(trapz(temp,4))./trapz(AIF(tidx));
    CBV_SE(~isfinite(CBV_SE)) = 0;
    CBV_SE(~ind) = 0;
end

function [CBF_map,CBFSE_map] = CBF_calc(dtemp_all,dtemp_allSE,DSC,U,S,V,maxS,S_orig,brain_img,threshold)
    % % %Calculate CBF*R(t)
    CBF_map = zeros(size(dtemp_all));
    CBFSE_map = zeros(size(dtemp_allSE));
    for i = 1:DSC.Parms.nx
        for j = 1:DSC.Parms.ny
            for k = 1:DSC.Parms.nz
                if brain_img(i,j,k) % much faster with this
                    S = S_orig;
                    S(S<(threshold(i,j,k)*maxS))=0;
                    %Inverts the diagonal of the S matrix produced through SVD
                    S = 1./S;
                    S(isinf(S)) = 0;
                    %Computes the inverted DeltaR2star values for the AIF
                    invDeltaR2starAIF=V*S*U';
                    CBF_map(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_all(i,j,k,:));
                    CBFSE_map(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_allSE(i,j,k,:));
                end
            end
        end
    end

    ind=find(brain_img);
    CBF_map(~ind)=0;
    CBFSE_map(~ind)=0;
    CBF_map(isnan(CBF_map))=0;
    CBFSE_map(isnan(CBFSE_map))=0;
    CBF_map(isinf(CBF_map))=0;
    CBFSE_map(isinf(CBFSE_map))=0;


    CBF_map = max(CBF_map,[],4);
    CBF_map(~isfinite(CBF_map)) = 0;
    CBFSE_map = max(CBFSE_map,[],4);
    CBFSE_map(~isfinite(CBFSE_map)) = 0;
end


function [vol] = loader(path)
    info = niftiinfo(path);
    vol = double(niftiread(info));
end
