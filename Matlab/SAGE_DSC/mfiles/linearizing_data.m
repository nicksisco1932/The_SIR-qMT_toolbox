function [dR2star_all,dR2star_SE,CTC_all,CTC_SE] = linearizing_data(DSC,filtered_image,brain_img)

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
    
    mask4d = repmat(brain_img,[1 1 1 nt]);
    dR2star_all(isnan(dR2star_all))=0;
    dR2star_SE(isnan(dR2star_SE))=0;
    dR2star_all(isinf(dR2star_all))=0;
    dR2star_SE(isinf(dR2star_SE))=0;
    dR2star_all(dR2star_all<0)=0; 
    dR2star_SE(dR2star_SE<0)=0;
    dR2star_all = dR2star_all.*mask4d;
    dR2star_SE = dR2star_SE.*mask4d;

    
end