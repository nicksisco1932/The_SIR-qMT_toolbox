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