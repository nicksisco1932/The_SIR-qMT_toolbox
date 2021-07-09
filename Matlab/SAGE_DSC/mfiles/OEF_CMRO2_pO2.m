function [OEF,CMRO2,pO2] = OEF_CMRO2_pO2(R2s,R2,CBF,CBV,DSC,brain_img)
close all
    dR2 = R2;
    R2prime=R2s-dR2;
    R2prime=nanmean(squeeze(R2prime(:,:,:,DSC.Parms.pk_tp)),4);
    R2prime(find(logical(brain_img)==0))=0;
    R2prime(R2prime<0)=0;
%     for jj = 1
%         i = 9;
%         slice_num = i;
%         % subplot(1,2,1);
%         img=permute(brain_img.*R2prime(:,:,slice_num,jj),[2 1 3]);
%         clim = round([prctile(nonzeros(img),10,'all') prctile(nonzeros(img),90,'all')],2,'decimals');
%         imagesc(img(:,:,slice_num),[0 25]); axis image; colorbar;title("R2 Prime")        
% %         imagesc(img(:,:,slice_num)); axis image; colorbar;title("R2 Prime")        
%         pause(0.1)
% 
%     end
    lamda_H=(267.522E6); %rad/s/T
    deltaChi = 0.0924E-7; % Susceptibility
%     deltaChi = 0.264E-6; % Susceptibility
    Hct=0.42*0.85;       % Hematocrit
    B0=3;                % T 
    OEF = R2prime./(4/3*pi()*lamda_H*deltaChi*Hct*B0.*CBV);
    OEF(find(logical(brain_img)==0))=0;
    % OEF(OEF>100)=0;
%     for slice_num = 1:nt
%         img=permute(brain_img.*OEF(:,:,slice_num,jj),[2 1 3]);
%         imagesc(img(:,:,slice_num)); axis image; colorbar;title("OEF")
%         pause(0.5)
%     end
    Ca=8.68; % mmol/mL arterial blood oxygen content
    
    CBF(CBF>1000)=0;
    
    CMRO2 = OEF.*CBF.*Ca;
%     for jj = 1
%         img=permute(brain_img.*CMR02(:,:,slice_num,jj),[2 1 3]);
%         imagesc(img(:,:,slice_num)); axis image; colorbar;title("CMR02")
%         pause(0.05)
%     end
    p50=27; % mmHg
%     h=0.27; % Hill coefficient of Hg binding O2
    h=2.7; % Hill coefficient of Hg binding O2
    L=4.4;  % mmol/Hg/min
    pO2 = p50.* power(2./OEF-1,1/h)-CMRO2/L;
    pO2(isinf(pO2))=0;
    pO2(isnan(pO2))=0;
    pO2(pO2<0)=0;
    pO2(pO2>200)=0;
%     for slice_num = 1:nt
%         img=permute(brain_img.*pO2(:,:,slice_num,jj),[2 1 3]);
%         imagesc(img(:,:,slice_num),[0 80]); axis image; colorbar;title("pO2")
%         pause(0.5)
%     end
end