function [OEF,CMRO2,pO2] = OEF_CMRO2_pO2(R2s,R2,CBF,CBV,DSC,brain_img)
     %% This is what I had to do to get this to work.
    % OEF is off by a factor of 25 somewhere.
    % 8.8 mmol/mL is probably the right number and unit, but it should be
    % in umol/mL so 8.8/1000. These two changes put the maps in realistic
    % units. There is so much conflicting information on how this is
    % calculated, I'm not actually sure it is worth the hassle.
    
    % The units are all over the place from all literature. What a mess.
    
    mask4d = repmat(brain_img,[1 1 1 size(R2s,4)]);
    R2s = R2s.*mask4d;
    R2 = R2.*mask4d;
    R2prime=R2s-R2;
    R2prime=R2prime(:,:,:,DSC.Parms.pk_tp);
%     R2prime=nanmean(squeeze(R2prime(:,:,:,DSC.Parms.pk_tp)),4);
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
    gammaH = 267.522E6;     %rad/s/T
    deltaChi = 0.264e-6; % 0.264 ppm Spees et al. for RBC 
%     deltaChi = 0.18E-6;     % (ppm) cgs units https://pubmed.ncbi.nlm.nih.gov/1569876/ & https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5833262/
    Hct=0.42*0.85;          % Hematocrit volume fraction (L of RBC/L of solute)
    B0 = 3;                 % Tesla (T)
    delta_omegaS = (4*pi()/3)*gammaH*deltaChi*Hct*B0;

    OEF = R2prime./(delta_omegaS.*CBV); % CBV in mL/100g 
    OEF(find(logical(brain_img)==0))=0; 
    % OEF(OEF>100)=0;
%     for slice_num = 1:nt
%         img=permute(brain_img.*OEF(:,:,slice_num,jj),[2 1 3]);
%         imagesc(img(:,:,slice_num)); axis image; colorbar;title("OEF")
%         pause(0.5)
%     end
    
%% OEF is a factor of 25 off
%     OEF=25.*OEF.*100; % as a percentage in a lot of publication figures 
    % Is this b/c I'm converting CBV to mL/100g with 0.733*100?
    OEF=OEF.*100; % is it supposed to be a percentage?

    Ca = 8.8; % mmol/mL arterial blood oxygen content (wrongly cited) umol/mL https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4605053/  https://onlinelibrary.wiley.com/doi/10.1002/mrm.23283  
    CBF(CBF>1000)=0;

%     CMRO2 = OEF.*CBF.*Ca/1000; % units should be umol/100g/min according to  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4605053/
    CMRO2 = OEF.*CBF.*Ca/1000;
%     for jj = 1
%         img=permute(brain_img.*CMR02(:,:,slice_num,jj),[2 1 3]);
%         imagesc(img(:,:,slice_num)); axis image; colorbar;title("CMR02")
%         pause(0.05)
%     end
    p50=26; % mmHg, O2 affinity for Hg
    h=2.55; % Hill coefficient of Hg binding O2, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3949039/
    L=4.4;  % mmol/Hg/min % now it's back to mmol
%     pO2 = p50.* power(2./OEF-1,1/h)-CMRO2/L;
    pO2 = p50.* (2./(OEF)-1).^(1/h)-CMRO2/L;
    pO2(isinf(pO2))=0;
    pO2(isnan(pO2))=0;
    pO2(pO2<0)=0;
    pO2(pO2>2000)=0;
    pO2=real(pO2);
    
    warning('There is something wrong with the pO2 calculation')
%     for slice_num = 1:nt
%         img=permute(brain_img.*pO2(:,:,slice_num,jj),[2 1 3]);
%         imagesc(img(:,:,slice_num),[0 80]); axis image; colorbar;title("pO2")
%         pause(0.5)
%     end
end