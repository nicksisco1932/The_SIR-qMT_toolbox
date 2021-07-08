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
    
    CBV_all = brain_img.*CBV_all;
    CBV_SE = brain_img.*CBV_SE;
end