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