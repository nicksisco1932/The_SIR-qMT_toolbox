function [CBF_map,CBFSE_map] = CBF_calc_volterra(DSC,brain_img,threshold)
        
    %--------------------
    DeltaT = DSC.Parms.TR; % sampling interval
    
    options.nT=DSC.Parms.nt;
    disp('Doing Volterra Method.')
    nt = DSC.Parms.nt;
    nTpad=2*nt;
    columnG=zeros(nTpad,1);
    AIF=DSC.AIF(2,:)';
    a = zeros(nTpad,1);
    a(1:nt) = AIF; % For continuity with literature and formal matrix notation
    G=zeros(nTpad);
    for ii = 2:nt-1
        for jj = ii:nt-1
            t=ii;
            G(jj,ii) = (2*a(t)*DeltaT/6)+(2*a(t)*DeltaT/6); % lower triangle matrix
        end
    end
    for ii = 2:nt-1
        jj = ii;
        n=ii-1;
        G(ii,jj) = diag(2*a(n)-a(n)*DeltaT/6); % a-  diag (aii)
    end
    G(1:nt-1,1) = 2*a(1:nt-1)-a(1:nt-1)*DeltaT/6; %a- column 1 (ai0)
    
    [U,S,V]=svd(G);
    maxS = max(diag(S));
    S_orig = S;
    
    %%%Calculate CBF*R(t)
    % CBFR_all = zeros(size(dtemp_all)); % don't need these, the max is taken
    % in the loop to eliminate a step
    % CBFR_allSE = zeros(size(dtemp_allSE));
    CBF_map=zeros(DSC.Parms.nx,DSC.Parms.ny,DSC.Parms.nz);
    CBFSE_map=zeros(DSC.Parms.nx,DSC.Parms.ny,DSC.Parms.nz);
    nTpad=size(G,1); % renamed here
    residual_cbf=zeros(nx,ny,nz,nTpad);
    residual_cbfse=zeros(nx,ny,nz,nTpad);
    mask=brain_img; % Adding the brain index here is much much faster
    for i = 1:DSC.Parms.nx
        for j = 1:DSC.Parms.ny
            for k = 1:DSC.Parms.nz
                if mask(i,j,k) % don't do anything outside the mask
                    % Apply adaptive threshold
                    S = S_orig;
                    adapt_thresh = threshold(i,j,k)*maxS;
                    %Computes the inverted DeltaR2star values for the AIF
                    eigenV=diag(S);
                    eigenV(eigenV<adapt_thresh)=0;
                    %                 newEigen=eigenV;
                    newEigen=eigenV.^-1;
                    newEigen(isinf(newEigen))=0;
                    Ginv=V*diag(newEigen)*(U'); % the second diagonal puts the eiganvector back into a mxn matrix
                    vettConc=zeros(nTpad,1);
                    vettConc(1:nt)=reshape(CTC_all(i,j,k,:),nt,1);
                    vettRes=(1/TR)*Ginv*vettConc;
                    CBF_map(i,j,k)=max(abs(vettRes));
                    residual_cbf(i,j,k,:)=vettRes;
                    vettConc=zeros(nTpad,1);
                    vettConc(1:nt)=reshape(CTC_SE(i,j,k,:),nt,1);
                    vettRes=(1/TR)*Ginv*vettConc;
                    CBFSE_map(i,j,k)=max(abs(vettRes));
                    residual_cbfse(i,j,k,:)=vettRes;
                end
            end
        end
    end
    
    CBF_map(CBF_map>1000)=0;
    CBFSE_map(CBFSE_map>1000)=0;
    
    CBF_map(isnan(CBF_map))=0;
    CBF_map(isinf(CBF_map))=0;
    CBFSE_map(~ind)=0;
    CBFSE_map(isnan(CBFSE_map))=0;
    CBFSE_map(isinf(CBFSE_map))=0;
    
    %         CBF_map = removeOutlierVolume(CBF_map,3,[3 3 3]);
    %         CBFSE_map = removeOutlierVolume(CBFSE_map,3,[3 3 3]);
    
    CBF_map(CBF_map>100)=0;
    CBFSE_map(CBFSE_map>100)=0;
    
    
end