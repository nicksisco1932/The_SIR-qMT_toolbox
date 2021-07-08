function [CBF_map,CBFSE_map] = CBF_calc_volterra(DSC,brain_img,threshold,CTC_all,CTC_SE)
    disp('Doing Standard cSVD AKA: Block Circulant, depending where it is cited.')

    % Block Circulant
    % The code to make this block circulant comes from DSC toolbox in a GitHub repository.

    ind=find(brain_img);
    nx=DSC.Parms.nx;
    ny=DSC.Parms.ny;
    nz=DSC.Parms.nz;
    ne=DSC.Parms.ne;
    nt=DSC.Parms.nt;
    TR = DSC.Parms.TR;
    nt = DSC.Parms.nt;
    options.nT=DSC.Parms.nt;
    nTpad=2*nt;
    %     nTpad=2*length(DSC.AIF(2,tidx)');
    columnG=zeros(nTpad,1);
    AIF = zeros(nTpad,1);
    AIF(1:nt)=DSC.AIF(2,:)';
    %     a = aif; % For continuity with literature and formal matrix notation
    a = AIF(1:nTpad); % For continuity with literature and formal matrix notation
    columnG(1)=a(1);
    columnG(nt)=DSC.Parms.TR/6*(AIF(nt-1)+4*AIF(nt)); % Equation 6 10.1002/mrm.10522, Wu, 2003
    
%--------------------------------------------------------------------------    
    % 1. Created G matrix
    aif = DSC.AIF(2,:)';
    nTpad=2*nt;
    columnG=zeros(nTpad,1);
    columnG(1)=aif(1);
    columnG(options.nT)=(aif(options.nT-1)+4*aif(options.nT))/6;
    columnG(options.nT+1)=aif(options.nT)/6;
    for k=2:(options.nT-1)
        columnG(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
    end
    rowG=zeros(1,nTpad);
    rowG(1)=columnG(1);
    for k=2:nTpad
        rowG(k)=columnG(nTpad+2-k);
    end
    
    G=toeplitz(columnG,rowG);
%--------------------------------------------------------------------------
    % 2. Singular value deconvolution
    [U,S,V]=svd(G);
    S_orig=S;
    maxS = max(diag(S));
%     eigenV=diag(S);
%     threshold=options.deconv.cSVD.threshold*max(eigenV);
  
    
%     newEigen=zeros(size(eigenV));
%     for k=1:length(eigenV);
%         if eigenV(k)>=threshold;
%             newEigen(k)=1/eigenV(k);
%         end
%     end
% 
%     Ginv=V*diag(newEigen)*(U');
%--------------------------------------------------------------------------
    % 3. Calculate CBF
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
    
%--------------------------------------------------------------------------
    

    CBF_map(CBF_map>1000)=0;
    CBFSE_map(CBFSE_map>1000)=0;

    CBF_map(isnan(CBF_map))=0;
    CBF_map(isinf(CBF_map))=0;
    
    CBFSE_map(isnan(CBFSE_map))=0;
    CBFSE_map(isinf(CBFSE_map))=0;
    
end