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



    %     columnG(nt+1)=a(nt)*TR/6;
    columnG(nt+1)=a(nt)*TR/6;
    %     for i=2:(nt-1)
    for i=2:(nTpad-1)
        a_minus = 2*a(i)+a(i-1)*TR/6;
        a_plus = 2*a(i)+a(i+1)*TR/6;
        a_plusminus = a_plus+a_minus;
        columnG(i)=a_plusminus;
    end

    %     for i=nt
    for i=nTpad
        a_minus = 2*a(i)+a(i-1)*TR/6;
        a_plus = 2*a(i)*TR/6;
        a_plusminus = a_plus+a_minus;
        columnG(i)=a_plusminus;
    end

    rowG=zeros(1,nTpad);
    rowG(1)=columnG(1);
    for j=2:nTpad
        rowG(j)=columnG(nTpad+2-j);
    end

    G=toeplitz(columnG,rowG);
    [U,S,V]=svd(G);
    maxS = max(diag(S));
    S_orig = S;

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


                    eigenV=diag(S);

                    newEigen(eigenV<adapt_thresh) = 0;

                    Ginv=V*diag(newEigen)*(U'); % the second diagonal puts the eiganvector back into a mxn matrix

                    vettConc=zeros(nTpad,1);

                    vettConc(1:nt)=squeeze(CTC_all(i,j,k,1:nt));
                    vettConc(nTpad+1:nTpad)=0;
                    vettRes=(1/TR)*Ginv*vettConc;
                    CBF_map(i,j,k)=max(abs(vettRes));
                    residual_cbf(i,j,k,:)=vettRes;

                    vettConc=zeros(nTpad,1);

                    vettConc(1:nt)=squeeze(CTC_SE(i,j,k,1:nt));
                    vettConc(nTpad+1:nTpad)=0;
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
    
end