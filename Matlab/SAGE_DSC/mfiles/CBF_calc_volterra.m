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