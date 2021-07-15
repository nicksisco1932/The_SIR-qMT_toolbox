function [M,J] = signalSIR(ti,td,kmf,R1,M0,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function M = signalSIR(ti,td,kmf,M0,S)
%
% Calculates signal for SIR-TSE sequences assuming 2 pools undergoing MT;
% f = free water, m = macromolecular.
%
% Input:
% ti - array of inversion times (N x 1)
% td - array of delay times (N x 1)
% kmf - exchange rate from m to f
% M0 - M0 for each pool [M0f Mom]'
% R1 - R1 for each pool [R1f R1m]'
% S -  sat/inv effeciency for each pool [Sf Sm]'
%
% Output:
% M - z-mag [Mf Mm]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No longer using my implentation based upon the notation of:
% Dortch et al. NeuroImage 64 (2013) 640?649
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Make sure vectors is the correct size
% if size(M0,1) == 1
%     M0 = M0';
% end
% if size(S,2) == 1
%     S = S';
% end
%
% % Setup exchange matrix
% kfm = kmf*M0(2)/M0(1);
% K = [-kfm kmf ; kfm -kmf];
%
% % Calculate L1
% L1 = -diag(R1) + K;
%
% % Eigen-expansion of L1
% [U1,S1] = eig(L1);
% V1 = eye(size(U1))/U1;
%
% % Diagonalize S
% S = diag(S);
%
% % Set identitiy matrix
% I = eye(2);
%
% % Calculate signal
% Ei = zeros(size(V1));
% Ed = zeros(size(V1));
% M  = zeros(length(ti),2);
% for ii = 1:length(ti)
%     Ei =  U1*diag(exp(diag(S1*ti(ii))))*V1;
%     Ed =  U1*diag(exp(diag(S1*td(ii))))*V1;
%     M(ii,:) = (Ei*S*(I - Ed) + (I - Ei))*M0;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ke's implementation is faster
% Notation from Li et al. MRM 64:491-500 (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (ti,td,kmf,R1,M0,S);
% M0 = [1 p(1)]*p(4);
% R1 = [p(2) p(2)];
% S  = [p(3) Sm];
% 
% % Calculate signal
% Mm = signalSIR(ti,td,kmf,R1,M0,S);

% Get pmf and Mfinf
pmf = M0(2)/M0(1);
Mfinf = M0(1);

% Get R1f/R1m and Sf/Sm
R1f = R1(1);  
R1m = R1(2);
Sf = S(1);    
Sm = S(2);

% Calculate R1+/- in Eq. 4
R1diff  = sqrt((R1f - R1m + (pmf - 1) * kmf)^ 2 + 4 * pmf * kmf^2);
R1plus  = (R1f + R1m + (1 + pmf) * kmf + R1diff) / 2;
R1minus = R1plus - R1diff;

% Component amplitude terms for td terms (Eq. 5)
bftdplus = -(R1f - R1minus) / R1diff;
bftdminus = (R1f - R1plus) / R1diff;
bmtdplus  = -(R1m - R1minus) / R1diff;
bmtdminus = (R1m - R1plus) / R1diff;

% Signal recovery during td (Eq. 5)
Mftd = bftdplus * exp(-R1plus * td) + bftdminus * exp(-R1minus * td) + 1 ;
Mmtd = bmtdplus * exp(-R1plus * td) + bmtdminus * exp(-R1minus * td) + 1 ;

% Component amplitude terms for ti terms (Eq. 5)
bfplus= ((Sf * Mftd -1) * (R1f - R1minus)  +  (Sf * Mftd - Sm * Mmtd) * pmf * kmf) / R1diff;
bfminus= -((Sf * Mftd -1) * (R1f - R1plus) + (Sf * Mftd - Sm * Mmtd) * pmf * kmf) / R1diff;

% Signal equation (Eq. 3)
M = (bfplus .* exp (-R1plus * ti) + bfminus .* exp (-R1minus * ti) + 1) * Mfinf;


