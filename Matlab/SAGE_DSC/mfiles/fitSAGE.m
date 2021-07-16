function [p,res,J]= fitSAGE(echos,Yy,p0,plotFlag)

options = optimoptions('lsqnonlin',"Display","off","FunctionTolerance",1E-10,...
    "UseParallel",false,"Algorithm","levenberg-marquardt","StepTolerance",1E-10);

[p,~,res,~,~,~,J] = lsqnonlin(@optfitSIR,p0,[],[],options,echos,Yy,plotFlag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function er = optfitSIR(p,echos,M,plotFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract parameters
% S0I=p(1);
% S0II=p(2);
% R2s=p(3);
% R2=p(4);

% Calculate signal
Mm = SAGE_biexp4p(p,echos);

% Misfit
er = Mm - M;


% Plot data % save some time by not evaluating these if loops.
if plotFlag == 'y'
    plot(ti,Mm,'b-',ti,M,'rs','LineWidth',2)
    
    xlabel('t_i (s)'), ylabel('Signal (a.u.)')
    grid on; drawnow
end

function [SI_sage]=SAGE_biexp4p(p,echos)

%x should be [SI_I SI_II R2s R2] where d = SI_I/SI_II; 
tn=echos;
TE=echos(end);
% % TE = 0.09;
% % TE = 0.087;

SI_sage = zeros(size(tn));
%Create piece-wise function based on Schmiedeskamp H et al. MRM 2012
%67:378-388
SI_sage(tn<(TE/2)) = p(1).*exp(-tn(tn<(TE/2)).*p(3));
SI_sage(tn>(TE/2)) = p(2)*exp(-TE*(p(3)-p(4))).*exp(-tn(tn>(TE/2))*(2*p(4)-p(3)));