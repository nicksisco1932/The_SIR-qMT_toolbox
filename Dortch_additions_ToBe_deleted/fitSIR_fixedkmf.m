function [p,res,J]= fitSIR_fixedkmf(ti,td,M,Sm,p0,lb,ub,magFlag,plotFlag,kmf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [p,res]= fitSIR_fixedkmf(ti,td,M,Sm,p0,lb,ub,magFlag,plotFlag,kmf)
%
% Fits signal for SIR-TSE sequences two a 2 pool model of MT (with fixed kmf);  
% f = free water, m = macromolecular.
% 
% Input:
% ti - array of inversion times (N x 1)
% td - array of delay times (N x 1)
% M  - measured signal (N x 1)
% Sm - numerically estimaged saturation fraction of m pool due to inv pulse
% p0 - initial parameter guess (see p for description)
% lb,ub - lower and upper bounds for p
% magFlag - magnitude data? ('y' or 'n')
% plotFlag - plot data? ('y' or 'n')
% kmf - fixed exchange rate
% 
% Output:
% p - fitted parameters 
%      p(1) = pmf
%      p(2) = R1f = R1m 
%      p(3) = Sf
%      p(4) = M0f
% ci - corresponding 95% CIs
% res - residuals
% J = Jacobian of FUN at p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Make sure ti/td/M are column vectors
if size(ti,1) == 1
    ti = ti';
end
if size(td,1) == 1
    td = td';
end
if size(M,1) == 1
    M = M';
end


% Perform nonlinear least square fit
%options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
options = optimset('Display','off','Algorithm','levenberg-marquardt');


[p,~,res,~,~,~,J] = lsqnonlin(@optfitSIR,p0,lb,ub,options,ti,td,M,Sm,magFlag,plotFlag,kmf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function er = optfitSIR(p,ti,td,M,Sm,magFlag,plotFlag,kmf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract parameters
M0 = [1 p(1)]*p(4);
R1 = [p(2) p(2)];
S  = [p(3) Sm];

% Calculate signal
Mm = signalSIR(ti,td,kmf,R1,M0,S);
if magFlag == 'y'
    Mm = abs(Mm);
end

% Misfit
er = Mm - M;

% Plot data
if plotFlag == 'y'
    plot(ti,Mm,'b-',ti,M,'rs','LineWidth',2)
    xlabel('t_i (s)'), ylabel('Signal (a.u.)')
    grid on, drawnow
end


