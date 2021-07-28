function dCs = twomodel(t,Cs,~,AIFfull,tspan,Ftm,PStm,vptm,vetm)

% % global AIFfull;
% % global Ftm;  
% % global PStm; 
% % % global vptm;
% % % global vetm;
% % global tspan;

%AIF = AIFfull(round(t+1));

AIF = interp1(tspan, AIFfull, t);

dCs = zeros(size(Cs));
dCs(1) = (Ftm/vptm)*(AIF-Cs(1))-(PStm/vptm)*(Cs(1)-Cs(2));  %Cp
dCs(2) = (PStm/vetm)*(Cs(1)-Cs(2));  %Ce



