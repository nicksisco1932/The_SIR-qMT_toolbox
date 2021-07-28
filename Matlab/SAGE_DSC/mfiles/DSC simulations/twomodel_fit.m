function [S] = twomodel_fit(B,tspan,AIFfull,~)

x0 = [0 0];

[t,Cs]=ode45(@twomodel,tspan,x0,[],[],AIFfull,tspan);

    function dCs = twomodel(t,Cs,~,AIFfull,tspan)

    AIF = interp1(tspan, AIFfull, t);

    dCs = zeros(size(Cs));
    dCs(1) = (B(1)/B(3))*(AIF-Cs(1))-(B(2)/B(3))*(Cs(1)-Cs(2));  %Cp
    dCs(2) = (B(2)/B(4))*(Cs(1)-Cs(2));  %Ce
    
    end

% % S = B(3)/100*Cs(:,1) + B(4)/100*Cs(:,2); %if fitting Ct curves
S = Cs; %if fitting Cp and Ce curves

end