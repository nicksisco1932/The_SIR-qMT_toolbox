%Clears screen and all variables in the workspace.
clearvars
close all
clc

idj = 1; idk = 1;
%% 
% % clear all;
tic
% % global AIFfull;
% % global Ftm;  
% % global PStm; 
% % global vptm;
% % global vetm;
% % global tspan;

% % Ktrans_all = zeros(30,21);
% % Vp_all = zeros(30,21);
% % CBV_meas = zeros(30,21);
% % CBV_true = zeros(30,21);
% % CBV_corr = zeros(30,21);
% % 
% % CBVmeas_perror = zeros(30,21);
% % CBVcorr_perror = zeros(30,21);

Cs0=[0 0];

%Npts  = Number of seconds the simulation should run.
Npts = 300;
dT = 1; %delta T

%t = array for time points in study
%The first 60 time points are zero because the imaging is being carried out
%60 seconds prior to the first injection. 
t = 0:dT:Npts;
%Ns = (Npts+dT)/dT - 1; %total number of time pts/samples
Ns = size(t,2); %total number of time pts/samples
tspan = 0:dT:Npts;

%%%%%NEW AIF BASED ON MGH DATA - JULY 2012
to = 60; %actual injection time
toS = to/dT; %to index (in Ns)
Cp3(1:Ns)=0.0;
Cp4(1:Ns)=0.0;
r = 3.95;%orig = 3.5
b = 2.95;
Co=1;
for y = (toS+1):1:Ns, Cp3(y) = Co*((tspan(y)-to)^r)*exp(-(tspan(y)-to)/b);end
tor = 2;
Rc(1:Ns)=0.0;
taur = 625;
for y = 1:1:Ns, Rc(y) = (1/taur)*exp(-(tspan(y)-tor)/taur);end;
CpR = dT*conv(Cp3,Rc);
Cp4(toS:Ns) = Cp3(toS:Ns)+8*CpR(toS-tor:Ns-tor);
AIFt = Cp4/(max(Cp4)/5.1845);  %Scales the AIF so that the maximum concentration in a major artery (e.g. carotid) is 6 mM (REF multiple sources)
%This section adds dispersion to the AIF (probably a calamante REF was used
%to do this) It was adjusted until the rise/fall time of signals created
%for conditions found in gray matter matched clinical data
AIFt2 = [AIFt(1:Ns)];
beta = 1/2.0;
HP = beta*exp(-tspan*beta);
ConvH = dT*conv(AIFt2,HP);
AIFfull = ConvH(1:Ns);
AIFfull = AIFfull/(max(AIFfull)/12.6845);%5.6845


% % for idk = 1:21
% % for idj = 1:30
% % clearvars -except Ftm PSvptm PSvetm PStm AIFfull Cs0 Ns tspan idj idk CBVmeas* CBVcorr* Ktrans_all Vp_all CBV_*


%% -------------------PHYSIOLOGICAL PARAMETER DEFINITION---------------------
%This section defines the physiological parameters. Only one physiological
%parameter's graphing is scripted, thus only one parameter is given a range
%of values.


%Cerebral Blood Flow
%Cerebral Blood Volume
%Permeability, Surface Area Product
%Fractional volume of CA in the EES


% % CBFmat = [1.5*60]; %in C6 tumor
% % CBVmat = [1.4*4]; %in C6 tumor
% % PSmat = [0, 0.03/1.04/(1.4*4/100)]; %PS/Vp in C6  = Ktrans/rho (from Tofts PSrho = Ktrans)
% % Vemat = [0.28]; %ve in C6 tumors
% % idx = 3;

% % CBFmat = [1.2*60]; %in 9L tumor
% % CBVmat = [4*1.6]; %in 9L tumor
% % PSmat = [0, 0.046/1.04/(1.4*4/100)]; %PS/Vp in 9L = Ktrans/rho (from Tofts PSrho = Ktrans)
% % Vemat = [0.23]; %in 9L tumors
% % idx = 3;
Hct_small = 0.25; 

CBFmat = [110]; %in avg tumor 1.35*60
% % if idj == 1; CBVmat = [0.5]; %in C6 tumor
% % else CBVmat = [100*Vp_all(idj-1,idk)/(1-Hct_small)/1.04+0.5]; %in C6 tumor
% % end
CBVmat = 6  %if not using loop for CBV and Ktrans
% % if idk == 1; Ktrans = 0; %in C6 tumor
% % else Ktrans = [Ktrans_all(idj,idk-1)+0.01]; %in C6 tumor
% % end
Ktrans = [0 0.1 0.25 0.45];  %if not using loop for CBV and Ktrans
PSmat = [100*Ktrans/1.04]; %PS/Vp in avg  = Ktrans/rho (from Tofts PSrho = Ktrans)
Vemat = [0.25]; %ve in avg tumors   0.25    

Aincrementer = 0;
for Findex = 1:length(CBFmat);
    for Vpind = 1:length(CBVmat);
            for Veind = 1:length(Vemat);
                for PSindex = 1:length(PSmat);
                    Aincrementer = Aincrementer + 1;

                    PS = PSmat(PSindex)/60; %unit conversion 

                    CBV = CBVmat(Vpind); %unit conversion

                    CBF = CBFmat(Findex)/60;

                    Ve = Vemat(Veind);

                    BFval = num2str(CBF*60);
                    BVval = num2str(CBV);
                    PSval = num2str(PS*60);
                    Veval = num2str(Ve);

                    %Vp = Volume fraction of the plasma
                    Vp = CBV/100 * 1.04;

                    %Vi = Volume fraction of the extravascular intracellular space
                    Vi= 1 - Ve - Vp;

                    Ftm = CBF*(1-Hct_small)*1.04;

                    %KPSVe = PS/Vp * (vp / ve) = KPS / Ve (like in Brix paper) - here PS is
                    %PS/VP
                    % % KPSVe = PS*(Vp/Ve);
                    % % PStm = KPSVe;
                    % % PSvptm = PS;
                    vptm = Vp*100;
                    vetm = Ve*100;
                    PStm = PS;

                    [t,Cs]=ode45('twomodel',tspan,Cs0,[],AIFfull,tspan,Ftm,PStm,vptm,vetm);

                    %VeCe
                    Ce = Cs(:,2);
                    VeCe = Ce*Ve;

                    %beta0 = [0.06/60, 0.3];
                    %[beta] = nlinfit(tspan,VeCe',@Toftsmodel,beta0);
                    %beta = lsqcurvefit(@Toftsmodel,beta0,tspan,VeCe',[0 0], [0.01 1.0]);

                    %VpCp
                    Cp = Cs(:,1);
                    VpCp=Cp*Vp;


                    A(Aincrementer,1:Ns) = VeCe;
                    C(Aincrementer,1:Ns) = VpCp;

                end;
            end;
    end;
end;

r2sp = 90;
r2se = 30;
DR2_meas = r2sp*C + r2se*A;DR2_meas = DR2_meas';


Eq1 = r2sp*Vp*((Ve*abs(C/Vp-A/Ve)) + (Vi*C/Vp)) + r2se*Vi*A;Eq1 = Eq1';
DR2s_vasc = r2sp*C;
DR2_corr = Eq1 - (r2sp*Vp + r2se*Ve)*(C+A)';
DR2s_true = r2sp*Vp*(C/Vp+A/Ve);


% % CBV_true(idj,idk) = trapz(DR2s_true(1:end))/trapz(AIFfull(1:end))
% % CBV_meas(idj,idk) = trapz(DR2s_meas(1:120))/trapz(AIFfull(1:120))
% % CBV_vasc(idj,idk) = trapz(DR2s_vasc(1:120))/trapz(AIFfull(1:120))
% % CBV_corr(idj,idk) = trapz(DR2s_corr(1:120))/trapz(AIFfull(1:120))

% % CBVmeas_perror(idj,idk) = abs(CBV_true(idj,idk) - CBV_meas(idj,idk))./CBV_true(idj,idk);
% % CBVcorr_perror(idj,idk) = abs(CBV_true(idj,idk) - CBV_corr(idj,idk))./CBV_true(idj,idk);

% % end
% % idk
toc
% % end

% % figure;plot(1:601,Eq1(:,1),'k',1:601,Eq1(:,2),1:601,Eq1(:,3),1:601,DR2_corr(:,2),'--',1:601,DR2_corr(:,3),'--','linewidth',2)

%% fit Brix model - input from Brix model -- proof of principle


C_aif = AIFfull;

B0 = [85/60*(1-Hct_small)*1.04 0.1*100/1.04/60 4.18 38]
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',10000,'Display','off');
[B] = lsqcurvefit('twomodel_fit',B0,tspan',[C(2,:)'/Vp A(2,:)'/Ve],[0 0 0 0],[150/60*(1-Hct_small)*1.04 1*100/1.04/60 20 80],options,C_aif');

disp(B)

FitData = twomodel_fit(B,tspan',C_aif'); 

Cp = FitData(:,1);
Ce = FitData(:,2); 

figure;plot(tspan,C(2,:)'/Vp,tspan, A(2,:)'/Ve)
figure;plot(tspan,Cp,tspan,Ce)

FitData = twomodel_fit([Ftm PStm vptm vetm],tspan',C_aif'); 

Cp = FitData(:,1);
Ce = FitData(:,2); 

figure;plot(tspan,C(2,:)'/Vp,tspan, A(2,:)'/Ve)
figure;plot(tspan,Cp,tspan,Ce)

%% fit Brix model - input average DeltaR2* curves in tumor and normal tissue, plus AIF

% % clearvars -except tumorROI normalROI aif
% % 
% % C_aif = aif./60;
% % C_tumor = tumorROI./87;
% % C_normal = normalROI./87;

Hct_small = 0.25;
tspan = 1:600;

B0 = [85/60*(1-Hct_small)*1.04 0.03*100/1.04/60 4.18 28]
% % B0 = [75/60*(1-Hct_small)*1.04 0.08*100/1.04/60 6 34]
options = optimset('TolFun',1e-15,'Tolx',1e-15,'MaxIter',10000,'Display','off');

[B,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@twomodel_fit,B0,tspan',C_tumor',[0 0 0 0],[150/60*(1-Hct_small)*1.04 1*100/1.04/60 20 80],options,C_aif);

disp(B)

FitData = twomodel_fit(B,tspan',C_aif); 

figure;plot(tspan,C_tumor(1,:),tspan,FitData,':','linewidth',2);set(gca,'LineWidth',1,'FontSize',18);legend('Input C_t', 'Fit C_t')


%%
FitData = twomodel_fit([90/60*(1-Hct_small)*1.04 0.11*100/1.04/60 9 45],tspan',C_aif); 

figure;plot(tspan,C_tumor(1,:),tspan,FitData,':','linewidth',2);set(gca,'LineWidth',1,'FontSize',18);legend('Input C_t', 'Fit C_t',2)


FitData = twomodel_fit([Ftm PStm vptm vetm],tspan',AIFfull(1:600)'); 

figure;plot(tspan,C_tumor(1,:),tspan,FitData,':','linewidth',2);set(gca,'LineWidth',1,'FontSize',18);legend('Input C_t', 'Fit C_t',2)


%% BRIX
clear H* temp
CBFmat = [90];
% % CBFmat = [5 30 60 90 120 150]; 

% % CBVmat = [0.5 3 5.5 8 10.5 13];
CBVmat = [5.5];

Ktrans = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5];
% % Ktrans = [0 0.25];
% % Ktrans = [0 0.1 0.4 0.8];
PSmat = [100*Ktrans/1.04]; %PS = Ktrans/rho (from Tofts PSrho = Ktrans)

Vemat = [0.45];
% % Vemat = [0 0.1 0.2 0.3 0.45 0.6];   
Hct_small = 0.25; 

Incrementor = 0;
for PSindex = 1:length(PSmat);
    for Vpind = 1:length(CBVmat);
            for Veind = 1:length(Vemat);
                for Findex = 1:length(CBFmat);
                    Incrementor = Incrementor + 1;
                    PS = PSmat(PSindex)/60;
    
                    CBV = CBVmat(Vpind); %unit conversion

                    CBF = CBFmat(Findex)/60;

                    Ve = Vemat(Veind);

                    %Vp = Volume fraction of the plasma
                    Vp = CBV/100 * 1.04;

                    Ftm = CBF*(1-Hct_small)*1.04;

                    vptm = Vp*100;
                    vetm = Ve*100;
    
                    T_P = vptm./(PS+Ftm);
                    T_E = vetm./PS;
                    T_B = vptm./Ftm;
                    T = (vptm + vetm)./Ftm;
                    alpha(Incrementor) = 1./(T_B*T_E);
                    beta(Incrementor) = (1./T_B) + (T)./(T_B*T_E);
                    gamma(Incrementor) = (Ftm*T)./(T_B*T_E);

                    K_plus = 0.5*(1/T_P + 1/T_E + sqrt(((1/T_P+1/T_E)^2)-4*1/T_E*1/T_B));KPVal(Incrementor) = K_plus;
                    K_minus = 0.5*(1/T_P + 1/T_E - sqrt(((1/T_P+1/T_E)^2)-4*1/T_E*1/T_B));KMVal(Incrementor) = K_minus;

                    E_minus = (K_plus - 1/T_B)./(K_plus - K_minus);EMVal(Incrementor) = E_minus;

                    H_p = (1-E_minus)*K_plus*exp(-tspan.*K_plus) + E_minus*K_minus*exp(-tspan.*K_minus);%
                    H_e = (exp(-tspan.*K_minus) - exp(-tspan.*K_plus))./((1/K_minus) - (1/K_plus));

                    H_p_all(Incrementor,:) = H_p;
                    H_e_all(Incrementor,:) = H_e;

% % H_p2 = ((1-E_minus)*K_plus*exp(-tspan.*K_plus) + E_minus*K_minus*exp(-tspan.*K_minus)) .* (exp(tspan.*0.0427));%
% % figure; plot(tspan(1:20),H_p_nL(1:20),tspan(1:20),H_p_wL(1:20),tspan(1:20),H_p2(1:20),'--','linew',2)
%%
% % H_pAIF = conv(H_p,AIFfull);
% % H_pAIF = H_pAIF(1:size(tspan,2));
% %     
H_pAIF2(1) = 0;
clear expo* crp* t2 int_t dummy_t
for k = 2:size(tspan,2)
int_t = tspan(k);
for j = 1:k
dummy_t = tspan(j);
expo1(j) = exp(-(int_t-dummy_t).*K_plus);
crpexp1(j) = AIFfull(j)*expo1(j);
expo2(j) = exp(-(int_t-dummy_t).*K_minus);
crpexp2(j) = AIFfull(j)*expo2(j);

end
t2 = tspan(1:k);
crpexp_integral1(k,Incrementor) = trapz(t2,crpexp1);crpexp_integral2(k,Incrementor) = trapz(t2,crpexp2);
H_pAIF2(k) = (1-E_minus)*K_plus*crpexp_integral1(k,Incrementor) + E_minus*K_minus*crpexp_integral2(k,Incrementor);
end
tempP(:,Incrementor,1) = crpexp_integral1(:,Incrementor);tempP(:,Incrementor,2) = crpexp_integral2(:,Incrementor);
H_pAIF_all(Incrementor,:) = H_pAIF2;

% % figure;plot(tspan,H_pAIF2,tspan,H_pAIF,'--','linew',2)
% % figure;plot(tspan,H_pAIF2,tspan,C(2,:)./Vp,'linew',2)

LC_CpAIF2(1) = 0;
clear expo* crp* t2 int_t dummy_t
for k = 2:size(tspan,2)
int_t = tspan(k);
for j = 1:k
dummy_t = tspan(j);
expo1(j) =exp(-(int_t-dummy_t).*Ftm./vptm);
crpexp1(j) = AIFfull(j)*expo1(j);
end
t2 = tspan(1:k);
crpexp_integral1 = trapz(t2,crpexp1);
LC_CpAIF2(k) = crpexp_integral1;
end

%%



H_eAIF2(1) = 0;
clear crp*
for k = 2:size(tspan,2)
int_t = tspan(k);
for j = 1:k
dummy_t = tspan(j);
expo1(j) =exp(-(int_t-dummy_t).*K_minus);
crpexp1(j) = AIFfull(j)*expo1(j);
expo2(j) = exp(-(int_t-dummy_t).*K_plus);
crpexp2(j) = AIFfull(j)*expo2(j);
end
t2 = tspan(1:k);
crpexp_integral1(k,Incrementor) = trapz(t2,crpexp1);crpexp_integral2(k,Incrementor) = trapz(t2,crpexp2);
H_eAIF2(k) = (crpexp_integral1(k,Incrementor) - crpexp_integral2(k,Incrementor))./((1/K_minus) - (1/K_plus));   %extended tofts model from Yankeelov MRM 2003
end
tempsE(:,Incrementor,1) = crpexp_integral1(:,Incrementor);tempsE(:,Incrementor,2) = -crpexp_integral2(:,Incrementor);
H_eAIF_all(Incrementor,:) = H_eAIF2;



H_eAIF = conv(H_e,AIFfull);
H_eAIF = H_eAIF(1:size(tspan,2));

% % figure;plot(tspan,H_eAIF2,tspan,H_eAIF,'--',tspan,A(2,:)./Ve,':','linew',2)

%%

%%

%% difference in the propagators for H_p_nL and H_p_wL

LC_CpAIF(1) = 0;
clear expo* crp* t2 int_t dummy_t
for k = 2:size(tspan,2)
int_t = tspan(k);
for j = 1:k
dummy_t = tspan(j);   
expo1(j) = exp(-(int_t-dummy_t)*((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm)))...
    *(((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 - Ftm/vptm + PS/(2*vetm))/(((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2) - 1)...
    *((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm));
crpexp1(j) = AIFfull(j)*expo1(j);
expo2(j) = (Ftm*exp(-(Ftm*(int_t-dummy_t))/vptm))/vptm;
crpexp2(j) = AIFfull(j)*expo2(j);
expo3(j) = (exp(-(int_t-dummy_t)*((Ftm + PS)/(2*vptm) - (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm)))*((Ftm + PS)/(2*vptm) - (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm))*((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 - Ftm/vptm + PS/(2*vetm)))...
    /(((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2);
crpexp3(j) = AIFfull(j)*expo3(j);
end
t2 = tspan(1:k);
crpexp_integral1(k,Incrementor) = trapz(t2,crpexp1);crpexp_integral2(k,Incrementor) = trapz(t2,crpexp2);crpexp_integral3(k,Incrementor) = trapz(t2,crpexp3);
LC_CpAIF(k) = crpexp_integral1(k,Incrementor) + crpexp_integral2(k,Incrementor) - crpexp_integral3(k,Incrementor); 
end
temps(:,Incrementor,1) = crpexp_integral1(:,Incrementor);temps(:,Incrementor,2) = crpexp_integral2(:,Incrementor);temps(:,Incrementor,3) = crpexp_integral3(:,Incrementor);
H_LCAIF_all(Incrementor,:) = LC_CpAIF;

                end
            end
    end
end

% % diff = exp(-t*((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve)))...
% %     *(((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 - F/Vp + PS/(2*Ve))/(((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2) - 1)...
% %     *((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve))...
% %     + (F*exp(-(F*t)/Vp))/Vp...
% %     - (exp(-t*((F + PS)/(2*Vp) - (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve)))*((F + PS)/(2*Vp) - (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve))*((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 - F/Vp + PS/(2*Ve)))...
% %     /(((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2);
C_t = H_pAIF_all.*vptm+H_eAIF_all.*vetm;

figure;plot(tspan,H_pAIF_all,tspan,H_eAIF_all)
figure;plot(tspan,Ftm.*cumsum(AIFfull,2) - alpha(1)*cumsum(cumsum(C_t(1,:),2),2)  - beta(1).*cumsum(C_t(1,:),2) + gamma(1).*cumsum(cumsum(AIFfull,2),2),tspan,Ftm.*cumsum(AIFfull,2) - alpha(9)*cumsum(cumsum(C_t(9,:),2),2) - beta(9)*cumsum(C_t(9,:),2) + gamma(9).*cumsum(cumsum(AIFfull,2),2),'--',tspan,Ftm.*cumsum(AIFfull,2)  - alpha(9)*cumsum(cumsum(C_t(9,:),2),2) - beta(9)*cumsum(C_t(9,:),2)  + gamma(9).*cumsum(cumsum(AIFfull,2),2) - (-alpha(9)*cumsum(cumsum(C_t(9,:),2),2) - beta(9)*cumsum(C_t(9,:),2) + gamma(9).*cumsum(cumsum(AIFfull,2),2) + Ftm.*Ftm./vptm.*cumsum(LC_CpAIF2(1,:),2)),'--','linew',2)

%% TM and ETM
clear H* temp
CBFmat = [90];
% % CBFmat = [5 30 60 90 120 150]; 

% % CBVmat = [0.5 3 5.5 8 10.5 13];
CBVmat = [5.5];

Ktrans = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5];
% % Ktrans = [0 0.25];
% % Ktrans = [0 0.1 0.4 0.8];
PSmat = [100*Ktrans/1.04]; %PS = Ktrans/rho (from Tofts PSrho = Ktrans)

Vemat = [0.45]
% % Vemat = [0 0.1 0.2 0.3 0.45 0.6];   
Hct_small = 0.25; 

Incrementor = 0;
for PSindex = 1:length(PSmat);
    for Vpind = 1:length(CBVmat);
            for Veind = 1:length(Vemat);
                for Findex = 1:length(CBFmat);
                    Incrementor = Incrementor + 1;
                    PS = PSmat(PSindex)/60;
    
                    CBV = CBVmat(Vpind); %unit conversion

                    CBF = CBFmat(Findex)/60;

                    Ve = Vemat(Veind);

                    %Vp = Volume fraction of the plasma
                    Vp = CBV/100 * 1.04;

                    Ftm = CBF*(1-Hct_small)*1.04;

                    vptm = Vp*100;
                    vetm = Ve*100;
    
                    T_P = vptm./(PS+Ftm);
                    T_E = vetm./PS;
                    T_B = vptm./Ftm;

                    K_plus = 0.5*(1/T_P + 1/T_E + sqrt(((1/T_P+1/T_E)^2)-4*1/T_E*1/T_B));KPVal(Incrementor) = K_plus;
                    K_minus = 0.5*(1/T_P + 1/T_E - sqrt(((1/T_P+1/T_E)^2)-4*1/T_E*1/T_B));KMVal(Incrementor) = K_minus;

                    E_minus = (K_plus - 1/T_B)./(K_plus - K_minus);EMVal(Incrementor) = E_minus;

                    H_p = zeros(size(tspan));H_p(1) = 1;%
                    H_e = PS/vetm*(exp(-tspan.*PS/vetm));

                    H_p_all(Incrementor,:) = H_p;
                    H_e_all(Incrementor,:) = H_e;

%%
H_pAIF = conv(H_p,AIFfull);
H_pAIF = H_pAIF(1:size(tspan,2));
H_pAIF_all(Incrementor,:) = H_pAIF;

%%

H_eAIF2(1) = 0;
clear crp*
for k = 2:size(tspan,2)
int_t = tspan(k);
for j = 1:k
dummy_t = tspan(j);
expo1(j) =exp(-(int_t-dummy_t).*PS/vetm);
crpexp1(j) = AIFfull(j)*expo1(j);
end
t2 = tspan(1:k);
crpexp_integral1(k,Incrementor) = trapz(t2,crpexp1);
H_eAIF2(k) = PS*crpexp_integral1(k,Incrementor)/vetm;
end
tempsE(:,Incrementor,1) = crpexp_integral1(:,Incrementor);
H_eAIF_all(Incrementor,:) = H_eAIF2;

%%

%%

%% difference in the propagators for H_p_nL and H_p_wL

LC_CpAIF(1) = 0;
clear expo* crp* t2 int_t dummy_t
for k = 2:size(tspan,2)
int_t = tspan(k);
for j = 1:k
dummy_t = tspan(j);   
expo1(j) = exp(-(int_t-dummy_t)*((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm)))...
    *(((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 - Ftm/vptm + PS/(2*vetm))/(((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2) - 1)...
    *((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm));
crpexp1(j) = AIFfull(j)*expo1(j);
expo2(j) = (Ftm*exp(-(Ftm*(int_t-dummy_t))/vptm))/vptm;
crpexp2(j) = AIFfull(j)*expo2(j);
expo3(j) = (exp(-(int_t-dummy_t)*((Ftm + PS)/(2*vptm) - (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm)))*((Ftm + PS)/(2*vptm) - (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 + PS/(2*vetm))*((Ftm + PS)/(2*vptm) + (((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2)/2 - Ftm/vptm + PS/(2*vetm)))...
    /(((Ftm + PS)/vptm + PS/vetm)^2 - (4*Ftm*PS)/(vetm*vptm))^(1/2);
crpexp3(j) = AIFfull(j)*expo3(j);
end
t2 = tspan(1:k);
crpexp_integral1(k,Incrementor) = trapz(t2,crpexp1);crpexp_integral2(k,Incrementor) = trapz(t2,crpexp2);crpexp_integral3(k,Incrementor) = trapz(t2,crpexp3);
LC_CpAIF(k) = crpexp_integral1(k,Incrementor) + crpexp_integral2(k,Incrementor) - crpexp_integral3(k,Incrementor); 
end
temps(:,Incrementor,1) = crpexp_integral1(:,Incrementor);temps(:,Incrementor,2) = crpexp_integral2(:,Incrementor);temps(:,Incrementor,3) = crpexp_integral3(:,Incrementor);
H_LCAIF_all(Incrementor,:) = LC_CpAIF;

                end
            end
    end
end


figure;plot(tspan,H_pAIF_all,tspan,H_eAIF_all)


% % diff = exp(-t*((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve)))...
% %     *(((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 - F/Vp + PS/(2*Ve))/(((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2) - 1)...
% %     *((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve))...
% %     + (F*exp(-(F*t)/Vp))/Vp...
% %     - (exp(-t*((F + PS)/(2*Vp) - (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve)))*((F + PS)/(2*Vp) - (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 + PS/(2*Ve))*((F + PS)/(2*Vp) + (((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2)/2 - F/Vp + PS/(2*Ve)))...
% %     /(((F + PS)/Vp + PS/Ve)^2 - (4*F*PS)/(Ve*Vp))^(1/2);



