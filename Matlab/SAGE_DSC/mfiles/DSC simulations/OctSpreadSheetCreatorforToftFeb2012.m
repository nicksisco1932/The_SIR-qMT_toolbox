%Clears screen and all variables in the workspace.
clear all
close all
clc

idj = 1;
idk = 1;
%% 
tic
global Ftm;
global PStm;
global AIFfull;
global tspan;

Ktrans_all = zeros(30,21);
Vp_all = zeros(30,21);
CBV_meas = zeros(30,21);
CBV_true = zeros(30,21);
CBV_corr = zeros(30,21);

CBVmeas_perror = zeros(30,21);
CBVcorr_perror = zeros(30,21);

for idk = 1:21
for idj = 1:30
% % clearvars -except idj idk CBVmeas* CBVcorr* Ktrans_all Vp_all CBV_*


Cs0=[0 0];

%Npts  = Number of seconds the simulation should run.
Npts = 1000;
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
AIFfull = AIFfull/(max(AIFfull)/5.6845);


%-------------------PHYSIOLOGICAL PARAMETER DEFINITION---------------------
%This section defines the physiological parameters. Only one physiological
%parameter's graphing is scripted, thus only one parameter is given a range
%of values.


%Cerebral Blood Flow
%Cerebral Blood Volume
%Permeability, Surface Area Product
%Fractional volume of CA in the EES

% % CBFmat = [80]; %in brain
% % CBVmat = [8]; %in brain
% % PSmat = [0, 9.6]; %PS/Vp in brain (ref?)
% % Vemat = [0.35];
% % idx = 1;
% % 
% % CBFmat = [10.9]; %in muscle (from Goh et al. AJR 2006)
% % CBVmat = [1.8]; %in muscle (from Goh et al. AJR 2006)
% % PSmat = [0, 1.6]; %PS/Vp in muscle (from Brix et al. Radiology 1999)
% % Vemat = [0.1]; %in muscle (from Brix et al. Radiology 1999)
% % idx = 2;

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

CBFmat = [1.35*60]; %in avg tumor
if idj == 1; CBVmat = [0.5]; %in C6 tumor
else CBVmat = [100*Vp_all(idj-1,idk)+0.5]; %in C6 tumor
end
if idk == 1; Ktrans = 0; %in C6 tumor
else Ktrans = [Ktrans_all(idj,idk-1)+0.01]; %in C6 tumor
end
PSmat = [0, Ktrans/1.04/(CBVmat/100)]; %PS/Vp in avg  = Ktrans/rho (from Tofts PSrho = Ktrans)
Vemat = [0.25]; %ve in avg tumors
idx = 3;

%*********************************************************
%*********************************************************
%CREATES AN INCREMENTING VARIABLE TO TOGGLE ROW IN A (MATRIX FOR EXCEL)
%*********************************************************
%*********************************************************
%NOTE: Starts at 1 because it will be increased at beginning of each loop
%and AIF is in row 1

Betaincrementer = 0;
Aincrementer = 0;
%-----------------------PHYSIOLOGICAL VARIANCE LOOP OPEN------------------
%This for loop varies the parameter that is being graphed. This loop allows
%for curves with the associated parameter value to be placed on one figure.
                        %Vary Kees
%                         for KpandKeesindex = 1:length(Keesmat);
%                             Kees = Keesmat(KpandKeesindex);
%                             Kp = Kpmat(KpandKeesindex);
                       
                 
                        %Loop for PS
                        for PSindex = 1:length(PSmat);
                            PS = PSmat(PSindex)/60; %unit conversion
                     
                                         
                        %Loop for CBV
                        for CBVindex = 1:length(CBVmat);
                            CBV = CBVmat(CBVindex); %unit conversion
                            
                        %Loop for CBF
                        for CBFindex = 1:length(CBFmat);
                            CBF = CBFmat(CBFindex)/60;
                            
                        
                        %Loop for Ve
                        for Veindex = 1:length(Vemat);
                            Ve = Vemat(Veindex);
                        
                        BFval = num2str(CBF*60);
                        BVval = num2str(CBV);
                        PSval = num2str(PS*60);
                        Veval = num2str(Ve);
%*****************************************************
%******************************************************
%INCREASES THE ROW INCREMENTER FOR A (MATRIX FOR EXCEL)
%******************************************************
%******************************************************
%NOTE: Increments by 2 because adding 2 rows per loop (VeCe and VpCp). If
%you want to add rows (Possibly Ce and Cp) then increase the integer in
%this statement accordingly
Aincrementer = Aincrementer + 1;
Betaincrementer = Betaincrementer + 1;

%Vp = Volume fraction of the plasma
Vp = CBV/100;

%Vi = Volume fration of the extravascular intracellular space
Vi= 1 - Ve - Vp;

FVp = CBF/CBV;
Ftm = FVp;
%Ftm = CBF*(1-KH)*1.04;

%KPSVe = PS/Vp * (vp / ve) = KPS / Ve (like in Brix paper) - here PS is
%PS/VP
KPSVe = PS*(Vp/Ve);
PStm = KPSVe;
%PStm = PS;

[t,Cs]=ode45('twomodel',tspan,Cs0);

%VeCe
Ce = Cs(:,2);
VeCe = Ce*Ve;

%beta0 = [0.06/60, 0.3];
%[beta] = nlinfit(tspan,VeCe',@Toftsmodel,beta0);
%beta = lsqcurvefit(@Toftsmodel,beta0,tspan,VeCe',[0 0], [0.01 1.0]);

%VpCp
Cp = Cs(:,1);
VpCp=Cp*Vp;


%*******************************************************************
%*******************************************************************
%PLACES VeCe AND VpCp INTO THE NEXT TWO ROWS IN A (MATRIX FOR EXCEL)
%*******************************************************************
%*******************************************************************
%NOTE: Row is specified relative to Aincrementer. If you decide to add data
%(possibly Ce and Cp) then you will need to change these accordingly.
A(Aincrementer,1:Ns) = VeCe;
C(Aincrementer,1:Ns) = VpCp;
%beta(Aincrementer,1:2) = Betas(Betaincrementer,1:2);


                        end;
                        %B=[Vemat', A];
                        %B2 = B';
                        %D=[Vemat', C];
                        %D2 = D';
                        %sheetname1 = strcat('VpCp-BF=',BFval,', BV=',BVval,', kPS=',PSval);
                        %sheetname2 = strcat('VeCe-BF=',BFval,', BV=',BVval,', kPS=',PSval);
                        %xlswrite('Lym-Data',B2)%,sheetname2);
                        %xlswrite('Lym-Data',D2)%,sheetname1);
                        %xlswrite('Lym-Data',beta,sheetname2,'R1:S15');
                        end;
                        end;                                      
                        end;
                        
A2 = A(2,45:600)';
C2 = C(:,45:600)';

%plot(tspan/60,VeCe/Ve,tspan/60,VpCp/Vp)
                        
%beta(1)*60
%beta(2)
%plot(tspan/60,VeCe,tspan/60,Toftsmodel(beta,tspan))

%plot(tspan,A/Vemat)
%xlswrite('VeCe-Oct',B);
%xlswrite('VpCp-Oct',D);
%xlswrite('VeCe-BF60-BV4-Ve1',B);
%xlswrite('VpCp-BF60-BV4-Ve1',D);
if idx == 1
    Ce_brain = Ce;
    Cp_brain = Cp;
    BVval_brain = BVval;
    BFval_brain = BFval;
    PSval_brain = PSval;
    Veval_brain = Veval;
elseif idx == 2
    Ce_muscle = Ce;
    Cp_muscle = Cp;
    BVval_muscle = BVval;
    BFval_muscle = BFval;
    PSval_muscle = PSval;
    Veval_muscle = Veval;
elseif idx == 3
    Ce_tumor = Ce;
    Cp_tumor = Cp;
    BVval_tumor = BVval;
    BFval_tumor = BFval;
    PSval_tumor = PSval;
    Veval_tumor = Veval;
end

% % if exist('Ce_muscle','var') && exist('Ce_brain','var') && exist('Ce_tumor','var')
% % figure;set(gca,'LineWidth',1,'FontSize',18);set(gcf,'Color','w');
% % subplot(2,4,1);plot(AIFfull,'k','LineWidth',2);ylim([0 6]);title('AIF','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % subplot(2,4,2);plot(Ce_muscle,'k','LineWidth',2);ylim([0 2.25]);title('Muscle C_e','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % text(1100,0.85,{['BV = ' BVval_muscle],['BF = ' BFval_muscle],['PS = ' PSval_muscle],['V_e = ' Veval_muscle]},'FontSize',18);
% % subplot(2,4,3);plot(Ce_brain,'k','LineWidth',2);ylim([0 2.25]);title('Brain C_e','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % text(1100,0.85,{['BV = ' BVval_brain],['BF = ' BFval_brain],['PS = ' PSval_brain],['V_e = ' Veval_brain]},'FontSize',18);
% % subplot(2,4,4);plot(Ce_tumor,'k','LineWidth',2);ylim([0 2.25]);title('Tumor C_e','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % text(1100,0.85,{['BV = ' BVval_tumor],['BF = ' BFval_tumor],['PS = ' PSval_tumor],['V_e = ' Veval_tumor]},'FontSize',18);
% % subplot(2,4,6);plot(Cp_muscle,'k','LineWidth',2);ylim([0 4.65]);title('Muscle C_p','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % subplot(2,4,7);plot(Cp_brain,'k','LineWidth',2);ylim([0 4.65]);title('Brain C_p','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % subplot(2,4,8);plot(Cp_tumor,'k','LineWidth',2);ylim([0 4.65]);title('Tumor C_p','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % end

% % if exist('Ce_muscle','var') && exist('Ce_brain','var')
% % figure;set(gca,'LineWidth',1,'FontSize',18);set(gcf,'Color','w');
% % subplot(2,4,1);plot(AIFfull,'k','LineWidth',2);ylim([0 6]);title('AIF','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % subplot(2,4,2);plot(Ce_muscle,'k','LineWidth',2);ylim([0 2.25]);title('Muscle C_e','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % text(1100,0.85,{['BV = ' BVval_muscle],['BF = ' BFval_muscle],['PS = ' PSval_muscle],['V_e = ' Veval_muscle]},'FontSize',18);
% % subplot(2,4,3);plot(Ce_brain,'k','LineWidth',2);ylim([0 2.25]);title('Brain C_e','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % text(1100,0.85,{['BV = ' BVval_brain],['BF = ' BFval_brain],['PS = ' PSval_brain],['V_e = ' Veval_brain]},'FontSize',18);
% % subplot(2,4,4);plot(VeCe,'k','LineWidth',2);ylim([0 0.7]);title('Brain V_eC_e','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % % % text(1100,0.85,{['BV = ' BVval_tumor],['BF = ' BFval_tumor],['PS = ' PSval_tumor],['V_e = ' Veval_tumor]},'FontSize',18);
% % subplot(2,4,6);plot(Cp_muscle,'k','LineWidth',2);ylim([0 4.95]);title('Muscle C_p','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % subplot(2,4,7);plot(Cp_brain,'k','LineWidth',2);ylim([0 4.95]);title('Brain C_p','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % subplot(2,4,8);plot(VpCp,'k','LineWidth',2);ylim([0 0.7]);title('Brain V_pC_p','FontSize',18);set(gca,'LineWidth',1,'FontSize',18);
% % end


% % figure(1);plot(tspan,AIFfull,'k',tspan, Cp_tumor, 'r',tspan, Ce_tumor,'b', 'linewidth',2)
% % set(gca,'LineWidth',1,'FontSize',18);legend('AIF','Tumor Cp','Tumor Ce');
% % 
% % figure(2);plot(tspan, VpCp, 'r',tspan, VeCe,'b', 'linewidth',2)
% % set(gca,'LineWidth',1,'FontSize',18);legend('Tumor VpCp','Tumor VeCe');

%
r2sp = 87.7;
r2se = 30;
DR2s_meas = r2sp*VpCp + r2se*VeCe;
DR2s_true = r2sp*VpCp;
DR2s_corr = DR2s_meas - (r2sp*Vp + r2se*Ve)*(VpCp+VeCe);

% % figure(3);plot(tspan, DR2s_true, 'k',tspan, DR2s_meas,'r', tspan, DR2s_corr,'b','linewidth',2)
% % set(gca,'LineWidth',1,'FontSize',18);legend('True','Measured','Corrected');

Ktrans_all(idj,idk) = Ktrans;
Vp_all(idj,idk) = Vp;
CBV_meas(idj,idk) = trapz(DR2s_meas(1:120));
CBV_true(idj,idk) = trapz(DR2s_true(1:120));
CBV_corr(idj,idk) = trapz(DR2s_corr(1:120));

CBVmeas_perror(idj,idk) = abs(CBV_true(idj,idk) - CBV_meas(idj,idk))./CBV_true(idj,idk);
CBVcorr_perror(idj,idk) = abs(CBV_true(idj,idk) - CBV_corr(idj,idk))./CBV_true(idj,idk);
idj
toc
end
idk
toc
end
