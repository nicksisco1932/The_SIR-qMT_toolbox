%% Simulations based on Quarles et al. "A theoretical framework to model DSC-MRI data acquired in the 
%  presence of contrast agent extravasation" 2009 Phys. Med. Biol. 54 5749

clear all
clc
close all

%% normal tissue
%Initialize all variables
CBV = 0.04; 
CBF = 60/6000; %normal tissue
flip = 90; 
R10 = 1/1.44; 
R20 = 1/0.035;%normal tissue
PS = 0;%normal tissue


%% tumor tissue
CBV = 0.04; 
CBF = 80/6000; %tumor tissue
flip = 90; 
R10 = 1/1.74; 
R20 = 1/0.045;%tumor tissue
PS = 30/6000;%tumor tissue
% PS = 0;

%% calculations
A = 200; 
B = 2; 
C = 1;
tAIF = 1:300;
t = 1:360; 
tp = 2;
rho = 1; kH = 0.733;

E = 1 - exp(-PS/(CBF*(1-kH)));
Ktrans = E*CBF*rho*(1-kH);
ve = 0.25; 
TE1 = 8.6e-3; 
TE2 = 35e-3; 
TR = 1;
r1 = 3.9; 
r2 = 5.3;
Kp = 0.55; 
Ke = 0.5*Kp;
vp = CBV; 
vi = 1 - vp - ve; 
S0 = 1;

%Define arterial input function
CAIFt = zeros(1,360);
CAIFt(61:end) = A*(tAIF/tp^2).*exp(-tAIF/tp) + B*(1-exp(-C*tAIF/tp));

figure; plot(CAIFt,'Linewidth',2); xlim([0 200]); ylim([-5 1.1*max(CAIFt)]); title('AIF');

%%


%Define vpCpt (product of blood volume fraction and blood plasma CA concentration)
vpCpt = (rho/kH) * CBF * conv(CAIFt, exp(-t/(CBV/CBF)));
vpCpt = vpCpt(1:max(t));
Cpt = vpCpt/vp;
figure; plot(vpCpt,'Linewidth',2); xlim([0 200]); ylim([-0.1*max(vpCpt) 1.1*max(vpCpt)]); title('plasma vp*Cp(t)');

%Define veCet (from Kety-Tofts model)
veCet = (Ktrans) * conv(CAIFt, exp(-(Ktrans/ve)*t));
veCet = veCet(1:max(t));
Cet = veCet/ve;
figure; plot(veCet,'Linewidth',2); xlim([0 200]); %ylim([-0.1*max(veCet) 1.1*max(veCet)]); title('EES ve*Ce(t)');

%Define signal intensity
R1micro = (r1*(veCet + vpCpt));
R1comp = S0*sin(flip*pi/180)*(1-exp(-TR*(R1micro+R10))) ./ (1-cos(flip*pi/180)*exp(-TR*(R1micro+R10)));
R2micro = (r2*(veCet + vpCpt));
R2meso = (Kp*vp*(ve*abs(Cpt - Cet)+vi*Cpt) + Ke*vi*veCet);
R2comp_TE1 = exp(-TE1 * (R2micro + R2meso + R20));
R2comp_TE2 = exp(-TE2 * (R2micro + R2meso + R20));
figure; plot(1:360,R1comp,'Linewidth',2); xlim([0 300]); ylim([0.8*min(R1comp) 1.1*max(R1comp)]); title('R1 component of signal');
figure; plot(1:360,R2comp_TE1,1:360,R2comp_TE2,'Linewidth',2); xlim([0 300]); title('R2 component of signal');legend('TE1','TE2')

St_TE1 = R1comp .* R2comp_TE1; St_TE2 = R1comp .* R2comp_TE2; 
figure; 
plot(1:360,St_TE1,1:360,St_TE2,'Linewidth',2); 
xlim([0 300]); ylim([-0.1*max(St_TE2) 1.1*max(St_TE1)]); 
title('Signal');legend('TE1','TE2')

%Define delta_R2star curves
delta_R2star = (-1/TE2)*log(St_TE2/mean(St_TE2(30:50)));
figure; 
plot(delta_R2star,'r','Linewidth',2); 
xlim([0 300]); title('Delta R2*'); %ylim([-0.1*max(delta_R2star) 2*max(delta_R2star)]);

delta_R2stard = (1/(TE2-TE1))*log((St_TE1/mean(St_TE1(30:50)))./(St_TE2/mean(St_TE2(30:50))));
figure(100); plot(delta_R2stard,'r','Linewidth',2); xlim([0 300]); title('Dual-Echo Delta R2*'); %ylim([-0.1*max(delta_R2star) 2*max(delta_R2star)]);

corr = delta_R2stard - Ke*veCet;
temp = delta_R2stard - 3.9*Ke*veCet;

