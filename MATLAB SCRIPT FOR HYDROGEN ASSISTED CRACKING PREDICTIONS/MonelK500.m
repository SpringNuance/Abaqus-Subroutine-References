% Matlab script to compute the variation with the applied potential of
% (a) the cracking threshold, and (b) the stage II crack growth rate.
% Corresponds to Fig. 3 of E. Martínez-Pañeda, C.F Niordson, R.P. Gangloff. 
% Strain gradient plasticity-based modeling of hydrogen environment 
% assisted cracking. Acta Materialia 117 (2016) 321-332
%
% E. Martínez-Pañeda
% mail@empaneda.com

clear all
clc

% Flag variable,
% 0 - Plot model predictions for a given Sh and Alpha
% 1 - Compute Alpha for a given Kth, Sh and Ea
% 2 - Compute Ccrit for a given Alpha and da/dt
flag=0;

if flag==0 %Model predictions
% Experimental data
 x=[-1100 -1005 -1000 -995 -950]; % ATI Allvac (intergranular HEAC)
 y=[16 19 14 17 30];
 yA=[0.17 0.023 0.018 0.019 0.004];
 x1=-934; % Special Metals (intergranular HEAC)
 y1=32;
 y1A=0.016;
 x2=[-900 -800]; % ATI Allvac (grain boundary segregation)
 y2=[18.5 44];
 y2A=[0.01 0.003];
 x3=[-850 -850 -850 -850 -750 -700]; % ATI Allvac
 y3=[73.9 70.6 84.0 120.0 90 77];
 y3A=[2e-5 2.25e-5 2e-5 2.5e-5 2e-5 2e-5];
 x4=-690; % Special Metals
 y4=77;
 y4A=2e-5;

% Assign manually: 
 alpha=6.36; %Alpha
 Sy=767 ;%Sy - yield stress
 Beta=0.2; % Beta' [see Gerberich et al. (1991)]
 Alpha2=0.0002; %Alpha'' [see Gerberich et al. (1991)]
 Kig=0.880; %Kig
 Ccrit=407; %Ccrit
 Sh=8.15*Sy; % Hydrostatic stress at r=0.1 \mu m
% Please visit www.empaneda.com/codes to download a finite element SGP code 
% to compute the hydrostatic stress accounting for the role of GNDs

 Ea=-1150:10:-760; % Applied potential (mVsce)

% Kth vs Ea
 Ch=(-52.5-0.0687.*Ea)*exp((Sh*1.73)/(8.3142*298));
% Conversion
 Ch1=((Ch./10000)./1.00794)./((Ch./10000)./1.00794+(100-(Ch./10000))./58.6934);
 Kth=(1./Beta).*(exp(((Kig-alpha.*Ch1).^2)./(Sy*Alpha2)));

 plot(Ea,Kth,'k')
 xlabel('$E_{APP}$ (m$V_{SCE}$)','Interpreter','LaTex','FontSize',16) 
 ylabel('$K_{TH}$ (MPa$\sqrt{m}$)','Interpreter','LaTex','FontSize',16)
 xlim([-1150 -670])
 ylim([0 160])
 xlim manual
 hold on
 plot(x,y,'ms','MarkerFaceColor','m')
 plot(x1,y1,'y^','MarkerFaceColor','y')
 plot(x2,y2,'gs','MarkerFaceColor','g')
 plot(x3,y3,'ks',x4,y4,'k^')

% da/dtII vs Ea
 Ea=-1150:0.5:-840;
 Ch=(-52.5-0.0687.*Ea)*exp((Sh*1.73)/(8.3142*298));
 dadt=0.04.*(erfinv(1-Ccrit./Ch)).^2;
 figure
 semilogy(Ea,dadt,'k')
 xlabel('$E_{APP}$ (m$V_{SCE}$)','Interpreter','LaTex','FontSize',16) 
 ylabel('$da/dt_{II}$ ($\mu m/s)$','Interpreter','LaTex','FontSize',16)
 xlim([-1150 -670])
 ylim([10e-8 10])
 xlim manual
 hold on
 plot(x,yA,'ms','MarkerFaceColor','m')
 plot(x1,y1A,'y^','MarkerFaceColor','y')
 plot(x2,y2A,'gs','MarkerFaceColor','g')
 plot(x3,y3A,'ks',x4,y4A,'k^')

elseif flag==1 %Computation of alpha

% Assign manually: 
 Ea=-1000; % Applied potential (mVsce)
 Kth=17.33; % Threshold Stress Intensity
 Sy=773; %Sy - yield stress
 Sh=8.1*Sy; % Hydrostatic stress at r=0.1 \mu m
 Beta=0.2; % Beta'
 Alpha2=0.0002; %Alpha''
 Kig=0.880; %Kig

% Compute Ch in wppm
 Ch=(-52.5-0.0687*Ea)*exp((Sh*1.73)/(8.3142*298));
% Conversion
 Ch1=(Ch/10000)/((Ch/10000)/1.00794+(100-(Ch/10000))/58.6934);

 syms alpha
 Eqn=Kth==(1/Beta)*(exp(((Kig-alpha*Ch1)^2)/(Sy*Alpha2)));
 Alpha=double(solve(Eqn,alpha))

else % Computation of Ccrit
% Assign manually
 alpha=6.36;
 dadt=0.018;
 Kig=0.880;
 alpha2=0.0002;
 Sy=773;
 Beta=0.2;
 Kth=17.33;

 Ccrit=(1/alpha)*(Kig-sqrt(alpha2*Sy*log(Kth*Beta)))*...
     (1-erf(sqrt(dadt/0.04)));

% Conversion
 syms Ch
 Eqn=Ccrit==(Ch/10000)/((Ch/10000)/1.00794+(100-(Ch/10000))/58.6934);
 Ccrit=double(solve(Eqn,Ch)) 
end
