%% This is the main code for parameter optimization

clear all
close all
clc

% Allow enough decimal values
format long g 

% Load the experimental data. tEXP is the desorption time in [seconds] and
% yEXP is the hydrogen desorption rate in [mol/m^3/s]
load('expData','tEXP','yEXP');
time = tEXP';

% Define vector with initial D and C0 values for iterating
D = 5; %*1e-11 [m^2/s] (diffusivity is already multiplied by 1e-11 in funCOMSOL)
C0 = 3; %[mol/m^3]
DuCu0 = [D C0]; 

% Define vectors with lower (lb) and upper bounds (ub) for C0 and D as iteration limits,
% the order of limits corresponding to D and C0 in vectors lb and ub is [D C0]   
lb = [0.1 1];
ub = [20 5];

% Fitting algorithm 
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt'...
,'FunctionTolerance',1e-15,'Display','iter-detailed','OptimalityTolerance',1e-15...
,'StepTolerance', 1e-20,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-1);

[DuCu,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(@funCOMSOL,DuCu0,time,yEXP,lb,ub,options);

% Plot experimental and fitted curve
tPlot = 0:50:6000;
yPlot=funCOMSOL(DuCu,tPlot);
figure
plot(tEXP,yEXP,'b*')
hold on
plot(tPlot,yPlot,'r','linewidth',2)
grid on

% Determine the R2 of the fit
yfit = funCOMSOL(DuCu,tEXP);
meanyEXP = sum(yEXP,'all')/size(yEXP,1); % mean of experimental data
SSres = sum((yEXP - yfit).^2); % the sum of squares of residuals
SStot = sum((yEXP - meanyEXP).^2); % the total sum of squares 
Rsquared = 1 - SSres/SStot; 
R2 = round(Rsquared,5) % R^2 
