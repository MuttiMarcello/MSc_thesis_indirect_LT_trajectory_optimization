clear all; close all; clc;
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% Time/fuel/energy optimal problem set-up

% cspice_furnsh('kernels\naif0012.tls'); %LSK leapseconds
% cspice_furnsh('kernels\gm_de440.tpc'); %PCK gravitational parameters
% cspice_furnsh('kernels\de440.bsp');    %SPK solar system eph
% 
% cspice_furnsh('kernels\3702319_2014YD.bsp');    % 2014 YD ephemerides
% cspice_furnsh('kernels\3550232_2010UE51.bsp');  % 2010 UE51 ephemerides
% cspice_furnsh('kernels\3568303_2011MD.bsp')     % 2011 MD ephemerides
% cspice_furnsh('kernels\3054374_2000SG344.bsp'); % 2000 SG334 ephemerides
% cspice_furnsh('kernels\20478784_2012UV136.bsp');% 2012 UV136 ephemerides

%-shooting-problem-set-up--------------------------------------------------

% launch window (redundant)
str_wo='2022-12-31 12:00 UTC';
str_wc='2024-12-31 12:00 UTC';

% mission parameters
m0=1500;  % [kg]
Tmax=10;  % [N]
Isp=3000; % [s]

% boundary constraints
x0=[-0.019488511458668 -0.016033479812051 0 8.918881923678198 -4.081793688818725 0].';
xf=[0.823385182067467 0 -0.022277556273235 0 0.134184170262437 0].';

%--------------------------------------------------------------------------

cspice_furnsh('..\kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

t_wo=cspice_str2et(str_wo); % [8400.5 mjd2000]
t_wc=cspice_str2et(str_wc); % [9131.5 mjd2000]

mu=1.21506683e-2;   % [-]

% dinemsional units
LU=3.84405000e5;    % [km]
TU=3.75676967e5;    % [s]
MU=m0;              % [kg]
VU=LU/TU;           % [km/s]
FU=MU*LU/TU^2*1e3;  % [N]

% exhaust velocity
g0=9.80665e-3;                              % [km/s^2]
c=Isp*g0;

% TO struct assembly
isFO=0;
TO_prob=struct_assembly([t_wo t_wc],m0,[Tmax/FU c/VU mu],x0,xf,isFO);

% EO/FO struct assembly
isFO=1;
FO_prob=struct_assembly([t_wo t_wc],m0,[Tmax/FU c/VU mu],x0,xf,isFO);

% system dynamics creation from symbolic setup
% CR3BP_sym_dyn(); % to KEEP
addpath('CR3BP_Dyn')

% TO solution search
TO_prob=TO_sol(TO_prob);

% scaling parameter lower bound
ll_TO=TO_prob.y0(8:14);
gL=-1/max(TO_prob.S);

% FO solution through solution scaling
FO_prob.t0=TO_prob.t0;
FO_prob.tf=TO_prob.tf;
FO_prob.tf_ad=TO_prob.tf_ad;
FO_prob.epsilon=0;

ll_FO=gL*ll_TO;
[~,~,FO_prob]=FO_ZFP(ll_FO,FO_prob);
FO_prob=DispRes(FO_prob);

% EO solution through solution scaling
EO_prob=FO_prob;
EO_prob.epsilon=1;

ll_EO=2*gL*ll_TO;
[~,~,EO_prob]=FO_ZFP(ll_EO,EO_prob);
EO_prob=DispRes(EO_prob);

%% Set default font size for axes labels, tick labels and legend
fs=15;
set(0,'DefaultAxesFontSize',fs);
set(0,'DefaultLegendFontSize',fs);

% Set the default linewidth
set(0,'DefaultLineLineWidth',2);

% Set the default text and legend interpreter
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f=figure;
f.Position=[1,1 800,400];
plot(TO_prob.tt,TO_prob.S)
hold on
plot(FO_prob.tt,FO_prob.S)
plot(EO_prob.tt,EO_prob.S)
hold off
grid on
grid minor
axis tight
xlabel('$t$')
ylabel('$S$')
legend('$S_{TO}$','$S_{FO}$','$S_{EO}$','Location','southwest')

axes('position',[0.608928571428571,0.201666666666667,0.268214285714286,0.249285714285715])
box on
plot(TO_prob.tt,TO_prob.S)
hold on
plot(FO_prob.tt,FO_prob.S)
plot(EO_prob.tt,EO_prob.S)
hold off
grid on
grid minor
axis tight
xlim([5.9 6.5])
ylim([-2,1])

%% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');