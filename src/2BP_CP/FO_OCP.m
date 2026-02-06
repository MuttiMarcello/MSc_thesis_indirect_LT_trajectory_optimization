clear all; close all; clc;
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP constant parameters

% Fuel optimal exploration space problem set-up

% Loaded kernels
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

% launch window
str_wo='2022-12-31 12:00 UTC';
str_wc='2024-12-31 12:00 UTC';

% target & mission parameters
targ='3054374';
m0=22.6;  % [kg]
Pmax=120; % [W]
Pmin=20;  % [W]
isFO=1;

% 3054374  2000 SG344 [538]
% 3550232  2010 UE51  [392]
% 3568303  2011 MD    [504]
% 20478784 2012 UV136 [428]
% 3702319  2014 YD    [682]

%--------------------------------------------------------------------------

cspice_furnsh('..\kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

t_wo=cspice_str2et(str_wo); % [8400.5 mjd2000]
t_wc=cspice_str2et(str_wc); % [9131.5 mjd2000]

% struct assembly
FO_prob=struct_assembly(targ,[t_wo t_wc],m0,[Pmin Pmax],isFO);

% system dynamics creation from symbolic setup
% FO_SEP_sym_dyn(FO_prob); % to KEEP
addpath('FO_dyn')

load('temp.mat');

FO_prob=FO_CONT(FO_prob,TO_prob);