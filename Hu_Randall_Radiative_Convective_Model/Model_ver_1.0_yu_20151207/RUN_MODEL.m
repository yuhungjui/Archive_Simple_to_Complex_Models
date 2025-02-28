clear; close all; clc;
% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% =========================================================================
%% Load data:

run ./Model/DATA_LOADING;

% =========================================================================
%% Settings.

run ./Model/SET_CONSTANT;

run ./Model/SET_NAMELIST;

% =========================================================================
%% Time Integration and Calculation:

run ./Model/TIME_INTEGRATION;

% =========================================================================
%% Output

output_name = 'Test_1';

clear ARM_PBL_Gan Gan_ARM_MET_3hr Gan_Sonde_Full Revelle_MET_Full

save(['./Output/',output_name],'*');

copyfile('./Model/NAMELIST.dat',['./Output/',output_name,'_Namelist.dat'])

% =========================================================================
%% ========================================================================