clear all;close all;clc;
addpath(genpath('../../../../../Work/MatLAB_TOOLS/'));
tic;
% =========================================================================
% 
% Calculate the parameters for Hu and Randall Radiative-COnvective Systems
% Model (Hu and Randall 1995) from DYNAMO sounding data at given site,
% first of all, Gan Island.
% 
% The parameters derived from the observational datasets are needed for
% identifying the climatological equilibrium state and also as the input
% for model integration.
% 
% Those parameters are as follows,
% theta     =   potential temperature (K)
% theta_e   =   equivalent potentail temperature (K)
% PB        =   pressure at top of mixed-layer (hPa)
% ZB        =   height at top of mixed-layer (m)
% S         =   dry static energy (J/kg=m^2/s^2)
% q         =   water vapor mixing ratio (kg/kg)
% 
% =========================================================================
%% Load data:
% Load sounding:
load ../Data/DATA_Sonde_Gan_incMissingT.mat

% Load Revelle MET data:
load ../Data/DATA_Revelle_MET.mat

% Set time table:
Set_TimeTable_daily_3hrly;

% =========================================================================
%% Settings
% Set constants:
g       =       9.80665; % gravity accelation (m/s^2)
R0      =           287; % gas constant for dry air (J/(K*kg))
cp      =          1004; % specific heat of dry air at constant pressure (J/(K*kg))
Gamma_d =          g/cp; % dry adiabatic temperature lapse rate (K/m) 
Gamma_m =        6.8e-3; % moist adiabatic tempreature lapse rate (K/m)
Lc      =         2.5e6; % latent heat of condensation (J/kg)
CD      =        1.5e-3; % surface turbulent transfer coefficient over oceans

% Set functions:
% Potential temperature, P in [Pa], T in [K].
theta   = @(P,T) T.*(1e5./P).^(R0/cp);

% Equivalent Potential temperature, theta in [K], qs in [kg/kg], T in [K].
theta_e = @(theta,qs,T) theta.*exp((Lc.*qs)./(cp.*T));

% Dry Static Energy, T in [K], z in [m].
dse     = @(T,z) cp.*T+g.*z;

% Moist Static Energy, dse in [J/kg], q in [kg/kg].
mse     = @(dse,q) dse+Lc.*q;

% =========================================================================
%% Organize data:
Gan_Sonde.Time      = TT;
Gan_Sonde.P         = repmat(SOD.D_PZ,[1,736]).*100;
Gan_Sonde.Z         = SOD.D_ALT;
Gan_Sonde.T         = SOD.D_T+273.15;
Gan_Sonde.Td        = SOD.D_Td+273.15;
Gan_Sonde.RH        = SOD.D_RH;
Gan_Sonde.RHice     = SOD.D_RHI;
Gan_Sonde.U         = SOD.D_U;
Gan_Sonde.V         = SOD.D_V;
Gan_Sonde.Speed     = sqrt(SOD.D_U.^2 + SOD.D_V.^2);
Gan_Sonde.q         = SOD.D_MR;
Gan_Sonde.qq        = SOD.D_Q;
Gan_Sonde.P_freeze  = SOD.P_zero;
Gan_Sonde.Z_freeze  = SOD.A_zero;
clear SOD

% =========================================================================
%% Calculation:
% Potential temperarure [K].
Gan_Sonde.Theta     = theta(Gan_Sonde.P,Gan_Sonde.T);

% Equivalent Potential temperarure [K].
Gan_Sonde.Theta_e   = theta_e(Gan_Sonde.Theta,Gan_Sonde.q,Gan_Sonde.T);

% Dry Static Energy [J/kg].
Gan_Sonde.DSE       = dse(Gan_Sonde.T,Gan_Sonde.Z);

% Moist Static Energy [J/kg].
Gan_Sonde.MSE       = mse(Gan_Sonde.DSE,Gan_Sonde.q);

% =========================================================================
%% Save the data:
save ./DATA_Gan_SondeParameter.mat Gan_Sonde
% =========================================================================
%% Display runnung time.
time_cost = toc;
disp(time_cost);


