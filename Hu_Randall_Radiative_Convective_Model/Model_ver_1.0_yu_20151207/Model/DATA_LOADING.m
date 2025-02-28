% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% Procedure 1. Data Loading
% 
% Load data from DYNAMO Field Campaign for initial inputs.
% For now, we focus on the data from Gan Island (Addu Atoll).
% 
% =========================================================================

% Radiosonde at Gan:
load ../../Data/DATA_Gan_SondeParameter_Full;

% ARM MET Surface Data at Gan:
load ../../Data/DATA_Gan_ARM_MET_3hr;

% ARM PBL at Gan:
load ../../Data/DATA_ARM_PBL_Gan;

% R/V Revelle SST Data:
load ../../Data/DATA_Revelle_MET_Full;

% =========================================================================