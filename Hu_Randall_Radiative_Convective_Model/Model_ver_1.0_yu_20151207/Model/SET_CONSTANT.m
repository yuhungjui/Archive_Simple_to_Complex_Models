% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% Procedure 2. Set Constants.
% 
% Set the constants / functions / parameters for model integration.
% 
% =========================================================================
% 
% Set Constants.
% 
% =========================================================================

g       =       9.80665; % gravity acceleration (m/s^2)
R0      =           287; % gas constant for dry air (J/(K*kg))
cp      =          1004; % specific heat of dry air at constant pressure 
                         % (J/(K*kg))
Gamma_d =          g/cp; % dry adiabatic temperature lapse rate (K/m) 
Gamma_m =        6.8e-3; % moist adiabatic tempreature lapse rate (K/m)
Lc      =         2.5e6; % latent heat of condensation (J/kg)
CD      =        1.5e-3; % surface turbulent transfer coefficient over 
                         % oceans

% =========================================================================
% 
% Set Scientific Functions.
% 
% =========================================================================

% Saturated Mixing Ratio.
% Functions 1:

Pws1    = @(T) 6.116441*10^((7.591386*(T-273.15))/((T-273.15)+240.7263));
        % Calculate saturated water vapor pressure (Pa) based on Vaisala 
        % HUMIDITY CONVERSION FORMULAS, T in unit of (K)

qs1     = @(P,Pws) 621.97*((Pws)/(P-Pws))*(1/1000); 
        % Saturated mixing ratio according to water vapor saturation 
        % pressure (kg/kg)

% Function 2:

Pws2    = @(T) calculate_saturation_vapor_pressure_liquid(T); 
        % Calculate saturated water vapor pressure (Pa) based on Murphy 
        % and Koop 2005 method, T in unit of (K)

qs2     = @(P,T,in,type_in,type_out) convert_humidity(P,T,in,type_in,type_out);
        % Calculate saturated mixing ratio according to water vapor 
        % saturation pressure (kg/kg)

% =========================================================================

