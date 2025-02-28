clear all;close all;clc;
addpath(genpath('../../../../../Work/MatLAB_TOOLS/'));
tic;
% =========================================================================
% 
% Radiative-Convective Systems Idealized Model - Equilibrium.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% This formulation calculates the equilibrium solutions of the idealized
% model from Hu and Randall 1995. The main purpose is to reproduce its Fig.
% 3 (a), (b), (c) and test the linear solutions with DYNAMO data input.
% 
% =========================================================================
% %% Load data:
% load ../Data/DATA_Gan_ARM_MET_3hr.mat
% load ../Data/DATA_Gan_PBL.mat
% load ../Data/DATA_Gan_SondeParameter.mat

% =========================================================================
%% Settings.

% *Set constants:
g       =       9.80665; % gravity acceleration (m/s^2)
R0      =           287; % gas constant for dry air (J/(K*kg))
cp      =          1004; % specific heat of dry air at constant pressure (J/(K*kg))
Gamma_d =          g/cp; % dry adiabatic temperature lapse rate (K/m) 
Gamma_m =        6.8e-3; % moist adiabatic tempreature lapse rate (K/m)
Lc      =         2.5e6; % latent heat of condensation (J/kg)
CD      =        1.5e-3; % surface turbulent transfer coefficient over oceans

% *Set functions:
Hs      = @(Tm) R0*Tm/g;                % layer mean scale height for layer mean temperature Tm(K) (m)
Zh      = @(Hs,p,p0) -Hs*log(p/p0);     % Hypsometric equation

% *Set functions for saturated mixing ratio 1: 
Pws1    = @(T) 6.116441*10^((7.591386*(T-273.15))/((T-273.15)+240.7263));
                                        % Calculate saturated water vapor
                                        % pressure (Pa) based on 
                                        % Vaisala HUMIDITY CONVERSION FORMULAS,
                                        % T in unit of (K)
qs1     = @(P,Pws) 621.97*((Pws)/(P-Pws))*(1/1000); 
                                        % saturated mixing ratio according
                                        % to water vapor saturation pressure
                                        % (kg/kg)

% *Set functions for saturated mixing ratio 2:
% Pws2    = @(T) calculate_saturation_vapor_pressure_liquid(T); 
                                        % Calculate saturated water vapor
                                        % pressure (Pa) based on 
                                        % Murphy and Koop 2005 method,
                                        % T in unit of (K)
% qs2     = @(P,T,in,type_in,type_out) convert_humidity(P,T,in,type_in,type_out);
                                        % Calculate saturated mixing ratio
                                        % according to water vapor
                                        % saturation pressure
                                        % (kg/kg)                                 
                                       
% =========================================================================
% *Repeatedly setting parameters and calculating the linear solutions with
% varying SST and Surface Wind (Vs):
 
xx = linspace(298,304,100);     % SST range (K)
yy = linspace(1.0,10.0,100);    % Surface WInd (Vs) range (m/s)

for xi = 1:length(xx);
    for yi = 1:length(yy);
        
        %% Settings.
        % Set Parameters:
        P0      =     1000.*100; % pressure at surface (Pa)
        PB      =      800.*100; % pressure at mixed-layer top (Pa)
        Pmid    =      500.*100; % pressure at mid-level (Pa)
        PT      =      200.*100; % pressure at top of troposphere (Pa)
        del_PM  =         P0-PB; % mixed layer depth (Pa)
        
        % Varying parameter:
        SST     =        xx(xi); % sea surface temperature (K)
        
        TRM     =           288; % radiative equilibrium temperature in the mixed-layer (K)
        TRB     =           278; % radiative equilibrium temperature at mixed-layer top (K)
        Z0      =             0; % surface height (m)
        ZB      = Zh(Hs(mean([SST,TRB])),PB,P0); % mixed-layer top height (m)
        Zmid    =          5000; % mid-level height (m)
        S00     =        cp*SST; % dry static energy at the surface (J/kg=m^2/s^2)
        SRM     =        cp*TRM; % dry static energy reference in mixed layer (J/kg=m^2/s^2)
        SRB     =   cp*TRB+g*ZB; % dry static energy at mixed-layer top (J/kg=m^2/s^2)
        
        q001    = qs1(P0./100,Pws1(SST)); q00 = q001;
        % saturated mixing ratio at the surface (kg/kg) based on method 1
        % q002    = qs2(P0,SST,Pws2(SST),'partial pressure','mixing ratio'); q00 = q002;
        % saturated mixing ratio at the surface (kg/kg) based on method 2
        
        c       =           0.6; % precipitating efficiency factor
        
        % Varying parameter:
        Vs      =        yy(yi); % surface wind (m/s)
        
        rho0    =         1.225; % standard sea level air density (kg/m^3)
        tau     =      20*86400; % radiation relaxation time scale (s)
        v       =  rho0.*CD.*Vs; % symbol for rho0.*CD.*Vs (kg/m^2/s))
        
        R       = (Pmid-PB)/(PT-PB); % For calculating G in Hu and Randall 1995
        W       = cp*Zmid*(Gamma_d-Gamma_m)/q00; % For calculating G in Hu and Randall 1995
        
        % *Set convective mass flux (Mc) of linear solution,
        % in which coeff_B is needed for calculation:
        coeff_B = ( g.*v.*S00./del_PM)+(SRM./tau ); 
        
        Mc(xi,yi) = ...
            ( ...
            ( ( SRB./tau )-((q00.*(W-Lc.*R))./(tau.*(1-R))) ) ...
            .* ( (g.*v./del_PM)+(1./tau) ) ...
            -  ( coeff_B./tau ) ...
            ) ...
            ./( ...
              ( c.*coeff_B./(v.*tau) ) ...
            - ( (c.*SRB./(v.*tau)).*((g.*v./del_PM)+(1./tau)) ) ...
            + ( (q00.*(W-Lc.*R).*g.*v)./(tau.*(1-R).*del_PM) ) ...
            + ( (g.*q00./(PB-PT)).*((W-Lc.*R)./(1-R)-c).*(g.*v./del_PM+(1./tau)) ) ...
            );
        
        if ( Mc(xi,yi) <= 0 )
            Mc(xi,yi) = 0;
        end
          
        % *Set coefficients for linear solutions:
        coeff_A = (g.*v./del_PM) + (1./tau) + (g.*Mc(xi,yi)./del_PM);
        coeff_B = (g.*v.*S00./del_PM) + (SRM./tau);
        coeff_C = (g.*Mc(xi,yi)./del_PM);
        coeff_D = (g.*Mc(xi,yi)./(PB-PT));
        coeff_F = (g.*(1-c).*Mc(xi,yi)./(PB-PT));
        coeff_H = (1+c.*Mc(xi,yi)./v);
        
        %% Calculation:
        % Calculate the linear solutions:
        
        % 1. dry static energy at top of mixed-layer (J/kg=m^2/s^2)
        SBt(xi,yi) = ((coeff_B.*coeff_D./coeff_A)+(SRB./tau)+(Lc.*q00.*(coeff_D-coeff_F)./coeff_H)) ...
                   ./((1./tau)+(coeff_D)-(coeff_C.*coeff_D./coeff_A));
        
        % 2. dry static energy in mixed-layer (J/kg=m^2/s^2)
        SM(xi,yi) = (coeff_B+coeff_C.*SBt(xi,yi))./(coeff_A);
        
        % 3. water vapor mixing ratio at top of mixed-layer (kg/kg)
        qBt(xi,yi) = (q00.*(1-c))./(1+(c.*Mc(xi,yi)./v));
        
        % 4. water vapor mixing ratio in mixed-layer (kg/kg)
        qM(xi,yi) = (q00)./(1+(c.*Mc(xi,yi)./v));
        
    end
end

% =========================================================================
%% Plotting:
% contourf(xx,yy,SM./cp);
contourf(xx,yy,qM.*1000);
% contourf(xx,yy,c.*qM.*Mc);
colormap;
colorbar;
axis square;
% =========================================================================
%% Display runnung time.
time_cost = toc;
disp(time_cost);