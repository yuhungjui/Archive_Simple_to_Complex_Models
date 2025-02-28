clear all;close all;clc;
addpath(genpath('../../../../../Work/MatLAB_TOOLS/'));
tic;
% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% =========================================================================
%% Load data:


% =========================================================================
%% Settings.
% Set constants:
g       =       9.80665; % gravity acceleration (m/s^2)
R0      =           287; % gas constant for dry air (J/(K*kg))
cp      =          1004; % specific heat of dry air at constant pressure (J/(K*kg))
Gamma_d =          g/cp; % dry adiabatic temperature lapse rate (K/m) 
Gamma_m =        6.8e-3; % moist adiabatic tempreature lapse rate (K/m)
Lc      =         2.5e6; % latent heat of condensation (J/kg)
CD      =        1.5e-3; % surface turbulent transfer coefficient over oceans

% Set functions:
Hs      = @(Tm) R0*Tm/g;                % layer mean scale height for layer mean temperature Tm(K) (m)
Zh      = @(Hs,p,p0) -Hs*log(p/p0);     % Hypsometric equation

% Set functions for saturated mixing ratio 1: 
Pws1    = @(T) 6.116441*10^((7.591386*(T-273.15))/((T-273.15)+240.7263));
                                        % Calculate saturated water vapor
                                        % pressure (Pa) based on 
                                        % Vaisala HUMIDITY CONVERSION FORMULAS,
                                        % T in unit of (K)
qs1     = @(P,Pws) 621.97*((Pws)/(P-Pws))*(1/1000); 
                                        % saturated mixing ratio according
                                        % to water vapor saturation pressure
                                        % (kg/kg)

% Set functions for saturated mixing ratio 2:
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
%% Settings.
% Set Parameters:
P0      =      1000*100; % pressure at surface (Pa)
PB      =       800*100; % pressure at mixed-layer top (Pa)
Pmid    =       500*100; % pressure at mid-level (Pa)
PT      =       100*100; % pressure at top of troposphere (Pa)
del_PM  =         P0-PB; % mixed layer depth (hPa)
SST     =           300; % sea surface temperature (K)
TRM     =           300; % radiative equilibrium temperature in the mixed-layer (K)
TRB     =           285; % radiative equilibrium temperature at mixed-layer top (K)             
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

Pq      =      100.*100; % pressure scale for moisture in free atmosphere (Pa)
c       =           0.6; % precipitating efficiency factor

Vs      =             2; % surface wind (m/s)
rho0    =         1.225; % standard sea level air density (kg/m^3)
tao     =      10*86400; % radiation relaxation time scale (s)

R       = (Pmid-PB)/(PT-PB); % For calculating G in Hu and Randall 1995
W       = cp*Zmid*(Gamma_d-Gamma_m)/q00; % For calculating G in Hu and Randall 1995

dt      =      0.5*3600; % time step (s)
tmax    =     100*86400; % total integration time (s)
smax    =  fix(tmax/dt); % total integration steps

tt      = [0:smax].*dt./86400; % time-axis for plotting (daily)

% Setting Variables (initial conditions):
TM(1)   =           295; % initial mean temperature in mixed-layer (K)
SM(1)   =        cp.*TM; % initial dry static energy in mixed-layer (J/kg=m^2/s^2)
TBt(1)  =           290; % initial mean temperature at top of mixed-layer (K)
SBt(1)  =   cp*TBt+g*ZB; % initial dry static energy at top of mixed-layer (J/kg=m^2/s^2)
qM(1)   =         12e-3; % initial water vapor mixing ratio in mixed-layer (kg/kg)
qBt(1)  =         16e-3; % initial water vapor mixing ratio at top of mixed-layer (kg/kg)

% =========================================================================
%% Time Integration:
for i = 1:smax;
    
    Fso         = rho0.*CD.*abs(Vs).*(S00-SM(i));
    Fqo         = rho0.*CD.*abs(Vs).*(q00-qM(i));
    a           = (SM(i)+Lc.*qM(i)-SBt(i))./(PT-PB);
    LHS         = (R-1).*(g.*(SBt(i)-SM(i))./del_PM+g.*a+g.*(1-c)./(PB-PT).*Lc.*qM(i))+(Lc.*R-W).*(g./del_PM.*(qM(i)-qBt(i)));
    RHS         = (1-R).*(g.*Fso./del_PM-(SM(i)-SRM)./tao+(SBt(i)-SRB)./tao)-(W-Lc.*R).*(g.*Fqo./del_PM);
    Mc(i)       = RHS./LHS;
    
    if ( Mc(i) <= 0 )
        Mc(i) = 0;
    end
        
    dSM         = (g.*Fso./del_PM-(SM(i)-SRM)./tao+g.*Mc(i).*(SBt(i)-SM(i))./del_PM).*dt;
    dSBt        = (-(SBt(i)-SRB)./tao-g.*Mc(i).*a-g.*(1-c)./(PB-PT).*Lc.*qM(i).*Mc(i)).*dt;
    dqM         = (g.*Fqo/del_PM-g.*Mc(i)./del_PM.*(qM(i)-qBt(i))).*dt;
    dqBt        = ((1-c).*qM(i)-qBt(i)).*g.*Mc(i)./Pq.*dt;
    
    SM(i+1)     = SM(i) + dSM;
    SBt(i+1)    = SBt(i) + dSBt;
    qM(i+1)     = qM(i) + dqM;
    qBt(i+1)    = qM(i) + dqBt;

end
% =========================================================================
%% Calculation

PR = c.*qM(1:smax).*Mc(1:smax).*3600./1000; % precipitation assumed according to convective mass flux (mm/hr)
PR(PR<=0) = 0;

% =========================================================================
%% Plotting:
subplot(2,1,1)
plot(tt(1:smax),SBt(1:smax))
subplot(2,1,2)
plot(tt(1:smax),PR(1:smax))
% =========================================================================
%% Display runnung time.
time_cost = toc;
disp(time_cost);