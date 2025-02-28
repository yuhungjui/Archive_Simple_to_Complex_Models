function F=balanced(X,Vs,SST)


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

Pws1    = @(T) 6.116441*10^((7.591386*(T-273.15))/((T-273.15)+240.7263));
                                        % Calculate saturated water vapor
                                        % pressure (Pa) based on 
                                        % Vaisala HUMIDITY CONVERSION FORMULAS,
                                        % T in unit of (K)
qs1     = @(P,Pws) 621.97*((Pws)/(P-Pws))*(1/1000); 
                                        % saturated mixing ratio according
                                        % to water vapor saturation pressure
                                        % (kg/kg)
%% Settings.
% Set Parameters:
P0      =      1000*100; % pressure at surface (Pa)
PB      =       800*100; % pressure at mixed-layer top (Pa)
Pmid    =       500*100; % pressure at mid-level (Pa)
PT      =       100*100; % pressure at top of troposphere (Pa)
del_PM  =         P0-PB; % mixed layer depth (hPa)
TRM     =           295; % radiative equilibrium temperature in the mixed-layer (K)
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
rho0    =         1.225; % standard sea level air density (kg/m^3)
tao     =      20*86400; % radiation relaxation time scale (s)
R       = (Pmid-PB)/(PT-PB); % For calculating G in Hu and Randall 1995
W       = cp*Zmid*(Gamma_d-Gamma_m)/q00; % For calculating G in Hu and Randall 1995

%%


F(1)= (g*rho0*CD*abs(Vs)*(S00-X(1)))/del_PM    -(X(1)-SRM)/tao                           +g*X(5)*(X(2)-X(1))./del_PM;
F(2)= -(X(2)-SRB)/tao                          -g*X(5)*((X(1)+Lc*X(3)-X(2))/(PT-PB))     -g*(1-c)/(PB-PT)*Lc*X(3)*X(5);
F(3)= g*(rho0*CD*abs(Vs)*(q00-X(3)))/del_PM    -g*X(5)/del_PM*(X(3)-X(4));
F(4)= ((1-c)*X(3)-X(4))*g*X(5)/Pq;
F(5)= (1-R)*(X(1)-X(2))-(W-Lc*R)*X(3);


end