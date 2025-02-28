clear all
% Hu and Randall 1994- Radiative Convective Feedback

SST=300;                % sea surface temperature
g=9.8;                  % gravity accelation (m/s^2)
cp=1004;                % J/K
delta_PM=200;           % mixed layer depth (hPa)
SRM=cp*300;             % dry static reference value in mixed layer (J/m^2)
tao=6*86400;            % Radiation-Relaxation (s)
dt=0.5*3600;            % time step (s)
SRB=cp*290+g*1000;      % dry static energy value in the base of free tropos (J/m^2)
c=0.6;                  % precipitating fraction 
L=2.5*10^6;             % Latent heat release (J/kg)
Pq=100;                 % scale-height of moisture (m)
Vs=2;                   % surface wind (m/s)
rao=1.25;               % kg/m^3 (surface air density)
PT=200;                 % top of troposphere
PB=1000;                % surface pressure
qoo=24*10^(-3);         % kg/kg      
soo=cp*SST;             % surface dry static energy 
Gamma_d=9.8;            % dry lapse rate (K/m) 
Gamma_m=6.8;            % moist lapse rate (K/m)
R=(500-800)/(200-800);  
W=cp*5000*(Gamma_d-Gamma_m)/qoo;
CD=2.5*10^-3;
tend=200;               % days

steps=fix(tend*86400/dt);
xx=[1:1:steps+1]/48;

% initial condition
SM(1)=1004*288;
SBt(1)=1004*279+g*2000;
qM(1)=6*10^(-3); 
qBt(1)=16*10^(-3);




for i=1:steps

    
Fso=rao*CD*abs(Vs)*(soo-SM(i));
Fqo=rao*CD*abs(Vs)*(qoo-qM(i));    
a=(SM(i)+2.5*10^6*qM(i)-SBt(i))/(PT-PB);
RR=(R-1)*(g*(SBt(i)-SM(i))/delta_PM+g*a+g*(1-c)/(PB-PT)*L*qM(i))+(L*R-W)*(g/delta_PM*(qM(i)-qBt(i)));
BB=(1-R)*(g*Fso/delta_PM-(SM(i)-SRM)/tao+(SBt(i)-SRB)/tao)-(W-L*R)*(g*Fqo/delta_PM);
Mc=BB/RR;

if Mc<=0
Mc=0;
end


dSM=(g*Fso/delta_PM-(SM(i)-SRM)/tao+g*Mc*(SBt(i)-SM(i))/delta_PM)*dt;
dSBt=(-(SBt(i)-SRB)/tao-g*Mc*a-g*(1-c)/(PB-PT)*L*qM(i)*Mc)*dt;
dqM=(g*Fqo/delta_PM-g*Mc/delta_PM*(qM(i)-qBt(i)))*dt;
dqBt=((1-c)*qM(i)-qBt(i))*g*Mc/Pq*dt;

SM(i+1)=SM(i)+dSM;
SBt(i+1)=SBt(i)+dSBt;
qM(i+1)=qM(i)+dqM;
qBt(i+1)=qM(i)+dqBt;

end

plot(xx(1:end),SBt(1:end));
xlabel('days')
set(gca,'fontsize',15)
