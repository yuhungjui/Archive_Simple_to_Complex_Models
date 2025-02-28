% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% Procedure 3. Set Namelist.
% 
% Set the parameters and initial values for model integration.
% While setting, write out the model namelist to a txt.
% 
% 
% =========================================================================
% 
% Prepare saving namelist to NAMELIST.dat.
% 
% =========================================================================

fid = fopen('./NAMELIST.dat','w');

% =========================================================================
% 
% Find Interested Time Duration.
% 
% =========================================================================

Time_ini = '100100';  % The interested start time.
Time_end = '123121';  % The interested end time.

TimeID_ini = find( (Gan_Sonde_Full.Time(:,1) == str2double(Time_ini(1:2))) ... 
                 & (Gan_Sonde_Full.Time(:,2) == str2double(Time_ini(3:4))) ...
                 & (Gan_Sonde_Full.Time(:,3) == str2double(Time_ini(5:6))) );
             
TimeID_end = find( (Gan_Sonde_Full.Time(:,1) == str2double(Time_end(1:2))) ...
                 & (Gan_Sonde_Full.Time(:,2) == str2double(Time_end(3:4))) ...
                 & (Gan_Sonde_Full.Time(:,3) == str2double(Time_end(5:6))) );

% Write out
fprintf(fid, '============================================================= \n\n');
fprintf(fid, '  Parameters \n\n');
fprintf(fid, '============================================================= \n\n');
fprintf(fid, 'Start Time = %s \n', Time_ini);
fprintf(fid, 'End Time   = %s \n\n', Time_end);

% =========================================================================
% 
% Set/Write-out Parameters: Model parameters.
% 
% =========================================================================

Z0      = 0; 
        % surface height (m)
        
ZB      = mean(ARM_PBL_Gan.PBLH_LL(TimeID_ini:TimeID_end));    
        % mean mixed-layer top height (m)
      
Zmid    = interp1(Gan_Sonde_Full.P(:,1),mean(Gan_Sonde_Full.Z(:,TimeID_ini:TimeID_end),2),500*100); 
        % mid-level height (500-hPa height) (m)

P0      = mean(Gan_ARM_MET_3hr.P(TimeID_ini:TimeID_end))*100;               
        % mean pressure at surface (Pa)

PB      = interp1(mean(Gan_Sonde_Full.Z(:,TimeID_ini:TimeID_end),2),Gan_Sonde_Full.P(:,1),ZB); 
        % pressure at mixed-layer top (Pa)

Pmid    = 500*100; 
        % pressure at mid-level (Pa)
        
PT      = 100*100; 
        % pressure at top of troposphere (Pa)
        
del_PM  = P0-PB; 
        % mixed layer depth (Pa)

SST     = mean(Revelle_MET_Full.SST)+273.15; 
        % mean sea surface temperature (K)
        
TRM     = mean(mean(Gan_Sonde_Full.T(Gan_Sonde_Full.P(:,1) >= PB,TimeID_ini:TimeID_end))); 
        % radiative equilibrium temperature in the mixed-layer (K)
        % set to the mean temperature in the mixed-layer (K) for now

TRB     = interp1(mean(Gan_Sonde_Full.Z(:,TimeID_ini:TimeID_end),2), ...
                  mean(Gan_Sonde_Full.T(:,TimeID_ini:TimeID_end),2), ZB);
        % radiative equilibrium temperature at mixed-layer top (K)
        % set to the mean temperature at top of the mixed-layer (K) for now
        
S001    = cp*SST; 
S002    = cp.*mean(Gan_ARM_MET_3hr.T(TimeID_ini:TimeID_end)+273.15);
S003    = mean(Gan_Sonde_Full.DSE(1,TimeID_ini:TimeID_end));
S00     = S001;
        % 3 approaches to dry static energy at the surface (J/kg=m^2/s^2)
        
SRM     = mean(mean(Gan_Sonde_Full.DSE(Gan_Sonde_Full.P(:,1) >= PB,TimeID_ini:TimeID_end))); 
        % dry static energy reference in mixed layer (J/kg=m^2/s^2)

SRB     = interp1(mean(Gan_Sonde_Full.Z(:,TimeID_ini:TimeID_end),2), ...
                  mean(Gan_Sonde_Full.DSE(:,TimeID_ini:TimeID_end),2), ...
                  ZB);
        % dry static energy at mixed-layer top (J/kg=m^2/s^2)

q001    = qs1(P0./100,Pws1(SST));  
q002    = mean(Revelle_MET_Full.qsea(TimeID_ini:TimeID_end))./1000; 
q00     = q001;
        % 2 approaches for saturated mixing ratio at sea surface (kg/kg)

Pq      = 100.*100; 
        % pressure scale for moisture in free atmosphere (Pa)

c       = 0.6; 
        % precipitating efficiency factor

Vs1     = mean(Gan_Sonde_Full.Speed(1,TimeID_ini:TimeID_end));
Vs2     = mean(ARM_PBL_Gan.WSPD(1,TimeID_ini:TimeID_end)); 
Vs      = Vs1;
        % 2 approaches for mean surface wind (m/s)
        
rho0    = 1.225; 
        % standard sea level air density (kg/m^3)

tao     = 20*86400; 
        % radiation relaxation time scale (s)

R       = (Pmid-PB)/(PT-PB); 
        % For calculating G in the model, defied in Hu and Randall 1995

W       = cp*Zmid*(Gamma_d-Gamma_m)/q00; 
        % For calculating G in the model, defied in Hu and Randall 1995


% Write out:
fprintf(fid, 'Surface Height (m) = %7.2f \n', Z0);
fprintf(fid, 'PBL Top Height (m) = %7.2f \n', ZB);  
fprintf(fid, '500-hPa Height (m) = %7.2f \n\n', Zmid);        
fprintf(fid, 'Surface Pressure (hPa) = %7.1f \n', P0./100);        
fprintf(fid, 'PBL Top Pressure (hPa) = %7.1f \n', PB./100);        
fprintf(fid, '500-hPa Pressure (hPa) = %7.1f \n', Pmid./100);        
fprintf(fid, 'Top Pressure (hPa)     = %7.1f \n\n', PT./100);
fprintf(fid, 'PBL Depth (hPa) = %5.1f \n\n', del_PM./100);        
fprintf(fid, 'SST (degC) = %5.2f \n\n', SST-273.15);
fprintf(fid, 'PBL Equilibrium Temp. (degC) 	   = %5.2f \n', TRM-273.15);
fprintf(fid, 'PBL Top Equilibrium Temp. (degC) = %5.2f \n\n', TRB-273.15);
fprintf(fid, 'Surface DSE. (J/kg) = %6.3e \n', S00);
fprintf(fid, 'PBL DSE. (J/kg)     = %6.3e \n', SRM);
fprintf(fid, 'PBL Top DSE. (J/kg) = %6.3e \n\n', SRB);
fprintf(fid, 'Surface Saturated q. (g/kg) = %5.2f \n\n', q00.*1000);
fprintf(fid, 'Pressure Scale for Moisture. (hPa) = %7.1f \n\n', Pq./100);
fprintf(fid, 'Precipitation Efficiency Factor = %3.1f \n\n', c);
fprintf(fid, 'Mean Surface Wind (m/s) = %4.1f \n\n', Vs);
fprintf(fid, 'Surface Air Density (kg/m^3) = %5.3f \n\n', rho0);
fprintf(fid, 'Radiation Relaxation Time (day) = %4.1f \n\n', tao./86400);
        
% =========================================================================
% 
% Set/Write-out Parameters: Integration Time.
% 
% =========================================================================

dt      = 60*1; 
        % time step (s)

tmax    = length(TimeID_ini:TimeID_end).*3.*3600; 
        % total integration time (s)

smax    = fix(tmax/dt); 
        % total integration steps

tday    = 0:dt/86400:tmax/86400; 
        % time-axis for plotting (daily)
        
% Write out:
fprintf(fid, 'Time Step (s)          = %5.0f \n', dt);
fprintf(fid, 'Time Integration (day) = %5.1f \n\n', tmax./86400);

% =========================================================================
% 
% Set/Write-out Initial Values (Initial Condition).
% 
% =========================================================================

TM(1)   = mean(Gan_Sonde_Full.T(Gan_Sonde_Full.P(:,1) >= PB,TimeID_ini)); 
        % initial mean temperature in mixed-layer (K)

TBt(1)  = interp1(Gan_Sonde_Full.Z(:,TimeID_ini), ...
                  Gan_Sonde_Full.T(:,TimeID_ini), ZB); 
        % initial temperature at top of mixed-layer (K)

SM(1)   = mean(Gan_Sonde_Full.DSE(Gan_Sonde_Full.P(:,1) >= PB,TimeID_ini));
        % initial mean dry static energy in mixed-layer (J/kg=m^2/s^2)

SBt(1)  = interp1(Gan_Sonde_Full.Z(:,TimeID_ini), ...
                  Gan_Sonde_Full.DSE(:,TimeID_ini), ZB); 
        % initial dry static energy at top of mixed-layer (J/kg=m^2/s^2)

qM(1)   = mean(Gan_Sonde_Full.q(Gan_Sonde_Full.P(:,1) >= PB,TimeID_ini)); 
        % initial water vapor mixing ratio in mixed-layer (kg/kg)

qBt(1)  = interp1(Gan_Sonde_Full.Z(:,TimeID_ini), ...
                  Gan_Sonde_Full.q(:,TimeID_ini), ZB); 
        % initial water vapor mixing ratio at top of mixed-layer (kg/kg)

% Write out:      
fprintf(fid, '============================================================= \n\n');
fprintf(fid, '  Initial Conditions \n\n');
fprintf(fid, '============================================================= \n\n');
fprintf(fid, 'Initial PBL Temp. (degC)     = %5.2f \n', TM-273.15);
fprintf(fid, 'Initial PBL Top Temp. (degC) = %5.2f \n\n', TBt-273.15);
fprintf(fid, 'Initail PBL DSE. (J/kg)      = %6.3e \n', SM);
fprintf(fid, 'Initial PBL Top DSE. (J/kg)  = %6.3e \n\n', SBt);
fprintf(fid, 'Initial PBL q. (g/kg)        = %5.2f \n\n', qM.*1000);
fprintf(fid, 'Initial PBL Top q. (g/kg)    = %5.2f \n\n', qBt.*1000);

fclose(fid);

% =========================================================================


