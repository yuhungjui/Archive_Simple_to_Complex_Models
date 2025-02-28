clear;close all;clc;
addpath(genpath('../../../../../Work/MatLAB_TOOLS/'));
tic;
% =========================================================================
% 
% Reading in the high temporal resolution Revelle MET data (includes SSST),
% average it into 3-hourly, then output.
% 
% =========================================================================
%% Load data:
% Set starting time:
yd0 = datenum(datestr('01-Jan-2011 00:00:00'))-1;

% Load data 1:
load ../Data/ESRL-PSD_Revelle_Met/Revelle10minutesLeg1_r3.mat

yday_1      = datestr(yday + yd0,30); % UTC
U10_1       = U10;      % Wind speed (m/s) relative to earth adjusted to 10 m
Pair10_1    = Pair10;   % Pressure (mb) adjusted to 10 m
RH10_1      = RH10;     % Relative humidity(%) adjusted to 10 m
T10_1       = T10;      % Temperature (C) adjusted to 10 m
Tsea_1      = Tsea;     % Near surface sea temperature (C) from Sea snake with warm layer correction
SST_1       = SST;      % Sea surface (skin) temperature (C) from Tsea minus cool skin
Q10_1       = Q10;      % Specific humidity (g/Kg) adjusted to 10 m
Qsea_1      = Qsea;     % Specific humidity (g/Kg) 'near' ocean surface from sea snake
SSQ_1       = SSQ;      % Sea surface specific humidity (g/Kg) from Qsea minus cool skin
shf_1       = shf;      % Sensible heat flux (W/m2)
lhf_1       = lhf;      % Latent heat flux (W/m2)

% Load data 2:
load ../Data/ESRL-PSD_Revelle_Met/Revelle10minutesLeg2_r3.mat

yday_2      = datestr(yday + yd0,30); % UTC
U10_2       = U10;      % Wind speed (m/s) relative to earth adjusted to 10 m
Pair10_2    = Pair10;   % Pressure (mb) adjusted to 10 m
RH10_2      = RH10;     % Relative humidity(%) adjusted to 10 m
T10_2       = T10;      % Temperature (C) adjusted to 10 m
Tsea_2      = Tsea;     % Near surface sea temperature (C) from Sea snake with warm layer correction
SST_2       = SST;      % Sea surface (skin) temperature (C) from Tsea minus cool skin
Q10_2       = Q10;      % Specific humidity (g/Kg) adjusted to 10 m
Qsea_2      = Qsea;     % Specific humidity (g/Kg) 'near' ocean surface from sea snake
SSQ_2       = SSQ;      % Sea surface specific humidity (g/Kg) from Qsea minus cool skin
shf_2       = shf;      % Sensible heat flux (W/m2)
lhf_2       = lhf;      % Latent heat flux (W/m2)

% Load data 3:
load ../Data/ESRL-PSD_Revelle_Met/Revelle10minutesLeg3_r3.mat

yday_3      = datestr(yday + yd0,30); % UTC
U10_3       = U10;      % Wind speed (m/s) relative to earth adjusted to 10 m
Pair10_3    = Pair10;   % Pressure (mb) adjusted to 10 m
RH10_3      = RH10;     % Relative humidity(%) adjusted to 10 m
T10_3       = T10;      % Temperature (C) adjusted to 10 m
Tsea_3      = Tsea;     % Near surface sea temperature (C) from Sea snake with warm layer correction
SST_3       = SST;      % Sea surface (skin) temperature (C) from Tsea minus cool skin
Q10_3       = Q10;      % Specific humidity (g/Kg) adjusted to 10 m
Qsea_3      = Qsea;     % Specific humidity (g/Kg) 'near' ocean surface from sea snake
SSQ_3       = SSQ;      % Sea surface specific humidity (g/Kg) from Qsea minus cool skin
shf_3       = shf;      % Sensible heat flux (W/m2)
lhf_3       = lhf;      % Latent heat flux (W/m2)

% Load data 4:
load ../Data/ESRL-PSD_Revelle_Met/Revelle10minutesLeg4_r3.mat

yday_4      = datestr(yday + yd0,30); % UTC
U10_4       = U10;      % Wind speed (m/s) relative to earth adjusted to 10 m
Pair10_4    = Pair10;   % Pressure (mb) adjusted to 10 m
RH10_4      = RH10;     % Relative humidity(%) adjusted to 10 m
T10_4       = T10;      % Temperature (C) adjusted to 10 m
Tsea_4      = Tsea;     % Near surface sea temperature (C) from Sea snake with warm layer correction
SST_4       = SST;      % Sea surface (skin) temperature (C) from Tsea minus cool skin
Q10_4       = Q10;      % Specific humidity (g/Kg) adjusted to 10 m
Qsea_4      = Qsea;     % Specific humidity (g/Kg) 'near' ocean surface from sea snake
SSQ_4       = SSQ;      % Sea surface specific humidity (g/Kg) from Qsea minus cool skin
shf_4       = shf;      % Sensible heat flux (W/m2)
lhf_4       = lhf;      % Latent heat flux (W/m2)

% Organizing loaded data:
clearvars -EXCEPT *_1 *_2 *_3 *_4
datatime    = [yday_1;yday_2;yday_3;yday_4];            % UTC
U10         = [U10_1,U10_2,U10_3,U10_4];                % Wind speed (m/s) relative to earth adjusted to 10 m
P10         = [Pair10_1,Pair10_2,Pair10_3,Pair10_4];    % Pressure (mb) adjusted to 10 m
RH10        = [RH10_1,RH10_2,RH10_3,RH10_4];            % Relative humidity(%) adjusted to 10 m
T10         = [T10_1,T10_2,T10_3,T10_4];                % Temperature (C) adjusted to 10 m
Tsea        = [Tsea_1,Tsea_2,Tsea_3,Tsea_4];            % Near surface sea temperature (C) from Sea snake with warm layer correction
SST         = [SST_1,SST_2,SST_3,SST_4];                % Sea surface (skin) temperature (C) from Tsea minus cool skin
Q10         = [Q10_1,Q10_2,Q10_3,Q10_4];                % Specific humidity (g/Kg) adjusted to 10 m
Qsea        = [Qsea_1,Qsea_2,Qsea_3,Qsea_4];            % Specific humidity (g/Kg) 'near' ocean surface from sea snake
SSQ         = [SSQ_1,SSQ_2,SSQ_3,SSQ_4];                % Sea surface specific humidity (g/Kg) from Qsea minus cool skin
shf         = [shf_1,shf_2,shf_3,shf_4];                % Sensible heat flux (W/m2)
lhf         = [lhf_1,lhf_2,lhf_3,lhf_4];                % Latent heat flux (W/m2)

clear *_1 *_2 *_3 *_4

% Set time table:
Set_TimeTable_daily_3hrly;

% =========================================================================
%% Calculate the 3-hourly averaged data:
Revelle_MET.Time = TT;

for ti = 1:length(TT)
    ti
    
    % Find the id for every 3 hour time:
    TTid_3hr = find( ( str2num(datatime(:,5:6)) == TT(ti,1) ) ...
                   & ( str2num(datatime(:,7:8)) == TT(ti,2) ) ...
                   & ( str2num(datatime(:,10:11)) == TT(ti,3) ) ...
                   & ( str2num(datatime(:,12:13)) == 0 ) );
    
    % Calculate the 3-hourly averages of variables and save:
    Revelle_MET.U10(ti)     = mean(U10(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.P10(ti)     = mean(P10(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.RH10(ti)    = mean(RH10(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.T10(ti)     = mean(T10(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.Tsea(ti)    = mean(Tsea(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.SST(ti)     = mean(SST(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.q10(ti)     = mean(Q10(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.qsea(ti)    = mean(Qsea(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.SSQ(ti)     = mean(SSQ(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.shf(ti)     = mean(shf(TTid_3hr-9:TTid_3hr+9));
    Revelle_MET.lhf(ti)     = mean(lhf(TTid_3hr-9:TTid_3hr+9));
               
end

% =========================================================================
%% Save the data:
save ./DATA_Revelle_MET.mat Revelle_MET 

% =========================================================================
%% Display runnung time.
time_cost = toc;
disp(time_cost);
