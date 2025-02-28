Yday       Decimal yearday (UTC)
           To convert to matlab time format use:
           yd0=datenum(datestr('01-Jan-2011 00:00:00'))-1;
           tmatlab=yday+yd0;

 Lat       Latitude (deg)
 Lon       Longitude (deg)
 SOG       Speed over ground (m/s)
 COG       Course over ground (deg)
 Heading   Ship's heading (deg)
 Cspd      Current speed (m/s)
 Cdir      Current direction (deg) from
 U10       Wind speed (m/s) relative to earth adjusted to 10 m
 Wdir      Wind direction (deg) from relative to earth
 Ur10      Wind speed (m/s) relative to water adjusted to 10 m
 WdirR     Wind direction (deg) from relative to water
 Pair10    Pressure (mb) adjusted to 10 m
 RH10      Relative humidity(%) adjusted to 10 m
 T10       Temperature (C) adjusted to 10 m
 Tsea      Near surface sea temperature (C) from Sea snake with warm layer correction     
 SST       Sea surface (skin) temperature (C) from Tsea minus cool skin
 Q10       Specific humidity (g/Kg) adjusted to 10 m
 Qsea      Specific humidity (g/Kg) 'near' ocean surface from sea snake
 SSQ       Sea surface specific humidity (g/Kg) from Qsea minus cool skin
 stress    Surface stress (N/m2) measured relative to water
 shf       Sensible heat flux (W/m2)
 lhf       Latent heat flux (W/m2)
 rhf       Sensible heat flux from rain (W/m2)
 Solarup   Reflected solar (W/m2) estimated from Payne (1972)
 Solardn   Measured downwelling solar (W/m2)
 IRup      Upwelling IR (W/m2) computed from SST 
 IRdn      Measured downwelling IR (W/m2)
 E         Evaporation rate (mm/hr)
 P         Precipitation rate (mm/hr)
 Evap      Accumulated evaporation for Leg (mm)
 Precip    Accumulated precipitation for Leg (mm)
 Interped  1=data interpolated due to poor relative winds 
           0=no interpolation
 TseaTSG   Sea temperature from the Thermosalinograph at 5 m depth (C)
 SalTSG    Salinity from the Thermosalinograph at 5 m depth (PSU)
 T02       Temperature (C) adjusted to 2 m
 Q02       Specific humidity (g/Kg) adjusted to 2 m
 RH02      Relative humidity(%) adjusted to 2 m
 sigH      Significant wave height (m)
 cp        Phase speed of dominant waves (m/s)
 sigDir    Direction of dominant wave (deg) - placeholder
 dSolar    Measured diffuse radiation (W/m2)
 zenith    Solar zenith angle (deg)
 Smax      Solar radiation at surface with no atmosphere (W/m2)
 Sclr      Modeled clear sky radiation at ocean surface (W/m2)
   
12/10/11

Notes:  Fluxes are defined as negative downward and positive upwards.  For
example, the net heat flux is defined as:

Qnet = Solarup+Solardn+IRup+IRdn+lhf+shf+rhf

Qnet<0 is heating ocean

The wind and current directions are in meteorological convention (i.e., 
direction from).

Tair is taken from the calibrated PSD and UConn aspirated air temperature 
sensors on the bow mast.  These were least affected by solar heating. Qair 
and Pair are computed the calibrated UConn RH/T/P sensors on the on the bow 
mast.  Q is less sensitive to solar heating as long as the temperature and
RH are measured simultaneously.  RH is reconstructed from the Q, aspirated 
Tair and P measurements to remove the effects of solar heating.  The sonic
anemometers on the bow mast are used to measure the wind speed and 
direction.  Relative wind speed is taken into consideration to minimize
flow distortion.  Tsea is primarily measured by the sea snake with a few 
values provide by the IR radiometer deployed by LDEO.  SST is estimated 
after correction for cool skin and this accounts for the difference 
between Tsea and SST.  Similar corrections are applied to SSQ from Qsea.  
Solardn is provided by the ship's pyranometer on the top of the forward 
mast. IRdn represents an average of the gyrostabilized PSD purgeometer on 
top of its van and the ship's purgeometer on the top of the forward mast.  
The ship's purgeometer was first corrected for the effects of solar 
heating. Solarup is taken from a commonly used parameterization for 
surface albedo of the ocean (Payne, 1972).  IRup was derived from the
SST measurements using the COARE 3.0 algorithm.  The bulk fluxes of stress
(momentum), sensible heat, latent heat and sensible heat due to rain were
provided by the COARE 3.0 algorithm.  The COARE 3.0 algorithm was also 
used to compute the 10-m values of wind speed, temperature and humidity. 
SOG, COG and Gyro were taken from the PSD GPS compass.  These were used 
to compute the wind speed relative to earth.  Surface currents are measured 
by the ship's ADCP and have been QCed by OSU Ocean Mixing group.  These 
were used to compute the wind speed relative to water. The wind speed 
relative to water are used to compute the fluxes. 

10/08/12

New calibrations were applied to the RH measurements from the UConn and
ETL sensors.  The calibration raised the RH values 1-2%, which modifies 
slightly the flux estimates and any variabiles adjusted using MO similarity.

Values of temperature, specific humidity and relative humidity adjusted to 
2-m were added to the matlab and ASCII files. 

07/23/13

New calibrations were applied to the RH/Tair measurements from the UConn and
ETL sensors using Vaisala factory calibration and a RH/T calibration chamber
at UConn.  This results in small changes to the previous revision.  

The sea-snake measurements were reduced by 0.058 based on a comparison with 
the OSU T-chain SBE device.

The TOGA-COARE warm layer correction is applied to the seasnake (causing a 
small increase during the day.  The depth of the sea-snake was set to 5 cm.
This value is then used to compute the fluxes using the TOGA-COARE 3.5 
algorithm described by Edson et al., 2013: On the exchange of momentum over 
the open ocean,” J. Phys. Oceanogr., 43, 1589-1610.  Therefore, this value of
Tsea can be used in the COARE algorithm without the need to run the warm 
layer code.  However, the warm layer correction is saved in the matlab files 
for investigators that want to convert Tsea back to it measured value.  This 
variable is called dsea in the matlab files and the measured value equals: 

Tsea_measured=Tsea-dsea.

Measurements of the sea temperature and salinity from the thermosalinograph
(TSG) are provided.   The intake for these measurements is reported as 5 m.
The TSG temperature measurements were reduced by 0.05 C based on a comparison 
with the calibrated sea-snake.

Estimates of the significant wave height and phase speed of the dominant wave
derived from laser altimeter measurementsare provided.   The 1 and 10 minutes 
values are determined from 1 hour averages, so the variability is
representative of that time scale.  Values that did not contain enough points
resolve 2 second waves (on average) were removed and interpolated through.
The phase speed of the dominant wave was determined from the frequency at 
the spectral peak.  The search for the peak was limited to frequencies 
between 0.01 and 0.5 Hz or 100 down to 2 second waves.   These values are 
preliminary and further refinements are expected with the inclusion of 
WaMOS data from the ship's radar.   A placeholder has been included for
the direction of the dominant waves once that becomes available.  For now,
this value is given by NaN.  Also note that the laser altimeter was not
operational during Leg 1 and those values are also given by NaN. 

Measurements of the diffuse solar radiation from a sensor deployed on the top
of the bow mast by NCAR are given.

The solar zenith angle is provided, which takes into account the hour angle 
for the local solar time, inclination and latitude. 

Model estimates of the solar radiation at the surface in the absence of an 
atmosphere is provided using a solar constant of 1367 W/m2 and allowing for 
changes due to variations in the distance between Earth and Sun.

Model estimates of the clear sky solar radiation (i.e. and atmosphere with 
no clouds) are provided using a parameterization technique from M. Iqbal 
in "Physical Climatology for Solar and Wind Energy", 1988, World Scientific, 
pp. 196-242.  The model accounts for absortion and scattering by gases, 
aerosols, ozone and water vapor.  Model coefficients were tuned using 
observations on the clearest days.

8/20/13

The 10-m wind speed was incorrectly given as the 2-m wind speed in revision 2.
This error is corrected in revision 3, i.e., the 10-m values are correctly 
reported in Ur10 and U10.

 