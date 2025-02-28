; IDL WAVEACTIVITY.PRO
; Main Program to compute wave activity and flux
; according to Guinn and Schubert (1993)
; Written: 05/29/08
; Eric A. Hendricks

dir = '../run/e7o20hres/'
var = ''
temp = ''
outfile = ''
infilene = 'e7o20'
infile2 = ''
varstring = ''
xmin = 0.0d0
xmax = 0.0d0
fcor = 0.0d0
beta = 0.0d0
timesec = 0.0d0
ymin = 0.0d0
ymax = 0.0d0

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
READ, lbegin, PROMPT="Enter beginning file number: "
READ, lend, PROMPT="Enter ending file number: "
lbegin = FIX(lbegin)
lend = FIX(lend)

nxd = 1024
nyd = 1024
cval = 205.0    ; m/s
grav = 9.81

;
; BEGIN FILE LOOP
;

FOR l = lbegin, lend DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

if (l lt 100) then begin

if (l lt 10) then begin
     tstamp = 't0' + lstring
     ostamp = '00' + lstring
     outfile =  dir+infilene+'2dwa.'+ostamp+'.ps'
     outfile1 = dir+infilene+'dAdt.'+ostamp+'.ps'
  endif else begin
     tstamp = 't' + lstring
     ostamp = '0' + lstring
     outfile = dir+infilene+'2dwa.'+ostamp+'.ps'
     outfile1 = dir+infilene+'dAdt.'+ostamp+'.ps'
  endelse
endif else begin
  tstamp = lstring
  outfile = dir+infilene+'2dwa.'+tstamp+'.ps'
  outfile1 = dir+infilene+'dAdt.'+tstamp+'.ps'
endelse

;
; OPEN AND READ THE OUTPUT FILE
; 

fieldh = DBLARR(nxd+1,nyd+1)
fieldu = DBLARR(nxd+1,nyd+1)
fieldv = DBLARR(nxd+1,nyd+1)
fieldz = DBLARR(nxd+1,nyd+1)

FOR m = 0, 3 DO BEGIN

IF m EQ 0 THEN varstring = 'p'
IF m EQ 1 THEN varstring = 'u'
IF m EQ 2 THEN varstring = 'v'
IF m EQ 3 THEN varstring = 'z'
infile = dir + infilene + varstring + '.' + tstamp 
PRINT, "** waveactivity.pro: opening and reading " + infile + "..."
OPENR, 1, infile, ERROR = err 
if (err NE 0) then begin 
   PRINT, "** waveactivity.pro: file not found, exiting..."  
   break
endif
READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
READF, 1, xmin, xmax, ymin, ymax
READF, 1, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
field = DBLARR(nx+1,ny+1)
y = DBLARR(ny+1)
x = DBLARR(nx+1)
xc = fix(nx / 2.0)
yc = fix(ny / 2.0)
timehour = STRCOMPRESS(STRING(timesec/3600.), /REMOVE_ALL)
deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR j = 0, ny DO BEGIN
y(j) = j*deltay - 0.5d*(ymax-ymin)
FOR i = 0, nx DO BEGIN
x(i) = i*deltax - 0.5d*(xmax-xmin)
	READF, 1, dum1
	field(i,j) = dum1  
ENDFOR
ENDFOR
CLOSE, 1

if m EQ 0 then fieldh = field
if m EQ 1 then fieldu = field
if m EQ 2 then fieldv = field
if m EQ 3 then fieldz = field

ENDFOR

;
; CONVERT P to H FOR DIAGNOSTICS
;

fieldh = (fieldh*cval + cval^2) / grav

;
; DETERMINE CENTER FROM VORTICITY
;

PRINT, "** waveactivity.pro: Determining center..."

boxsize = 40
zeta_high = 0.0                ; max zeta 
zeta_mn = FLTARR(nx+1,ny+1)    ; vorticity averaged                          
xCenter_z = 0                  ; integer index of x-center
yCenter_z = 0                  ; integer index of y-center

FOR j = 0, ny-1-boxsize DO BEGIN
FOR i = 0, nx-1-boxsize DO BEGIN

     zeta_mn(i+boxsize/2,j+boxsize/2) = $
                      MEAN(fieldz(i:i+boxsize,j:j+boxsize))

ENDFOR
ENDFOR

tmpvar = 0
zeta_high = MAX(zeta_mn,tmpvar,/NAN)
yCenter_z = (tmpvar/FIX(nx+1))
xCenter_z = tmpvar - ((tmpvar/FIX(nx+1))*FIX(nx+1))  ;Note: IDL drops decimals
xc = xCenter_z
yc = yCenter_z

PRINT, "** waveactivity.pro: xCenter, yCenter indeces = ", xCenter_z, yCenter_z

;
; AZIMUTHAL MEAN, PERTURBATION QUANTITIES ABOUT CENTER
;

PRINT, "** waveactivity.pro: Determining azimuthal mean diagnostics..."

nazm = 50
nr = FIX(nx/2.0)
deltar = deltax
r = FINDGEN(nr) * deltar

fieldu_azmean = DBLARR(nr)
fieldu_prime = DBLARR(nr,nazm)
fieldu_polar = DBLARR(nr,nazm)

fieldv_azmean = DBLARR(nr)
fieldv_prime = DBLARR(nr,nazm)
fieldv_polar = DBLARR(nr,nazm)

fieldh_azmean = DBLARR(nr)
fieldh_prime = DBLARR(nr,nazm)
fieldh_polar = DBLARR(nr,nazm)

fieldz_azmean = DBLARR(nr)
fieldz_prime = DBLARR(nr,nazm)
fieldz_polar = DBLARR(nr,nazm)

vr = DBLARR(nx+1,ny+1)
vt = DBLARR(nx+1,ny+1)

CARTESIANTOCYLIND2D, fieldu, fieldv, $
                   vr, vt, $
                   deltax, deltay, nazm, $
		       xc, $
		       yc

INTERPTOCYLIND2D, vr, fieldu_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, vt, fieldv_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldh, fieldh_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldz, fieldz_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

fieldu_azmean = AZIMUTHALMEAN2D(fieldu_polar, $
                        perturbation = fieldu_prime)
fieldv_azmean = AZIMUTHALMEAN2D(fieldv_polar, $
                        perturbation = fieldv_prime)
fieldh_azmean = AZIMUTHALMEAN2D(fieldh_polar, $
                        perturbation = fieldh_prime)
fieldz_azmean = AZIMUTHALMEAN2D(fieldz_polar, $
                        perturbation = fieldz_prime)

;
; WAVE ACTIVITY DIAGNOSTICS
;

wavea_polar = DBLARR(nr,nazm)
mhbsr = DBLARR(nr)
mhtor = DBLARR(nr,nazm)

term1 = DBLARR(nr,nazm)
term2 = DBLARR(nr,nazm)

rhbar = DBLARR(nr)
rhbar = fieldh_azmean * r

mhbsr(1) = rhbar(1)*deltar
FOR i = 2, nr-1 DO BEGIN 
   integral1d, rhbar, $
                0, i, $
                deltar, $
                var_integrated
   mhbsr(i) = var_integrated
ENDFOR

FOR j = 0, nazm-1 DO BEGIN
fieldh_polar_one = fieldh_polar(*,j)
fieldh_polar_tmp = REFORM(fieldh)
rhbar = DBLARR(nr)
rhbar = fieldh_polar_tmp * r

mhtor(1,*) = rhbar(1)*deltar
FOR i = 2, nr-1 DO BEGIN   
   integral1d, rhbar, $
                0, i, $
                deltar, $
                var_integrated
   mhtor(i,j) = var_integrated
ENDFOR
ENDFOR

; COMPUTE THE POTENTIAL VORTICITY  

fieldp_azmean = DBLARR(nr)
fieldp_prime = DBLARR(nr,nazm)
fieldp_polar = DBLARR(nr,nazm)
fieldp_polar = ( fcor + fieldz_polar ) / ( fieldh_polar )
fieldp_azmean = AZIMUTHALMEAN2D(fieldp_polar, $
                        perturbation = fieldp_prime)

FOR j = 0, nazm-1 DO BEGIN
FOR i = 0, nr-1 DO BEGIN   
   term2(i,j) = ( mhtor(i,j) - mhbsr(i) ) * fieldh_polar(i,j) * $
                 fieldp_prime(i,j)  
ENDFOR
ENDFOR

FOR j = 0, nazm-1 DO BEGIN
FOR i = 0, nr-1 DO BEGIN   
     term1(i,j) = r(i) * fieldv_prime(i,j) * fieldh_prime(i,j)
     wavea_polar(i,j) = term1(i,j) + term2(i,j) 
ENDFOR
ENDFOR

wavea_azmean = DBLARR(nr)
wavea_prime  = DBLARR(nr,nazm)
wavea_azmean = AZIMUTHALMEAN2D(wavea_polar, perturbation = wavea_prime)

wavea_t   = DBLARR(nr)       ; local time rate of change of wave activity
wavea_rad = DBLARR(nr)

u_wavea_azmean = DBLARR(nr)
u_wavea_prime  = DBLARR(nr,nazm)
u_wavea_azmean = AZIMUTHALMEAN2D(wavea_polar, perturbation = wavea_prime)

ue_ve_azmean  = DBLARR(nr)
ue_ve_azmean  = DBLARR(nr,nazm)
ue_ve_azmean  = AZIMUTHALMEAN2D(wavea_polar, perturbation = wavea_prime)

FOR i = 0, nr-1 DO BEGIN
  wavea_rad(i) = r(i) * ( u_wavea_azmean(i) + $
                 fieldh_azmean(i) * r(i) * ue_ve_azmean(i)  )

ENDFOR

; differentiate with radius

FOR i = 1, nr-2 DO BEGIN
  wavea_t(i) = ( -1.0d0 / r(i) ) * ( wavea_rad(i+1) - wavea_rad(i-1) ) $
               / (2.0d0 * deltar)
ENDFOR
wavea_t(0) = 2.0d0 * wavea_t(1) - wavea_t(2)
wavea_t(nr-1) = 2.0d0 * wavea_t(nr-2) - wavea_t(nr-3)

PCONTOUR_RTHETA1, wavea_polar(0:50,*), 51, nazm, deltar, $
      8, outfile

xaxis = "r (km)"
yaxis = "dAdt"
titles = ["dA/dt (m^3/s^2)",xaxis,yaxis]
r = r / 1000.0d0
bounds = [0,100.0,min(wavea_t),max(wavea_t),timesec/3600.0d0] 
onedplotter, wavea_t, r, TITLES, BOUNDS, outfile1
r = r * 1000.0d0

ENDFOR

END


