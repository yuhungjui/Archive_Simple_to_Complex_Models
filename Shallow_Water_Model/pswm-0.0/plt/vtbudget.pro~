; IDL VTBUDGET.PRO
; Tangential momentum budget for shallow water equations
; Written: 06/29/08
; Eric A. Hendricks

dir = '/home/eric/scratch3/modelruns/pswm-0.0/e7o20hres/'
infilene = 'e7o20'
var = ''
temp = ''
outfile = STRARR(5)
outfile1 = STRARR(5)
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
deltat = 180.0 ; s

;
; DEFINE THE CONTOUR LEVELS
;

fac = 1000.0
clevels_per_h = 0.1*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_h  = "h' (m)"
clevels_per_z = 0.0001*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_z  = "!7f!3' (x 10!E-3!N s!E-1!N)"
clevels_per_d = 0.00001*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_d  = "!7d!3' (x 10!E-4!N s!E-1!N)"
clevels_per_u = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_u  = "u' (m s!E-1!N))"
clevels_per_v = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_v  = "v' (m s!E-1!N))"

;
; AZIMUTHAL MEAN ARRAYS
;

count = 0
nrd = FIX(FLOAT(nxd) / 2.0)
nt = lend - lbegin + 1
bs = STRCOMPRESS(STRING(lbegin), /REMOVE_ALL)
es = STRCOMPRESS(STRING(lend),   /REMOVE_ALL)
outfile = dir + infilene + 'vtbud.' + bs + '.' + es + '.ps'
uz_p_tend = DBLARR(nrd,nt)
un_m_tend = DBLARR(nrd,nt)
fieldv_azm_beg = DBLARR(nrd)
fieldv_azm_end = DBLARR(nrd)

;
; BEGIN FILE LOOP
;
FOR l = lbegin, lend DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

;
; OPEN AND READ THE OUTPUT FILE
; 

fieldu = DBLARR(nxd+1,nyd+1)
fieldv = DBLARR(nxd+1,nyd+1)
fieldz = DBLARR(nxd+1,nyd+1)

FOR m = 0, 2 DO BEGIN

IF m EQ 0 THEN varstring = 'u'
IF m EQ 1 THEN varstring = 'v'
IF m EQ 2 THEN varstring = 'z'

if (l lt 100) then begin
  if (l lt 10) then begin
     tstamp = 't0' + lstring
     ostamp = '00' + lstring
  endif else begin
     tstamp = 't' + lstring
     ostamp = '0' + lstring
  endelse
endif else begin
  tstamp = lstring
endelse

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

if m EQ 0 then fieldu = field
if m EQ 1 then fieldv = field
if m EQ 2 then fieldz = field

ENDFOR
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
fieldz_azmean = AZIMUTHALMEAN2D(fieldz_polar, $
                        perturbation = fieldz_prime)
;
; COMPUTE u'z'_bar
;

fielduz_prime2 = DBLARR(nr,nazm)
fielduz_azmean = DBLARR(nr)
fielduz_prime = fieldu_prime * fieldz_prime
fielduz_azmean = AZIMUTHALMEAN2D(fielduz_prime, $
                        perturbation = fielduz_prime2)

;
; STORE TENDENCIES IN ARRAYS
;
uz_p_tend(*,count) = fielduz_azmean
un_m_tend(*,count) = fieldu_azmean * (fieldz_azmean + fcor)
IF l EQ lbegin THEN fieldv_azm_beg = fieldv_azmean
IF l EQ lend THEN fieldv_azm_end = fieldv_azmean
count = count + 1

ENDFOR

;
; COMPUTE THE TANGENTIAL MOMENTUM BUDGET
;
eddyterm = DBLARR(nr)
FOR j = 0, nr-1 DO BEGIN
ref1 = uz_p_tend(j,*)
uz_p_tend_ref = REFORM(ref1)
eddyterm(1) = uz_p_tend_ref(1)*deltat
integral1d, uz_p_tend_ref, $
                2, nt-1, $
                deltat, $
                var_integrated
eddyterm(j) = var_integrated
ENDFOR

meanterm = DBLARR(nr)
FOR j = 0, nr-1 DO BEGIN
ref1 = un_m_tend(j,*)
un_m_tend_ref = REFORM(ref1)
meanterm(1) = un_m_tend_ref(1)*deltat
integral1d, un_m_tend_ref, $
                2, nt-1, $
                deltat, $
                var_integrated
meanterm(j) = var_integrated
ENDFOR

allterms = DBLARR(nr,4)
allterms(*,0) = meanterm
allterms(*,1) = eddyterm
allterms(*,2) = meanterm + eddyterm
allterms(*,3) = fieldv_azm_end - fieldv_azm_beg 

xaxis = "r (km)"
yaxis = "vtbudget"
titles = ['',xaxis,yaxis]
r = r / 1000.0d0
bounds = [0,100.0,min(allterms),max(allterms),timesec/3600.0d0] 
onedoverplotter, allterms, r, TITLES, BOUNDS, outfile
r = r * 1000.0d0



END


