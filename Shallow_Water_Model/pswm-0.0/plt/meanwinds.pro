; IDL MEANWINDS.PRO
; Azimuthal mean wind profiles plotted at two different times
; Written: 07/05/08
; Eric A. Hendricks

CLOSE,1,2,3
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
nmult = 0
dum1 = 0.0d0
dum2 = 0.0d0
dum3 = 0.0d0

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
READ, var, PROMPT="Enter variable to plot (p,u,v,z,d,q,h): "
READ, lone, PROMPT="Enter beginning file number: "
READ, ltwo, PROMPT="Enter ending file number: "
if var EQ 'u' then begin 
    nmult = 1
    compstring = STRCOMPRESS('v', /REMOVE_ALL)
endif
if var EQ 'v' then begin 
    nmult = 1
    compstring = STRCOMPRESS('u', /REMOVE_ALL)
endif        
lone = FIX(lone)
ltwo = FIX(ltwo)
; TEMP FIX
lthr = 5
varstring = STRCOMPRESS(var, /REMOVE_ALL)

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
bs = STRCOMPRESS(STRING(lone), /REMOVE_ALL)
es = STRCOMPRESS(STRING(ltwo),   /REMOVE_ALL)
outfile = dir + infilene + 'azmean.' + bs + '.' + es + '.ps'
lstring1 = STRCOMPRESS(STRING(lone), /REMOVE_ALL)
lstring2 = STRCOMPRESS(STRING(ltwo), /REMOVE_ALL)

;
; OPEN AND READ THE OUTPUT FILE
; 

fieldone = DBLARR(nxd+1,nyd+1)
fieldtwo = DBLARR(nxd+1,nyd+1)
fieldonetmp = DBLARR(nxd+1,nyd+1)
fieldtwotmp = DBLARR(nxd+1,nyd+1)
field    = DBLARR(nxd+1,nyd+1)
field2   = DBLARR(nxd+1,nyd+1)
fieldz   = DBLARR(nxd+1,nyd+1)   ; for centering

; TEMP CODE
fieldthr = DBLARR(nxd+1,nyd+1)
fieldthrtmp = DBLARR(nxd+1,nyd+1)

FOR m = 0,2 DO BEGIN

IF m EQ 0 THEN l = lone
IF m EQ 1 THEN l = ltwo

IF m EQ 2 THEN l = lthr

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

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
PRINT, "** meanwinds.pro: opening and reading " + infile + "..."
OPENR, 1, infile, ERROR = err 
READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
READF, 1, xmin, xmax, ymin, ymax
READF, 1, fcor, beta, cval, timesec
if (var EQ 'u' OR var EQ 'v') then begin
   infile2 = dir + infilene + compstring + '.' + tstamp 
   OPENR, 2, infile2, ERROR = err 
   READF, 2, temp
   READF, 2, temp
   READF, 2, nx, ny
   READF, 2, xmin, xmax, ymin, ymax
   READF, 2, fcor, beta, cval, timesec
endif
infilez = dir + infilene + 'z' + '.' + tstamp 
OPENR, 3, infilez, ERROR = err 
READF, 3, temp
READF, 3, temp
READF, 3, nx, ny
READF, 3, xmin, xmax, ymin, ymax
READF, 3, fcor, beta, cval, timesec
if (err NE 0) then begin 
   PRINT, "** meanwinds.pro: file not found, exiting..."  
   break
endif

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
        if (var EQ 'u' OR var EQ 'v') then begin
           READF, 2, dum2
           field2(i,j) = dum2
        endif  
        READF, 3, dum3
        fieldz(i,j) = dum3
ENDFOR
ENDFOR
CLOSE, 1, 2, 3

if m EQ 0 then begin
   fieldone = field
   fieldonetmp = field2
endif
if m EQ 1 then begin
   fieldtwo = field
   fieldtwotmp = field2
endif
if m EQ 2 then begin
   fieldthr = field
   fieldthrtmp = field2
endif

;
; DETERMINE CENTER FROM VORTICITY
;

PRINT, "** meanwinds.pro: Determining center..."

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
IF m EQ 0 THEN BEGIN
  xc1 = xCenter_z
  yc1 = yCenter_z
ENDIF
IF m EQ 1 THEN BEGIN
  xc2 = xCenter_z
  yc2 = yCenter_z
ENDIF
; temp code
IF m EQ 2 THEN BEGIN
  xc3 = xCenter_z
  yc3 = yCenter_z
ENDIF
PRINT, "** meanwinds.pro: xCenter, yCenter indeces = ", xCenter_z, yCenter_z

ENDFOR

;
; AZIMUTHAL MEAN, PERTURBATION QUANTITIES ABOUT CENTER
;

PRINT, "** meanwinds.pro: Determining azimuthal mean diagnostics..."

nazm = 50
nr = FIX(nx/2.0)
deltar = deltax
r = FINDGEN(nr) * deltar

fieldone_azmean = DBLARR(nr)
fieldone_prime = DBLARR(nr,nazm)
fieldone_polar = DBLARR(nr,nazm)

fieldtwo_azmean = DBLARR(nr)
fieldtwo_prime = DBLARR(nr,nazm)
fieldtwo_polar = DBLARR(nr,nazm)

; TEMP CODE
fieldthr_azmean = DBLARR(nr)
fieldthr_prime = DBLARR(nr,nazm)
fieldthr_polar = DBLARR(nr,nazm)

vr = DBLARR(nx+1,ny+1)
vt = DBLARR(nx+1,ny+1)

if (var EQ 'u') then begin
CARTESIANTOCYLIND2D, fieldone, fieldonetmp, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc1, $
		       yc1

fieldone = vr
CARTESIANTOCYLIND2D, fieldtwo, fieldtwotmp, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc2, $
		       yc2
		     
fieldtwo = vr
; TEMP CODE
CARTESIANTOCYLIND2D, fieldthr, fieldthrtmp, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc3, $
		       yc3
		     
fieldthr = vr
endif
if (var EQ 'v') then begin
CARTESIANTOCYLIND2D, fieldonetmp, fieldone, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc1, $
		       yc1

fieldone = vt	
CARTESIANTOCYLIND2D, fieldtwotmp, fieldtwo, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc2, $
		       yc2     
fieldtwo = vt
; TEMP CODE
CARTESIANTOCYLIND2D, fieldthrtmp, fieldthr, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc3, $
		       yc3     
fieldthr = vt
endif


INTERPTOCYLIND2D, fieldone, fieldone_polar, $
                nr, nazm, $
                xc1, $
                yc1, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldtwo, fieldtwo_polar, $
                nr, nazm, $
                xc2, $
                yc2, $
                xPositions, $
                yPositions
; TEMP CODE
INTERPTOCYLIND2D, fieldthr, fieldthr_polar, $
                nr, nazm, $
                xc3, $
                yc3, $
                xPositions, $
                yPositions


fieldone_azmean = AZIMUTHALMEAN2D(fieldone_polar, $
                        perturbation = fieldone_prime)
fieldtwo_azmean = AZIMUTHALMEAN2D(fieldtwo_polar, $
                        perturbation = fieldtwo_prime)
; TEMP CODE
fieldthr_azmean = AZIMUTHALMEAN2D(fieldthr_polar, $
                        perturbation = fieldthr_prime)



;
;  PLOTTING
;

fieldall_azmean = DBLARR(nr,3)
fieldall_azmean(*,0) = fieldone_azmean
fieldall_azmean(*,1) = fieldtwo_azmean
; TEMP CODE
fieldall_azmean(*,2) = fieldthr_azmean
xaxis = "r (km)"
yaxis = "v (m s!E-1!N)"
titles = ['',xaxis,yaxis]
bounds = [0,300,min(fieldall_azmean),max(fieldall_azmean),timesec/3600.0d0] 

; temp code for angular velocity
omega = DBLARR(nr,3)
FOR i = 0, 2 DO BEGIN
omega(*,i) = fieldall_azmean(*,i) / r
ENDFOR
omega(0:9,0) = omega(10,0)
omega(0:4,1) = omega(5,1)
omega(0:4,2) = omega(5,2)

r = r / 1000.0d0
onedoverplotter, fieldall_azmean, r, TITLES, BOUNDS, outfile
titles = ['',xaxis, '!7x!3 (s!E-1!N)']
bounds = [0,300,0,1.0d-2,timesec/3600.0d0] 
;onedoverplotter, omega, r, TITLES, BOUNDS, outfile
r = r * 1000.0d0


; kinetic energy
keone = 0.0
ketwo = 0.0
kethr = 0.0
FOR i = 0, nr-2 DO BEGIN
   dm = 3.14159 * (r(i+1)*r(i+1) - r(i) * r(i)) * 1.13 * 4284.0
   keone = keone + dm * (fieldone_azmean(i) * fieldone_azmean(i) / 2.0)
   ketwo = ketwo + dm * (fieldtwo_azmean(i) * fieldtwo_azmean(i) / 2.0)
   kethr = kethr + dm * (fieldthr_azmean(i) * fieldthr_azmean(i) / 2.0)
ENDFOR
PRINT, "KE (J) 1,2,3 = ", keone, ketwo, kethr

END


