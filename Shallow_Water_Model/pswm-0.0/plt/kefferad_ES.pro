PRINT, 'Running IDL kefferad_ES.pro'
PRINT, 'Effective Diffusivity Routine'
PRINT, 'for Pseudospectral Model Output'
PRINT, 'Equivalent Radius Version'

;  This program computes it the way Emily Shuckburgh suggested
;  Slightly different than my way.

dir = ''
var = ''
outfile = ''
default = ''
READ, dir, PROMPT="Enter output directory path:"
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)
READ, default, PROMPT="Test case (1-yes, 0-no):"
READ, nerad, PROMPT="Enter number of equivalent radius points:"
READ, nu, PROMPT="Enter the model diffusion coefficient:"
READ, lbegin, PROMPT="Enter the begining file index:"
nerad = FIX(nerad)
nu = double(nu)
lbegin = FIX(lbegin)

CLOSE, 1, 2, 3, 15

temp = ''
nfiles = 500
pi = 3.14159d
g = 9.81d
len = 600000.0d
half = 0.5d0
two = 2.0d0
one = 1.0d0
; regular viscosity for run
charzeta = 0.001982d
tracer = 1

ediff_erad_oneplot = FLTARR(nerad,nfiles)
lcount = 0

FOR l = lbegin, nfiles DO BEGIN

lcount = lcount + 1

;
; Set the time stamp
;
lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)
if (l lt 100) then begin
   if (l lt 10) then begin
      tstamp = '00' + lstring
   endif else begin
      tstamp = '0' + lstring
   endelse
endif else begin
   tstamp = lstring
endelse

;
; Open the correct input and output files
;
if (tracer EQ 0) then begin
   OPENR, 1, dir + infile + tstamp, ERROR = err   
   outfile = dir+infile+tstamp+'.ps' 
endif else begin
   OPENR, 1, dir + infile + 'tr-' + tstamp, ERROR = err
   outfile = dir+infile+'keff.tr-'+tstamp+'.ps'
endelse
if (err NE 0) then break

READF, 1, temp
READF, 1, temp
READF, 1, time
READF, 1, nx, ny

nx = fix(nx)
ny = fix(ny)
zeta = DBLARR(nx,ny)
y = DBLARR(ny)
x = DBLARR(nx)

timestr = STRCOMPRESS(STRING(FIX(time/3600.)), /REMOVE_ALL)

deltax = len/double(nx)
deltay = len/double(ny)

FOR j = 0, ny-1 DO BEGIN
y(j) = j*deltay
FOR i = 0, nx-1 DO BEGIN
x(i) = i*deltax
	READF, 1, dum
	zeta(i,j) = dum
ENDFOR
ENDFOR
if (tracer EQ 0) then begin
    zeta = zeta * charzeta
endif
x = x / 1000.
y = y / 1000.

if (default EQ '1') then begin
inc = 4
   FOR j = 0, ny-1 DO BEGIN
   FOR i = 0, nx-1 DO BEGIN
	zeta(i,j) = (double(i)/double((nx-1)))*inc
   ENDFOR
   ENDFOR
endif

CLOSE, 1

; Begin routine to calculate the effective diffusivity
; in equivalent radius coordinates

; number of equiv. radius points
narea = nerad
; the array of equiv. radii
erad = DBLARR(nerad+1)
area = DBLARR(narea+1)
; maximum area defining domain
eradmax = len / two
maxarea = pi * eradmax^two
; vortex center indices
ic = fix (nx / two - one)
jc = fix (ny / two - one) 

; increment of zeta versus eradius and area
zeta_inc = (MAX(zeta) - MIN(zeta)) / double(nerad)
; the array of C(r_e)
; this is specified, and area under C(r_e) is calculated
zeta_erad = MIN(zeta) + FINDGEN(nerad+1)*zeta_inc
dzeta_derad = DBLARR(nerad)

; Now calculate grad zeta in rectangular coordinates across the domain
PRINT, "** kefferad: calculating tracer gradient array..."
gradzeta = DBLARR(nx,ny)
dzetadx = DBLARR(nx,ny)
dzetady = DBLARR(nx,ny)
FIRST_DERIVATIVE2D, zeta, deltax, deltay, 1, dzetadx
FIRST_DERIVATIVE2D, zeta, deltax, deltay, 2, dzetady
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
      gradzeta(i,j) = (dzetadx(i,j)^two + dzetady(i,j)^two)^half
ENDFOR
ENDFOR

; this routine calculates the area coordinates
; and equivalent radius coordinates
gradzeta_area = DBLARR(narea+1)
dgradzeta_derad = DBLARR(nerad)
PRINT, "** kefferad: determining area and equivalent radius coordinates..."
FOR k = 0, nerad DO BEGIN
areatmp = 0.0d0
gradzetatmp = 0.0d0
print, 'zeta_erad = ', zeta_erad(k)
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
    rad = (((i-ic)*deltax)^two + ((j-jc)*deltay)^two)^half
    if (rad GT eradmax) then continue
    if (zeta(i,j) LE zeta_erad(k)) then begin
          areatmp = areatmp + deltax*deltay 
          gradzetatmp = gradzetatmp + gradzeta(i,j)*gradzeta(i,j)* $
                        deltax*deltay  
    endif      
ENDFOR
ENDFOR
if (areatmp EQ 0.0) then begin
    print, 'Warning: Zero area detected between contours!'
endif
if (k EQ 0) then areatmp = 0.0d0
if (k EQ 0) then gradzetatmp = 0.0d0
area(k) = areatmp
erad(k) = ((area(k) - maxarea) / (-one * pi))^(half)	
gradzeta_area(k) = gradzetatmp
ENDFOR

; The array of delta equiv. radius
derad = DBLARR(nerad+1)
; dC / dr_e or dZeta / dr_e
; Staggered grid for centered differencing
; C --- dC/dr --- C --- dC/dr --- C
PRINT, "** kefferad: calculating tracer gradient with equivalent radius..."
FOR k = 0, nerad-1 DO BEGIN
     derad(k) = erad(k)-erad(k+1)
     dzeta_derad(k) = (zeta_erad(k+1) - zeta_erad(k)) / $
                      (erad(k+1) - erad(k))
ENDFOR
derad(nerad) = derad(nerad-1)

eradc = DBLARR(nerad)
FOR k = 0, narea-1 DO BEGIN
     eradc(k) = (erad(k) + erad(k+1)) / two
ENDFOR

; differentiate gradzeta_area with erad
FOR k = 0, nerad-1 DO BEGIN 
;  Note the -one is there because erad is monotonically decreasing
     dgradzeta_derad(k) = -one * (gradzeta_area(k+1) - gradzeta_area(k)) / $
                          (two*pi*eradc(k)*(erad(k+1) - erad(k)))
     print, 'dgradzeta_derad', dgradzeta_derad(k)
ENDFOR

PRINT, "** kefferad: calculating effective diffusivity as a function of equivalent radius..."
ediff_erad = DBLARR(nerad)
elength = DBLARR(nerad)
FOR k = 0, narea-1 DO BEGIN
     ediff_erad(k) = nu*dgradzeta_derad(k) / $
                    (dzeta_derad(k)^two)
     ; equivalent length in meters  
     elength(k) = ((ediff_erad(k) / nu) * $
                  (two * pi * eradc(k))^(two))^(half)                     
     ; for one plot at the end, store in array
     ediff_erad_oneplot(k,l) = nu*dgradzeta_derad(k) / $
                    (dzeta_derad(k)^two)
ENDFOR

;
;  Now loop over all points and assign effective diffusivity to 
;  the points (per email from E. Shuckburgh)
;

;  the effective diffusivity as a function of x and y
PRINT, "** kefferad: mapping effective diffusivity to x and y..."
ediff_xy = DBLARR(nx,ny)
elength_xy = DBLARR(nx,ny)
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN 
    FOR k = 0, nerad-1 DO BEGIN
         if (zeta(i,j) GE zeta_erad(k) AND $
             zeta(i,j) LE zeta_erad(k+1)) then begin
                ediff_xy(i,j) = ediff_erad(k)
                elength_xy(i,j) = elength(k)
         endif
    ENDFOR
ENDFOR
ENDFOR
PRINT, "** kefferad: writing data files..."

datafile = dir + infile + 'keff.tr-'+ tstamp + '.dat'
openw, 15, datafile

printf, 15, nerad, nx, ny
printf, 15, deltax, deltay, nu
printf, 15, "Effective Diffusivity & Equiv. Length vs. Equiv. Radius..."
FOR k = 0, narea-1 DO BEGIN
     printf, 15, eradc(k), ediff_erad(k), elength(k)
ENDFOR
printf, 15, "Effective Diffusivity vs. X and Y..."
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN 
     printf, 15, x(i), y(j), ediff_xy(i,j), elength_xy(i,j)
ENDFOR
ENDFOR
close, 15

ENDFOR

PRINT, "** kefferad: DONE."
END
