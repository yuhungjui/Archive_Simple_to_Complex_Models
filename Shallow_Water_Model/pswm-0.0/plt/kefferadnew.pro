PRINT, 'Running IDL kefferadnew.pro'
PRINT, 'Effective Diffusivity Routine'
PRINT, 'for Pseudospectral Model Output'
PRINT, 'Equivalent Radius Version'
PRINT, 'The ED calculation is exactly like Hendricks and Schubert (2008)'

dir = ''
var = ''
outfile = ''
default = ''
READ, dir, PROMPT="Enter output directory path:"
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)
READ, narea, PROMPT="Enter number of area points:"
READ, lbegin, PROMPT="Enter the begining file index:"
narea = FIX(narea)
lbegin = FIX(lbegin)

CLOSE, 1, 2, 3, 15

;
; Some initial parameters
;
temp = ''
nfiles = 500
pi = 3.14159d
g = 9.81d
half = 0.5d0
two = 2.0d0
one = 1.0d0

;
; First open the .in.out file to get domain size and tracer diffusion
;
OPENR, 2, dir+infile + 'out'
FOR i = 0, 55 DO BEGIN
    readf, 2, temp
    print, i, temp
;
;   Note may need to change to 22 and 52 for some .in.out files (evortex and rturb)
;   Generally set to 21 and 51
;
    if (i EQ 22) then len_string = STRMID(temp,31,14)
    if (i EQ 52) then trdiff_string = STRMID(temp,44,14)
ENDFOR
CLOSE, 2
; domain length in meters
len   = double(len_string) * 1000.0d0
; tracer diffusivity in m^2 / s
kappa = double(trdiff_string)
PRINT, "**kefferadnew: domain length (m) and tracer diffusivity (m^2 / s)"
PRINT, len, kappa


ediff_erad_oneplot = FLTARR(narea,nfiles)
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
OPENR, 1, dir + infile + 'tr-' + tstamp, ERROR = err
if (err NE 0) then break

READF, 1, temp
READF, 1, temp
READF, 1, time
READF, 1, nx, ny

nx = fix(nx)
ny = fix(ny)
conc = DBLARR(nx,ny)     ; array of tracer concentrations
y = DBLARR(ny)
x = DBLARR(nx)

deltax = len/double(nx)
deltay = len/double(ny)

;
;  read in the tracer concentration 2D array from
;  pseudospectral model .tr-000 files
;
FOR j = 0, ny-1 DO BEGIN
y(j) = j*deltay
FOR i = 0, nx-1 DO BEGIN
x(i) = i*deltax
	READF, 1, dum
	conc(i,j) = dum
ENDFOR
ENDFOR
x = x / 1000.
y = y / 1000.
CLOSE, 1

;
; Begin routine to calculate the effective diffusivity
; in the area coordinate (Hendricks and Schubert 2008)
; 

; the array of areas
area = DBLARR(narea+1)
; maximum area defining domain
maxarea = len * len
; vortex center indices
ic = fix (nx / two - one)
jc = fix (ny / two - one) 

; increment of zeta versus eradius and area
conc_inc = (MAX(conc) - MIN(conc)) / double(narea)
; the array of C(A)
conc_area = MIN(conc) + FINDGEN(narea+1)*conc_inc
; allocate dC/dA array
dconc_darea = DBLARR(narea)

; Calculate grad zeta in rectangular coordinates across the domain
PRINT, "** kefferadnew: calculating tracer gradient array..."
gradconc = DBLARR(nx,ny)
dconcdx = DBLARR(nx,ny)
dconcdy = DBLARR(nx,ny)
FIRST_DERIVATIVE2D, conc, deltax, deltay, 1, dconcdx
FIRST_DERIVATIVE2D, conc, deltax, deltay, 2, dconcdy
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
      gradconc(i,j) = (dconcdx(i,j)^two + dconcdy(i,j)^two)^half
ENDFOR
ENDFOR

; Calculate the A(C) function
gradconc_area = DBLARR(narea+1)
PRINT, "** kefferadnew: determining area and equivalent radius coordinates..."
area(0) = 0.0d0
gradconc_area(0) = 0.0d0
FOR k = 1, narea DO BEGIN
areatmp = 0.0d0
gradconctmp = 0.0d0
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
    if (conc(i,j) LE conc_area(k) AND conc(i,j) GT conc_area(k-1)) then begin
          areatmp = areatmp + deltax*deltay 
          gradconctmp = gradconctmp + gradconc(i,j)*gradconc(i,j)* $
                        deltax*deltay  
    endif      
ENDFOR
ENDFOR
if (areatmp EQ 0.0) then begin
    print, '** kefferadnew: WARNING: Zero area detected between contours!'
endif
if (k EQ 0) then areatmp = 0.0d0
if (k EQ 0) then gradzetatmp = 0.0d0
area(k) = area(k-1) + areatmp
gradconc_area(k) = gradconc_area(k-1) + gradconctmp	
print, "A, C(A) = ", area(k), conc_area(k)
ENDFOR

;
; dC / dA on staggered grid for centered differencing
; C --- dC/dA --- C --- dC/dA --- C
;
PRINT, "** kefferadnew: calculating tracer gradient with equivalent radius..."
FOR k = 0, narea-1 DO BEGIN
     dconc_darea(k) = (conc_area(k+1) - conc_area(k)) / $
                      (area(k+1) - area(k))
ENDFOR

;
; interpolate the area points to the center point
;
areac = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     areac(k) = (area(k) + area(k+1)) / two
ENDFOR

;
; differentiate gradconc_area with area
;
dgradconc_darea = DBLARR(narea)
PRINT, "** kefferadnew: calculating dgradconc_darea..."
FOR k = 0, narea-1 DO BEGIN 
     dgradconc_darea(k) = (gradconc_area(k+1) - gradconc_area(k)) / $
                          (area(k+1) - area(k))
ENDFOR

;
; calculate the effective diffusivity
;
PRINT, "** kefferadnew: calculating effective diffusivity as a function of area..."
ediff_area = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     ediff_area(k) = kappa*dgradconc_darea(k) / $
                     (dconc_darea(k)^two)
ENDFOR

;
; BEGIN CODE TO CALCULATE EQUIVALENT RADIUS DIAGNOSTICS
; ALL THAT IS NEEDED IS ediff_area FROM ABOVE
; HENDRICKS AND SCHUBERT (2008)
;

PRINT, "** kefferadnew: calculating equivalent radius diagnostics..."

;
; calculate equivalent radius array 
;
re = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     re(k) = ((areac(k) - maxarea) / (-one * pi))^(half)
ENDFOR

;
; calculate ediff_re (m^2 / s)
;
ediff_re = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
    ; ediff_re(k) = ediff_area(k) / (two * two * pi * areac(k))
      ediff_re(k) = ediff_area(k) / ((two * pi * re(k))^(two))
ENDFOR

;
; calculate Le_re (m) 
;
elength_re = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     elength_re(k) = ( ediff_re(k) / kappa ) * ( two * pi * re(k) )
ENDFOR

;
; calculate tau_re (s)
;
tau_re = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     tau_re(k) = ( two * pi * re(k) )^(two) / ( ediff_re(k) )
ENDFOR

;
; calculate lambda_re 
;
lambda_re = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     lambda_re(k) = ( ediff_re(k) / kappa )
ENDFOR

;
; calculate conc_re 
;
conc_re = DBLARR(narea+1)
FOR k = 0, narea-1 DO BEGIN
     conc_re(k) = conc_area(k)
ENDFOR

;
;  Now loop over all points and assign effective diffusivity to 
;  the points (per email from E. Shuckburgh)
;

;  the effective diffusivity as a function of x and y
PRINT, "** kefferadnew: mapping effective diffusivity to x and y..."
ediff_xy   = DBLARR(nx,ny)
elength_xy = DBLARR(nx,ny)
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN 
    FOR k = 0, narea-1 DO BEGIN
         if (conc(i,j) GT conc_re(k) AND $
             conc(i,j) LE conc_re(k+1)) then begin
                ediff_xy(i,j)   = ediff_re(k)
                elength_xy(i,j) = elength_re(k)
         endif
    ENDFOR
ENDFOR
ENDFOR

;
; Change units for output
;

re         = re / 1000.0d0          ; m to km
elength_re = elength_re / 1000.0d0  ; m to km
tau_re     = tau_re / 3600.0d0      ; s to h

;
; WRITE THE DATA TO OUTPUT FILES
;

PRINT, "** kefferadnew: writing data files..."

datafile = dir + infile + 'keff.tr-'+ tstamp + '.dat'

openw, 15, datafile
printf, 15, "Number of area points, nx, ny: "
printf, 15, narea, nx, ny
printf, 15, "dx (m), dy (m), kappa (m^2/s), time (s), domain size (m): "
printf, 15, deltax, deltay, kappa, time, len
printf, 15, "re (km), C(re), ediff_re (m^2/s), elength_re (km), tau_re (h), lambda_re..."
FOR k = 0, narea-1 DO BEGIN
     printf, 15, re(k), conc_re(k), ediff_re(k), $
                 elength_re(k), tau_re(k),lambda_re(k), $
                 FORMAT='(F7.2,F9.2,F15.5,F10.2,F20.5,F8.2)'
ENDFOR
printf, 15, "Effective Diffusivity vs. X and Y..."
; do not print out x, y arrays to save space
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN 
      printf, 15, ediff_xy(i,j), elength_xy(i,j)
ENDFOR
ENDFOR
close, 15

ENDFOR

PRINT, "** kefferadnew: DONE."
END
