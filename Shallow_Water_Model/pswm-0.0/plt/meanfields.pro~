dir = ''
var = ''
outfile = ''
PRINT, 'Running IDL meanfields.pro'
PRINT, 'Computes the azimuthal mean vorticity'
PRINT, 'and tangential velocity for a psndbp run' 
PRINT, 'and plots it.'
READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)
READ, vorzeta, PROMPT="Vorticity (0) or Tang. Velocity (1): "

CLOSE, 1, 2, 3
temp = ''
len = 600000.0d0   ; m
pi = 3.14159d0
g = 9.81d0

nfiles = 500
nrmax = 500
ntheta = 50

; First open .in.out file to get characteristic vorticity
OPENR, 2, dir+infile + 'out'
FOR i = 0, 30 DO BEGIN
    readf, 2, temp
ENDFOR
CLOSE, 2
zetachar_string = STRMID(temp,22,14)
zetachar = float(zetachar_string)
print, "Charateristic vorticity: ", zetachar

FOR l = 2, nfiles DO BEGIN

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
PRINT, "** meanfields: opening input file: " +dir+infile+tstamp+ " ..."

OPENR, 1, dir + infile + tstamp, ERROR = err
;if (err NE 0) then break  

READF, 1, temp
READF, 1, temp
READF, 1, time
READF, 1, nx, ny

nx = fix(nx)
ny = fix(ny)
zeta = FLTARR(nx,ny)
y = FLTARR(ny)
x = FLTARR(nx)

timestr = STRCOMPRESS(STRING(time/3600.), /REMOVE_ALL)

deltax = len/nx
deltay = len/ny

FOR j = 0, ny-1 DO BEGIN
y(j) = j*deltay
FOR i = 0, nx-1 DO BEGIN
x(i) = i*deltax
	READF, 1, dum
	zeta(i,j) = dum
ENDFOR
ENDFOR
zeta = zeta * zetachar
CLOSE, 1

nr = fix(0.5d0 * nx)
xcen = 384
ycen = 384
zetamean = FLTARR(nr)
zetapolar = FLTARR(nr,ntheta)
zetaprime = FLTARR(nr,ntheta)

INTERPTOCYLIND2D, zeta, zetapolar, nr, ntheta, $
		  xcen, ycen, xPositions, yPositions

zetamean = AZIMUTHALMEAN2D(zetapolar, PERTURBATION = zetaprime)

deltar = deltax
radius = FINDGEN(nr) * deltar
namestr = "zm."

if (vorzeta GT 0.5) then begin
rzetamean = FLTARR(nr)
vmean = FLTARR(nr)
; compute the tangential velocity from mean vorticity
rzetamean = radius * zetamean
vmean(0) = 0.0
FOR i = 2, nr-1 DO BEGIN
   integral1d, rzetamean, 0, i, deltar, integrand
   vmean(i) = integrand / radius(i)
ENDFOR
vmean(1) = (vmean(2) + vmean(0)) / 2.0
namestr = "vm."
endif

!P.MULTI=[0,1,1]
SET_PLOT,'ps'
;SET_PLOT, 'X'

zetafile = dir + infile + namestr + tstamp + ".ps"

DEVICE,bits=8,filename=zetafile,/portrait,xoffset=1,yoffset=4,xsize=6, $
       ysize=6,/inches, /COLOR

if (vorzeta LT 0.5) then begin
   PLOT, radius/1000.0, zetamean, XTITLE = "Radius (km)", YTITLE = "Zeta (1/s)", $
      TITLE = "Mean Vorticity (1/s), t = " + timestr + " h." , $, 
      YRANGE = [-0.0001,0.004], XRANGE = [0,150] 
endif else begin
   PLOT, radius/1000.0, vmean, XTITLE = "Radius (km)", YTITLE = "Tangential Velocity (m/s)", $
      TITLE = "Mean Tangential Velocity (m/s), t = " + timestr + " h.", $, 
      YRANGE = [-2,40], XRANGE = [0,150]
endelse

DEVICE, /CLOSE

PRINT, "** meanfields: completed diagnostics for t = " + timestr + " h..."

ENDFOR

PRINT, "** meanfields: DONE."

END
