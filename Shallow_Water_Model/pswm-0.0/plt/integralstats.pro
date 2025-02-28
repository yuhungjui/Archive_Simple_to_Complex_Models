dir = ''
var = ''
outfile = ''
PRINT, 'Running IDL integralstats.pro'
PRINT, 'Computes the enstrophy and palinstrophy
PRINT, 'for a psndbp run.'
READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)

CLOSE, 1, 2, 3
temp = ''
len = 600000.0   ; m
nfiles = 500
pi = 3.14159
g = 9.81

enstrophy = FLTARR(nfiles)
palinstrophy = FLTARR(nfiles)
timearray = FLTARR(nfiles)

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

OPENR, 1, dir + infile + tstamp, ERROR = err
if (err NE 0) then break  

READF, 1, temp
READF, 1, temp
READF, 1, time
READF, 1, nx, ny

timearray(l) = time

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

tempnum = 0.0
tempnum2 = 0.0
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
	tempnum2 = tempnum + 0.5 * zeta(i,j)*zeta(i,j)*deltax*deltay
        tempnum = tempnum2
ENDFOR
ENDFOR
enstrophy(l) = tempnum

dzetadx = FLTARR(nx,ny)
dzetady = FLTARR(nx,ny)
first_derivative2d, zeta, deltax, deltay, 1, dzetadx
first_derivative2d, zeta, deltax, deltay, 2, dzetady

tempnum = 0.0
tempnum2 = 0.0
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
	tempnum2 = tempnum + 0.5*deltax*deltay*( $
              dzetadx(i,j)*dzetadx(i,j) + dzetady(i,j)*dzetady(i,j))
        tempnum = tempnum2 
ENDFOR
ENDFOR
palinstrophy(l) = tempnum

PRINT, 'Enstrophy: ', enstrophy(l), ', Palinstrophy: ', palinstrophy(l)

x = x / 1000.
y = y / 1000.

CLOSE, 1

zetamean = total(zeta) / (double(nx)*double(ny))
PRINT, max(zeta), min(zeta), zetamean

ENDFOR

!P.MULTI=[0,1,1]
SET_PLOT,'ps'
;SET_PLOT, 'X'

DEVICE,bits=8,filename=dir+'integralstats.ps',/portrait,xoffset=1,yoffset=4,xsize=6, $
       ysize=5,/inches, /COLOR

OPENW, 2, dir + 'integralstats.dat'
FOR l = 2, nfiles DO BEGIN
   printf, 2, time(l), enstrophy(l), palinstrophy(l)
ENDFOR
CLOSE, 2

PLOT, time, enstrophy, XAXIS = "Time (sec)", YAXIS = "Enstrophy (1/s^2)"
PLOT, time, palinstrophy, XAXIS = "Time (sec)", YAXIS = "Palinstrophy (1/s^2)"

DEVICE, /CLOSE

END
