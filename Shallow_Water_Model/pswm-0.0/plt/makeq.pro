dir = ''
infilene = ''
infile = ''
temp1 = ''
temp2 = ''

READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
READ, infilene, PROMPT="Enter infile (no extension): "
READ, fcor, PROMPT="Enter the latitude (degrees): "
READ, cval, PROMPT="Enter the cval: "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)

omega = 7.292d-05
pi = 3.14159d
fcor = 2.0d*omega*sin(fcor*pi/180.0d)
nfiles = 1000 

PRINT, "fcor = ", fcor

FOR l = 0, nfiles DO BEGIN

CLOSE, 1, 2, 3
lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

if (l lt 100) then begin

if (l lt 10) then begin
     tstamp = 't0' + lstring
  endif else begin
     tstamp = 't' + lstring
  endelse
endif else begin
  tstamp = lstring
endelse

;
; Open and read zeta and pressure files
; 

infilez  = dir + infilene + 'z.' + tstamp 
infilep  = dir + infilene + 'p.' + tstamp
outfileq = dir + infilene + 'q.' + tstamp

PRINT, "** makeq.pro: processsing all files for " + tstamp + "..."

OPENR, 1, infilez, ERROR = err 
OPENR, 2, infilep, ERROR = err
OPENW, 3, outfileq
if (err NE 0) then begin 
   PRINT, "** makeq.pro: file not found, exiting..."  
   break
endif

for m = 1, 2 do begin
   READF, m, temp1
   READF, m, temp2
   READF, m, nx, ny
   READF, m, xmin, xmax, ymin, ymax
   READF, m, fcor, beta, cval, timesec
endfor

nx = FIX(nx)
ny = FIX(ny)

PRINTF, 3, temp1
PRINTF, 3, "field q for run " + infilene
PRINTF, 3, nx, ny
PRINTF, 3, xmin, xmax, ymin, ymax
PRINTF, 3, fcor, beta, cval, timesec

FOR j = 0, ny DO BEGIN
FOR i = 0, nx DO BEGIN
	READF,  1, dum1
        READF,  2, dum2
	PRINTF, 3, (fcor + dum1) / (1.0d + dum2 / cval)          
ENDFOR
ENDFOR

CLOSE, 1, 2, 3

ENDFOR

END
