PRINT, 'Running IDL keffarea.pro'
PRINT, 'Effective Diffusivity Routine'
PRINT, 'for Pseudospectral model output'
PRINT, 'Area-based version'

dir = ''
var = ''
outfile = ''
default = ''
READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)
READ, default, PROMPT="Test case (1-yes, 0-no):"

CLOSE, 1, 2, 3
temp = ''
nfiles = 500
pi = 3.14159d
g = 9.81d
len = 600000.0d
half = 0.5d0
two = 2.0d0
; regular viscosity for run
nu = 25.0d
charzeta = 0.001982d
tracer = 1

FOR l = 1, nfiles DO BEGIN

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
   OPENR, 1, dir + infile + tstamp   
   outfile = dir+infile+tstamp+'.ps' 
endif else begin
   OPENR, 1, dir + infile + 'tr-' + tstamp
   outfile = dir+infile+'keff.tr-'+tstamp+'.ps'
endelse

!P.MULTI=[0,1,1]
SET_PLOT,'ps'
;SET_PLOT, 'X'

DEVICE,bits=8,filename=outfile,/portrait,xoffset=1,yoffset=4,xsize=6, $
       ysize=5,/inches, /COLOR

READF, 1, temp
READF, 1, temp
READF, 1, time
READF, 1, nx, ny

nx = fix(nx)
ny = fix(ny)
zeta = DBLARR(nx,ny)
y = DBLARR(ny)
x = DBLARR(nx)

timestr = STRCOMPRESS(STRING(time/3600.), /REMOVE_ALL)

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
; in area coordinates

; number of area points
narea = 10
; the array of areas
area = DBLARR(narea+1)
; maximum area defining domain
maxarea = len * len 
; vortex center indices
ic = nx / two
jc = ny / two

; increment of zeta versus eradius and area
zeta_inc = (MAX(zeta) - MIN(zeta)) / double(narea)
; the array of zeta(area)
zeta_area = MIN(zeta) + FINDGEN(narea+1)*zeta_inc
dzeta_darea = DBLARR(narea)

; comment in for exponential area distribution
;zeta_area = DBLARR(narea+1)
;coeff = 7.6 / (narea+1)
;coeff2 = 4.0d6 / (narea+1)
;FOR k = 0, narea DO BEGIN
;zeta_area(k) = (0.000001d0)*exp(coeff*(k+1))
;zeta_area(k) = (0.000001d0)*(coeff2*(k+1))^half
;ENDFOR

; zeta_area is C(A) or zeta(A)
; area is A(C) or A(zeta)
FOR k = 0, narea DO BEGIN
areatmp = 0.0
print, 'zeta_area = ', zeta_area(k)
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
    if (zeta(i,j) LE zeta_area(k)) then begin
          areatmp = areatmp + deltax*deltay
    endif      
ENDFOR
ENDFOR
if (areatmp EQ 0.0) then begin
    print, 'Warning: Zero area detected between contours!'
endif
if (k EQ 0) then areatmp = 0
area(k) = areatmp	
ENDFOR

; dC / dA or dZeta / dA
FOR k = 0, narea-1 DO BEGIN
     dzeta_darea(k) = (zeta_area(k+1) - zeta_area(k)) / $
                      (area(k+1) - area(k))
ENDFOR

; Now calculate grad zeta in rectangular coordinates across the domain
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

; Now finally calculate the effective diffusivity
gradzeta_area = DBLARR(narea+1)
dgradzeta_darea = DBLARR(narea)
FOR k = 0, narea DO BEGIN
gradzetatmp = 0.0d
gcnt = 0
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
      if (zeta(i,j) LE zeta_area(k)) then begin  
          gradzetatmp = gradzetatmp + gradzeta(i,j)*gradzeta(i,j)* $
                        deltax*deltay     
      endif    
ENDFOR
ENDFOR
if (k EQ 0) then gradzetatmp = 0.0d
gradzeta_area(k) = gradzetatmp 
ENDFOR

; differentiate gradzeta_area with area
FOR k = 0, narea-1 DO BEGIN 
     dgradzeta_darea(k) = (gradzeta_area(k+1) - gradzeta_area(k)) / $
                          (area(k+1) - area(k))
     print, 'dgradzeta_darea', dgradzeta_darea(k)
ENDFOR

ediff_area = DBLARR(narea)
areac = DBLARR(narea)
FOR k = 0, narea-1 DO BEGIN
     areac(k) = (area(k) + area(k+1)) / two
     ediff_area(k) = nu*dgradzeta_darea(k) / $
                    (dzeta_darea(k)^two)                    
ENDFOR

; Plot scaled zeta versus area
PLOT, areac/maxarea, alog(ediff_area), XTITLE = 'area/maxarea', $
    YTITLE = 'LN effective diffusivity', TITLE = 'ediff vs area', $
    YRANGE = [0.0, max(alog(ediff_area))], $
    XRANGE = [0.98, 1.0]
     
; Plot equivalent length
PLOT, areac/maxarea, (ediff_area/nu)^0.5, XTITLE = 'area/maxarea', $
    YTITLE = 'Equivalent Length (m)', TITLE = 'elength vs area', $
    YRANGE = [0.0, max((ediff_area/nu)^0.5)], $
    XRANGE = [0.98, 1.0]

; Plot scaled zeta versus area
PLOT, area/maxarea, zeta_area, XTITLE = 'area/max area', $
    YTITLE = 'zeta', TITLE = 'zeta vs area', $
    YRANGE = [min(zeta_area), max(zeta_area)], $
    XRANGE = [0.0, 1.0]

; Plot scaled zeta versus area
PLOT, zeta_area, area/maxarea, XTITLE = 'zeta (1/s)', $
    YTITLE = 'area/areamax', TITLE = 'area vs zeta', $
    XRANGE = [min(zeta_area), max(zeta_area)], $
    YRANGE = [0.0, 1.0]

DEVICE, /CLOSE


ENDFOR

END
