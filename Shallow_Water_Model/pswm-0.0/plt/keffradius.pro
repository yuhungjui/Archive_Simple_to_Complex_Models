PRINT, 'Running IDL keffradius.pro'
PRINT, 'Effective Diffusivity Routine'
PRINT, 'for Pseudospectral model output'
PRINT, 'Equivalent radius based version'

dir = ''
var = ''
outfile = ''
READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)

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

FOR l = 1, nfiles DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

if (l lt 100) then begin

if (l lt 10) then begin
   tstamp = '00' + lstring
   outfile = dir+infile+tstamp+'.ps'
endif else begin
   tstamp = '0' + lstring
   outfile = dir+infile+tstamp+'.ps'
endelse

endif else begin

tstamp = lstring
outfile = dir+infile+tstamp+'.ps'

endelse

OPENR, 1, dir + infile + tstamp    

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
zeta = FLTARR(nx,ny)
y = FLTARR(ny)
x = FLTARR(nx)

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
zeta = zeta * charzeta
x = x / 1000.
y = y / 1000.

CLOSE, 1

; Begin routine to calculate the effective diffusivity
; in equivalent radius coordinates
; do this in both eradius and area coordinates

; number of area points
narea = 6
; the array of areas
area = FLTARR(narea+1)
radmax = len / two
; maximum area defining domain
maxarea = pi * radmax * radmax
; vortex center indices
ic = nx / two
jc = ny / two

; find min zeta
zetamin = 100.0d
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
    rad = (((i-ic)*deltax)^two + ((j-jc)*deltay)^two)^half
    if (rad GT radmax) then continue
    if (zeta(i,j) LE zetamin) then begin
         zetamin = zeta(i,j)
    endif      
ENDFOR
ENDFOR

; increment of zeta versus eradius and area
zeta_inc = (MAX(zeta) - zetamin) / narea

; the array of zeta(area)
zeta_eradius = zetamin + FINDGEN(narea+1)*zeta_inc
zeta_area = zetamin + FINDGEN(narea+1)*zeta_inc
dzeta_darea = FLTARR(narea)

eradius = FLTARR(narea+1)
; specify A = pi*rmax**2 - pi*eradius*2
dzeta_deradius = FLTARR(narea)

FOR k = 0, narea DO BEGIN
areatmp = 0.0
print, 'zeta_area = ', zeta_area(k)
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
    rad = (((i-ic)*deltax)^two + ((j-jc)*deltay)^two)^half
    if (rad GT radmax) then continue
    if (zeta(i,j) LE zeta_area(k)) then begin
         areatmp = areatmp + deltax*deltay
    endif      
ENDFOR
ENDFOR
area(k) = areatmp
;eradius(k) = (area(k) / pi)^half
eradius(k) = ((maxarea - area(k))/pi)^half	
ENDFOR

FOR k = 0, narea-1 DO BEGIN
     dzeta_darea(k) = (zeta_area(k+1) - zeta_area(k)) / $
                      (area(k+1) - area(k))
     dzeta_deradius(k) = (zeta_eradius(k+1) - zeta_eradius(k)) / $
                      (eradius(k+1) - eradius(k))
ENDFOR 


; Now calculate grad zeta in rectangular coordiantes across the domain
gradzeta = FLTARR(nx,ny)
FOR j = 1, ny-2 DO BEGIN
FOR i = 1, nx-2 DO BEGIN
      rad = (((i-ic)*deltax)^two + ((j-jc)*deltay)^two)^half
      if (rad GT radmax) then continue
      dzetadx = (zeta(i+1,j) - zeta(i-1,j)) / (2.0d * deltax)
      dzetady = (zeta(i,j+1) - zeta(i,j-1)) / (2.0d * deltay)	
      gradzeta(i,j) = (dzetadx^two + dzetady^two)^half
ENDFOR
ENDFOR

; Now finally calculate the effective diffusivity
gradzeta_area = FLTARR(narea+1)
dgradzeta_deradius = FLTARR(narea)
FOR k = 0, narea DO BEGIN
gradzetatmp = 0.0d
gcnt = 0
FOR j = 0, ny-1 DO BEGIN
FOR i = 0, nx-1 DO BEGIN
      rad = (((i-ic)*deltax)^two + ((j-jc)*deltay)^two)^half
      if (rad GT radmax) then continue
      if (zeta(i,j) LT zeta_eradius(k)) then begin  
          gradzetatmp = gradzetatmp + gradzeta(i,j)*gradzeta(i,j)* $
                        deltax*deltay     
      endif    
ENDFOR
ENDFOR
print, gradzetatmp
gradzeta_area(k) = gradzetatmp 
ENDFOR

; differentiate gradzeta_area with eradius
FOR k = 0, narea-1 DO BEGIN
     dgradzeta_deradius(k) = (gradzeta_area(k+1) $ 
                             - gradzeta_area(k)) / $
                             (eradius(k+1) - eradius(k))
ENDFOR

ediff_eradius = FLTARR(narea)
; integer level eradius points
eradius_center = FLTARR(narea)
FOR k = 0, narea-1 DO BEGIN
    fac = - (two * pi * eradius(k)) 
    ediff_eradius(k) = fac*nu*dgradzeta_deradius(k) / $
                       (dzeta_deradius(k)^two)
    eradius_center(k) = half*(eradius(k) + eradius(k+1))                
ENDFOR

; Plot scaled zeta versus area
PLOT, eradius_center(1:narea-2)/1000.0d, alog(ediff_eradius(1:narea-2)), $
    XTITLE = 'Equivalent Radius (km)', $
    YTITLE = 'LN Effective Diffusivity', TITLE = 'EDIFF vs ERADIUS', $
    YRANGE = [min(alog(ediff_eradius(1:narea-2))), $
              max(alog(ediff_eradius(1:narea-2)))], $
    XRANGE = [min(eradius_center(1:narea-2)/1000.0d), $
              max(eradius_center(1:narea-2)/1000.0d)]

PLOT, eradius_center/1000.0d, alog(ediff_eradius), $
    XTITLE = 'Equivalent Radius (km)', $
    YTITLE = 'LN Effective Diffusivity', TITLE = 'EDIFF vs ERADIUS', $
    YRANGE = [min(alog(ediff_eradius)), $
              max(alog(ediff_eradius))], $
    XRANGE = [min(eradius_center/1000.0d), $
              max(eradius_center/1000.0d)]
 
PLOT, eradius/1000.0d, zeta_eradius, $
    XTITLE = 'Equivalent Radius (km)', $
    YTITLE = 'Zeta (1/s)', TITLE = 'ZETA vs ERADIUS', $
    YRANGE = [min(zeta_eradius), $
              max(zeta_eradius)], $
    XRANGE = [min(eradius/1000.0d), $
              max(eradius/1000.0d)]

DEVICE, /CLOSE


ENDFOR

END
