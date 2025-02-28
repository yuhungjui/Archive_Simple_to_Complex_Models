dir = ''
READ, dir, PROMPT="Enter output directory path: "

var = STRCOMPRESS(STRING(var), /REMOVE_ALL)
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)

CLOSE, 1, 2, 3
temp = ''
nfiles = 100
pi = 3.14159
g = 9.81

OPENR, 1, dir + 'cstest_u.t00'
READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
CLOSE, 1

nx = fix(nx+1)
ny = fix(ny+1)
nyd2 = fix(ny/2.)
nxd2 = fix(nx/2.)

u = FLTARR(nx,ny)
v = FLTARR(nx,ny)

y = FLTARR(ny)
x = FLTARR(nx)

FOR l = 0, nfiles DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)
if (l lt 10) then begin
   tstamp = '0' + lstring
endif else begin
   tstamp = lstring
endelse

!P.MULTI=[0,1,1]
SET_PLOT,'ps'
;SET_PLOT, 'X'
DEVICE,bits=8,filename=dir+'cheb'+tstamp+'.ps',/portrait,xoffset=1,yoffset=4,xsize=6, $
       ysize=5,/inches, /COLOR

OPENR, 1, dir + 'cstest_u.t' +tstamp
READF, 1, temp
READF, 1, temp
READF, 1, dum, dum
READF, 1, xmin, xmax, ymin, ymax
READF, 1, f, beta, cval, time

timestr = STRCOMPRESS(STRING(time/3600.), /REMOVE_ALL)

xs2 = (xmax + xmin)/2
xd2 = (xmax - xmin)/2
ys2 = (ymax + ymin)/2
yd2 = (ymax - ymin)/2

delta = (xmax - xmin)/nx

mx = fix(3*nx/2)
my = fix(3*ny/2)

FOR j = 0, ny-1 DO BEGIN
	y(j) = j*delta-0.5*(ymax-ymin)
	FOR i = 0, nx-1 DO BEGIN
		x(i) = i*delta-0.5*(xmax-xmin)
		READF, 1, dum
		u(i,j) = dum
	ENDFOR
ENDFOR
CLOSE,1

OPENR, 2, dir + 'cstest_v.t' +tstamp
READF, 2, temp
READF, 2, temp
READF, 2, dum, dum
READF, 2, xmin, xmax, ymin, ymax
READF, 2, f, beta, cval, time

timestr = STRCOMPRESS(STRING(time/3600.), /REMOVE_ALL)

xs2 = (xmax + xmin)/2
xd2 = (xmax - xmin)/2
ys2 = (ymax + ymin)/2
yd2 = (ymax - ymin)/2

delta = (xmax - xmin)/nx

mx = fix(3*nx/2)
my = fix(3*ny/2)

FOR j = 0, ny-1 DO BEGIN
	y(j) = j*delta-0.5*(ymax-ymin)
	FOR i = 0, nx-1 DO BEGIN
		x(i) = i*delta-0.5*(xmax-xmin)
		READF, 2, dum
		v(i,j) = dum
	ENDFOR
ENDFOR

if (var eq 'p') then begin
    ; converts from scaled geopotential to height
    q = q*cval/g   
endif

if (l eq 0) then begin
      q0 = q
endif

x = x / 1000.
y = y / 1000.

CLOSE, 2

d = 40
contour_num = 8
plot_title = var

VELOVECT, u, v

DEVICE, /CLOSE

ENDFOR

END
