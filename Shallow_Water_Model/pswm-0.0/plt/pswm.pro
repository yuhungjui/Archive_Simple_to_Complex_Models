dir = ''
var = ''
temp = ''
outfile = ''
infilene = ''
infile2 = ''

dir = '/usbhdd1/ldata/pswm/ellipse2/'
infilene = 'ellipse'

;dir = '/usbhdd1/ldata/pswm/ring2/'
;infilene = 'ring2'

;dir = '/ldata0/hendric/pswm/pubsim2/f6u/'
;infilene = 'f6u'

dir='/ldata0/hendric/pswm/stratpolar/s2/'
infilene='s2'

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
READ, var, PROMPT="Enter variable to plot (p,u,v,z,d,q,h): "
varstring = STRCOMPRESS(var, /REMOVE_ALL)
READ, plotcolor, PROMPT="Color (0) or Black & White (1) plots: "

CLOSE, 1, 2, 3

;######################################################################
; Set up plotting colors
;######################################################################

new_color_table

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Rainbow Colors.
;rbw_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
rbw_red_color = [153, 125, 077, 000, 070, 070, 070, 070, 050, 200, 255, $
238, 238, 220, 238, 220, 255]
rbw_gre_color = [050, 095, 000, 000, 070, 150, 235, 255, 220, 255, 255, $
205, 154, 105, 069, 075, 000]
rbw_blu_color = [204, 186, 186, 205, 255, 255, 255, 160, 050, 050, 000, $
000, 000, 000, 000, 075, 000]

; THIS IS MODIFIED RAINBOW SO LAST COLOR IS WHITE (FOR NEGATIVE VALUES)
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Rainbow Colors.
;rbw_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
rbw_red_color = [254, 125, 077, 000, 070, 070, 070, 070, 050, 200, 255, $
                 238, 238, 220, 238, 220, 255]
rbw_gre_color = [254, 095, 000, 000, 070, 150, 235, 255, 220, 255, 255, $
                 205, 154, 105, 069, 075, 000]
rbw_blu_color = [254, 186, 186, 205, 255, 255, 255, 160, 050, 050, 000, $
                 000, 000, 000, 000, 075, 000]

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Cold to warm
;c2w_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
;c2w_red_color = [000, 032, 064, 096, 128, 160, 191, 223, 254, 255, 255, $
;255, 255, 255, 255, 255, 255]
;c2w_gre_color = [000, 032, 064, 096, 128, 160, 191, 223, 254, 213, 181, $
;150, 118, 086, 054, 022, 000]
;c2w_blu_color = [255, 255, 255, 255, 255, 255, 255, 255, 254, 223, 191, $
;160, 128, 096, 064, 032, 000]

;
; NOTE: THIS C2W HAS TWO WHITES IN THE CENTER
; DONE FOR THE DIVERGENCE PLOTS FOR THE RADIATIVE VORTEX
;

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Cold to warm
;c2w_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
c2w_red_color = [000, 032, 064, 096, 128, 160, 191, 254, 254, 255, 255, $
255, 255, 255, 255, 255, 255]
c2w_gre_color = [000, 032, 064, 096, 128, 160, 191, 254, 254, 213, 181, $
150, 118, 086, 054, 022, 000]
c2w_blu_color = [255, 255, 255, 255, 255, 255, 255, 254, 254, 223, 191, $
160, 128, 096, 064, 032, 000]

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Black-White Colors (low value - black, high value - white)
;rbw_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
bwh_red_color = [000, 015, 030, 045, 060, 075, 090, 105, 120, 135, 150, $
165, 180, 195, 210, 225, 240]
bwh_gre_color = [000, 015, 030, 045, 060, 075, 090, 105, 120, 135, 150, $
165, 180, 195, 210, 225, 240]
bwh_blu_color = [000, 015, 030, 045, 060, 075, 090, 105, 120, 135, 150, $
165, 180, 195, 210, 225, 240]

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Black-White Colors. (low value - white, high value - black)
;rbw_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
bwl_red_color = [250, 225, 210, 195, 180, 165, 150, 135, 120, 105, 090, $
075, 060, 045, 030, 015, 000]
bwl_gre_color = [250, 225, 210, 195, 180, 165, 150, 135, 120, 105, 090, $
075, 060, 045, 030, 015, 000]
bwl_blu_color = [250, 225, 210, 195, 180, 165, 150, 135, 120, 105, 090, $
075, 060, 045, 030, 015, 000]

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Black-White Colors HIGH CONTRAST. (low value - white, high value - black)
;rbw_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
;12   13   14   15   16   17]
bwc_red_color = [250, 170, 160, 150, 140, 130, 120, 110, 100, 090, 080, $
070, 060, 045, 030, 015, 000]
bwc_gre_color = [250, 170, 160, 150, 140, 130, 120, 110, 100, 090, 080, $
070, 060, 045, 030, 015, 000]
bwc_blu_color = [250, 170, 160, 150, 140, 130, 120, 110, 100, 090, 080, $
070, 060, 045, 030, 015, 000]

c2w_plot_colors = c32(c2w_red_color, c2w_gre_color, c2w_blu_color, 1)
rbw_plot_colors = c32(rbw_red_color, rbw_gre_color, rbw_blu_color, 1)
bwh_plot_colors = c32(bwh_red_color, bwh_gre_color, bwh_blu_color, 1)
bwl_plot_colors = c32(bwl_red_color, bwl_gre_color, bwl_blu_color, 1)
bwc_plot_colors = c32(bwc_red_color, bwc_gre_color, bwc_blu_color, 1)
NLEV1 = FIX(SIZE(c2w_plot_colors, /DIMENSIONS))
NLEV  = NLEV1(0)
black_cons = c32(intarr(NLEV), intarr(NLEV), intarr(NLEV), 1)

;
; DEFINE THE CONTOUR LEVELS
;

fac = 1.0
if (var EQ 'p') then begin
   clevels = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
   plot_title = "p (m s!E-1!N)"
endif else if (var EQ 'u' OR var EQ 'v') then begin
   clevels = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
   plot_title = "v (m s!E-1!N)"
endif else if (var EQ 'z') then begin
   clevels = 0.01*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
;   clevels = 0.1*[-10,2,3,4,5,6,8,10,12,15,20,25,30,40,50,75,100]
   plot_title = "!7f!3 (x 10!E-3!N s!E-1!N)"
   fac = 1000.0
; comment in for perturbation
;   clevels = 0.0000001*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
; comment in for same as ringbatch psndbp
;    clevels = [-3,0.75,1,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7]
endif else if (var EQ 'd') then begin
;   clevels = 0.001*[-30,-20,-10,-7,-5,-2,-0.5,-0.3,0,0.3,0.5,2,5,7,10,20,30]
   clevels = 0.03*[-30,-20,-10,-7,-5,-2,-0.5,-0.05,0,0.05,0.5,2,5,7,10,20,30]
   plot_title = "!7d!3 (x 10!E-4!N s!E-1!N)"
   fac = 10000.0
endif else if (var EQ 'q') then begin
 ;  clevels = 0.00004*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
   clevels = 0.1*[-30,7,7.2,7.4,7.6,8,9,10,12,15,20,25,30,40,50,75,100]
   fac=1000.0
endif else if (var EQ 'h') then begin
   clevels = 10.*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
   plot_title = "h (m)"
   fac = 1
endif

if (plotcolor LT 0.5) then begin
;    pcolors = c2w_plot_colors
   pcolors = rbw_plot_colors
endif else begin
   pcolors = bwl_plot_colors
endelse

omega = 7.292d-05
nfiles = 1000
pi = 3.14159d0
g = 9.81d0
ostamp = ''

;
; OPEN MAX WIND FILE
;

OPENW, 10, dir+infilene+'minpressure.dat' 
PRINTF, 10, "TIME (h), MINIMUM PRESSURE (m/s), MIN FLUID DEPTH (m)"

;
; BEGIN FILE LOOP
;

FOR l = 0, nfiles DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

if (l lt 100) then begin

if (l lt 10) then begin
     tstamp = 't0' + lstring
     ostamp = '00' + lstring
     outfile = dir+infilene+varstring+'.'+ostamp+'.ps'
  endif else begin
     tstamp = 't' + lstring
     ostamp = '0' + lstring
     outfile = dir+infilene+varstring+'.'+ostamp+'.ps'
  endelse
endif else begin
  tstamp = lstring
  outfile = dir+infilene+varstring+'.'+tstamp+'.ps'
endelse

;
; OPEN THE OUTPUT FILE
; 
if var EQ 'h' then varstring = 'p'
if (var NE 'q') then begin
   infile = dir + infilene + varstring + '.' + tstamp 
   OPENR, 1, infile, ERROR = err 
endif else begin
   infile  = dir + infilene + 'z.' + tstamp
   infile2 = dir + infilene + 'p.' + tstamp
   OPENR, 1, infile, ERROR = err 
   OPENR, 2, infile2, ERROR = err
endelse
if (err NE 0) then begin 
   PRINT, "** pswm.pro: file not found, exiting..."  
   break
endif
PRINT, "** pswm.pro: opening and reading " + infile + "..."
if (var EQ 'q') then PRINT, "** pswm.pro: opening and reading " + infile2 + "..."

READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
READF, 1, xmin, xmax, ymin, ymax
READF, 1, fcor, beta, cval, timesec

if (var EQ 'q') then begin
   READF, 2, temp
   READF, 2, temp
   READF, 2, nx, ny
   READF, 2, xmin, xmax, ymin, ymax
   READF, 2, fcor, beta, cval, timesec
endif

nx = FIX(nx)
ny = FIX(ny)
field = DBLARR(nx+1,ny+1)
field2 = DBLARR(nx+1,ny+1)
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
        if (var EQ 'q') then begin
 	   READF, 2, dum2
           field2(i,j) = dum2
        endif         
ENDFOR
ENDFOR

;
; COMPUTE POTENTIAL VORTICITY
;
grav = g
if (var EQ 'q') then begin
   field3 = DBLARR(nx+1,ny+1)
   field4 = DBLARR(nx+1,ny+1)
   field4 = (field2*cval + cval^2) / grav
   field3 = (cval^2/grav)*(field + fcor) / field4
   field = field3
endif

;
; COMPUTE FLUID DEPTH
;
if (var EQ 'h') then begin
   field3 = DBLARR(nx+1,ny+1)
   field3 = ( cval^2 + cval * field )  / g
   tmp1 = field3(nx-300,ny/2)
   field = field3 - tmp1
;   field = field3
endif
;
; MULTIPLY FIELD BY FACTOR
;

field = field * fac

x = x / 1000.0d0
y = y / 1000.0d0

CLOSE, 1, 2, 3

PRINT, "** pswm.pro: plotting " + infile + "..."

plot_title = var
xaxis = "x (km)"
yaxis = "y (km)"
titles = [plot_title,xaxis,yaxis]
bounds = [-50,50,-50,50,timesec/3600.0d0] 
bounds = [-6000,6000,-6000,6000,timesec/3600.0d0] 


;stop
twodplotter, field, x, y, titles, pcolors, $
             clevels, bounds, outfile

hour_run = STRMID(STRCOMPRESS(STRING(bounds(4)),/REMOVE_ALL),0,5)
;contourxy, field, x, y, 9, hour_run, titles(0), $
;                titles(1), titles(2), outfile, bounds

fieldmean = total(field) / (double(nx)*double(ny))
PRINT, "** pswm.pro: max, min, mean: ", max(field), min(field), fieldmean

if (var EQ 'h') then begin
   PRINT, " *** time, pressure: ", timesec/3600.0d0, min(field*9.81*1.13), min(field)
;   PRINTF, 10, timesec/3600.0d0, min(field*9.81*1.13), min(field)
    PRINTF, 10, timesec/3600.0d0, field(256,256)*9.81*1.13, field(256, 256)
endif

ENDFOR
CLOSE,10
PRINT, "** pswm.pro: DONE."

END
