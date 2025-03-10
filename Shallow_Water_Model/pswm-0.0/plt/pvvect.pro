dir = ''
var = ''
temp = ''
outfile = ''
infilene = ''
infile2 = ''

;dir = '/usbhdd1/ldata/pswm/ellipse2/'
;infilene = 'ellipse'

;dir = '/usbhdd1/ldata/pswm/ring2/'
;infilene = 'ring2'

dir = '/ldata0/hendric/pswm/pubsim2/ftmp2/'
infilene = 'ftmp2'

;dir = '/ldata0/hendric/pswm/stratpolar/s2/'
;infilene = 's2'

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
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
; DEFINE THE PV CONTOUR LEVELS
;

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

OPENW, 10, dir+infilene+'maxwinds.dat' 
PRINTF, 10, "TIME (h), MAX WIND SPEED (m/s)"

;
; BEGIN FILE LOOP
;

FOR l = 0, nfiles DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

if (l lt 100) then begin

if (l lt 10) then begin
     tstamp = 't0' + lstring
     ostamp = '00' + lstring
     outfile1 = dir+infilene+'qvect.'+ostamp+'.ps'
     outfile2 = dir+infilene+'dvect.'+ostamp+'.ps'
  endif else begin
     tstamp = 't' + lstring
     ostamp = '0' + lstring
     outfile1 = dir+infilene+'qvect.'+ostamp+'.ps'
     outfile2 = dir+infilene+'dvect.'+ostamp+'.ps'
  endelse
endif else begin
  tstamp = lstring
  outfile1 = dir+infilene+'qvect.'+tstamp+'.ps'
  outfile2 = dir+infilene+'dvect.'+tstamp+'.ps'
endelse

;
; OPEN THE OUTPUT FILES
; 

infile1 = dir + infilene + 'q.' + tstamp 
OPENR, 1, infile1, ERROR = err
infile2= dir + infilene + 'u.' + tstamp 
OPENR, 2, infile2, ERROR = err 
infile3= dir + infilene + 'v.' + tstamp 
OPENR, 3, infile3, ERROR = err 
infile4= dir + infilene + 'd.' + tstamp 
OPENR, 4, infile4, ERROR = err
if (err NE 0) then begin 
   PRINT, "** pswm.pro: file not found, exiting..."  
   break
endif

PRINT, "** pswm.pro: opening and reading " + infile1 + "..."
READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
READF, 1, xmin, xmax, ymin, ymax
READF, 1, fcor, beta, cval, timesec
READF, 2, temp
READF, 2, temp
READF, 2, nx, ny
READF, 2, xmin, xmax, ymin, ymax
READF, 2, fcor, beta, cval, timesec
READF, 3, temp
READF, 3, temp
READF, 3, nx, ny
READF, 3, xmin, xmax, ymin, ymax
READF, 3, fcor, beta, cval, timesec
READF, 4, temp
READF, 4, temp
READF, 4, nx, ny
READF, 4, xmin, xmax, ymin, ymax
READF, 4, fcor, beta, cval, timesec

nx = FIX(nx)
ny = FIX(ny)
field1 = DBLARR(nx+1,ny+1)
field2 = DBLARR(nx+1,ny+1)
field3 = DBLARR(nx+1,ny+1)
field4 = DBLARR(nx+1,ny+1)
speed = DBLARR(nx+1,ny+1)
y = DBLARR(ny+1)
x = DBLARR(nx+1)
xv = DBLARR(nx+1)
yv = DBLARR(ny+1)
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
	field1(i,j) = dum1
        READF, 2, dum1
	field2(i,j) = dum1
        READF, 3, dum1
	field3(i,j) = dum1
        READF, 4, dum1
	field4(i,j) = dum1
ENDFOR
ENDFOR
CLOSE,1,2,3,4

;
; MULTIPLY FIELD BY FACTOR
;
fac4=10000.0
fac2=1.0
fac3=1.0
fac1=1000.0
field1 = field1 * fac1
field2 = field2 * fac2
field3 = field3 * fac3
field4 = field4 * fac4

; REDUCE VECTORS FOR PLOTTING
field2_plot = DBLARR(nx+1,ny+1)
field3_plot = DBLARR(nx+1,ny+1)
FOR j = 0, ny DO BEGIN
FOR i = 0, nx DO BEGIN
   IF i MOD 5 EQ 0 AND j MOD 5 EQ 0 THEN BEGIN
       field2_plot(i,j) = field2(i,j)
       field3_plot(i,j) = field3(i,j)
   ENDIF ELSE BEGIN
       field2_plot(i,j) = 0.0
       field3_plot(i,j) = 0.0
   ENDELSE
ENDFOR
ENDFOR

x = x / 1000.0d0
y = y / 1000.0d0
yv = y
xv = x

CLOSE, 1, 2, 3

PRINT, "** pswm.pro: plotting " + infile1 + "..."

xaxis = "x (km)"
yaxis = "y (km)"
plot_title = "P (x 10!E-3!N s!E-1!N)"
titles = [plot_title,xaxis,yaxis]
bounds = [-50,50,-50,50,timesec/3600.0d0] 
clevels = 0.01*[-30,5,8,10,12,15,20,25,30,35,40,50,60,70,80,90,100]
;clevels =0.01*[58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90]
pcolors = rbw_plot_colors
;pcolors = bwl_plot_colors
twodplottervec, field1, field2_plot, field3_plot, x, y, xv, yv, titles, pcolors, $
             clevels, bounds, outfile1

xaxis = "x (km)"
yaxis = "y (km)"
plot_title = "DIV (x 10!E-4!N s!E-1!N)"
titles = [plot_title,xaxis,yaxis]
bounds = [-50,50,-50,50,timesec/3600.0d0] 
clevels = 0.05*[-150,-100,-50,-40,-30,-20,-10,-5,5,10,20,30,40,50,100,150,200]
pcolors = c2w_plot_colors
pcolors = bwl_plot_colors 
twodplottervec, field4, field2_plot, field3_plot, x, y, xv, yv, titles, pcolors, $
             clevels, bounds, outfile2

hour_run = STRMID(STRCOMPRESS(STRING(bounds(4)),/REMOVE_ALL),0,5)

;contourxy, field, x, y, 9, hour_run, titles(0), $
;                titles(1), titles(2), outfile, bounds

fieldmean1 = total(field1) / (double(nx)*double(ny))
fieldmean4 = total(field4) / (double(nx)*double(ny))
PRINT, "** pswm.pro: max, min, mean (PV) : ", max(field1), min(field1), fieldmean1
PRINT, "** pswm.pro: max, min, mean (DIV): ", max(field4), min(field4), fieldmean4
PRINT, "** pswm.pro: max, min (V): ", max(field2), min(field2)

speed = (field2^2.0d0 + field3^2.0d0)^0.5d0

PRINT, " *** time, speed: ", timesec/3600.0d0, max(speed)
PRINTF, 10, timesec/3600.0d0, max(speed)

ENDFOR
CLOSE, 10

PRINT, "** pswm.pro: DONE."

END
