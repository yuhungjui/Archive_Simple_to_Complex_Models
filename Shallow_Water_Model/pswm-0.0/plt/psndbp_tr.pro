dir = ''
var = ''
outfile = ''
READ, dir, PROMPT="Enter output directory path: "
READ, tracer, PROMPT="Vorticity (0) or Tracer (1): "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile_noext = 'psndbp.in.'
READ, infile_noext, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile_noext+'.in.'), /REMOVE_ALL)
READ, plotcolor, PROMPT="Color plots (0-yes, 1-no): "
 
CLOSE, 1, 2, 3
temp = ''
nfiles = 500
pi = 3.14159
g = 9.81

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
c2w_red_color = [000, 032, 064, 096, 128, 160, 191, 223, 254, 255, 255, $
255, 255, 255, 255, 255, 255]
c2w_gre_color = [000, 032, 064, 096, 128, 160, 191, 223, 254, 213, 181, $
150, 118, 086, 054, 022, 000]
c2w_blu_color = [255, 255, 255, 255, 255, 255, 255, 255, 254, 223, 191, $
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

if (plotcolor LT 0.5) then begin
   pcolors = rbw_plot_colors
endif else begin
   pcolors = bwl_plot_colors
endelse

lbegin = 1
; COMMENT IN FOR UNSTABLE RINGS
;if (tracer EQ 1) then begin
;    lbegin = 2   
;endif

; First open .in.out file to get characteristic vorticity
OPENR, 2, dir+infile + 'out'
;FOR i = 0, 31 DO BEGIN
; COMMENT IN FOR READIN VORTICITY INPUT
FOR i = 0, 32 DO BEGIN
    readf, 2, temp
ENDFOR
CLOSE, 2
zetachar_string = STRMID(temp,22,14)
zetachar = float(zetachar_string)

FOR l = lbegin, nfiles DO BEGIN

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
   outfile = dir+infile+'tr-'+tstamp+'.ps'
endelse

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

deltax = 600000./nx
deltay = 600000./ny

FOR j = 0, ny-1 DO BEGIN
y(j) = j*deltay
FOR i = 0, nx-1 DO BEGIN
x(i) = i*deltax
	READF, 1, dum
	zeta(i,j) = dum
ENDFOR
ENDFOR
if (tracer LT 0.5) then zeta = zeta * zetachar * 1000.0e0

x = x / 1000.
y = y / 1000.

CLOSE, 1

if (tracer LT 0.5) then begin
  plot_title = infile_noext + " !7f!3 (x 10!E-3!N s!E-1!N)"
  plot_title = "!7f!3 (x 10!E-3!N s!E-1!N)" 

;  FOR RANKINE VORTEX IN TURBULENCE
  clevels = [-3,0.05,0.10,0.75,1.5,1.75,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7]

;  FOR BINARY VORTEX INTERACTION
;  clevels = [-2,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,3,3.5,4,4.5,5,5.5,6]
endif else begin
  plot_title = infile_noext + " c"
;  clevels = [1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700, $
;             1750,1800,1850,1900,1950,2000]
  clevels = [200,250,300,350,400,450,500,550,600,650,700, $
             750,800,850,900,950,1000]
endelse
xaxis = "x (km)"
yaxis = "y (km)"
titles = [plot_title,xaxis,yaxis]
bounds = [200,400,200,400,time/3600.0d0] 

; Total field
;CONTOURXY, zeta, x, y, contour_num, timestr, plot_title, xaxis, yaxis
twodplotter, zeta, x, y, titles, pcolors, $
           clevels, bounds, outfile

;Anomaly
;CONTOURXY, q-q0, x, y, contour_num, timestr, plot_title
;CONTOURXY_COLOR, q-q0, x, y, contour_num, timestr, plot_title
;CONTOURXY, q(d:nx-d,d:ny-d)-q0(d:nx-d,d:ny-d), $
;        x(d:nx-d), y(d:ny-d), contour_num, timestr, plot_title


zetamean = total(zeta) / (double(nx)*double(ny))
PRINT, max(zeta), min(zeta), zetamean

ENDFOR

END
