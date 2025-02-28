PRINT, 'Running IDL kefferad_plot.pro'
PRINT, 'Effective Diffusivity PLOTTING Routine'
PRINT, 'for Pseudospectral Model Output'
PRINT, 'Equivalent Radius Version'
PRINT, 'Reads data format of KEFFERADNEW.PRO only'

lbegin = 2
tstamp = ''
close, 15
dir = ''
lstring = ''
timestr = ''
lcount = 0
len = 600000.0d0   ; default value, this is reset
neradmax = 200       ; default value, must be >= narea

;
;  This programs will dump out PS files for all
;  effective diffusivity diagnostics of re
;  To save time the 2D arrays are only output
;  if requested
;

READ, dir, PROMPT="Enter output directory path:"
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)
READ, lbegin, PROMPT="Enter beginning file number: "
lbegin = FIX(lbegin)
READ, lend, PROMPT="Enter end file number: "
lend = FIX(lend)
READ, oneortwod, PROMPT="oned or twod plots (0-1D, 1-2D): "
oneortwod = FIX(oneortwod)
READ, vartoplt, PROMPT="Plot type (0-ed, 1-elength, 2-C, 3-tau, 4-lambda): "
vartoplt = FIX(vartoplt)
READ, plotcolor, PROMPT="Color plots? (0-yes, 1-no)"
plotcolor = FIX(plotcolor)

nfiles = lend - lbegin + 1

;
; 2D radius-time arrays (hovmoller and mean plots)
;
ediff_erad_hov   = DBLARR(neradmax,nfiles)
elength_erad_hov = DBLARR(neradmax,nfiles)
conc_erad_hov    = DBLARR(neradmax,nfiles)
tau_erad_hov     = DBLARR(neradmax,nfiles)
lambda_erad_hov  = DBLARR(neradmax,nfiles)

;
; Array of times
;
time_array       = DBLARR(nfiles)
eradc_avg = DBLARR(neradmax)
eradc_avg = 500.0d0 + FINDGEN(neradmax) * 1500.0d0
eradc_avg = REVERSE(eradc_avg)
eradc_avg = eradc_avg / 1000.0d0    ;m to km
temp_edifferad = DBLARR(neradmax)
temp_elengtherad = DBLARR(neradmax)

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
bwl_red_color = [240, 225, 210, 195, 180, 165, 150, 135, 120, 105, 090, $
075, 060, 045, 030, 015, 000]
bwl_gre_color = [240, 225, 210, 195, 180, 165, 150, 135, 120, 105, 090, $
075, 060, 045, 030, 015, 000]
bwl_blu_color = [240, 225, 210, 195, 180, 165, 150, 135, 120, 105, 090, $
075, 060, 045, 030, 015, 000]

c2w_plot_colors = c32(c2w_red_color, c2w_gre_color, c2w_blu_color, 1)
rbw_plot_colors = c32(rbw_red_color, rbw_gre_color, rbw_blu_color, 1)
bwh_plot_colors = c32(bwh_red_color, bwh_gre_color, bwh_blu_color, 1)
bwl_plot_colors = c32(bwl_red_color, bwl_gre_color, bwl_blu_color, 1)
NLEV1 = FIX(SIZE(c2w_plot_colors, /DIMENSIONS))
NLEV  = NLEV1(0)
black_cons = c32(intarr(NLEV), intarr(NLEV), intarr(NLEV), 1)

if (plotcolor LT 0.5) then begin
   pcolors = rbw_plot_colors
endif else begin
   pcolors = bwh_plot_colors
endelse

;
; SET THE CONTOUR LEVELS
;

; for midrange ED
;clevels_ediff = [0,100,200,300,400,500,600,700,800,900,1000,1100, $
;           1300,1500,1700,1900,2100]

; for lowrange ED
clevels_ediff = [0,100,250,500,750,1000,1500,2000,2500,3000,3500, $
           4000,5000,6000,7000,8000,9000]

; for highrange ED
;clevels_ediff = [0,100,500,1000,1500,2000,3000,4000,5000, $
;           6000,7500,9000,11000,13000,15000,17500,20000]

; elength
clevels_elength = [0,100,200,300,400,500,600,700,800,900,1000,1100, $
           1300,1500,1700,1900,2100]

; C
clevels_conc = [0,50,100,150,200,250,300,350,400,450,500,550, $
           600,650,700,800,850]
clevels_conc = [500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900]

; tau
clevels_tau = 1000.0*[0.01,.5,1,10,30,50,70,100,150,200,300,400,500,600,700,800,900]

; lambda
clevels_lambda = 0.01 * [0,100,200,300,400,500,600,700,800,900,1000,1100, $
           1300,1500,1700,1900,2100]
clevels_lambda = [1,3,5,7,10,15,20,25,30,35,40,45,50,70,100,150,200]

IF vartoplt EQ 0 THEN BEGIN

    clevels = clevels_ediff

ENDIF ELSE IF vartoplt EQ 1 THEN BEGIN
  
    clevels = clevels_elength

ENDIF ELSE IF vartoplt EQ 2 THEN BEGIN
  
    clevels = clevels_conc

ENDIF ELSE IF vartoplt EQ 3 THEN BEGIN
  
    clevels = clevels_tau

ENDIF ELSE IF vartoplt EQ 4 THEN BEGIN
  
    clevels = clevels_lambda

ENDIF

;
; BEGIN THE FILE LOOP
;

FOR l = lbegin, lend DO BEGIN

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
; Open the correct output file
;
outfile = dir+infile+'keff.tr-'+tstamp+'.ps'

;
; Open and read the data from .dat file
;
datafile = dir + infile + 'keff.tr-'+ tstamp + '.dat'
openr, 15, datafile, ERROR=err
if (err NE 0) then break
temp = ''
readf, 15, temp
readf, 15, nerad, nx, ny

nerad = FIX(nerad)
nx = FIX(nx)
ny = FIX(ny)

readf, 15, temp
readf, 15, deltax, deltay, nu, time, len
readf, 15, temp
timestr = STRCOMPRESS(STRING(time/3600.0d0), /REMOVE_ALL)
time_array(lcount) = time/3600.0d0

; 1D ARRAYS
eradc = DBLARR(nerad)
conc_erad = DBLARR(nerad)
ediff_erad = DBLARR(nerad)
elength_erad = DBLARR(nerad)
tau_erad = DBLARR(nerad)
lambda_erad = DBLARR(nerad)

; 2D ARRAYS
ediff_xy = DBLARR(nx,ny)
elength_xy = DBLARR(nx,ny)

x = DBLARR(nx)
y = DBLARR(ny)

FOR k = 0, nerad-1 DO BEGIN
     readf, 15, t1, t2, t3, t4, t5, t6
     eradc(k) = t1
     conc_erad(k) = t2 
     ediff_erad(k) = t3
     elength_erad(k) = t4
     tau_erad(k) = t5
     lambda_erad(k) = t6
ENDFOR

;
;  THIS INTERPOLATES THE IRREGULAR ERAD TO 
;  LINEAR ERAD ARRAY
;

temp_edifferad   = INTERPOL(  ediff_erad, eradc, eradc_avg)
temp_elengtherad = INTERPOL(elength_erad, eradc, eradc_avg)
temp_concerad    = INTERPOL(   conc_erad, eradc, eradc_avg)
temp_tauerad     = INTERPOL(    tau_erad, eradc, eradc_avg)
temp_lambdaerad  = INTERPOL( lambda_erad, eradc, eradc_avg)

ediff_erad_hov(*,lcount)   = temp_edifferad
elength_erad_hov(*,lcount) = temp_elengtherad
conc_erad_hov(*,lcount)    = temp_concerad
tau_erad_hov(*,lcount)     = temp_tauerad
lambda_erad_hov(*,lcount)  = temp_lambdaerad

FOR i = 0, neradmax-1 DO BEGIN
   IF ediff_erad_hov(i,lcount) LT 0.0d0 THEN ediff_erad_hov(i,lcount) = 0.0d0
   IF elength_erad_hov(i,lcount) LT 0.0d0 THEN elength_erad_hov(i,lcount) = 0.0d0
   IF conc_erad_hov(i,lcount) LT 0.0d0 THEN conc_erad_hov(i,lcount) = 0.0d0
   IF tau_erad_hov(i,lcount) LT 0.0d0 THEN tau_erad_hov(i,lcount) = 0.0d0
   IF lambda_erad_hov(i,lcount) LT 0.0d0 THEN lambda_erad_hov(i,lcount) = 0.0d0
ENDFOR

;
; READ IN THE 2D DATA IF REQUESTED 
;

IF oneortwod EQ 1 THEN BEGIN
   readf, 15, temp
   FOR j = 0, ny-1 DO BEGIN
   y(j) = j*deltay
   FOR i = 0, nx-1 DO BEGIN 
   x(i) = i*deltax
      readf, 15, t1, t2
      ediff_xy(i,j) = t1
      elength_xy(i,j) = t2
   ENDFOR
   ENDFOR
ENDIF

close, 15

PRINT, "** kefferad_plot: creating all plots for t = " + timestr + " h..."

d = 250
hour_run = timestr
xax = 'x (km)'
yax = 'y (km)'
x = x / 1000.0d0
y = y / 1000.0d0
bounds = [200,400,200,400,time/3600.0d0]

IF oneortwod EQ 0 THEN BEGIN

     IF vartoplt EQ 0 THEN BEGIN

          onedbounds = [0,150,0,1000,time/3600.0d0]
          titles = ['ED vs Equiv. Radius','r!De!N', $
                    '!7j!3!Deff!N (m!E2!N s!E-1!N)']
          onedplotter, ediff_erad, eradc, $
                      titles, onedbounds, outfile

     ENDIF ELSE IF vartoplt EQ 1 THEN BEGIN

          onedbounds = [0,150,0,10000,time/3600.0d0]
          titles = ['ELENGTH vs Equiv. Radius','r!De!N', $
                    'L!De!N (km)']
          onedplotter, elength_erad, eradc, $
                      titles, onedbounds, outfile

     ENDIF ELSE IF vartoplt EQ 2 THEN BEGIN

          onedbounds = [0,150,0,1000,time/3600.0d0]
          titles = ['CONC vs Equiv. Radius','r!De!N', $
                    'C ']
          onedplotter, conc_erad, eradc, $
                      titles, onedbounds, outfile
       
     ENDIF ELSE IF vartoplt EQ 3 THEN BEGIN

          onedbounds = [0,150,0,10000,time/3600.0d0]
          titles = ['TAU vs Equiv. Radius','r!De!N', $
                    '!7s!3!De!N (h)']
          onedplotter, tau_erad, eradc, $
                      titles, onedbounds, outfile
 
     ENDIF ELSE IF vartoplt EQ 4 THEN BEGIN

          onedbounds = [0,150,0,100,time/3600.0d0]
          titles = ['LAMBDA vs Equiv. Radius','r!De!N', $
                    '!7K!3!De!N ']
          onedplotter, lambda_erad, eradc, $
                      titles, onedbounds, outfile       

     ENDIF
 
ENDIF ELSE BEGIN

     IF vartoplt EQ 0 THEN BEGIN
       
          print, clevels   
          titles = ['!7j!3!Deff!N (m!E2!N s!E-1!N)',xax,yax] 
          twodplotter, ediff_xy, x, y, titles, pcolors, $
                        clevels, bounds, outfile

     ENDIF ELSE IF vartoplt EQ 1 THEN BEGIN

          titles = ['L!De!N (km)',xax,yax] 
          twodplotter, elength_xy, x, y, titles, pcolors, $
                        clevels, bounds, outfile
      
     ENDIF

ENDELSE

lcount = lcount + 1

ENDFOR

; compute the average effective diffusivity
ediff_erad_avg   = DBLARR(neradmax)
elength_erad_avg = DBLARR(neradmax)
ediff_erad_avg   = TOTAL(  ediff_erad_hov, 2) / (double(lcount-1))
elength_erad_avg = TOTAL(elength_erad_hov, 2) / (double(lcount-1))

OPENW, 5, dir + infile + "meankeff.dat"
FOR i = 0, neradmax-1 DO BEGIN
   PRINTF, 5, eradc_avg(i), ediff_erad_avg(i), elength_erad_avg(i)
ENDFOR
CLOSE, 5 

PRINT, "** kefferad_plot: creating mean+multiple plot for all times..."

;
; Plot the mean profile for all times
;

;
; Open the correct output file
;

onedbounds = [0,150,0,max(clevels),-999]
outfile = dir+infile+'kefferadavg.ps'
titles = ['MEAN ED vs Equiv. Radius','r!De!N', '!7j!3!Deff!N (m!E2!N s!E-1!N)']
onedplotter, ediff_erad_avg, eradc_avg, $
             titles, onedbounds, outfile

onedbounds = [0,150,0,max(clevels),-999]
outfile = dir+infile+'elengtheradavg.ps'
titles = ['MEAN ELENGTH vs Equiv. Radius','r!De!N', 'L!De!N (km)']
onedplotter, elength_erad_avg, eradc_avg, $
             titles, onedbounds, outfile

bounds = [0,100,0,48,-999]
outfile = dir+infile+'keffhov.ps'
titles = ['!7j!3!Deff!N (m!E2!N s!E-1!N)','r!De!N (km)', 'time (hrs.)']
hovplotter, ediff_erad_hov, eradc_avg, time_array, $
             titles, pcolors, clevels_ediff, bounds, outfile

outfile = dir+infile+'elengthhov.ps'
titles = ['L!De!N (km)','r!De!N (km)', 'time (hrs.)']
hovplotter, elength_erad_hov, eradc_avg, time_array, $
             titles, pcolors, clevels_elength, bounds, outfile

outfile = dir+infile+'conchov.ps'
titles = ['C ','r!De!N (km)', 'time (hrs.)']
hovplotter, conc_erad_hov, eradc_avg, time_array, $
             titles, pcolors, clevels_conc, bounds, outfile

outfile = dir+infile+'tauhov.ps'
titles = ['!7s!3!De!N (h)','r!De!N (km)', 'time (hrs.)']
hovplotter, tau_erad_hov, eradc_avg, time_array, $
             titles, pcolors, clevels_tau, bounds, outfile

outfile = dir+infile+'lambdahov.ps'
titles = ['!7K!3!De!N ','r!De!N (km)', 'time (hrs.)']
hovplotter, lambda_erad_hov, eradc_avg, time_array, $
             titles, pcolors, clevels_lambda, bounds, outfile


PRINT, "** kefferad_plot: DONE. successful exit."
END
