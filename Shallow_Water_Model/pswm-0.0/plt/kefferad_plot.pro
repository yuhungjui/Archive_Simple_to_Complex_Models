PRINT, 'Running IDL kefferad_plot.pro'
PRINT, 'Effective Diffusivity PLOTTING Routine'
PRINT, 'for Pseudospectral Model Output'
PRINT, 'Equivalent Radius Version'

lbegin = 2
tstamp = ''
close, 15
dir = ''
lstring = ''
timestr = ''
lcount = 0
len = 600000.0d0   ;default value, this is reset

neradmax = 200     ; should be set to nerad in .dat files

READ, dir, PROMPT="Enter output directory path:"
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
infile = 'psndbp.in.'
READ, infile, PROMPT="Enter infile: "
infile = STRCOMPRESS(STRING(infile+'.in.'), /REMOVE_ALL)
READ, lbegin, PROMPT="Enter beginning file number: "
lbegin = FIX(lbegin)
READ, lend, PROMPT="Enter end file number: "
lend = FIX(lend)
READ, logplot, PROMPT="Log Plot (0-yes, 1-no): "
logplot = FIX(logplot)
READ, oneortwod, PROMPT="oned or twod plots (0-1D, 1-2D): "
oneortwod = FIX(oneortwod)
READ, lengthored, PROMPT="twod plot type (0-ed, 1-elength): "
lengthored = FIX(lengthored)

nfiles = lend - lbegin + 1

; 2D radius-time arrays (hovmoller and mean plots)
ediff_erad_hov   = DBLARR(neradmax,nfiles)
elength_erad_hov = DBLARR(neradmax,nfiles)
time_array       = DBLARR(nfiles)
eradc_avg = DBLARR(neradmax)
eradc_avg = 500.0d0 + FINDGEN(neradmax) * 1000.0d0
eradc_avg = REVERSE(eradc_avg)
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

c2w_plot_colors = c32(c2w_red_color, c2w_gre_color, c2w_blu_color, 1)
rbw_plot_colors = c32(rbw_red_color, rbw_gre_color, rbw_blu_color, 1)
NLEV1 = FIX(SIZE(c2w_plot_colors, /DIMENSIONS))
NLEV  = NLEV1(0)
black_cons = c32(intarr(NLEV), intarr(NLEV), intarr(NLEV), 1)

; for midrange ED
clevels = [0,100,200,300,400,500,600,700,800,900,1000,1100, $
           1300,1500,1700,1900,2100]

; for lowrange ED
;clevels = [0,100,250,500,750,1000,1500,2000,2500,3000,3500, $
;           4000,4500,5000,5500,6000,6500]

; for highrange ED
;clevels = [0,100,500,1000,1500,2000,3000,4000,5000, $
;           6000,7500,9000,11000,13000,15000,17500,20000]

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
readf, 15, nerad, nx, ny

nerad = FIX(nerad)
nx = FIX(nx)
ny = FIX(ny)

readf, 15, deltax, deltay, nu, time, len
readf, 15, temp
timestr = STRCOMPRESS(STRING(time/3600.0d0), /REMOVE_ALL)
time_array(lcount) = time/3600.0d0

eradc = DBLARR(nerad)
ediff_erad = DBLARR(nerad)
elength_erad = DBLARR(nerad)
ediff_xy = DBLARR(nx,ny)
elength_xy = DBLARR(nx,ny)
x = DBLARR(nx)
y = DBLARR(ny)

FOR k = 0, nerad-1 DO BEGIN
     readf, 15, t1, t2, t3
     eradc(k) = t1 
     ediff_erad(k) = t2
     elength_erad(k) = t3
ENDFOR

;
;  THIS INTERPOLATES THE IRREGULAR ERAD TO 
;  LINEAR ERAD ARRAY
;

temp_edifferad   = INTERPOL(  ediff_erad, eradc, eradc_avg)
temp_elengtherad = INTERPOL(elength_erad, eradc, eradc_avg)

ediff_erad_hov(*,lcount)   = temp_edifferad
elength_erad_hov(*,lcount) = temp_elengtherad
 
deltax = len/double(nx)
deltay = len/double(ny)

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

IF (logplot EQ 0) THEN BEGIN

    IF oneortwod EQ 0 THEN BEGIN

       IF lengthored EQ 0 THEN BEGIN

          onedbounds = [0,150,0,10,time/3600.0d0]
          titles = ['ln ED vs Equiv. Radius','r!De!N', $
                    'ln !7j!3!Deff!N (m!E2!N s!E-1!N)']
          onedplotter, alog(ediff_erad), eradc/1000.0d0, $
                      titles, onedbounds, outfile

       ENDIF ELSE IF lengthored EQ 1 THEN BEGIN

          onedbounds = [0,150,0,10,time/3600.0d0]
          titles = ['ln ELENGTH vs Equiv. Radius','r!De!N', $
                    'ln L!De!N (km)']
          onedplotter, alog(elength_erad), eradc/1000.0d0, $
                      titles, onedbounds, outfile
       
       ENDIF
 
    ENDIF ELSE BEGIN

       IF lengthored EQ 0 THEN BEGIN
          
          titles = ['ln !7j!3!Deff!N (m!E2!N s!E-1!N)',xax,yax] 
          twodplotter, alog(ediff_xy), x, y, titles, rbw_plot_colors, $
                        clevels, bounds, outfile

       ENDIF ELSE IF lengthored EQ 1 THEN BEGIN

          titles = ['ln L!De!N (km)',xax,yax] 
          twodplotter, alog(elength_xy), x, y, titles, rbw_plot_colors, $
                        clevels, bounds, outfile
      
       ENDIF

    ENDELSE

ENDIF ELSE BEGIN

    IF oneortwod EQ 0 THEN BEGIN

       IF lengthored EQ 0 THEN BEGIN

          onedbounds = [0,150,0,max(clevels),time/3600.0d0]
          titles = ['ED vs Equiv. Radius','r!De!N', $
                    '!7j!3!Deff!N (m!E2!N s!E-1!N)']
          onedplotter, ediff_erad, eradc/1000.0d0, $
                      titles, onedbounds, outfile

       ENDIF ELSE IF lengthored EQ 1 THEN BEGIN

          onedbounds = [0,150,0,max(clevels),time/3600.0d0]
          titles = ['ELENGTH vs Equiv. Radius','r!De!N', $
                    'ln L!De!N (km)']
          onedplotter, elength_erad, eradc/1000.0d0, $
                      titles, onedbounds, outfile
       
       ENDIF
 
    ENDIF ELSE BEGIN

       IF lengthored EQ 0 THEN BEGIN

           titles = ['!7j!3!Deff!N (m!E2!N s!E-1!N)',xax,yax] 
           twodplotter, ediff_xy, x, y, titles, rbw_plot_colors, $
                        clevels, bounds, outfile

       ENDIF ELSE IF lengthored EQ 1 THEN BEGIN

           titles = ['L!De!N (km)',xax,yax] 
           twodplotter, elength_xy, x, y, titles, rbw_plot_colors, $
                        clevels, bounds, outfile
      
       ENDIF

    ENDELSE

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
onedplotter, ediff_erad_avg, eradc_avg/1000.0d0, $
             titles, onedbounds, outfile

onedbounds = [0,150,0,max(clevels),-999]
outfile = dir+infile+'elengtheradavg.ps'
titles = ['MEAN ELENGTH vs Equiv. Radius','r!De!N', 'L!De!N (km)']
onedplotter, elength_erad_avg, eradc_avg/1000.0d0, $
             titles, onedbounds, outfile

bounds = [0,100,0,48,-999]
outfile = dir+infile+'keffhov.ps'
titles = ['!7j!3!Deff!N (m!E2!N s!E-1!N)','r!De!N (km)', 'time (hrs.)']
hovplotter, ediff_erad_hov, eradc_avg/1000.0d0, time_array, $
             titles, rbw_plot_colors, clevels, bounds, outfile

outfile = dir+infile+'elengthhov.ps'
titles = ['L!De!N (km)','r!De!N (km)', 'time (hrs.)']
hovplotter, ediff_erad_hov, eradc_avg/1000.0d0, time_array, $
             titles, rbw_plot_colors, clevels, bounds, outfile

PRINT, "** kefferad_plot: DONE. successful exit."
END
