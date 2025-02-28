;
; EPFLUX.PRO
; Calculate EP flux vectors 
;

PRINT, "Running EPFLUX.PRO..."

dir1 = ''
dir2 = ''
var = ''
temp = ''
outfile = ''
infilene = ''
infile1 = ''
infile2 = ''
infile3 = ''
infile4 = ''
infile5 = ''
compstring = ''

nazm = 0.0
pltyp = 0.0
slice = 0.0
nmult = 0
nfiles = 192
;nfiles = 15
nazm = 50
nr = 192
nl = 19

epflux_r = FLTARR(nr,nl)
epflux_theta = FLTARR(nr,nl)

mean_int_last = FLTARR(nr,nl)
eddy_int_last = FLTARR(nr,nl)

mean_integral = FLTARR(nr,nl)
eddy_integral = FLTARR(nr,nl)

actual_mom1 = FLTARR(nr,nl)
actual_mom2 = FLTARR(nr,nl)

; for files before restart up to t68 (or 68*15/60=17 h)
dir1 = "/home/eric/scratch3/modelruns/lappem-isen-1.1/h00t15ac/"
; for files after restart (must reset time in these files)
dir2 = "/home/eric/scratch3/modelruns/lappem-isen-1.2/h00t15ac_rs/"

infilene = "h00t15ac"

dir1 = STRCOMPRESS(STRING(dir1), /REMOVE_ALL)
dir2 = STRCOMPRESS(STRING(dir2), /REMOVE_ALL)
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)

FOR l = 0, nfiles DO BEGIN

IF l LT 69 THEN BEGIN
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
ENDIF ELSE BEGIN
    lstring = STRCOMPRESS(STRING(l-68), /REMOVE_ALL)
    if (l-68 lt 100) then begin
       if (l-68 lt 10) then begin
         tstamp = 't0' + lstring
       endif else begin
         tstamp = 't' + lstring
       endelse
    endif else begin
         tstamp = lstring
    endelse
ENDELSE


;
;  OPEN THE INPUT FILES
; 

IF l LT 69 THEN BEGIN
   dir = dir1
ENDIF ELSE BEGIN
   dir = dir2
ENDELSE

infile1 = dir + infilene + 'u.' + tstamp 
PRINT, "** angbudget.pro: opening and reading " + infile1 + "..."
OPENR, 1, infile1, ERR = err
if (err NE 0) then begin 
   PRINT, "** lappem_isen.pro: file not found, exiting..."  
   break
endif

infile2 = dir + infilene + 'v.' + tstamp 
PRINT, "** angbudget.pro: opening and reading " + infile2 + "..."
OPENR, 2, infile2, ERROR = err
if (err NE 0) then begin 
   PRINT, "** lappem_isen.pro: file not found, exiting..."  
   break
endif

infile3 = dir + infilene + 'q.' + tstamp 
PRINT, "** angbudget.pro: opening and reading " + infile3 + "..."
OPENR, 3, infile3, ERROR = err
if (err NE 0) then begin 
   PRINT, "** lappem_isen.pro: file not found, exiting..."  
   break
endif

infile4 = dir + infilene + 'm.' + tstamp 
PRINT, "** angbudget.pro: opening and reading " + infile4 + "..."
OPENR, 4, infile4, ERROR = err
if (err NE 0) then begin 
   PRINT, "** lappem_isen.pro: file not found, exiting..."  
   break
endif

infile5 = dir + infilene + 'p.' + tstamp 
PRINT, "** angbudget.pro: opening and reading " + infile5 + "..."
OPENR, 5, infile5, ERROR = err
if (err NE 0) then begin 
   PRINT, "** lappem_isen.pro: file not found, exiting..."  
   break
endif

;
;  READ THE DATA INTO MEMORY
;

READF, 1, temp
READF, 1, temp
READF, 1, nx, ny, nl
nl = FIX(nl)
thl = DBLARR(nl)
READF, 1, xmin, xmax, ymin, ymax
FOR k = 0, nl-1 DO BEGIN
   READF, 1, dum1
   PRINT, dum1
   thl(k) = dum1
ENDFOR
READF, 1, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
ufld = DBLARR(nx+1,ny+1,nl)
y = DBLARR(ny+1)
x = DBLARR(nx+1)
xCenter = fix(nx / 2.0) 
yCenter = fix(ny / 2.0) 
timehour = timesec/3600.
deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, ny DO BEGIN
y(j) = j*deltay - 0.5d*(ymax-ymin)
FOR i = 0, nx DO BEGIN
x(i) = i*deltax - 0.5d*(xmax-xmin)
	READF, 1, dum1
	ufld(i,j,k) = dum1
 
ENDFOR
ENDFOR
ENDFOR

READF, 2, temp
READF, 2, temp
READF, 2, nx, ny, nl
nl = FIX(nl)
thl = DBLARR(nl)
READF, 2, xmin, xmax, ymin, ymax
FOR k = 0, nl-1 DO BEGIN
   READF, 2, dum1
   PRINT, dum1
   thl(k) = dum1
ENDFOR
READF, 2, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
vfld = DBLARR(nx+1,ny+1,nl)
timehour = timesec/3600.
deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, ny DO BEGIN
FOR i = 0, nx DO BEGIN
	READF, 2, dum1
	vfld(i,j,k) = dum1
ENDFOR
ENDFOR
ENDFOR

READF, 3, temp
READF, 3, temp
READF, 3, nx, ny, nl
nl = FIX(nl)
thl = DBLARR(nl)
READF, 3, xmin, xmax, ymin, ymax
FOR k = 0, nl-1 DO BEGIN
   READF, 3, dum1
   PRINT, dum1
   thl(k) = dum1
ENDFOR
READF, 3, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
qfld = DBLARR(nx+1,ny+1,nl)
timehour = timesec/3600.
deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, ny DO BEGIN
FOR i = 0, nx DO BEGIN
	READF, 3, dum1
	qfld(i,j,k) = dum1
ENDFOR
ENDFOR
ENDFOR

READF, 4, temp
READF, 4, temp
READF, 4, nx, ny, nl
nl = FIX(nl)
thl = DBLARR(nl)
READF, 4, xmin, xmax, ymin, ymax
FOR k = 0, nl-1 DO BEGIN
   READF, 4, dum1
   PRINT, dum1
   thl(k) = dum1
ENDFOR
READF, 4, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
mfld = DBLARR(nx+1,ny+1,nl)
;timehour = timesec/3600.
timehour = timesec/3600.0 + 17.0
deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, ny DO BEGIN
FOR i = 0, nx DO BEGIN
	READF, 4, dum1
	mfld(i,j,k) = dum1/9.81d
ENDFOR
ENDFOR
ENDFOR

READF, 5, temp
READF, 5, temp
READF, 5, nx, ny, nl
nl = FIX(nl)
thl = DBLARR(nl)
READF, 5, xmin, xmax, ymin, ymax
FOR k = 0, nl-1 DO BEGIN
   READF, 5, dum1
   PRINT, dum1
   thl(k) = dum1
ENDFOR
READF, 5, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
gfld = DBLARR(nx+1,ny+1,nl)
;timehour = timesec/3600.
timehour = timesec/3600.0 + 17.0
deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, ny DO BEGIN
FOR i = 0, nx DO BEGIN
	READF, 5, dum1
	gfld(i,j,k) = dum1
ENDFOR
ENDFOR
ENDFOR

; calculate pressure

nlh=nl+1
ptop=10600.0   ;Pa
pres=FLTARR(nx+1,ny+1,nlh)
presl = FLTARR(nx+1,ny+1,nl)
thh=FLTARR(nlh)
thh=[360,351,343,336,330,325,321,318,316, $
     314,312,310,308,306,304,302,301,300,299,298]

pres(0)=ptop
FOR i = 1, nlh-1 DO BEGIN
   pres(*,*,i) = pres(*,*,i-1) - mfld(*,*,i-1)* $
             (thh(i)-thh(i-1))*9.81d
ENDFOR
FOR i = 0, nl-1 DO BEGIN
   presl(*,*,i) = (pres(*,*,i+1)+pres(*,*,i))/2.0
ENDFOR
;stop


nr = FIX(FLOAT(nx)/2.0)
deltar = deltax

CLOSE,1,2,3, 4,5

; 
; CONVERT U & V TO CYLINDRICAL COMPONENTS
;

PRINT, '% Converting 3D winds/wind z-gradients to ' $
          + 'cylindrical coordinates (r,azm)'

CARTESIANTOCYLIND, ufld, vfld, vr, vt, $
                   deltax, deltay, nazm, $
                   xCenter, yCenter;
;  INTERPOLATE TO CYLINDRICAL GRID
;
INTERPTOCYLIND, qfld, q_polar, $
                nr, nazm, nl, $
       xCenter, $
       yCenter, $
       xPositions, $
       yPositions

INTERPTOCYLIND, vr, vr_polar, $
                nr, nazm, nl, $
       xCenter, $
       yCenter, $
       xPositions, $
       yPositions

INTERPTOCYLIND, vt, vt_polar, $
                nr, nazm, nl, $
       xCenter, $
       yCenter, $
       xPositions, $
       yPositions

INTERPTOCYLIND, mfld, m_polar, $
                nr, nazm, nl, $
       xCenter, $
       yCenter, $
       xPositions, $
       yPositions

INTERPTOCYLIND, gfld, g_polar, $
                nr, nazm, nl, $
       xCenter, $
       yCenter, $
       xPositions, $
       yPositions

INTERPTOCYLIND, presl, pres_polar, $
                nr, nazm, nl, $
       xCenter, $
       yCenter, $
       xPositions, $
       yPositions


;
;  SEPARATE OUT AZIMUTHAL MEAN AND PERTURBATIONS
;

q_mean = FLTARR(nr,nl)
vr_mean = FLTARR(nr,nl)
q_massmean = FLTARR(nr,nl)
vr_massmean = FLTARR(nr,nl)
vt_mean = FLTARR(nr,nl)
m_mean = FLTARR(nr,nl)
q_prime = FLTARR(nr,nazm,nl)
vr_prime = FLTARR(nr,nazm,nl)
vt_prime = FLTARR(nr,nazm,nl)
q_massmean = AZIMUTHALMEAN(m_polar*q_polar, PERTURBATION = q_prime)
vr_massmean = AZIMUTHALMEAN(m_polar*vr_polar, PERTURBATION = vr_prime)
q_mean = AZIMUTHALMEAN(q_polar, PERTURBATION = q_prime)
vr_mean = AZIMUTHALMEAN(vr_polar, PERTURBATION = vr_prime)
m_mean = AZIMUTHALMEAN(m_polar, PERTURBATION = m_prime)
q_massmean = q_massmean / m_mean
vr_massmean = vr_massmean / m_mean
vt_mean = AZIMUTHALMEAN(vt_polar, PERTURBATION = vt_prime)

;
;  COMPUTE MASS PRIME TERMS
;
FOR k = 0, nl-1 DO BEGIN
FOR j = 0, nazm-1 DO BEGIN
FOR i = 0, nr-1 DO BEGIN
    vr_prime(i,j,k) = vr_polar(i,j,k) - vr_massmean(i,k)
    q_prime(i,j,k) = q_polar(i,j,k) - q_massmean(i,k)
ENDFOR
ENDFOR
ENDFOR
print, 'max ur_prime, q_prime = ', max(vr_prime), max(q_prime)
print, 'max vr_massmean, q_massmean = ', max(vr_massmean), max(q_massmean)
print, 'max vr_mean, q_mean = ', max(vr_mean), max(q_mean) 

;
; COMPUTE BUDGET TERMS (TENDENCIES)
;

rad = FINDGEN(nr)*deltar

; -------------MEAN TERM

mean = DBLARR(nr,nl)
temp1 = DBLARR(nr,nl)

FOR j = 0, nl-1 DO BEGIN
FOR i = 0, nr-1 DO BEGIN
    temp1(i,j) = rad(i) * vt_mean(i,j) + 0.5 * fcor * rad(i)*rad(i)
ENDFOR
ENDFOR

dtemp1dr = FLTARR(nr,nl)
FIRST_DERIVATIVE2D, temp1, deltar, deltay, 1, dtemp1dr

mean = -1.0 * vr_massmean * dtemp1dr

; -------------EDDY TERM

eddy = DBLARR(nr,nl)
temp3 = DBLARR(nr,nl)
temp3_prime = DBLARR(nr,nl)
temp2 = DBLARR(nr,nazm,nl)

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, nr-1 DO BEGIN
FOR i = 0, nazm-1 DO BEGIN
    temp2(j,i,k) = m_polar(j,i,k) * q_prime(j,i,k) * vr_prime(j,i,k)
ENDFOR
ENDFOR
ENDFOR

temp3_mean = AZIMUTHALMEAN(temp2, PERTURBATION = temp3_prime)

FOR j = 0, nr-1 DO BEGIN
    eddy(j,*) = -1.0 * rad(j) * temp3_mean(j,*)
ENDFOR

;
;  TRAPEZOIDAL INTEGRATION (REWRITE)

deltat = 900.0    ;seconds
IF l GE 1 THEN BEGIN
  
    eddy_integral = eddy_integral + 0.5*(eddy_int_last + $
                                         eddy * deltat   )
    mean_integral = mean_integral + 0.5*(mean_int_last + $
                                         mean * deltat   )

ENDIF

eddy_int_last = eddy * deltat
mean_int_last = mean * deltat 


; ------ACTUAL
;
IF l EQ 0 THEN BEGIN
   actual_mom1 = temp1
ENDIF
IF l EQ nfiles THEN BEGIN
   actual_mom2 = temp1
ENDIF

; --------------------EP FLUX VECTOR TERMS


; calculate mean and perts
mu_mean = FLTARR(nr,nl)
mu_prime = FLTARR(nr,nazm,nl)
ang_mean = FLTARR(nr,nl)
ang_prime = FLTARR(nr,nazm,nl)
pres_mean = FLTARR(nr,nl)
pres_prime = FLTARR(nr,nazm,nl)
mont_mean = FLTARR(nr,nl)
mont_prime = FLTARR(nr,nazm,nl)
dmont_primedphi = FLTARR(nr,nazm,nl)

mu_mean = AZIMUTHALMEAN(m_polar*vr_polar, PERTURBATION = mu_prime)
mont_mean = AZIMUTHALMEAN(g_polar, PERTURBATION = mont_prime)
pres_mean = AZIMUTHALMEAN(pres_polar, PERTURBATION = pres_prime)

dazm = 0.125664   ; 360/50 in radians
first_derivative_cyl3d, mont_prime, $
                        nr, nazm, nl, $
                        deltar, dazm, dl, $
                        2, dmont_primedphi

epflux_rmean = FLTARR(nr,nl)
epflux_rprime = FLTARR(nr,nazm,nl)
epflux_tmean = FLTARR(nr,nl)
epflux_tprime = FLTARR(nr,nazm,nl)

FOR i = 0, nr-1 DO BEGIN
ang_prime(i,*,*) = rad(i) * vt_prime(i,*,*) + $
                    0.5 * fcor * rad(i) * rad(i)
ENDFOR


epflux_rmean = AZIMUTHALMEAN(mu_prime*ang_prime, $
               PERTURBATION = mu_prime)
epflux_tmean = AZIMUTHALMEAN(pres_prime*dmont_primedphi, $ 
               PERTURBATION = mu_prime)

epflux_r = epflux_r + epflux_rmean
epflux_theta = epflux_theta + epflux_tmean

print, "DONE WITH ONE TIME"

ENDFOR

epflux_r = epflux_r / float(nfiles)
epflux_theta = epflux_theta / float(nfiles)

; reduce vectors for plotting

epflux_r_plt=FLTARR(nr,nl)
epflux_t_plt=FLTARR(nr,nl)

FOR i = 0, nr-1 DO BEGIN
IF i MOD 3 EQ 0 THEN BEGIN
   epflux_r_plt(i,*) = epflux_r(i,*)
   epflux_r_plt(i,*) = epflux_r(i,*)
ENDIF
ENDFOR



;
;  PLOTTING
;

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
;rbw_red_color = [254, 125, 077, 000, 070, 070, 070, 070, 050, 200, 255, $
;238, 238, 220, 238, 220, 255]
;rbw_gre_color = [254, 095, 000, 000, 070, 150, 235, 255, 220, 255, 255, $
;205, 154, 105, 069, 075, 000]
;rbw_blu_color = [254, 186, 186, 205, 255, 255, 255, 160, 050, 050, 000, $
;000, 000, 000, 000, 075, 000]

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

rad = rad / 1000.0 

;  -------------------VR

   titles  = ["MEAN (m!U2!N s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = [-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
   bounds = [min(rad),100.0,min(thl),max(thl),timehour]
   plotname = './vr_term.ps'
   twodplotterposneg, vr_mean, rad, thl, titles, p_colors, $
             clevels, bounds, plotname



;  -------------------VT


   titles  = ["MEAN (m!U2!N s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = [-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
   bounds = [min(rad),100.0,min(thl),max(thl),timehour]
   plotname = './vt_term.ps'
   twodplotterposneg, vt_mean, rad, thl, titles, p_colors, $
             clevels, bounds, plotname


;  -------------------MEAN TERM

   
 
   titles  = ["MEAN (10!U5!N m!U2!N s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 0.5*[-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
   bounds = [min(rad),100.0,min(thl),max(thl),timehour]
   plotname = './mean_term.ps'
   twodplotterposneg, mean_integral/1.0e+5, rad, thl, titles, p_colors, $
             clevels, bounds, plotname


;  -------------------EDDY TERM

  
   titles  = ["EDDY (10!U5!N m!U2!N s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 0.5*[-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
   bounds = [min(rad),100.0,min(thl),max(thl),timehour]
   plotname = './eddy_term.ps'
;   twodplotterposneg, eddy_integral/1.0e+5, rad, thl, titles, p_colors, $
;             clevels, bounds, plotname
   twodplotterposneg_vec, eddy_integral/1.0e+5, epflux_r_plt, epflux_t_plt, $
              rad, thl, rad, thl,  $
              titles, p_colors, $
              clevels, bounds, plotname

;  -------------------MEAN+EDDY TERM

  
   titles  = ["MEAN+EDDY (10!U5!N m!U2!N s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 0.5*[-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
   bounds = [min(rad),100.0,min(thl),max(thl),timehour]
   plotname = './meanpeddy_term.ps'
   twodplotterposneg, (mean_integral+eddy_integral)/1.0e+5, rad, thl, $
              titles, p_colors, $
             clevels, bounds, plotname

;  -------------------ACTUAL TERM

  
   titles  = ["ACTUAL (10!U5!N m!U2!N s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 0.5*[-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
   bounds = [min(rad),100.0,min(thl),max(thl),timehour]
   plotname = './actual_term.ps'
   twodplotterposneg, (actual_mom2-actual_mom1)/1.0e+5, $
              rad, thl, titles, p_colors, $
             clevels, bounds, plotname

; --------------------EP FLUX VECTORS

;   VELOVECT, epflux_r, epflux_theta



END
