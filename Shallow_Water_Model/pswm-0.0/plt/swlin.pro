;
; SWLIN.PRO
; This is a program that plots the solution
; to the linearized SW solution for a resting
; basic state in cyclindrical coordinates
;

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

clevels_d = 0.1*[-30,-20,-10,-7,-5,-2,-0.5,-0.3,0,0.3,0.5,2,5,7,10,20,30]
plotcolor = 0
if (plotcolor LT 0.5) then begin
   pcolors = c2w_plot_colors
endif else begin
   pcolors = bwl_plot_colors
endelse

;
; INITIAL PARAMETERS
;

pi  = acos(-1)
kr  = 0.00004        ; radial wavenumber (m^-1)

omegabar = 0.0005     ; basic state angular velocity (s^-1)
m = 1                  ; azimuthal wavenumber (rad^-1)
f = 0.000037              ; Coriolis parameter (s^-1)
;f = 0.0
hbar = 4288.0           ; the mean fluid depth (m)
;hamp = 70.0             ; perturbation height amplitude (m)
hamp = 10.0
p = f * hamp / hbar  ; perturbation zeta amplitude (s^-1)
g = 9.81                ; gravity acceleration (m/s^2)

nr   = 100
nazm = 100
nt   = 20

tfac = 60.0       ; seconds in output interval
Lr   = 300000.0        ; domain length (m)
dr   = Lr / float(nr)
dt   = 3.0          ; minutes
dazm = 2.0 * pi / float(nazm-1)

;
; ALLOCATE ARRAYS
;

r     = FINDGEN(nr) * dr + dr
phi   = FINDGEN(nazm) * dazm
t     = FINDGEN(nt) * dt
hfld  = FINDGEN(nr,nazm,nt)
zfld  = FINDGEN(nr,nazm,nt)
dfld  = FINDGEN(nr,nazm,nt)
qfld  = FINDGEN(nr,nazm,nt)
ufld  = FINDGEN(nr,nazm,nt)    ; the radial velocity
vfld  = FINDGEN(nr,nazm,nt)    ; the tangential velocity

;
; CALCULATE THE REAL PART OF SOLUTION 
;

FOR k = 0, nt-1 DO BEGIN
FOR j = 0, nr-1 DO BEGIN
FOR i = 0, nazm-1 DO BEGIN

bj    = BESELJ(kr * r(j), float(m))
by    = BESELY(kr * r(j), float(m))
bj1   = BESELJ(kr * r(j), float(m-1))
by1   = BESELY(kr * r(j), float(m-1))
freq = ( f^2 + g * hbar * kr^2 )^(0.5)

hfld(j,i,k) = hamp * ( bj * cos( float(m) * phi(i) - freq * t(k) * tfac ) - $
                       by * sin( float(m) * phi(i) - freq * t(k) * tfac ) )

zfld(j,i,k) = (hamp/hbar) * ( f + 2.0 * omegabar) * ( bj * cos( float(m) * phi(i) - freq * t(k) * tfac ) - $
                       by * sin( float(m) * phi(i) - freq * t(k) * tfac ) )

dfld(j,i,k) = (-hamp / hbar) * (freq - omegabar * float(m)) * ( by * cos( float(m) * phi(i) - freq * t(k) * tfac ) + $
                       bj * sin( float(m) * phi(i) - freq * t(k) * tfac ) )

ufld(j,i,k) = (hamp / (kr^2 * hbar)) * ( -1.0*(float(m) /r(j)) * (freq + f) * bj * sin( float(m) * phi(i) - freq * t(k) * tfac ) - $
                 (float(m) /r(j)) * (freq + f) * by * cos( float(m) * phi(i) - freq * t(k) * tfac ) + $
                 freq * kr * bj1 * sin( float(m) * phi(i) - freq * t(k) * tfac ) + $
                 freq * kr * by1 * cos( float(m) * phi(i) - freq * t(k) * tfac )  )       

vfld(j,i,k) = (hamp / (kr^2 * hbar)) * ( (float(m) /r(j)) * (freq + f) * bj * cos( float(m) * phi(i) - freq * t(k) * tfac ) - $
                 (float(m) /r(j)) * (freq + f) * by * sin( float(m) * phi(i) - freq * t(k) * tfac ) - $
                 f * kr * bj1 * cos( float(m) * phi(i) - freq * t(k) * tfac ) + $
                 f * kr * by1 * sin( float(m) * phi(i) - freq * t(k) * tfac )  )

qfld(j,i,k) = hbar * ( f + zfld(j,i,k) ) / ( hbar + hfld(j,i,k) )

ENDFOR
ENDFOR
ENDFOR

;dfld = dfld * 10000.0    ;convert for plotting

;
; CONVERT TO CARTESIAN COMPONENTS FOR VECTOR PLOTS
;

ufldc = FLTARR(nr,nazm,nt)
vfldc = FLTARR(nr,nazm,nt)

FOR k = 0, nt-1 DO BEGIN
FOR j = 0, nr-1   DO BEGIN
FOR i = 0, nazm-1 DO BEGIN

ufldc(j,i,k) = ufld(j,i,k) * cos( phi(i) ) - vfld(j,i,k) * sin( phi(i) )
vfldc(j,i,k) = ufld(j,i,k) * sin( phi(i) ) + vfld(j,i,k) * cos( phi(i) )

ENDFOR
ENDFOR
ENDFOR


;
; NOW POLAR PLOT THE SOLUTION
;


dir = "./"
SET_PLOT, 'ps'
DEVICE, /INCHES, /COLOR, $
 XSIZE = 6.5, YSIZE = 6.0, $
; XSIZE = 8.0, YSIZE = 6.0, $
 XOFFSET = 0.1, YOFFSET = 2.5, $
 ;XSIZE = 7.5, YSIZE = 10.5, $
 ;XOFFSET  = 0.3, YOFFSET = 0.3, $
 FILENAME = dir + "swlin.ps"

!P.MULTI = [0,1,1]

;FOR k = 0, 10 DO BEGIN
k = 0

timestr = STRCOMPRESS(STRING(t(k)),/REMOVE_ALL)
ptitle = "Perturbation Height at t = " + timestr + " h."
ptitle = "!7d!3 (s!E-1!N)" 

hfldr = hfld(*,*,k)
hfld_reform = REFORM(hfldr)
zfldr = zfld(*,*,k)
zfld_reform = REFORM(zfldr)
dfldr = dfld(*,*,k)
dfld_reform = REFORM(dfldr)
ufldr = ufld(*,*,k)
ufld_reform = REFORM(ufldr)
vfldr = vfld(*,*,k)
vfld_reform = REFORM(vfldr)
qfldr = qfld(*,*,k)
qfld_reform = REFORM(qfldr)

ufldcr = ufldc(*,*,k)
vfldcr = vfldc(*,*,k)
ufldc_reform = REFORM(ufldcr)
vfldc_reform = REFORM(vfldcr)

;ptitle = "h' (m)" 
;PCONTOUR_RTHETA, hfld_reform(0:99,*), 100, nazm, dr, 10, 10, -10, ptitle
;ptitle = "!7f'!3 (s!E-1!N)" 
;PCONTOUR_RTHETA, zfld_reform(0:99,*), 100, nazm, dr, 10, 2e-6, -2e-6, ptitle
;ptitle = "!7d'!3 (s!E-1!N)"
;ptitle = "!7d'!3 (x 10!E-4!N s!E-1!N)" 
;PCONTOUR_RTHETA, dfld_reform(0:45,*), 46, nazm, dr, 9, 1e-4, -1e-4, ptitle
;PCONTOUR_RTHETA, dfld_reform(0:99,*), 100, nazm, dr, 10, 1e-5, -1e-5, ptitle
;ptitle = "u (m s!E-1!N)" 
;PCONTOUR_RTHETA, ufld_reform(0:99,*), 100, nazm, dr, 10, 2.0, -2.0, ptitle
;ptitle = "v (m s!E-1!N)" 
;PCONTOUR_RTHETA, vfld_reform(0:99,*), 100, nazm, dr, 10, 2.0, -2.0, ptitle
;ptitle = "P (s!E-1!N)" 
;PCONTOUR_RTHETA, qfld_reform(0:99,*), 100, nazm, dr, 100, 2.5e-4, -2.5e-4, ptitle
;ptitle = "!7d'!3 (s!E-1!N)" 
;PCONTOUR_RTHETA_COLOR, dfld_reform(0:99,*), 100, nazm, dr, clevels_d, 1e-4, -1e-4, ptitle, pcolors

; VECTOR OVERLAYS
;ptitle = "h' (m)" 
;PCONTOUR_RTHETA_VEC, hfld_reform, ufldc_reform, vfldc_reform, $
;                     100, nazm, dr, 10, 5, -5, 0.35, ptitle

;PCONTOUR_RTHETA_VEC, hfld_reform(0:99,*), ufldc_reform(0:99,*), vfldc_reform(0:99,*), $
;                     100, nazm, dr, 10, 5, -5, 1.0, ptitle

clevels = 0.1*[-500,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
titles = ["h' (m)","x (km)", "y (km)"]
outfile2='./swlin.ps'
bounds = [-200,200,-200,200,0.0/3600.0d0]
scalefac = 0.20
PCONTOUR_RTHETA_VEC2, hfld_reform, ufldc_reform, vfldc_reform, nr, nazm, dr, $
      clevels, titles, pcolors, scalefac, outfile2, bounds


;ONEDPLOTTER, DATA, XPLOT, TITLES, BOUNDS, PLOTNAME


;PCONTOUR_RTHETA, field_reform, nr, nazm, dr, 10, ptitle

;ENDFOR

DEVICE, /CLOSE


END
