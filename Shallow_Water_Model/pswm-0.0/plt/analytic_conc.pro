PRINT, "Running IDL analytic_conc.pro"
PRINT, "Creates a hovmoller plot of the"
PRINT, "analytic C(r,t) to compare to"
PRINT, "the numerical C(r,t) from the"
PRINT, "pseudospectral model."

;
; ANALYTIC SOLUTION TO DIFFUSION EQUATION
; GIVEN INITIAL GAUSSIAN PROFILE
;

;
; PARAMETERS
;
b     = 150000.0d0        ; m 
cmax  = 1000.0d0          ; dimensionless
kappa = 2500.0d0          ; m^2/s tracer diffusion

;
; SOLUTION ARRAY
;
nt   = 3600
nre  = 1000        
dre  = 500.0d0               ; m
dt   = 60.0d0                ; sec
re   = FINDGEN(nre) * dre 
time = FINDGEN(nt)  * dt
conc_re_analytic   = DBLARR(nre,nt)

FOR k = 0, nt-1  DO BEGIN
FOR j = 0, nre-1 DO BEGIN

conc_re_analytic(j,k) = cmax * ( b * b / (b*b + 4.0d0 * kappa * time(k)) ) $   
                        * exp( - re(j) * re(j) / (b*b + 4.0d0 * kappa * time(k)) )

ENDFOR
ENDFOR

;
; THE INITIAL CONDITION
;
c_re_gauss = DBLARR(nre)
FOR j = 0, nre-1 DO BEGIN
   c_re_gauss(j) = cmax * exp( - re(j) * re(j) / ( b * b ) )    
ENDFOR

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

re   = re / 1000.0d0    
time = time / 3600.0d0
;clevels = [0,50,100,150,200,250,300,350,400,450,500,550,600,700,$
;           800,900,1000,1100]
clevels = [500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900]
bounds = [0,100,0,48,-999]
outfile = 'analytic_conc.ps'
;titles = ['C(r!De!N)','r!De!N (km)', 'time (hrs.)']
titles = ['C','r!De!N (km)', 'time (hrs.)']
hovplotter, conc_re_analytic, re, time, $
            titles, bwh_plot_colors, clevels, bounds, outfile

;
; simple gaussian tracer field
;
bounds = [0,600,0,1200,-999]
;onedplotter, c_re_gauss, re, titles, bounds, outfile
;
PRINT, "** analytic_conc: DONE. successful exit."

END
