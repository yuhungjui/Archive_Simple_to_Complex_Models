; IDL PLOTASYMMETRIES.PRO
; Plot asymmetries h', etc
; Written: 05/29/08
; Eric A. Hendricks

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


;dir = '../run/e7o20hres/'
;dir='/home/eric/scratch3/modelruns/pswm-0.0/e7o20hres/'
;infilene = 'e7o20'

dir = '/usbhdd1/ldata/pswm/ellipse2/'
infilene = 'ellipse'

;dir = '/usbhdd1/ldata/pswm/ring1/'
;infilene = 'ring1'

var = ''
temp = ''
outfile = ''
infile2 = ''
varstring = ''
xmin = 0.0d0
xmax = 0.0d0
fcor = 0.0d0
beta = 0.0d0
timesec = 0.0d0
ymin = 0.0d0
ymax = 0.0d0
dum1 = 0.0d0

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
READ, lbegin, PROMPT="Enter beginning file number: "
READ, lend, PROMPT="Enter ending file number: "
lbegin = FIX(lbegin)
lend = FIX(lend)

nxd = 512
nyd = 512
cval = 205.0    ; m/s
grav = 9.81

;
; BEGIN FILE LOOP
;

FOR l = lbegin, lend DO BEGIN

lstring = STRCOMPRESS(STRING(l), /REMOVE_ALL)

if (l lt 100) then begin

if (l lt 10) then begin
     tstamp = 't0' + lstring
     ostamp = '00' + lstring
     outfile =  'PRIME.'+ostamp+'.ps'
     outfile1 = dir+infilene+'dAdt.'+ostamp+'.ps'
  endif else begin
     tstamp = 't' + lstring
     ostamp = '0' + lstring
     outfile = 'PRIME.'+ostamp+'.ps'
     outfile1 = dir+infilene+'dAdt.'+ostamp+'.ps'
  endelse
endif else begin
  tstamp = lstring
  outfile = 'PRIME.'+tstamp+'.ps'
  outfile1 = dir+infilene+'dAdt.'+tstamp+'.ps'
endelse

;
; OPEN AND READ THE OUTPUT FILE
; 

fieldh = DBLARR(nxd+1,nyd+1)
fieldu = DBLARR(nxd+1,nyd+1)
fieldv = DBLARR(nxd+1,nyd+1)
fieldz = DBLARR(nxd+1,nyd+1)
fieldq = DBLARR(nxd+1,nyd+1)
fieldd = DBLARR(nxd+1,nyd+1)

FOR m = 0, 4 DO BEGIN

IF m EQ 0 THEN varstring = 'p'
IF m EQ 1 THEN varstring = 'u'
IF m EQ 2 THEN varstring = 'v'
IF m EQ 3 THEN varstring = 'z'
IF m EQ 4 THEN varstring = 'd'
infile = dir + infilene + varstring + '.' + tstamp 
PRINT, "** plotasymmetries.pro: opening and reading " + infile + "..."
OPENR, 1, infile, ERROR = err 
if (err NE 0) then begin 
   PRINT, "** plotasymmetries.pro: file not found, exiting..."  
   break
endif
READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
READF, 1, xmin, xmax, ymin, ymax
READF, 1, fcor, beta, cval, timesec
nx = FIX(nx)
ny = FIX(ny)
field = DBLARR(nx+1,ny+1)
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
ENDFOR
ENDFOR
CLOSE, 1

if m EQ 0 then fieldh = field
if m EQ 1 then fieldu = field
if m EQ 2 then fieldv = field
if m EQ 3 then fieldz = field
if m EQ 4 then fieldd = field

ENDFOR

;
; CONVERT P to H FOR DIAGNOSTICS
;

fieldh = (fieldh*cval + cval^2) / grav

;
; COMPUTE PV
;

fieldq = (cval^2/grav)*(fieldz + fcor) / fieldh

;
; DETERMINE CENTER FROM VORTICITY
;

PRINT, "** plotasymmetries.pro: Determining center..."

boxsize = 40
zeta_high = 0.0                ; max zeta 
zeta_mn = FLTARR(nx+1,ny+1)    ; vorticity averaged                          
xCenter_z = 0                  ; integer index of x-center
yCenter_z = 0                  ; integer index of y-center

FOR j = 0, ny-1-boxsize DO BEGIN
FOR i = 0, nx-1-boxsize DO BEGIN

     zeta_mn(i+boxsize/2,j+boxsize/2) = $
                      MEAN(fieldz(i:i+boxsize,j:j+boxsize))

ENDFOR
ENDFOR

tmpvar = 0
zeta_high = MAX(zeta_mn,tmpvar,/NAN)
yCenter_z = (tmpvar/FIX(nx+1))
xCenter_z = tmpvar - ((tmpvar/FIX(nx+1))*FIX(nx+1))  ;Note: IDL drops decimals
xc = xCenter_z
yc = yCenter_z
xc = 256
yc = 256

PRINT, "** plotasymmetries.pro: xCenter, yCenter indeces = ", xc, yc

;
; AZIMUTHAL MEAN, PERTURBATION QUANTITIES ABOUT CENTER
;

PRINT, "** plotasymmetries.pro: Determining azimuthal mean diagnostics..."

nazm = 50
nr = FIX(nx/2.0)
deltar = deltax
r = FINDGEN(nr) * deltar

fieldu_azmean = DBLARR(nr)
fieldu_prime = DBLARR(nr,nazm)
fieldu_polar = DBLARR(nr,nazm)

fieldv_azmean = DBLARR(nr)
fieldv_prime = DBLARR(nr,nazm)
fieldv_polar = DBLARR(nr,nazm)

fieldh_azmean = DBLARR(nr)
fieldh_prime = DBLARR(nr,nazm)
fieldh_polar = DBLARR(nr,nazm)

fieldz_azmean = DBLARR(nr)
fieldz_prime = DBLARR(nr,nazm)
fieldz_polar = DBLARR(nr,nazm)

fieldq_azmean = DBLARR(nr)
fieldq_prime = DBLARR(nr,nazm)
fieldq_polar = DBLARR(nr,nazm)

fieldd_azmean = DBLARR(nr)
fieldd_prime = DBLARR(nr,nazm)
fieldd_polar = DBLARR(nr,nazm)

vr = DBLARR(nx+1,ny+1)
vt = DBLARR(nx+1,ny+1)

CARTESIANTOCYLIND2D, fieldu, fieldv, $
                   vr, vt, $
                   deltax, deltay, nazm, $
		       xc, $
		       yc

INTERPTOCYLIND2D, vr, fieldu_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, vt, fieldv_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldh, fieldh_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldz, fieldz_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldq, fieldq_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

INTERPTOCYLIND2D, fieldd, fieldd_polar, $
                nr, nazm, $
                xc, $
                yc, $
                xPositions, $
                yPositions

fieldu_azmean = AZIMUTHALMEAN2D(fieldu_polar, $
                        perturbation = fieldu_prime)
fieldv_azmean = AZIMUTHALMEAN2D(fieldv_polar, $
                        perturbation = fieldv_prime)
fieldh_azmean = AZIMUTHALMEAN2D(fieldh_polar, $
                        perturbation = fieldh_prime)
fieldz_azmean = AZIMUTHALMEAN2D(fieldz_polar, $
                        perturbation = fieldz_prime)
fieldq_azmean = AZIMUTHALMEAN2D(fieldq_polar, $
                        perturbation = fieldq_prime)
fieldd_azmean = AZIMUTHALMEAN2D(fieldd_polar, $
                        perturbation = fieldd_prime)


;PCONTOUR_RTHETA1_COLOR, fieldv_polar(*,*), nr, nazm, deltar, $
;      clevels, titles, pcolors, outfile, bounds

;PCONTOUR_RTHETA1_COLOR, fieldq_prime(*,*), nr, nazm, deltar, $
;      clevels, titles, pcolors, outfile, bounds

scalefac = 0.015
bounds = [-300,300,-300,300,timesec/3600.0d0]
bounds = [-50,50,-50,50,timesec/3600.0d0]
clevels = 0.00003*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]

;
; convert to rectangular components, but defined on polar grid for plotting
;
fieldu_rect=DBLARR(nr,nazm)
fieldv_rect=DBLARR(nr,nazm)
dazm = 2.0 * 3.14159 / double(nazm-1)
phi   = FINDGEN(nazm) * dazm
FOR i = 0, nr-1 DO BEGIN
FOR j = 0, nazm-1 DO BEGIN
    fieldu_rect(i,j) = fieldu_prime(i,j) * cos( phi(j) ) - fieldv_prime(i,j) * sin( phi(j) ) 
    fieldv_rect(i,j) = fieldu_prime(i,j) * sin( phi(j) ) + fieldv_prime(i,j) * cos( phi(j) )
;    fieldu_rect(i,j) = fieldu_polar(i,j) * cos( phi(j) ) - fieldv_polar(i,j) * sin( phi(j) ) 
;    fieldv_rect(i,j) = fieldu_polar(i,j) * sin( phi(j) ) + fieldv_polar(i,j) * cos( phi(j) ) 
ENDFOR
ENDFOR

xaxis = "x (km)"
yaxis = "y (km)"
titles = ["P' (s!E-1!N)",xaxis,yaxis]
pcolors = c2w_plot_colors

clevels = 0.00003*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
;bounds = [-300,300,-300,300,timesec/3600.0d0]
bounds = [-50,50,-50,50,timesec/3600.0d0]
bounds = [-100,100,-100,100,timesec/3600.0d0]
outfile2=dir+infilene+'P_'+outfile
PCONTOUR_RTHETA_VEC2, fieldq_prime, fieldu_rect, fieldv_rect, nr, nazm, deltar, $
      clevels, titles, pcolors, scalefac, outfile2, bounds

clevels = 0.2*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
titles = ["h' (m)",xaxis,yaxis]
outfile2=dir+infilene+'h_'+outfile
PCONTOUR_RTHETA_VEC2, fieldh_prime, fieldu_rect, fieldv_rect, nr, nazm, deltar, $
      clevels, titles, pcolors, scalefac, outfile2, bounds

clevels = 0.00004*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
titles = ["h' (m)",xaxis,yaxis]
outfile2=dir+infilene+'PP_'+outfile
bounds = [-100,100,-100,100,timesec/3600.0d0]
scalefac = 0.01
titles = ["P (s!E-1!N)",xaxis,yaxis]
PCONTOUR_RTHETA_VEC2, fieldq_polar, fieldu_rect, fieldv_rect, nr, nazm, deltar, $
      clevels, titles, pcolors, scalefac, outfile2, bounds


titles = ["P' (s^-1)",xaxis,yaxis]
;twodplotter, fieldz, x, y, titles, pcolors, $
;             clevels, bounds, outfile


ENDFOR

END


