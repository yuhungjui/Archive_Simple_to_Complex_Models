; IDL MEANWINDS2.PRO
; Hovmoller plot of azimuthal mean wind and pressure
; Written: 08/11/08
; Eric A. Hendricks

CLOSE,1,2,3
dir = '/home/eric/scratch3/modelruns/pswm-0.0/e7o20hres/'
infilene = 'e7o20'
var = ''
temp = ''
outfile = STRARR(5)
outfile1 = STRARR(5)
infile2 = ''
varstring = ''
xmin = 0.0d0
xmax = 0.0d0
fcor = 0.0d0
beta = 0.0d0
timesec = 0.0d0
ymin = 0.0d0
ymax = 0.0d0
nmult = 0
dum1 = 0.0d0
dum2 = 0.0d0
dum3 = 0.0d0

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
READ, var, PROMPT="Enter variable to plot (p,u,v,z,d,q,h): "
READ, lone, PROMPT="Enter beginning file number: "
READ, ltwo, PROMPT="Enter ending file number: "
if var EQ 'u' then begin 
    nmult = 1
    compstring = STRCOMPRESS('v', /REMOVE_ALL)
endif
if var EQ 'v' then begin 
    nmult = 1
    compstring = STRCOMPRESS('u', /REMOVE_ALL)
endif        
lone = FIX(lone)
ltwo = FIX(ltwo)
varstring = STRCOMPRESS(var, /REMOVE_ALL)

nxd = 1024
nyd = 1024
cval = 205.0    ; m/s
grav = 9.81
deltat = 180.0 ; s
nrmax = 512
nt = ltwo - lone + 1

;
; DEFINE THE CONTOUR LEVELS
;

fac = 1000.0
clevels_per_h = 0.1*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_h  = "h' (m)"
clevels_per_z = 0.0001*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_z  = "!7f!3' (x 10!E-3!N s!E-1!N)"
clevels_per_d = 0.00001*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_d  = "!7d!3' (x 10!E-4!N s!E-1!N)"
clevels_per_u = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
plot_title_u  = "u' (m s!E-1!N))"
clevels_per_v = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
clevels_per_v = [0,5,10,20,30,40,50,55,60,65,70,75,80,85,90,95,100]
plot_title_v  = "v' (m s!E-1!N))"

;
; COLOR TABLES
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

plotcolor = 0.0
if (plotcolor LT 0.5) then begin
   pcolors = rbw_plot_colors
   ;pcolors = c2w_plot_colors
endif else begin
   pcolors = bwl_plot_colors
endelse


;
; AZIMUTHAL MEAN ARRAYS
;

count = 0
nrd = FIX(FLOAT(nxd) / 2.0)
bs = STRCOMPRESS(STRING(lone), /REMOVE_ALL)
es = STRCOMPRESS(STRING(ltwo),   /REMOVE_ALL)
outfile = dir + infilene + 'azmeanhov.' + bs + '.' + es + '.ps'
lstring1 = STRCOMPRESS(STRING(lone), /REMOVE_ALL)
lstring2 = STRCOMPRESS(STRING(ltwo), /REMOVE_ALL)

;
; OPEN AND READ THE OUTPUT FILE
; 

fieldone = DBLARR(nxd+1,nyd+1)
fieldonetmp = DBLARR(nxd+1,nyd+1)
field    = DBLARR(nxd+1,nyd+1)
field2   = DBLARR(nxd+1,nyd+1)
fieldz   = DBLARR(nxd+1,nyd+1)   ; for centering
time_array = DBLARR(nt)

fieldhov = DBLARR(nrmax,nt)

FOR m = lone,ltwo DO BEGIN

lstring = STRCOMPRESS(STRING(m), /REMOVE_ALL)

if (m lt 100) then begin
  if (m lt 10) then begin
     tstamp = 't0' + lstring
     ostamp = '00' + lstring
  endif else begin
     tstamp = 't' + lstring
     ostamp = '0' + lstring
  endelse
endif else begin
  tstamp = lstring
endelse

infile = dir + infilene + varstring + '.' + tstamp 
PRINT, "** meanwinds.pro: opening and reading " + infile + "..."
OPENR, 1, infile, ERROR = err 
READF, 1, temp
READF, 1, temp
READF, 1, nx, ny
READF, 1, xmin, xmax, ymin, ymax
READF, 1, fcor, beta, cval, timesec
if (var EQ 'u' OR var EQ 'v') then begin
   infile2 = dir + infilene + compstring + '.' + tstamp 
   OPENR, 2, infile2, ERROR = err 
   READF, 2, temp
   READF, 2, temp
   READF, 2, nx, ny
   READF, 2, xmin, xmax, ymin, ymax
   READF, 2, fcor, beta, cval, timesec
endif
infilez = dir + infilene + 'z' + '.' + tstamp 
OPENR, 3, infilez, ERROR = err 
READF, 3, temp
READF, 3, temp
READF, 3, nx, ny
READF, 3, xmin, xmax, ymin, ymax
READF, 3, fcor, beta, cval, timesec
if (err NE 0) then begin 
   PRINT, "** meanwinds.pro: file not found, exiting..."  
   break
endif
time_array(m) = timesec/3600.0d0

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
        if (var EQ 'u' OR var EQ 'v') then begin
           READF, 2, dum2
           field2(i,j) = dum2
        endif  
        READF, 3, dum3
        fieldz(i,j) = dum3
ENDFOR
ENDFOR
CLOSE, 1, 2, 3

fieldone = field
fieldonetmp = field2

;
; DETERMINE CENTER FROM VORTICITY
;

PRINT, "** meanwinds.pro: Determining center..."

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
xc1 = xCenter_z
yc1 = yCenter_z
PRINT, "** meanwinds.pro: xCenter, yCenter indeces = ", xCenter_z, yCenter_z

;
; AZIMUTHAL MEAN, PERTURBATION QUANTITIES ABOUT CENTER
;

PRINT, "** meanwinds.pro: Determining azimuthal mean diagnostics..."

nazm = 50
nr = FIX(nx/2.0)
deltar = deltax
r = FINDGEN(nr) * deltar

fieldone_azmean = DBLARR(nr)
fieldone_prime = DBLARR(nr,nazm)
fieldone_polar = DBLARR(nr,nazm)

vr = DBLARR(nx+1,ny+1)
vt = DBLARR(nx+1,ny+1)

if (var EQ 'u') then begin
CARTESIANTOCYLIND2D, fieldone, fieldonetmp, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc1, $
		       yc1

fieldone = vr
endif
if (var EQ 'v') then begin
CARTESIANTOCYLIND2D, fieldonetmp, fieldone, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc1, $
		       yc1

fieldone = vt	
endif

INTERPTOCYLIND2D, fieldone, fieldone_polar, $
                nr, nazm, $
                xc1, $
                yc1, $
                xPositions, $
                yPositions

fieldone_azmean = AZIMUTHALMEAN2D(fieldone_polar, $
                        perturbation = fieldone_prime)


fieldhov(*,m) = fieldone_azmean

ENDFOR

;
;  PLOTTING
;

bounds = [0,300,double(lone)*3.0/60.0,(double(ltwo)-1)*3.0/60.0, -999]
r = r / 1000.0d0
titles = ['v (m s!E-1!N)','r (km)', 'time (hrs.)']
hovplotter, fieldhov, r, time_array, $
             titles, pcolors, clevels_per_v, bounds, outfile
r = r * 1000.0d0



END


