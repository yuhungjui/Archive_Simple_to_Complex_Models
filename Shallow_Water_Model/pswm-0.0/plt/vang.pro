;
; VANG.PRO
; Program to plot angular momentum surfaces in r-z plane
; over the tangential velocity
;

PRINT, "Running VANG.PRO..."

dir = ''
var = ''
temp = ''
outfile = ''
infilene = ''
infile2 = ''
compstring = ''

lin=0
nazm = 0.0
pltyp = 0.0
slice = 0.0
nmult = 0

;dir = "/ldata/users/modelruns/isen/h00t15_restart/"
dir = "/home/eric/scratch3/modelruns/lappem-isen-1.1/h00t15ac/"
;dir= "/home/eric/scratch3/modelruns/lappem-isen-1.2/h00t15ac_rs/"
;infilene = "h00t15r"
infilene = "h00t15ac"

;READ, dir, PROMPT="Enter output directory path: "
dir = STRCOMPRESS(STRING(dir), /REMOVE_ALL)
;READ, infilene, PROMPT="Enter infile (no extension): "
READ, lin, PROMPT="Enter beginning time level: "
lin=FIX(lin)
infilene = STRCOMPRESS(STRING(infilene+'_'), /REMOVE_ALL)
READ, var, PROMPT="Enter variable to plot (m,u,v,z,d,q,p,e,pp): "
READ, azorslice, PROMPT="Azimuthal mean or slice? (0-azmean, 1-slice): "
if (azorslice LT 0.5) then begin
   READ, nazm, PROMPT="Enter the number of azimuthal points: " 
   if var EQ 'u' then begin 
      nmult = 1
      compstring = STRCOMPRESS('v', /REMOVE_ALL)
   endif
   if var EQ 'v' then begin 
      nmult = 1
      compstring = STRCOMPRESS('u', /REMOVE_ALL)
   endif        
endif else begin
   READ, pltyp, PROMPT="Enter the slice (0-xy, 1-xz, 2-yz): "
   READ, slice, PROMPT="Enter the index (0-max): "
   pltyp = FIX(pltyp)
endelse
varstring = STRCOMPRESS(var, /REMOVE_ALL)
READ, anomortotal, PROMPT="Total field (0), Anomaly (1): "

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

;
; DEFINE THE CONTOUR LEVELS AND TITLES
; THESE ARE INITIALLY SET FOR AZIMUTHAL MEAN PLOTS
;

if (var EQ 'p') then begin
   titles  = ["M(r,!7h!3) (J kg!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 1000.0*[288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320]
   p_colors = c2w_plot_colors
endif else if (var EQ 'pp') then begin
   titles  = ["M'(r,!7h!3) (m s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = [-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
   p_colors = c2w_plot_colors
endif else if (var EQ 'm') then begin
   titles  = ["m(r,!7h!3) (Pa K!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = [10,20,30,50,100,150,200,300,400,500,1000,1500,2000,2500,3000,4000,5000]
   p_colors = rbw_plot_colors
endif else if (var EQ 'u') then begin
   titles  = ["u(r,!7h!3) (m s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = [-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   p_colors = c2w_plot_colors
endif else if (var EQ 'v') then begin
;   titles  = ["v(r,!7h!3) (m s!U-1!N)", "r (km)", "!7h!3 (K)"]
   titles  = ["!S!A-!R!N!8v!N!3 (m s!U-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = [-20,-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15,20]
   clevels = [-10,0,1,2,4,6,8,10,15,20,25,30,35,40,45,50,55]
   p_colors = rbw_plot_colors
;   p_colors = bwl_plot_colors
endif else if (var EQ 'z') then begin
   titles  = ["!7f!3 (x 10!E-3!N s!E-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 0.1*[-2,0,1,3,5,10,15,20,25,30,40,50,60,70,80,90,100]   
   p_colors = rbw_plot_colors
endif else if (var EQ 'd') then begin
   titles  = ["!7d!3 (x 10!E-3!N s!E-1!N)", "r (km)", "!7h!3 (K)"]
   clevels = 0.1*[-100,-50,-30,-20,-10,-5,-3,-1,0,1,3,5,10,20,30,50,100]
   p_colors = c2w_plot_colors
endif else if (var EQ 'q') then begin
;   titles  = ["PV(r,!7h!3) (PVU)", "r (km)", "!7h!3 (K)"]
    titles  = ["PV (PVU)", "r (km)", "!7h!3 (K)"]
   clevels = 0.4*[-20,0,2.5,5,7.5,10,15,20,25,30,40,50,60,70,80,90,100]
   p_colors = rbw_plot_colors
endif

omega = 7.292d-05
nfiles = 1000
pi = 3.14159d0
g = 9.81d0
ostamp = ''

FOR l = lin, nfiles DO BEGIN

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
; Open the output file
; 
if (var EQ 'pp') then varstring = 'p'

infile = dir + infilene + varstring + '.' + tstamp 
OPENR, 1, infile, ERROR = err
if (azorslice LT 0.5 AND (var EQ 'u' OR var EQ 'v')) then begin
   infile2 = dir + infilene + compstring + '.' + tstamp 
   OPENR, 2, infile2, ERROR = err 
endif
if (err NE 0) then begin 
   PRINT, "** lappem_isen.pro: file not found, exiting..."  
   break
endif
PRINT, "** lappem_isen.pro: opening and reading " + infile + "..."

FOR i = 0, nmult DO BEGIN

READF, i+1, temp
READF, i+1, temp
READF, i+1, nx, ny, nl
nl = FIX(nl)
thl = DBLARR(nl)
READF, i+1, xmin, xmax, ymin, ymax
FOR k = 0, nl-1 DO BEGIN
   READF, i+1, dum1
   PRINT, dum1
   thl(k) = dum1
ENDFOR
READF, i+1, fcor, beta, cval, timesec

ENDFOR

nx = FIX(nx)
ny = FIX(ny)
field = DBLARR(nx+1,ny+1,nl)
field2 = DBLARR(nx+1,ny+1,nl)
y = DBLARR(ny+1)
x = DBLARR(nx+1)
xc = fix(nx / 2.0) 
yc = fix(ny / 2.0) 

timehour = timesec/3600.

deltax = (xmax-xmin)/nx
deltay = (ymax-ymin)/ny

FOR k = 0, nl-1 DO BEGIN
FOR j = 0, ny DO BEGIN
y(j) = j*deltay - 0.5d*(ymax-ymin)
FOR i = 0, nx DO BEGIN
x(i) = i*deltax - 0.5d*(xmax-xmin)
	READF, 1, dum1
	field(i,j,k) = dum1
        if (azorslice LT 0.5 AND (var EQ 'u' OR var EQ 'v')) then begin
           READF, 2, dum1
           field2(i,j,k) = dum1
        endif     
ENDFOR
ENDFOR
ENDFOR
CLOSE,1,2,3

if (anomortotal EQ 1.0 AND l EQ 0) then begin
    field_timezero = DBLARR(nx+1,ny+1,nl)
    field_timezero = field
    field_mean = DBLARR(nl)   ; set to zero
endif
if (anomortotal LT 0.5 AND l EQ 0) then begin
    field_timezero = DBLARR(nx+1,ny+1,nl)
    field_mean = DBLARR(nl)   ; set to zero
endif
if (var EQ 'pp' AND l EQ 0) then begin
    field_tzero = DBLARR(nx+1,ny+1,nl)
    field_tzero = field
endif


if (var EQ 'pp') then begin
     tmpfield = DBLARR(nx+1,ny+1,nl)
     tmpfield = field
     FOR i = 0, nl-1 DO BEGIN
        field(*,*,i) = tmpfield(*,*,i) - field_tzero(0,0,i)
     ENDFOR 
endif


;
; MULTIPLY BY 1000 FOR VORTICITY AND DIVERGENCE
;

if (var EQ 'z' OR var EQ 'd') then begin
   field = field * 1000.0
endif

;
; CONVERT TO PVU FOR POTENTIAL VORTICITY
;

if (var EQ 'q') then begin
   field = field * 1.0e+06
endif


x = x / 1000.0d0
y = y / 1000.0d0

CLOSE, 1, 2, 3

; tempfix
; timehour = timehour + 17.0

;
; COMPUTE THE AZIMUTHAL MEAN
;

if (azorslice LT 0.5) then begin

   nr = FIX(nx/2.0)
   r = FINDGEN(nr) * deltax
   r = r / 1000.0d0
   field_azmean = DBLARR(nr,nl)
   field_prime = DBLARR(nr,nazm,nl)
   field_polar = DBLARR(nr,nazm,nl)
   vr = DBLARR(nx+1,ny+1,nl)
   vt = DBLARR(nx+1,ny+1,nl)

   if (var EQ 'u') then begin
   CARTESIANTOCYLIND, field, field2, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc, $
		       yc
		     
   field = vr
   endif
   if (var EQ 'v') then begin
   CARTESIANTOCYLIND, field2, field, $
                       vr, vt, $
                       deltax, deltay, nazm, $
		       xc, $
		       yc
		     
   field = vt
   endif

   INTERPTOCYLIND, field, field_polar, $
                nr, nazm, nl, $
                xc, $
                yc, $
                xPositions, $
                yPositions

   field_azmean = AZIMUTHALMEAN(field_polar, $
                        perturbation = field_prime)

endif

PRINT, "** lappem_isen.pro: plotting " + infile + "..."

d = 0
plotname = outfile

if (azorslice LT 0.5) then begin
  
;  bounds = [min(r),max(r),min(thl),max(thl),timehour]
  bounds = [0,100.0,min(thl),max(thl),timehour]
;  bounds = [0,100.0,300,360,timehour]


; ANGULAR MOMENTUM ONLY IF FIELD IS V
  field2_azmean = DBLARR(nr,nl)
  FOR j = 0, nr-1 DO BEGIN
  FOR i = 0, nl-1 DO BEGIN
     field2_azmean(j,i) = ( 1000.0 * r(j) * field_azmean(j,i) + $
                          0.5 * fcor * (1.0e+6) * r(j)^2 ) / (1.0e+6)

  ENDFOR
  ENDFOR
;  stop

  twodplotter2fld, field_azmean, field2_azmean, r, thl, titles, p_colors, $
             clevels, bounds, plotname

;  twodstndplot, field_azmean, r, thl, titles, $
;             clevels, bounds, plotname

endif else begin

  if (pltyp eq 0) then begin

    qref = field(*,*,slice) - field_timezero(*,*,slice) 
    qred = REFORM(qref)
    print, 'Max/Min: ', max(qred), min(qred)

;    bounds = [min(x),max(x),min(y),max(y),timehour]
    bounds = [-100.0,100.0,-100.0,100.0,timehour]
;    clevels = min(qref) + FINDGEN(nlev) * (max(qref) - min(qref))/nlev
    titles = ["PV (PVU) on !7h!3=333 K", "x (km)", "y (km)"]

    twodplotter, qred, x, y, titles, p_colors, $
             clevels, bounds, plotname

  endif

  if (pltyp eq 1) then begin

    qref = field(*,slice,*) - field_timezero(*,slice,*) 
    qred = REFORM(qref)
    print, 'Max/Min: ', max(qred), min(qred)

;    bounds = [min(x),max(x),min(thl),max(thl),timehour]
    bounds = [-300.0,300.0,min(thl),max(thl),timehour]
    clevels = min(qref) + FINDGEN(nlev) * (max(qref) - min(qref))/nlev
    titles = [var, "x (km)", "Theta (K)"]
;    clevels = -4.25 + FINDGEN(nlev)*0.5 
;
    twodplotter, qred, x, thl, titles, p_colors, $
             clevels, bounds, plotname

  endif

  if (pltyp eq 2) then begin

    qref = field(slice,*,*) - field_timezero(slice,*,*) 
    qred = REFORM(qref)
    print, 'Max/Min: ', max(qred), min(qred)

 ;   bounds = [min(y),max(y),min(thl),max(thl),timehour]
    bounds = [-300.0,300.0,min(thl),max(thl),timehour]
    clevels = min(qref) + FINDGEN(nlev) * (max(qref) - min(qref))/nlev
    titles = [var, "y (km)", "Theta (K)"]

    twodplotter, qred, y, thl, titles, p_colors, $
             clevels, bounds, plotname

  endif

endelse

stop
ENDFOR


PRINT, "** lappem_isen.pro: DONE."

END
