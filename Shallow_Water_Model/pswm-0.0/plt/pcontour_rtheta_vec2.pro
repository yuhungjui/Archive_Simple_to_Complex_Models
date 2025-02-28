PRO PCONTOUR_RTHETA_VEC2, field, datau, datav, numradial, numtheta, deltar, $
      contourlevels, titles, colorscheme, scalefac, plot_title, bounds

;PCONTOUR_RTHETA_VEC2, fieldq_prime, fieldu_prime, fieldv_prime, nr, nazm, deltar, $
;      clevels, titles, pcolors, scalefac, outfile, bounds

; Plots contoured R-THETA plots of a 2D(number r points, number theta points)
; field with vector winds overlayed.  Datau and Datav are the cartesian components
; in polar space.  Scalefac is arbitrary for magnitude of vectors in ARROW.
; Inputs: field to be contoured, number radial pts, number theta points,
; delta r, number of contours, multiple file loop run, and generic
; title for all plots in series.  Program will automatically update hour
; using the loop_run info.  Note: loop_run must be the actualy file hour.


NLEV1 = FIX(SIZE(COLORSCHEME, /DIMENSIONS))
NLEV  = NLEV1(0)
black_cons = c32(intarr(NLEV), intarr(NLEV), intarr(NLEV), 1)

BLEF = 0.78
BRIG = 0.83
BBOT = 0.12
BTOP = 0.94
BHEI = (BTOP - BBOT) / FLOAT(NLEV)

TIMELEFT = 0.61
TIMEBOTT = 0.14

TEXTOFFR = 0.01
TEXTOFFB = -0.005

TIMESTRING = STRMID(STRCOMPRESS(STRING(BOUNDS(4)), /remove_all), 0, 5) + ' hrs.'

;!P.CHARTHICK = 5.0
;!P.CHARSIZE = 1.3
;!P.MULTI = [0, 1, 1]

set_plot, 'ps'
device, FILENAME = plot_title, /COLOR, /INCHES, XSIZE = 6.0, YSIZE = 6.0, $
  XOFFSET = 0.1, YOFFSET = 2.5

deltar = FLOAT(deltar)
tmp = SIZE(contourlevels)
contour_num = tmp(1)

rcoord = FLOAT(FINDGEN(numradial) * deltar) / 1000.0
thetacoord = FLOAT(FINDGEN(numtheta) * 2.0 * (!PI) / (numtheta-1))

deltar = FLOAT(deltar)
contour_num = FIX(contour_num)
cont_sty = FLTARR(contour_num)

flag_lt = WHERE(contourlevels lt 0)
flag_gt = WHERE(contourlevels gt 0)

IF (flag_lt[0] NE -1) THEN cont_sty[flag_lt] = 1
IF (flag_gt[0] NE -1) THEN cont_sty[flag_gt] = 0
maxstr = STRMID(STRCOMPRESS(MAX(field), /REMOVE_ALL),0,4)
minstr = STRMID(STRCOMPRESS(MIN(field), /REMOVE_ALL),0,4)

fieldplot = FLTARR(numtheta, numradial)
fieldplot = TRANSPOSE(field)

POLAR_CONTOUR, fieldplot, thetacoord, rcoord, $
    TITLE = TITLES(0), XTITLE = TITLES(1), YTITLE = TITLES(2), $
    LEVELS = contourlevels, C_COLORS = COLORSCHEME, /FILL, /FOLLOW, $
;    C_LINESTYLE = cont_sty, $
    /ISOTROPIC,$
;    C_CHARSIZE = 0.5, $ 
    XRANGE = [BOUNDS(0), BOUNDS(1)], $
    YRANGE = [BOUNDS(2), BOUNDS(3)], $
    XSTYLE = 1, CHARTHICK = 4.0, $
    YSTYLE = 1, CHARSIZE = 1.0, $
    POSITION = [0.1,0.1,0.9,0.9]


;XYOUTS, TIMELEFT, TIMEBOTT, TIMESTRING, /NORMAL, CHARSIZE = 1.5, CHARTHICK = 7.0, COLOR = black_cons(0)

;for i = 0, NLEV - 1 do begin;

;    CURRBOT = BBOT + i * BHEI
;    polyfill, [BLEF, BLEF, BRIG, BRIG], [CURRBOT, CURRBOT + BHEI, CURRBOT + BHEI, CURRBOT], $
;      /NORMAL, COLOR = COLORSCHEME(i)
;    PLOTS, [BLEF, BLEF], [CURRBOT, CURRBOT + BHEI], /NORMAL, COLOR = black_cons(i), THICK = 3.0
;    PLOTS, [BRIG, BRIG], [CURRBOT, CURRBOT + BHEI], /NORMAL, COLOR = black_cons(i), THICK = 3.0
;    PLOTS, [BLEF, BRIG], [CURRBOT, CURRBOT], /NORMAL, COLOR = black_cons(i), THICK = 3.0
;    PLOTS, [BLEF, BRIG], [CURRBOT + BHEI, CURRBOT + BHEI], /NORMAL, COLOR = black_cons(i), THICK = 3.0
;    CONSTRING = STRMID(STRCOMPRESS(STRING(CONTOURLEVELS(i)), /remove_all), 0, 5)
;   CONSTRING = STRCOMPRESS(STRING(CONTOURLEVELS(i)), /remove_all)
;    XYOUTS, BRIG + TEXTOFFR, CURRBOT + TEXTOFFB, CONSTRING, /NORMAL, $
;      CHARSIZE = 0.8, COLOR = black_cons(i);

;endfor

;
; OVERLAY THE VECTORS USING ARROW
;

halfd = 0.5  
tmparr = WHERE(rcoord GT bounds(1))
maxnr = tmparr(0)
rfac = 0.4/bounds(1)
           
FOR j = 0, maxnr-1 DO BEGIN
FOR i = 0, numtheta-1  DO BEGIN

   xa = halfd + rfac * rcoord(j) * COS( thetacoord(i) ) 
   ya = halfd + rfac * rcoord(j) * SIN( thetacoord(i) )
   xb = xa + scalefac * datau(j,i)
   yb = ya + scalefac * datav(j,i)
    
   PRINT, xa, ya, xb, yb
  
   if (j LT 10) then continue
   if (j MOD 5 EQ 0 AND i MOD 2 EQ 0) THEN BEGIN
;   if (j MOD 3 EQ 0) THEN BEGIN
      ARROW, xa, ya, xb, yb, /NORMALIZED, THICK=5.0, HSIZE=120.0
   endif
   
ENDFOR
;stop
ENDFOR


;ARROW, 0.78, 0.14, 0.81, 0.14, /NORMALIZED, HSIZE=120.0, THICK=3.0
; XYOUTS, 0.8, 0.15, '1 m s!E-1!N', /NORM for hamp = 70
;XYOUTS, 0.78, 0.12, '2 m s!E-1!N', /NORM , CHARTHICK=3.0

; FOR TESTING
;ARROW, 0.5, 0.5, 0.7, 0.7, /NORMALIZED, THICK=3.0
;ARROW, 0.96, 0.5, 0.7, 0.7, /NORMALIZED, THICK=3.0

;POLAR_CONTOUR, fieldplot, thetacoord, rcoord, /overplot,$
;    c_colors=255*(1+findgen(contour_num)-findgen(contour_num)), $
;    LEVELS = cont_lvls, /follow


device, /close

;!P.CHARTHICK = 1.0
;!P.CHARSIZE = 1.0
;!P.MULTI = [0, 1, 1]


RETURN

END
