PRO PCONTOUR_RTHETA_VEC, field, datau, datav, numradial, numtheta, deltar, $
      contour_num, maxfld, minfld, scalefac, plot_title

; Plots contoured R-THETA plots of a 2D(number r points, number theta points)
; field with vector winds overlayed.  Datau and Datav are the cartesian components
; in polar space.  Scalefac is arbitrary for magnitude of vectors in ARROW.
; Inputs: field to be contoured, number radial pts, number theta points,
; delta r, number of contours, multiple file loop run, and generic
; title for all plots in series.  Program will automatically update hour
; using the loop_run info.  Note: loop_run must be the actualy file hour.

deltar = FLOAT(deltar)
contour_num = FIX(contour_num)

rcoord = FLOAT(FINDGEN(numradial) * deltar) / 1000.0
thetacoord = FLOAT(FINDGEN(numtheta) * 2.0 * (!PI) / (numtheta-1))

;max_field=max(field)
;min_field=min(field)
max_field = maxfld
min_field = minfld
cont_incr=(max_field-min_field)/(contour_num*1.0)
cont_lvls=(findgen(contour_num)*cont_incr)+min_field
cont_lbl = FLTARR(contour_num)
cont_lbl = cont_lbl + 1
cont_sty = FLTARR(contour_num)
;cont_lvls = [-70,-50,-30,-15,-5,5,15,30,50,70]  for hamp = 70
cont_lvls = [-12,-8,-4,-2,-1,1,2,4,8,12]

flag_lt = WHERE(cont_lvls lt 0)
flag_gt = WHERE(cont_lvls gt 0)

IF (flag_lt[0] NE -1) THEN cont_sty[flag_lt] = 1
IF (flag_gt[0] NE -1) THEN cont_sty[flag_gt] = 0
maxstr = STRMID(STRCOMPRESS(MAX(field), /REMOVE_ALL),0,4)
minstr = STRMID(STRCOMPRESS(MIN(field), /REMOVE_ALL),0,4)

fieldplot = FLTARR(numtheta, numtheta)
fieldplot = TRANSPOSE(field)

POLAR_CONTOUR, fieldplot, thetacoord, rcoord, $
    TITLE = plot_title, $
    XTITLE = 'x (km)', $
;  $Max = ' +maxstr+ ' Min = ' +minstr, $
    YTITLE = 'y (km)', $
    LEVELS = cont_lvls, $
    C_LINESTYLE = cont_sty, $
    C_CHARSIZE = 1.2, $
    C_LABELS = cont_lbl, $
    XSTYLE = 1, CHARTHICK = 4.0, $
    YSTYLE = 1, CHARSIZE = 1.5, $    
    POSITION = [0.1,0.1,0.9,0.9]

;
; OVERLAY THE VECTORS USING ARROW
;

halfd = 0.5   ; meters
rfac = 0.001333333 ;set for 0.1-0.9 grid
FOR j = 35, numradial-1 DO BEGIN
FOR i = 0, numtheta-1  DO BEGIN

   xa = halfd + rfac * rcoord(j) * COS( thetacoord(i) ) 
   ya = halfd + rfac * rcoord(j) * SIN( thetacoord(i) )
   xb = xa + scalefac * datau(j,i)
   yb = ya + scalefac * datav(j,i) 
   PRINT, xa, ya, xb, yb
 
   if (j MOD 3 EQ 0 AND i MOD 2 EQ 0) THEN BEGIN
      ARROW, xa, ya, xb, yb, /NORMALIZED
   endif

ENDFOR
;stop
ENDFOR
ARROW, 0.75, 0.17, 0.8, 0.17, /NORMALIZED
; XYOUTS, 0.8, 0.15, '1 m s!E-1!N', /NORM for hamp = 70
XYOUTS, 0.75, 0.15, '0.15 m s!E-1!N', /NORM 

;ARROW, 0.5, 0.5, 0.7, 0.7, /NORMALIZED

;POLAR_CONTOUR, fieldplot, thetacoord, rcoord, /overplot,$
;    c_colors=255*(1+findgen(contour_num)-findgen(contour_num)), $
;    LEVELS = cont_lvls, /follow


RETURN

END
