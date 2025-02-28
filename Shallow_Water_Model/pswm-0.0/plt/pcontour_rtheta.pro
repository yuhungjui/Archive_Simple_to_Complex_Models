PRO PCONTOUR_RTHETA, field, numradial, numtheta, deltar, $
      contour_num, maxfld, minfld, plot_title

; Plots contoured R-THETA plots of a 2D(number r points, number theta points)
; field.
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

flag_lt = WHERE(cont_lvls lt 0)
flag_gt = WHERE(cont_lvls gt 0)

IF (flag_lt[0] NE -1) THEN cont_sty[flag_lt] = 2
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
    YSTYLE = 1, CHARSIZE = 1.5
; POSITION = [0.2,0.33,0.8,0.8], /FILL

;POLAR_CONTOUR, fieldplot, thetacoord, rcoord, /overplot,$
;    c_colors=255*(1+findgen(contour_num)-findgen(contour_num)), $
;    LEVELS = cont_lvls, /follow


RETURN

END
