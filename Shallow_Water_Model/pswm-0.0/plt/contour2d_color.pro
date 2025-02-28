PRO CONTOUR2D_COLOR, field, xcoord, ycoord, $
                hour_run, plot_title, $ 
                xaxistitle, yaxistitle, $
                red, green, blue, cont_levels 

; Makes contour plots of a 2D field
      
temp = SIZE(red)
contour_num = temp(1)
cont_int = (MAX(field) - MIN(field)) / contour_num
IF (cont_int EQ 0.0) THEN cont_int = 1.0
cont_lvls = FLTARR(contour_num)
cont_lvls = (MIN(field) + FINDGEN(contour_num) * cont_int)
cont_lbl = FLTARR(contour_num)
cont_lbl = cont_lbl + 1
cont_sty = FLTARR(contour_num)
cont_lvls = cont_levels

flag_lt = WHERE(cont_lvls lt 0)
flag_gt = WHERE(cont_lvls gt 0)

IF (flag_lt[0] NE -1) THEN cont_sty[flag_lt] = 1
IF (flag_gt[0] NE -1) THEN cont_sty[flag_gt] = 0

TVLCT, red, green, blue

CONTOUR, field, xcoord, ycoord, $
    TITLE = plot_title + ', t = ' + hour_run + 'h.', $
    XTITLE = xaxistitle, $
    YTITLE = yaxistitle, $
    LEVELS = cont_lvls, $
    C_LABELS = cont_lbl, $
    C_LINESTYLE = cont_sty, $
    C_CHARSIZE = 0.7, $
    XSTYLE = 1, $
    YSTYLE = 1, $
    CHARSIZE = 1.5, /FILL, $
    C_COLORS = [1,2,3,4,5,6,7,8]

RETURN

END
