PRO CONTOURXY, field, xcoord, ycoord, $
                contour_num, hour_run, plot_title, $
                xaxistitle, yaxistitle, plotname, bounds


; Plots contoured X-Y plots of a 2D(number x points, number y points)
; field.

cont_int = (MAX(field) - MIN(field)) / contour_num
IF (cont_int EQ 0.0) THEN cont_int = 1.0
cont_lvls = FLTARR(contour_num)
cont_lvls = (MIN(field) + FINDGEN(contour_num) * cont_int)
cont_lvls = [-0.5,-0.1,-0.05,-0.01,0,0.01,0.05,0.1,0.5]
cont_lbl = FLTARR(contour_num)
cont_lbl = cont_lbl + 1
cont_sty = FLTARR(contour_num)
flag_lt = WHERE(cont_lvls lt 0)
flag_gt = WHERE(cont_lvls gt 0)

IF (flag_lt[0] NE -1) THEN cont_sty[flag_lt] = 1
IF (flag_gt[0] NE -1) THEN cont_sty[flag_gt] = 0

!P.MULTI = [0, 1, 1]

set_plot, 'ps'
device, FILENAME = PLOTNAME, /COLOR, /INCHES, XSIZE = 6.5, YSIZE = 6.0, $
  XOFFSET = 0.1, YOFFSET = 2.5

CONTOUR, field, xcoord, ycoord, $
    TITLE = plot_title + ', t = ' + hour_run + ' h.', $
    XTITLE = xaxistitle, $
    YTITLE = yaxistitle, $
    LEVELS = cont_lvls, $
    C_LABELS = cont_lbl, $
    C_LINESTYLE = cont_sty, $
    C_CHARSIZE = 0.9, $
    XSTYLE = 1, $
    YSTYLE = 1, $
    CHARSIZE = 1.5, $
    XRANGE = [BOUNDS(0), BOUNDS(1)], $
    YRANGE = [BOUNDS(2), BOUNDS(3)] 

device, /close

RETURN

END
