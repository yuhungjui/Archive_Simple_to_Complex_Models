PRO CONTOURXY_COLOR, field, xcoord, ycoord, $
                contour_num, hour_run, plot_title, $ 
                xaxistitle, yaxistitle

; Plots contoured X-Y plots of a 2D(number x points, number y points)
; field.
      
cont_int = (MAX(field) - MIN(field)) / contour_num
IF (cont_int EQ 0.0) THEN cont_int = 1.0
cont_lvls = FLTARR(contour_num)
cont_lvls = (MIN(field) + FINDGEN(contour_num) * cont_int); activate for vorticity
cont_lvls = [0.00005,0.0001,0.0005,0.0010,0.0015,0.0020,0.0025,0.0030]
;cont_lvls = [0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.2]
;cont_lvls = [0.0001,0.0005,0.0010,0.0015,0.0020,0.0025,0.0030,0.0050]
;cont_lvls = [2,3,4,5,6,7,8,9]
; activate for tracer
;cont_lvls = [1600,1650,1700,1750,1800,1850,1900,1950]
; activate for ED
cont_lvls = [100,250,500,1000,2500,5000,10000,25000]
cont_lbl = FLTARR(contour_num)
cont_lbl = cont_lbl + 1
cont_sty = FLTARR(contour_num)

flag_lt = WHERE(cont_lvls lt 0)
flag_gt = WHERE(cont_lvls gt 0)

IF (flag_lt[0] NE -1) THEN cont_sty[flag_lt] = 1
IF (flag_gt[0] NE -1) THEN cont_sty[flag_gt] = 0

red = [0, 100,0,  0,  0  ,255,255,255,255]
green = [0, 0,  0,  100,255,255,170,100,0]
blue = [50,100,180,255  ,100,0  ,0  ,0  ,0]
TVLCT, red, green, blue

;cont_lvls = cont_lvls / 1.0d+10

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
