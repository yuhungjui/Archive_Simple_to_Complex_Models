PRO onedplotter, DATA, XPLOT, TITLES, BOUNDS, PLOTNAME

;************************************************************
; Creates a 1D plot.  DATA is a 1D array, XPLOT are ordinate
; values.  TITLES is a three element array (xtitle, ytitle, 
; and title).  XBOUNDS is the bounds of the ordinate, PLOTNAME
; if the filename. 
;************************************************************

;TIMESTRING = STRMID(STRCOMPRESS(STRING(BOUNDS(4)), /remove_all), 0, 5) + ' hrs.'
;IF BOUNDS(4) LT 0.0 THEN TIMESTRING = 'ALL TIMES'
;TIMELEFT = 0.75
;TIMEBOTT = 0.87

!P.CHARTHICK = 5.0
!P.CHARSIZE = 1.3
!P.MULTI = [0, 1, 1]

set_plot, 'ps'

device, FILENAME = PLOTNAME, /COLOR, /INCHES, XSIZE = 8.0, YSIZE = 6.0, $
  XOFFSET = 0.1, YOFFSET = 2.5

plot, XPLOT, DATA, TITLE = TITLES(0), XTITLE = TITLES(1), YTITLE = TITLES(2), $
  XRANGE = [BOUNDS(0), BOUNDS(1)], YRANGE = [BOUNDS(2), BOUNDS(3)], $ 
  XSTYLE = 1, YSTYLE=1, CHARSIZE = 1.5, CHARTHICK = 3.0, THICK=3.0

;XYOUTS, TIMELEFT, TIMEBOTT, TIMESTRING, /NORMAL, CHARSIZE = 1.5, CHARTHICK = 7.0

device, /close

!P.CHARTHICK = 1.0
!P.CHARSIZE = 1.0
!P.MULTI = [0, 1, 1]

set_plot, 'x'

;spawn, 'gv ' + PLOTNAME +  ' &'

end
