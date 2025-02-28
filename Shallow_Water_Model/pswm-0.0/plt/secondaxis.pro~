PRO secondaxis, DATA1, DATA2, XPLOT, TITLES, BOUNDS1, BOUNDS2, PLOTNAME

;************************************************************
; Creates a 1D plot.  DATA is a 2D array, XPLOT are ordinate
; values.  TITLES is a three element array (xtitle, ytitle, 
; and title).  XBOUNDS is the bounds of the ordinate, PLOTNAME
; if the filename
;************************************************************

!P.CHARTHICK = 5.0
!P.CHARSIZE = 1.3
;!P.MULTI = [0, 1, 1]

set_plot, 'ps'

device, FILENAME = PLOTNAME, /COLOR, /INCHES, XSIZE = 8.0, YSIZE = 6.0, $
  XOFFSET = 0.1, YOFFSET = 2.5

plot, XPLOT, DATA1, TITLE = TITLES(0), XTITLE = TITLES(1), YTITLE = TITLES(2), $
  XRANGE = [BOUNDS1(0), BOUNDS1(1)], YRANGE = [BOUNDS1(2), BOUNDS1(3)], $ 
  XSTYLE = 1, YSTYLE=1, CHARSIZE = 1.5, CHARTHICK = 3.0, THICK=5.0, $
  POSITION=[0.2,0.2,0.9,0.9], LINESTYLE=0
PLOT, XPLOT, DATA2, /noerase, xrange=[BOUNDS2(0), BOUNDS2(1)], ystyle=4, xstyle=1, $
  POSITION=[0.2,0.2,0.9,0.9], CHARSIZE = 1.5, CHARTHICK = 2.0, THICK=5.0, $
  LINESTYLE=1
AXIS, YAXIS=1, YRANGE=[BOUNDS2(2), BOUNDS2(3)], ytitle='Tangential Wind (m s!E-1!N)'

; TEMP FOR BUDGET
XYOUTS, 0.33, 0.86, 'Relative Vorticity', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
XYOUTS, 0.33, 0.82, 'Tangential Wind', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
PLOTS, [0.23, 0.32], [0.87, 0.87], LINESTYLE=0, THICK=5.0, /NORMAL
PLOTS, [0.23, 0.32], [0.83, 0.83], LINESTYLE=1, THICK=5.0, /NORMAL

device, /close

!P.CHARTHICK = 1.0
!P.CHARSIZE = 1.0
!P.MULTI = [0, 1, 1]

set_plot, 'x'

;spawn, 'gv ' + PLOTNAME +  ' &'

end
