PRO onedoverplotter3, DATA, XPLOT, TITLES, BOUNDS, PLOTNAME

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
DATA(*,0)=0.0
DATA(*,2)=0.0

device, FILENAME = PLOTNAME, /COLOR, /INCHES, XSIZE = 8.0, YSIZE = 6.0, $
  XOFFSET = 0.1, YOFFSET = 2.5

plot, XPLOT, DATA(*,0), TITLE = TITLES(0), XTITLE = TITLES(1), YTITLE = TITLES(2), $
  XRANGE = [BOUNDS(0), BOUNDS(1)], YRANGE = [BOUNDS(2), BOUNDS(3)], $ 
  XSTYLE = 1, YSTYLE=1, CHARSIZE = 1.5, CHARTHICK = 3.0, THICK=5.0, $
  POSITION=[0.2,0.2,0.9,0.9]
dsize = size(data)
;for i = 1, dsize(2)-1 do begin
 ;  oplot, XPLOT, DATA(*,i), LINESTYLE=2
;    oplot, XPLOT, DATA(*,i), THICK=0.5 + float(i)
;endfor
;
; Comment in for VTBUDGET
;
;stop
;for i = 1, dsize(2)-1 do begin
;   oplot, XPLOT, DATA(*,i), THICK=5.0, LINESTYLE=i
;endfor
oplot, XPLOT, DATA(*,1), THICK=5.0, LINESTYLE=1
oplot, XPLOT, DATA(*,2), THICK=5.0, LINESTYLE=3

;XYOUTS, TIMELEFT, TIMEBOTT, TIMESTRING, /NORMAL, CHARSIZE = 1.5, CHARTHICK = 7.0

; TEMP FOR BUDGET
;XYOUTS, 0.73, 0.86, 'THICK', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
XYOUTS, 0.73, 0.82, 'MIDDLE', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
;XYOUTS, 0.73, 0.78, 'MEAN+EDDY', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
;XYOUTS, 0.73, 0.74, 'ACTUAL', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
;XYOUTS, 0.73, 0.78, 'THIN', /NORMAL, CHARSIZE = 1.5, CHARTHICK = 5.0
;PLOTS, [0.63, 0.72], [0.87, 0.87], LINESTYLE=0, THICK=5.0, /NORMAL
PLOTS, [0.63, 0.72], [0.83, 0.83], LINESTYLE=1, THICK=5.0, /NORMAL
;PLOTS, [0.63, 0.72], [0.79, 0.79], LINESTYLE=2, THICK=5.0, /NORMAL
;PLOTS, [0.63, 0.72], [0.75, 0.75], LINESTYLE=3, THICK=5.0, /NORMAL
;PLOTS, [0.63, 0.72], [0.79, 0.79], LINESTYLE=3, THICK=5.0, /NORMAL

device, /close

!P.CHARTHICK = 1.0
!P.CHARSIZE = 1.0
!P.MULTI = [0, 1, 1]

set_plot, 'x'

;spawn, 'gv ' + PLOTNAME +  ' &'

end
