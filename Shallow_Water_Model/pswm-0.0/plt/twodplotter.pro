PRO twodplotter, DATA, XPLOT, YPLOT, TITLES, COLORSCHEME, CONTOURLEVELS, BOUNDS, PLOTNAME
;######################################################################
; Plots the DATA array in a PS file.  The output is put into
; PLOTNAME.  XPLOT and YPLOT are the x-axis and y-axis values,
; respectively.  TITLES is a three-element array of strings for the
; plot title, the x-axis title, and the y-axis title, respectively.
; COLORSCHEME is the color scheme (already put through c32) for the
; contour-filled plot.  BOUNDS is a four element array with the XRANGE
; and YRANGE values in that order, as well as the time (in hours).  
; CONTOURLEVELS is the array of the contour levels for the plot 
; (should be of same dimensionality as COLORSCHEME).
;------------------------------------------------------------
; Made especially for the RAMS runs that I am examining.  This
; function may be of little use to anyone else.
;######################################################################


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

!P.CHARTHICK = 5.0
!P.CHARSIZE = 1.3
!P.MULTI = [0, 1, 1]

set_plot, 'ps'

device, FILENAME = PLOTNAME, /COLOR, /INCHES, XSIZE = 8.0, YSIZE = 6.0, $
  XOFFSET = 0.1, YOFFSET = 2.5

contour, DATA, XPLOT, YPLOT, LEVELS = CONTOURLEVELS, C_COLORS = COLORSCHEME, $
  /FILL, /FOLLOW, /ISOTROPIC, TITLE = TITLES(0), XTITLE = TITLES(1), YTITLE = TITLES(2), $
  XRANGE = [BOUNDS(0), BOUNDS(1)], YRANGE = [BOUNDS(2), BOUNDS(3)], XSTYLE = 1, YSTYLE = 1

if ((CONTOURLEVELS(0) / CONTOURLEVELS(NLEV-1)) LT (-0.5)) then begin
;    contour, DATA, XPLOT, YPLOT, LEVELS = [0.0], THICK = 3.0, C_COLORS = 0, /FOLLOW, /OVERPLOT
endif

XYOUTS, TIMELEFT, TIMEBOTT, TIMESTRING, /NORMAL, CHARSIZE = 1.5, CHARTHICK = 7.0, COLOR = black_cons(0)

for i = 0, NLEV - 1 do begin

    CURRBOT = BBOT + i * BHEI
    polyfill, [BLEF, BLEF, BRIG, BRIG], [CURRBOT, CURRBOT + BHEI, CURRBOT + BHEI, CURRBOT], $
      /NORMAL, COLOR = COLORSCHEME(i)
    PLOTS, [BLEF, BLEF], [CURRBOT, CURRBOT + BHEI], /NORMAL, COLOR = black_cons(i), THICK = 3.0
    PLOTS, [BRIG, BRIG], [CURRBOT, CURRBOT + BHEI], /NORMAL, COLOR = black_cons(i), THICK = 3.0
    PLOTS, [BLEF, BRIG], [CURRBOT, CURRBOT], /NORMAL, COLOR = black_cons(i), THICK = 3.0
    PLOTS, [BLEF, BRIG], [CURRBOT + BHEI, CURRBOT + BHEI], /NORMAL, COLOR = black_cons(i), THICK = 3.0

    CONSTRING = STRMID(STRCOMPRESS(STRING(CONTOURLEVELS(i)), /remove_all), 0, 5)
;   CONSTRING = STRCOMPRESS(STRING(CONTOURLEVELS(i)), /remove_all)
    IF i NE 0 THEN BEGIN
    XYOUTS, BRIG + TEXTOFFR, CURRBOT + TEXTOFFB, CONSTRING, /NORMAL, $
      CHARSIZE = 0.8, COLOR = black_cons(i)
    ENDIF
endfor


device, /close

!P.CHARTHICK = 1.0
!P.CHARSIZE = 1.0
!P.MULTI = [0, 1, 1]

set_plot, 'x'

;spawn, 'gv ' + PLOTNAME +  ' &'

end
