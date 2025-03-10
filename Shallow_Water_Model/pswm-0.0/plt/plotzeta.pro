
FUNCTION shape_fun, s
shape_fun = (1-s)^2.0*(1+2*s)
RETURN, shape_fun
END 

;
; A PROGRAM TO PLOT DIFFERENT HEATING PROFILES
;
nvar=3
nr=1000

rad=FINDGEN(nr)*0.1

zeta=FLTARR(nr)
wind=FLTARR(nr)
both=FLTARR(nr,2)

r1=FLTARR(nvar)
r2=FLTARR(nvar)
r3=58
r4=62
r1=13
r2=17
zeta1 = 0.0007  
zeta2 = 0.0007

FOR i = 0, nr-1 DO BEGIN
   IF rad(i) LT r1 THEN BEGIN
      zeta(i) = zeta1
   ENDIF ELSE IF rad(i) LE r2 THEN BEGIN
      zeta(i) = zeta1*shape_fun((rad(i)-r1)/(r2-r1))+ $
               zeta2*shape_fun((r2-rad(i))/(r2-r1)) 
   ENDIF ELSE IF rad(i) LE r3 THEN BEGIN
      zeta(i) = zeta2
   ENDIF ELSE IF rad(i) LE r4 THEN BEGIN
      zeta(i) = zeta2*shape_fun((rad(i)-r3)/(r4-r3))
   ENDIF ELSE BEGIN
      zeta(i) = 0.0
   ENDELSE
   IF i GE 2 THEN BEGIN
      integral1d, (1000.0*rad)*zeta, 0, i, 0.1, tmp1
      wind(i) = (1.0/rad(i))*tmp1
   ENDIF
ENDFOR

outfile1 = './zetainit.ps'
xaxis = "Radius (km)"
yaxis = "Relative Vorticity (s!E-1!N)"
titles = ['Initial Vortex',xaxis,yaxis]
bounds1 = [0,80,0,0.001] 
bounds2 = [0,80,0,30] 
secondaxis, zeta, wind, rad, TITLES, BOUNDS1, BOUNDS2, outfile1

END
