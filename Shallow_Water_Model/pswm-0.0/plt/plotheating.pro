
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

Q=FLTARR(nr,nvar)

r1=FLTARR(nvar)
r2=FLTARR(nvar)
r3=58
r4=62
r1(0)=13
r2(0)=17
r1(1)=28
r2(1)=32
r1(2)=43
r2(2)=47
m2 = 125.0   ; K/day

FOR j = 0, nvar-1 DO BEGIN
heat = (1.0/1000000.0)*(m2 * 28935.2) / $
       ( ((r3+r4)/2.0)^2.0 - ((r1(j)+r2(j))/2.0)^2.0 ) 
heat=heat/30.0
print, heat
FOR i = 0, nr-1 DO BEGIN
   IF rad(i) LT r1(j) THEN BEGIN
      Q(i,j) = 0.0
   ENDIF ELSE IF rad(i) LE r2(j) THEN BEGIN
      Q(i,j) = heat*shape_fun((r2(j)-rad(i))/(r2(j)-r1(j))) 
   ENDIF ELSE IF rad(i) LE r3 THEN BEGIN
      Q(i,j) = heat
   ENDIF ELSE IF rad(i) LE r4 THEN BEGIN
      Q(i,j) = heat*shape_fun((rad(i)-r3)/(r4-r3))
   ENDIF ELSE BEGIN
      Q(i,j) = 0.0
   ENDELSE
ENDFOR
ENDFOR

outfile1 = './masssink.ps'
xaxis = "Radius (km)"
yaxis = "Mass Sink (s!E-1!N)"
titles = ['Mass Sink Profiles',xaxis,yaxis]
bounds = [0,80,0,0.0001] 
onedoverplotter2, Q, rad, TITLES, BOUNDS, outfile1

END
