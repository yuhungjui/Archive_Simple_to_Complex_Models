; Plots analytic solution to logistic equation

nt=120
dt=0.5    ;h
time=FINDGEN(nt)*dt
p=FLTARR(nt,3)
q=FLTARR(3)
pm=FLTARR(3)

pm(0)=0.003556    ;max PV
pm(1)=0.00444444    ;max PV
pm(2)=0.007619    ;max PV
po=0.00075       ;init PV
q(0)=7.65482e-05         ; thin heating mass sink inverse seconds
q(1)=4.46531e-05         ; middle
q(2)=3.57225e-05         ; thick

FOR k = 0, 2 DO BEGIN
FOR j = 0, nt-1 DO BEGIN
   p(j,k) = pm(k)*po/(po+(pm(k)-po)*exp(-1.0*q(k)*3600.0*time(j)))
ENDFOR
ENDFOR

outfile1 = './logistic.ps'
xaxis = "Radius (km)"
yaxis = "Logistic P (s!E-1!N)"
titles = ['Logistically-limited PV',xaxis,yaxis]
bounds = [0,60,0,0.01] 
onedoverplotter2, p, time, TITLES, BOUNDS, outfile1

END
   
