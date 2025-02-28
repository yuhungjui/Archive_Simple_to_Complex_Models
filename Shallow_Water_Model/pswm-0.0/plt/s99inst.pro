PRINT, 'Running IDL s99inst.pro'

; USER SPECIFIED DELTA AND GAMMA
;ndelta = 17
;ngamma = 11
;delta = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45, $
;         0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85]
;gamma = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

; HIGH DENSITY DELTA AND GAMMA
ndelta = 500
ngamma = 500
delta = FINDGEN(ndelta)*0.002
gamma = FINDGEN(ngamma)*0.002

nfreq = 15

freqp = COMPLEXARR(ndelta,ngamma,nfreq)
freqn = COMPLEXARR(ndelta,ngamma,nfreq)
freq  = COMPLEXARR(ndelta,ngamma,nfreq)

fastwn_arr = INTARR(ndelta,ngamma)

tempf1 = COMPLEX(1.0d0)
tempf2 = COMPLEX(1.0d0)
zetaavg = COMPLEX(2.0d-03)
one = COMPLEX(1.0d0)
two = COMPLEX(2.0d0)
half = COMPLEX(0.5d0)

FOR j = 0, ngamma-1 DO BEGIN
FOR i = 0, ndelta-1 DO BEGIN
FOR m = 0, nfreq-1 DO BEGIN

md = COMPLEX(double(m) + one + two)

parta = COMPLEX(md) 
partb = COMPLEX((md - one)*gamma(j))
cone = COMPLEX(md)
ctwo = COMPLEX((md - one)*gamma(j))
cthr = COMPLEX(two*(one - gamma(j)*delta(i)^two) / $
       (one - delta(i)^two))
cfiv = COMPLEX(two*two*(one - gamma(j)*delta(i)^two) / $
       (one - delta(i)^two))
csix = COMPLEX(gamma(j)) 
csev = COMPLEX((one - gamma(j)*delta(i)^two) / $
       (one - delta(i)^two)) 
ceig = COMPLEX(delta(i)^(two*md))
partc = COMPLEX(((cone - ctwo - cthr)^two + $
        cfiv*(csix - csev)*ceig)^half)

freqp(i,j,m) = half*half*(parta + partb + partc) 
freqn(i,j,m) = half*half*(parta + partb - partc) 

tempf1 = freqp(i,j,m)
tempf2 = freqn(i,j,m)

freq(i,j,m) = MAX(tempf2,tempf1)

ENDFOR
ENDFOR
ENDFOR

!P.MULTI=[0,1,1]
SET_PLOT,'ps'
;SET_PLOT, 'X'

DEVICE,bits=8,filename='s99inst.ps',/portrait,xoffset=1,yoffset=4,xsize=6, $
       ysize=5,/inches, /COLOR

openw, 3, 's99inst.out'
printf, 3, 'Gamma, Delta, fast WN, Max Dimless Growth Rate'
FOR j = 0, ngamma-1 DO BEGIN
FOR i = 0, ndelta-1 DO BEGIN

maxfreq = MAX(IMAGINARY(freqp(i,j,*)), fastwn)
fastwn = fastwn + 3
fastwn_arr(i,j) = fastwn
if (maxfreq LT 1.0d-06) then fastwn = 0

printf, 3, gamma(j), delta(i), fastwn

ENDFOR
ENDFOR
CLOSE, 3

FOR m = 0, nfreq-1 DO BEGIN
freqtmp = freqp(*,*,m)
freq2d = REFORM(freqtmp)
freq2dimag = IMAGINARY(freq2d)

hour_run = '00'
wnstr = STRCOMPRESS(STRING(m+3),/REMOVE_ALL)
plot_title = 'Wavenumber ' + wnstr
contour_num = 100
xtitle = 'Delta'
ytitle = 'Gamma'
CONTOURXY, freq2dimag, delta, gamma, $
           contour_num, hour_run, plot_title, $
           xtitle, ytitle

ENDFOR
plot_title = 'Most Unstable WN'
CONTOURXY_COLOR, fastwn_arr, delta, gamma, $
           contour_num, hour_run, plot_title, $
           xtitle, ytitle

DEVICE, /CLOSE
END
