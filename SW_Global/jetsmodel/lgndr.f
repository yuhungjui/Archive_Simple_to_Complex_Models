      subroutine lgndr (my2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
c
c  generate legendre polynomials and their derivatives on the
c  gaussian latitudes
c
c ***input***
c
c  my2:  number of gaussian latitudes from south pole and equator
c  jtrun:  zonal wavenumber truncation limit
c  mlmax: total number of triangular truncation spherical harmonics
c  mlsort: pointer array of 1-d indexs at functions of zonal and
c          total wavenumbers
c  sinl: sin of gaussian latitudes
c
c  ***output***
c
c  poly: associated legendre coefficients
c  dpoly: d(poly)/d(sinl)
c
c ******************************************************************
c
c ref= belousov, s. l., 1962= tables of normalized associated
c        legendre polynomials. pergamon press, new york
c
      dimension poly(mlmax,my2),dpoly(mlmax,my2),sinl(my2)
     *, mlsort(jtrun,jtrun)
c
cjh      parameter (jtrunx= 130)
cjh      dimension pnm(jtrunx+1,jtrunx),dpnm(jtrunx+1,jtrunx)
      dimension pnm(jtrun+1,jtrun),dpnm(jtrun+1,jtrun)
c
c sinl is sin(latitude) = cos(colatitude)
c pnm(np,mp) is legendre polynomial p(n,m) with np=n+1, mp=m+1
c pnm(mp,np+1) is x derivative of p(n,m) with np=n+1, mp=m+1
c
      jtrunp= jtrun+1
      do 1001 j=1,my2
      xx= sinl(j)
      sn= sqrt(1.0-xx*xx)
	sn2i = 1.0/(1.0 - xx*xx)
      rt2= sqrt(2.0)
	c1 = rt2
c
	pnm(1,1) = 1.0/rt2
      theta=-atan(xx/sqrt(1.0-xx*xx))+2.0*atan(1.0)
c
      do 20 n=1,jtrun
	np = n + 1
      fn=n
	fn2 = fn + fn
	fn2s = fn2*fn2
c eq 22
      c1= c1*sqrt(1.0-1.0/fn2s)
      c3= c1/sqrt(fn*(fn+1.0))
	ang = fn*theta
	s1 = 0.0
	s2 = 0.0
	c4 = 1.0
	c5 = fn
	a = -1.0
	b = 0.0
c
      do 27 kp=1,np,2
	k = kp - 1
      s2= s2+c5*sin(ang)*c4
      if (k.eq.n) c4 = 0.5*c4
      s1= s1+c4*cos(ang)
	a = a + 2.0
	b = b + 1.0
      fk=k
	ang = theta*(fn - fk - 2.0)
	c4 = (a*(fn - b + 1.0)/(b*(fn2 - a)))*c4
	c5 = c5 - 2.0
   27 continue
c eq 19
	pnm(np,1) = s1*c1
c eq 21
	pnm(np,2) = s2*c3
   20 continue
c
      do 4 mp=3,jtrunp
	m = mp - 1
      fm= m
	fm1 = fm - 1.0
	fm2 = fm - 2.0
	fm3 = fm - 3.0
      c6= sqrt(1.0+1.0/(fm+fm))
c eq 23
	pnm(mp,mp) = c6*sn*pnm(m,m)
      if (mp - jtrunp) 3,4,4
    3 continue
	nps = mp + 1
c
      do 41 np=nps,jtrunp
	n = np - 1
      fn= n
	fn2 = fn + fn
	c7 = (fn2 + 1.0)/(fn2 - 1.0)
	c8 = (fm1 + fn)/((fm + fn)*(fm2 + fn))
      c= sqrt((fn2+1.0)*c8*(fm3+fn)/(fn2-3.0))
      d= -sqrt(c7*c8*(fn-fm1))
      e= sqrt(c7*(fn-fm)/(fn+fm))
c eq 17
	pnm(np,mp) = c*pnm(np-2,mp-2)
     1            + xx*(d*pnm(np-1,mp-2) + e*pnm(np - 1,mp))
   41 continue
    4 continue
c
      do 50 mp=1,jtrun
      fm= mp-1.0
	fms = fm*fm
      do 50 np=mp,jtrun
      fnp= np
	fnp2 = fnp + fnp
	cf = (fnp*fnp - fms)*(fnp2 - 1.0)/(fnp2 + 1.0)
      cf= sqrt(cf)
c der
      dpnm(np,mp)   = -sn2i*(cf*pnm(np+1,mp) - fnp*xx*pnm(np,mp))
   50 continue
c
      do 71 m=1,jtrun
      do 71 l=m,jtrun
      ml= mlsort(m,l)
      poly(ml,j)= pnm(l,m)
      dpoly(ml,j)=dpnm(l,m)
   71 continue
      dpoly(1,j)= 0.0
 1001 continue
      return
      end
