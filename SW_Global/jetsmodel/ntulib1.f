      subroutine gausl3 (n,xa,xb,wt,ab)
c
c weights and abscissas for nth order gaussian quadrature on (xa,xb).
c input arguments
c
c n  -the order desired
c xa -the left endpoint of the interval of integration
c xb -the right endpoint of the interval of integration
c output arguments
c ab -the n calculated abscissas
c wt -the n calculated weights
c
      implicit double precision (a-h,o-z)
c
      real  ab(n) ,wt(n),xa,xb
c
c machine dependent constants---
c  tol - convergence criterion for double precision iteration
c  pi  - given to 15 significant digits
c  c1  -  1/8                     these are coefficients in mcmahon"s
c  c2  -  -31/(2*3*8**2)          expansions of the kth zero of the
c  c3  -  3779/(2*3*5*8**3)       bessel function j0(x) (cf. abramowitz,
c  c4  -  -6277237/(3*5*7*8**5)   handbook of mathematical functions).
c  u   -  (1-(2/pi)**2)/4
c
      data tol/1.d-14/,pi/3.14159265358979/,u/.148678816357662/
      data c1,c2,c3,c4/.125,-.080729166666667,.246028645833333,
     1                -1.82443876720609 /
c
c maximum number of iterations before giving up on convergence
c
      data maxit /5/
c
c arithmetic statement function for converting integer to double
c
      dbli(i) = dble(float(i))
c
      ddif = .5d0*(dble(xb)-dble(xa))
      dsum = .5d0*(dble(xb)+dble(xa))
      if (n .gt. 1) go to 101
      ab(1) = 0.
      wt(1) = 2.*ddif
      go to 107
  101 continue
      nnp1 = n*(n+1)
      cond = 1./sqrt((.5+float(n))**2+u)
      lim = n/2
c
      do 105 k=1,lim
	 b = (float(k)-.25)*pi
	 bisq = 1./(b*b)
c
c rootbf approximates the kth zero of the bessel function j0(x)
c
	 rootbf = b*(1.+bisq*(c1+bisq*(c2+bisq*(c3+bisq*c4))))
c
c      initial guess for kth root of legendre poly p-sub-n(x)
c
	 dzero = cos(rootbf*cond)
	 do 103 i=1,maxit
c
	    dpm2 = 1.d0
	    dpm1 = dzero
c
c       recursion relation for legendre polynomials
c
	    do 102 nn=2,n
		dp = (dbli(2*nn-1)*dzero*dpm1-dbli(nn-1)*dpm2)/dbli(nn)
		dpm2 = dpm1
		dpm1 = dp
  102       continue
	    dtmp = 1.d0/(1.d0-dzero*dzero)
	    dppr = dbli(n)*(dpm2-dzero*dp)*dtmp
	    dp2pri = (2.d0*dzero*dppr-dbli(nnp1)*dp)*dtmp
	    drat = dp/dppr
c
c       cubically-convergent iterative improvement of root
c
	    dzeri = dzero-drat*(1.d0+drat*dp2pri/(2.d0*dppr))
	    ddum= dabs(dzeri-dzero)
	 if (ddum .le. tol) go to 104
	    dzero = dzeri
  103    continue
	 print 504
  504    format(1x,' in gausl3, convergence failed')
  104    continue
	 ddifx = ddif*dzero
	 ab(k) = dsum-ddifx
	 wt(k) = 2.d0*(1.d0-dzero*dzero)/(dbli(n)*dpm2)**2*ddif
	 i = n-k+1
	 ab(i) = dsum+ddifx
	 wt(i) = wt(k)
  105 continue
c
      if (mod(n,2) .eq. 0) go to 107
      ab(lim+1) = dsum
      nm1 = n-1
      dprod = n
      do 106 k=1,nm1,2
	 dprod = dbli(nm1-k)*dprod/dbli(n-k)
  106 continue
      wt(lim+1) = 2.d0/dprod**2*ddif
  107 return
      end
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
      subroutine sortml (jtrun,mlmax,msort,lsort,mlsort)
c
c  sortml builds pointer arrays for functional dependency between
c  1-d spectral index and zonal and total wavenumber indices.  this
c  subroutine reflects the coefficient storage strategy used in the
c  model
c
c ***input***
c
c  jtrun:  zonal and total wavenumber limit
c  mlmax:  total number of 1-d spectral index for triangular trunc
c
c  ***output***
c
c  msort:  zonal wavenumber as function of 1-d spectral index
c  lsort:  total wavenumber as function of 1-d spectral index
c  mlsort: total wavenumber index as function of zonal and total
c          wavenumber
c
c *******************************************************************
c
      dimension msort(mlmax),lsort(mlmax),mlsort(jtrun,jtrun)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
      ml= 0
      do 1 k=1,jtrun-1,2
      do 1 m=1,jtrun-k
      ml= ml+1
      mlp= ml+mlx
      mlsort(m,m+k)= mlp
      mlsort(m,m+k-1)= ml
      msort(ml)= m
      lsort(ml)= m+k-1
      msort(mlp)= m
      lsort(mlp)= m+k
    1 continue
c
      ml= mlp
      do 2 m=2,jtrun,2
      ml= ml+1
      mlsort(m,jtrun)= ml
      msort(ml)= m
      lsort(ml)= jtrun
    2 continue
      return
      end
      subroutine tranrs (jtrun,mlmax,nx,my,ll,poly,w,r,s)
c
c  subroutine to transform a scalar grid point field to spectral
c  coefficients
c
c *** input ***
c
c  jtrun: zonal wavenumber truncation limit
c  mlmax: total number of spherical harmonic coeff  (horizontal field)
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  ll: number of vertical levels to transform
c  poly: legendre polynomials
c  w: gaussian quadrature weights
c  r: 3-dim input grid pt. field to be transformed
c
c *** output ***
c
c  s: spectral coefficient fields
c
c  **********************************
c
      dimension poly(mlmax,my/2),s(mlmax,2,ll),r(nx,ll,my),w(my)
      include '../include/fftcom.h'
csun  include '../include/paramt.h' .. change im,jm to nx,my
      dimension cc(nx+3,my),work(nx*my,2)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
c
      do 30 k=1,ll
c
c  put grid point fields into two dimensional horizontal array
c
      do 23 j=1,my
      do 23 i=1,nx
      cc(i,j)= r(i,k,j)
   23 continue
c
c  fft for each guassian latitude of 2-d field
c
      call rfftmlt(cc,work,trigs,ifax,1,nx+3,nx,my,-1)
c
c  to start quadrature integral we compute contribution without
c  adding to existing sum
c
      m1= 0
      do 62 l=jtrun-1,1,-2
!ocl novrec
cdir$ ivdep
         do 63 m = 1, l
         mm= 2*m-1
         mp= mm+1
         ml= m+m1
         mk= ml+mlx
            s(ml,1,k) = w(1)*poly(ml,1)*(cc(mm,1)+cc(mm,my))
            s(ml,2,k) = w(1)*poly(ml,1)*(cc(mp,1)+cc(mp,my))
            s(mk,1,k) = w(1)*poly(mk,1)*(cc(mm,1)-cc(mm,my))
            s(mk,2,k) = w(1)*poly(mk,1)*(cc(mp,1)-cc(mp,my))
   63    continue
      m1= m1+l
   62 continue
c
      ml= mlx*2
!ocl novrec
cdir$ ivdep
      do 64 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
      s(ml,1,k)= w(1)*poly(ml,1)*(cc(mm,1)+cc(mm,my))
      s(ml,2,k)= w(1)*poly(ml,1)*(cc(mp,1)+cc(mp,my))
   64 continue
c
c  for rest of quadrature integral we add contribution to
c  existing sum
c
      do 70 j=2,my/2
      jj= my-j+1
      m1= 0
      do 72 l=jtrun-1,1,-2
!ocl novrec
cdir$ ivdep
         do 73 m = 1, l
         mm= 2*m-1
         mp= mm+1
         ml= m+m1
         mk= ml+mlx
            s(ml,1,k) = s(ml,1,k)+w(j)*poly(ml,j)*(cc(mm,j)+cc(mm,jj))
            s(ml,2,k) = s(ml,2,k)+w(j)*poly(ml,j)*(cc(mp,j)+cc(mp,jj))
            s(mk,1,k) = s(mk,1,k)+w(j)*poly(mk,j)*(cc(mm,j)-cc(mm,jj))
            s(mk,2,k) = s(mk,2,k)+w(j)*poly(mk,j)*(cc(mp,j)-cc(mp,jj))
   73    continue
      m1= m1+l
   72 continue
c
      ml= mlx*2
!ocl novrec
cdir$ ivdep
      do 65 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
      s(ml,1,k)= s(ml,1,k)+w(j)*poly(ml,j)*(cc(mm,j)+cc(mm,jj))
      s(ml,2,k)= s(ml,2,k)+w(j)*poly(ml,j)*(cc(mp,j)+cc(mp,jj))
   65 continue
   70 continue
   30 continue
c
      return
      end
      subroutine transr (jtrun,mlmax,nx,my,ll,poly,s,r)
c
c  subroutine to transform a spectral coefficient field to
c  grid point form
c
c *** input ***
c
c  jtrun: zonal wavenumber resolution limit
c  mlmax: number of spectral coefficients (horizontal field)
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  ll: number of levels to transform
c  poly: legendre polynomials
c  s: spectral coefficient array to transform
c
c *** output ***
c
c  r: 3-d output grid point fields
c
c  **************************************
c
      dimension poly(mlmax,my/2),s(mlmax,2,ll),r(nx,ll,my)
      include '../include/fftcom.h'
csun  include '../include/paramt.h' .. change im,jm to nx,my
      dimension cc(nx+3,my),work(nx*my,2)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
c
ccc$omp  parallel do
ccc$omp1 private (cc,work,k,j,i,m1,l,m,jj,mm,mp,ml,mk)
ccc$omp2 shared  (trigs,ifax,nx,my,ll,poly,s,w,jtrun,r,mlx)
c
      do 20 k=1,ll
cjh      do 55 m=1,(nx+3)*my/2
cjh      cc(m,1)= 0.0
cjh      cc(m,my/2+1)= 0.0
cjh   55 continue
      do 55 j=1,my
      do 55 i=1,nx+3
      cc(i,j)= 0.0
   55 continue
c
      do 5 j=1,my/2
      jj= my+1-j
      ml= 2*mlx
!ocl novrec
cdir$ ivdep
      do 3 m=2,jtrun,2
      ml= ml+1
      mm= 2*m-1
      mp= mm+1
      cc(mm,j)= poly(ml,j)*s(ml,1,k)
      cc(mp,j)= poly(ml,j)*s(ml,2,k)
      cc(mm,jj)= cc(mm,j)
      cc(mp,jj)= cc(mp,j)
    3 continue
c
      m1= 0
      do 5 l=jtrun-1,1,-2
!ocl novrec
cdir$ ivdep
      do 6 m=1,l
      mm= 2*m-1
      mp= mm+1
      ml= m+m1
      mk= ml+mlx
      cc(mm,j)= cc(mm,j)+poly(ml,j)*s(ml,1,k)+poly(mk,j)*s(mk,1,k)
      cc(mm,jj)=cc(mm,jj)+poly(ml,j)*s(ml,1,k)-poly(mk,j)*s(mk,1,k)
      cc(mp,j)= cc(mp,j)+poly(ml,j)*s(ml,2,k)+poly(mk,j)*s(mk,2,k)
      cc(mp,jj)=cc(mp,jj)+poly(ml,j)*s(ml,2,k)-poly(mk,j)*s(mk,2,k)
    6 continue
      m1= m1+l
    5 continue
c
      call rfftmlt(cc,work,trigs,ifax,1,nx+3,nx,my,1)
c
      do 22 j=1,my
      do 22 i=1,nx
      r(i,k,j)= cc(i,j)
   22 continue
c
   20 continue
ccc$omp end parallel do
c
      return
      end
      subroutine tranuv (jtrun,mlmax,nx,my,ll,radsq,onocos,eps4,cim
     *, poly,dpoly,vor,div,ut,vt)
c
c  subroutine to transform vorticity and divergence to velocity
c  components
c
c *** input ***
c
c  jtrun: zonal wavenumber truncation limit
c  mlmax: total number of spectral coefficients (horizontal field)
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  ll: number of veritical levels to transform
c  radsq: (rad of earth)**2
c  onocos: 1.0/(cos(lat)**2)
c  eps4: spherical harmonic laplacian operator
c  cim: zonal wavenumber array
c  poly: legendre polynomials
c  dpoly: d(poly)/d(sin(lat))
c  vor: spectral vorticity
c  div: spectral divergence
c
c *** output ***
c
c  ut: e-w velocity component
c  vt: n-s velocity component
c
c  ****************************************
c
      dimension onocos(my),eps4(mlmax),cim(mlmax),poly(mlmax,my/2)
     *, dpoly(mlmax,my/2),vor(mlmax,2,ll)
     *, div(mlmax,2,ll),ut(nx,ll,my),vt(nx,ll,my)
c
      include '../include/fftcom.h'
csun  include '../include/paramt.h' .. change im,jm,mlm to nx,my,mlmax
      dimension cu(nx+3,my),cv(nx+3,my),work(nx*my,2)
     *, cfac(mlmax),dfac(mlmax)
c
!ocl novrec
cdir$ ivdep
      do 10 ml=2,mlmax
      dfac(ml)= 1.0/(radsq*eps4(ml))
      cfac(ml)= cim(ml)*dfac(ml)
   10 continue
      dfac(1)= 0.0
      cfac(1)= 0.0
      mlx= (jtrun/2)*((jtrun+1)/2)
c
ccc$omp  parallel do
ccc$omp1 private (cc,work,k,j,i,m1,l,m,jj,mm,mp,ml,mk,rcos,cu,cv)
ccc$omp2 shared  (trigs,ifax,nx,my,ll,dpoly,poly,jtrun,mlx,onocos)
ccc$omp3 shared  (cfac,dfac,vor,div,ut,vt)
c
      do 15 k=1,ll
      do 55 m=1,(nx+3)*my/2
      cu(m,1)= 0.0
      cu(m,my/2+1)= 0.0
      cv(m,1)= 0.0
      cv(m,my/2+1)= 0.0
   55 continue
c
      do 20 j=1,my/2
      rcos= 1.0/onocos(j)
      jj= my+1-j
c
      ml= 2*mlx
!ocl novrec
cdir$ ivdep
      do 50 m=2,jtrun,2
      ml= ml+1
      mm= 2*m-1
      mp= mm+1
      cu(mm,j)=         +cfac(ml)*poly(ml,j)*div(ml,2,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*vor(ml,1,k)
c
      cu(mm,jj)=          +cfac(ml)*poly(ml,j)*div(ml,2,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*vor(ml,1,k)
c
      cu(mp,j)=         -cfac(ml)*poly(ml,j)*div(ml,1,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*vor(ml,2,k)
c
      cu(mp,jj)=          -cfac(ml)*poly(ml,j)*div(ml,1,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*vor(ml,2,k)
c
      cv(mm,j)=         +cfac(ml)*poly(ml,j)*vor(ml,2,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*div(ml,1,k)
c
      cv(mm,jj)=          +cfac(ml)*poly(ml,j)*vor(ml,2,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*div(ml,1,k)
c
      cv(mp,j)=         -cfac(ml)*poly(ml,j)*vor(ml,1,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*div(ml,2,k)
c
      cv(mp,jj)=          -cfac(ml)*poly(ml,j)*vor(ml,1,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*div(ml,2,k)
c
   50 continue
c
      m1= 0
      do 30 l=jtrun-1,1,-2
!ocl novrec
cdir$ ivdep
      do 40 m=1,l
      mm= 2*m-1
      mp= mm+1
      ml= m+m1
      mk= ml+mlx
c
      cu(mm,j)= cu(mm,j)+cfac(ml)*poly(ml,j)*div(ml,2,k)
     *                  +cfac(mk)*poly(mk,j)*div(mk,2,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*vor(ml,1,k)
     *                  +dfac(mk)*dpoly(mk,j)*rcos*vor(mk,1,k)
c
      cu(mm,jj)= cu(mm,jj)+cfac(ml)*poly(ml,j)*div(ml,2,k)
     *                  -cfac(mk)*poly(mk,j)*div(mk,2,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*vor(ml,1,k)
     *                  +dfac(mk)*dpoly(mk,j)*rcos*vor(mk,1,k)
c
      cu(mp,j)= cu(mp,j)-cfac(ml)*poly(ml,j)*div(ml,1,k)
     *                  -cfac(mk)*poly(mk,j)*div(mk,1,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*vor(ml,2,k)
     *                  +dfac(mk)*dpoly(mk,j)*rcos*vor(mk,2,k)
c
      cu(mp,jj)= cu(mp,jj)-cfac(ml)*poly(ml,j)*div(ml,1,k)
     *                  +cfac(mk)*poly(mk,j)*div(mk,1,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*vor(ml,2,k)
     *                  +dfac(mk)*dpoly(mk,j)*rcos*vor(mk,2,k)
c
      cv(mm,j)= cv(mm,j)+cfac(ml)*poly(ml,j)*vor(ml,2,k)
     *                  +cfac(mk)*poly(mk,j)*vor(mk,2,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*div(ml,1,k)
     *                  -dfac(mk)*dpoly(mk,j)*rcos*div(mk,1,k)
c
      cv(mm,jj)= cv(mm,jj)+cfac(ml)*poly(ml,j)*vor(ml,2,k)
     *                  -cfac(mk)*poly(mk,j)*vor(mk,2,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*div(ml,1,k)
     *                  -dfac(mk)*dpoly(mk,j)*rcos*div(mk,1,k)
c
      cv(mp,j)= cv(mp,j)-cfac(ml)*poly(ml,j)*vor(ml,1,k)
     *                  -cfac(mk)*poly(mk,j)*vor(mk,1,k)
     *                  -dfac(ml)*dpoly(ml,j)*rcos*div(ml,2,k)
     *                  -dfac(mk)*dpoly(mk,j)*rcos*div(mk,2,k)
c
      cv(mp,jj)= cv(mp,jj)-cfac(ml)*poly(ml,j)*vor(ml,1,k)
     *                  +cfac(mk)*poly(mk,j)*vor(mk,1,k)
     *                  +dfac(ml)*dpoly(ml,j)*rcos*div(ml,2,k)
     *                  -dfac(mk)*dpoly(mk,j)*rcos*div(mk,2,k)
c
   40 continue
c
      m1= m1+l
   30 continue
   20 continue
c
      call rfftmlt(cu,work,trigs,ifax,1,nx+3,nx,my,1)
      call rfftmlt(cv,work,trigs,ifax,1,nx+3,nx,my,1)
c
      do 22 j=1,my
!ocl novrec
cdir$ ivdep
      do 22 i=1,nx
      ut(i,k,j)= cu(i,j)
      vt(i,k,j)= cv(i,j)
   22 continue
c
   15 continue
cccomp$ end parallel do
      return
      end
      subroutine trandv (jtrun,mlmax,nx,my,ll,ut,vt,w,cim
     *, onocos,poly,dpoly,vor,div)
c
c  subroutine to do grid point velocities to spectral vorticity
c  and divergence
c
c **** input ****
c
c  jtrun: zonal wavenumber truncation limit
c  mlmax: total number of spectral coefficients (horizontal field)
c  nx: e-w dimension no.
c  my: n-w dimension no.
c  ll: number of vertical levels to be transformed
c  ut: e-w velocity component
c  vt: n-s velocity component
c  w: gaussian quadrature weights
c  cim: zonal wavenumber array
c  onocos: 1.o/(cos(lat)**2)
c  poly: legendre polynomials
c  dpoly: d(poly)/d(sin(lat))
c
c *** output ***
c
c  vor: spectral vorticity
c  div: spectral divergence
c
c ********************************************
c
      dimension poly(mlmax,my/2),dpoly(mlmax,my/2),cim(mlmax)
     *, onocos(my),w(my),ut(nx,ll,my),vt(nx,ll,my),vor(mlmax,2,ll)
     *, div(mlmax,2,ll)
c
      include '../include/fftcom.h'
csun  include '../include/paramt.h' .. change im,jm,mlm to nx,my,mlmax
      dimension cc(nx+3,my),dd(nx+3,my)
      dimension work(nx*my,2),cfac(mlmax),dfac(mlmax)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
c
ccc$omp  parallel do
ccc$omp1 private (cc,dd,work,cfac,dfac,k,j,i,m1,l,m,jj,mm,mp,ml,mk)
ccc$omp2 shared  (trigs,ifax,nx,my,mlmax,ll,dpoly,onocos,cim,poly,vor)
ccc$omp3 shared  (div,jtrun,w,ut,vt,mlx)
c
      do 30 k=1,ll
      do 23 j=1,my
!ocl novrec
cdir$ ivdep
      do 23 i=1,nx
      cc(i,j)= ut(i,k,j)
      dd(i,j)= vt(i,k,j)
   23 continue
c
      call rfftmlt(cc,work,trigs,ifax,1,nx+3,nx,my,-1)
      call rfftmlt(dd,work,trigs,ifax,1,nx+3,nx,my,-1)
c
c  to begin quadrature integral we compute contribution without
c  adding to an existing sum
c
!ocl novrec
cdir$ ivdep
      do 82 ml=1,mlmax
      cfac(ml)= w(1)*dpoly(ml,1)
      dfac(ml)= w(1)*onocos(1)*cim(ml)*poly(ml,1)
   82 continue
c
      m1= 0
      do 62 l=jtrun-1,1,-2
!ocl novrec
cdir$ ivdep
         do 63 m = 1, l
         mm= 2*m-1
         mp= mm+1
         ml= m+m1
         mk= ml+mlx
            vor(ml,1,k)=  cfac(ml)*(cc(mm,1)
     *     -cc(mm,my))-dfac(ml)*(dd(mp,1)+dd(mp,my))
c
            vor(mk,1,k)=  cfac(mk)*(cc(mm,1)
     *     +cc(mm,my))-dfac(mk)*(dd(mp,1)-dd(mp,my))
c
            vor(ml,2,k)=  cfac(ml)*(cc(mp,1)
     *     -cc(mp,my))+dfac(ml)*(dd(mm,1)+dd(mm,my))
c
            vor(mk,2,k)=  cfac(mk)*(cc(mp,1)
     *     +cc(mp,my))+dfac(mk)*(dd(mm,1)-dd(mm,my))
c
            div(ml,1,k)=- cfac(ml)*(dd(mm,1)
     *     -dd(mm,my))-dfac(ml)*(cc(mp,1)+cc(mp,my))
c
            div(mk,1,k)=- cfac(mk)*(dd(mm,1)
     *     +dd(mm,my))-dfac(mk)*(cc(mp,1)-cc(mp,my))
c
            div(ml,2,k)=- cfac(ml)*(dd(mp,1)
     *     -dd(mp,my))+dfac(ml)*(cc(mm,1)+cc(mm,my))
c
            div(mk,2,k)=- cfac(mk)*(dd(mp,1)
     *     +dd(mp,my))+dfac(mk)*(cc(mm,1)-cc(mm,my))
c
   63    continue
      m1= m1+l
   62 continue
c
      ml= mlx*2
!ocl novrec
cdir$ ivdep
      do 64 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
            vor(ml,1,k)=  cfac(ml)*(cc(mm,1)
     *     -cc(mm,my))-dfac(ml)*(dd(mp,1)+dd(mp,my))
c
            vor(ml,2,k)=  cfac(ml)*(cc(mp,1)
     *     -cc(mp,my))+dfac(ml)*(dd(mm,1)+dd(mm,my))
c
            div(ml,1,k)=- cfac(ml)*(dd(mm,1)
     *     -dd(mm,my))-dfac(ml)*(cc(mp,1)+cc(mp,my))
c
            div(ml,2,k)=- cfac(ml)*(dd(mp,1)
     *     -dd(mp,my))+dfac(ml)*(cc(mm,1)+cc(mm,my))
c
   64 continue
c
c  now do rest of guassian latitudes by adding to accumulating
c  sum
c
      do 70 j=2,my/2
c
!ocl novrec
cdir$ ivdep
      do 84 ml=1,mlmax
      cfac(ml)= w(j)*dpoly(ml,j)
      dfac(ml)= w(j)*onocos(j)*cim(ml)*poly(ml,j)
   84 continue
c
      jj= my-j+1
      m1= 0
      do 72 l=jtrun-1,1,-2
!ocl novrec
cdir$ ivdep
         do 73 m = 1, l
         mm= 2*m-1
         mp= mm+1
         ml= m+m1
         mk= ml+mlx
            vor(ml,1,k)= vor(ml,1,k)+cfac(ml)*(cc(mm,j)
     *     -cc(mm,jj))-dfac(ml)*(dd(mp,j)+dd(mp,jj))
c
            vor(mk,1,k)= vor(mk,1,k)+cfac(mk)*(cc(mm,j)
     *     +cc(mm,jj))-dfac(mk)*(dd(mp,j)-dd(mp,jj))
c
            vor(ml,2,k)= vor(ml,2,k)+cfac(ml)*(cc(mp,j)
     *     -cc(mp,jj))+dfac(ml)*(dd(mm,j)+dd(mm,jj))
c
            vor(mk,2,k)= vor(mk,2,k)+cfac(mk)*(cc(mp,j)
     *     +cc(mp,jj))+dfac(mk)*(dd(mm,j)-dd(mm,jj))
c
            div(ml,1,k)= div(ml,1,k)-cfac(ml)*(dd(mm,j)
     *     -dd(mm,jj))-dfac(ml)*(cc(mp,j)+cc(mp,jj))
c
            div(mk,1,k)= div(mk,1,k)-cfac(mk)*(dd(mm,j)
     *     +dd(mm,jj))-dfac(mk)*(cc(mp,j)-cc(mp,jj))
c
            div(ml,2,k)= div(ml,2,k)-cfac(ml)*(dd(mp,j)
     *     -dd(mp,jj))+dfac(ml)*(cc(mm,j)+cc(mm,jj))
c
            div(mk,2,k)= div(mk,2,k)-cfac(mk)*(dd(mp,j)
     *     +dd(mp,jj))+dfac(mk)*(cc(mm,j)-cc(mm,jj))
c
   73    continue
      m1= m1+l
   72 continue
c
      ml= mlx*2
!ocl novrec
cdir$ ivdep
      do 65 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
            vor(ml,1,k)= vor(ml,1,k)+cfac(ml)*(cc(mm,j)
     *     -cc(mm,jj))-dfac(ml)*(dd(mp,j)+dd(mp,jj))
c
            vor(ml,2,k)= vor(ml,2,k)+cfac(ml)*(cc(mp,j)
     *     -cc(mp,jj))+dfac(ml)*(dd(mm,j)+dd(mm,jj))
c
            div(ml,1,k)= div(ml,1,k)-cfac(ml)*(dd(mm,j)
     *     -dd(mm,jj))-dfac(ml)*(cc(mp,j)+cc(mp,jj))
c
            div(ml,2,k)= div(ml,2,k)-cfac(ml)*(dd(mp,j)
     *     -dd(mp,jj))+dfac(ml)*(cc(mm,j)+cc(mm,jj))
c
   65 continue
   70 continue
   30 continue
ccc$omp end parallel do
c
      return
      end
      subroutine alpha2s(jtrun,mlmax,nx,my,ll,ut,vt,w,cim
     *, onocos,poly,dpoly,div)
c
c  subroutine to do grid point velocities to spectral vorticity
c  and divergence
c
c **** input ****
c
c  jtrun: zonal wavenumber truncation limit
c  mlmax: total number of spectral coefficients (horizontal field)
c  nx: e-w dimension no.
c  my: n-w dimension no.
c  ll: number of vertical levels to be transformed
c  ut: e-w velocity component
c  vt: n-s velocity component
c  w: gaussian quadrature weights
c  cim: zonal wavenumber array
c  onocos: 1.o/(cos(lat)**2)
c  poly: legendre polynomials
c  dpoly: d(poly)/d(sin(lat))
c
c *** output ***
c
c  div: spectral divergence
c
c ********************************************
c
      dimension poly(mlmax,my/2),dpoly(mlmax,my/2),cim(mlmax)
     *, onocos(my),w(my),ut(nx,ll,my),vt(nx,ll,my),div(mlmax,2,ll)
c
      include '../include/fftcom.h'
c      include '../include/paramt.h'
      dimension cc(nx+3,my),dd(nx+3,my)
      dimension work(nx*my,2),cfac(mlmax),dfac(mlmax)
c
c      jtrun= 1+(im-1)/3
      mlx= (jtrun/2)*((jtrun+1)/2)
c
      do 30 k=1,ll
      do 23 j=1,my
cdir$ ivdep
      do 23 i=1,nx
      cc(i,j)= ut(i,k,j)
      dd(i,j)= vt(i,k,j)
   23 continue
c
      call rfftmlt(cc,work,trigs,ifax,1,nx+3,nx,my,-1)
      call rfftmlt(dd,work,trigs,ifax,1,nx+3,nx,my,-1)
c  to begin quadrature integral we compute contribution without
c  adding to an existing sum
c
cdir$ ivdep
      do 82 ml=1,mlmax
      cfac(ml)= w(1)*dpoly(ml,1)
      dfac(ml)= w(1)*onocos(1)*cim(ml)*poly(ml,1)
   82 continue
c      
      m1= 0
      do 62 l=jtrun-1,1,-2
cdir$ ivdep 
         do 63 m = 1, l
         mm= 2*m-1
         mp= mm+1
         ml= m+m1
         mk= ml+mlx
c
            div(ml,1,k)=- cfac(ml)*(dd(mm,1)
     *     -dd(mm,my))-dfac(ml)*(cc(mp,1)+cc(mp,my))
c
            div(mk,1,k)=- cfac(mk)*(dd(mm,1)
     *     +dd(mm,my))-dfac(mk)*(cc(mp,1)-cc(mp,my))
c
            div(ml,2,k)=- cfac(ml)*(dd(mp,1)
     *     -dd(mp,my))+dfac(ml)*(cc(mm,1)+cc(mm,my))
c
            div(mk,2,k)=- cfac(mk)*(dd(mp,1)
     *     +dd(mp,my))+dfac(mk)*(cc(mm,1)-cc(mm,my))
c
   63    continue
      m1= m1+l
   62 continue
c
      ml= mlx*2
cdir$ ivdep
      do 64 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
c
            div(ml,1,k)=- cfac(ml)*(dd(mm,1)
     *     -dd(mm,my))-dfac(ml)*(cc(mp,1)+cc(mp,my))
c
            div(ml,2,k)=- cfac(ml)*(dd(mp,1)
     *     -dd(mp,my))+dfac(ml)*(cc(mm,1)+cc(mm,my))
c
   64 continue
c
c  now do rest of guassian latitudes by adding to accumulating
c  sum
c  
      do 70 j=2,my/2
c  
cdir$ ivdep
      do 84 ml=1,mlmax
      cfac(ml)= w(j)*dpoly(ml,j)
      dfac(ml)= w(j)*onocos(j)*cim(ml)*poly(ml,j)
   84 continue
c
      jj= my-j+1
      m1= 0
      do 72 l=jtrun-1,1,-2
cdir$ ivdep
         do 73 m = 1, l
         mm= 2*m-1
         mp= mm+1
         ml= m+m1
         mk= ml+mlx
c
            div(ml,1,k)= div(ml,1,k)-cfac(ml)*(dd(mm,j)
     *     -dd(mm,jj))-dfac(ml)*(cc(mp,j)+cc(mp,jj))
c
            div(mk,1,k)= div(mk,1,k)-cfac(mk)*(dd(mm,j)
     *     +dd(mm,jj))-dfac(mk)*(cc(mp,j)-cc(mp,jj))
c
            div(ml,2,k)= div(ml,2,k)-cfac(ml)*(dd(mp,j)
     *     -dd(mp,jj))+dfac(ml)*(cc(mm,j)+cc(mm,jj))
c
            div(mk,2,k)= div(mk,2,k)-cfac(mk)*(dd(mp,j)
     *     +dd(mp,jj))+dfac(mk)*(cc(mm,j)-cc(mm,jj))
c
   73    continue
      m1= m1+l
   72 continue
c
      ml= mlx*2
cdir$ ivdep
      do 65 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
c
            div(ml,1,k)= div(ml,1,k)-cfac(ml)*(dd(mm,j)
     *     -dd(mm,jj))-dfac(ml)*(cc(mp,j)+cc(mp,jj))
c
            div(ml,2,k)= div(ml,2,k)-cfac(ml)*(dd(mp,j)
     *     -dd(mp,jj))+dfac(ml)*(cc(mm,j)+cc(mm,jj))
c
   65 continue
   70 continue
   30 continue
c
      return
      end
