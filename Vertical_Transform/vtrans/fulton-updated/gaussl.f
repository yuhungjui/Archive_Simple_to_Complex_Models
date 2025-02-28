      subroutine gaussl( n, xa, xb, ab, wt, ierr )
      integer            n, ierr
      double precision   xa, xb, ab(n), wt(n)
c
c   abscissas and weights for  n point gaussian quadrature on (xa,xb)
c   input arguments:
c       n     the number of gauss-legendre points desired
c       xa    the left  endpoint of the interval of integration
c       xb    the right endpoint of the interval of integration
c   output arguments:
c       ab    the  n  calculated abscissas
c       wt    the  n  calculated weights
c       ierr  set to  9  if there is an error (zero on normal return)
c   source: National Center for Atmospheric Research
c
c   modified by Scott R. Fulton:
c   converted to double precision
c   set the convergence tolerance  tol  internally
c
      integer           i, k, m, maxit, nnp1
      double precision  b, bisq, cond, c1, c2, c3, c4, pi
      double precision  rootbf, tol, u
      double precision  ddif, ddifx, dp, dpm1, dpm2, dppr, dprod
      double precision  dp2pri, drat, dsum, dtmp, dzeri, dzero
c
c   maximum number of iterations before giving up on convergence
c
      parameter  ( maxit = 5)
c
c   machine dependent constants:
c       tol     convergence criterion for double precision iteration
c       c1  )
c       c2  )   coefficients in McMahon's expansions of the  kth  zero
c       c3  )   of the Bessel function j0(x) (cf. Abramowitz and Stegun)
c       c4  )
c
c     parameter  ( tol = 1.0d-22 )
      parameter  ( c1 = 1.0d0/8, c2 = -31.0d0/(2*3*8**2) )
      parameter  ( c3 = 3779.0d0/(2*3*5*8**3) )
      parameter  ( c4 = -6277237.0d0/(3*5*7*8**5) )
c
c   arithmetic statement function for converting integer to double
c
      integer integr
      double precision  dbli
      dbli(integr) = dble(float(integr))
c
c   set the convergence tolerance to the unit roundoff
c
      tol = 1.0
    1 tol = tol/2.0
      if (1.0+tol.gt.1.0 ) go to 1
      tol = 2.0*tol
c
      ierr = 0
      if ( n.lt.1 )  return
      ddif = 0.5d0*(dble(xb) - dble(xa))
      dsum = 0.5d0*(dble(xb) + dble(xa))
      if ( n.eq.1 )  then
          ab(1) = dsum
          wt(1) = 2.0d0*ddif
          return
      end if
      pi = acos( -1.0d0 )
      u = (1.0d0 - (2.0d0/pi)**2)/4.0d0
      cond = 1.0d0/sqrt( (0.5d0+float(n))**2 + u )
      nnp1 = n*(n+1)
c
      do 40 k=1,n/2
          b = (float(k)-0.25d0)*pi
          bisq = 1.0d0/(b*b)
c
c       rootbf approximates the kth zero of the bessel function j0(x)
c
          rootbf = b*(1.0d0+bisq*(c1+bisq*(c2+bisq*(c3+bisq*c4))))
c
c       initial guess for kth root of legendre poly p-sub-n(x)
c
          dzero = cos( rootbf*cond )
          do 20 i=1,maxit
c
c           recursion relation for legendre polynomials
c
              dpm2 = 1.0d0
              dpm1 = dzero
              do 10 m=2,n
                 dp = (dbli(2*m-1)*dzero*dpm1 - dbli(m-1)*dpm2)/dbli(m)
                 dpm2 = dpm1
                 dpm1 = dp
   10         continue
              dtmp = 1.0d0/(1.0d0 - dzero*dzero)
              dppr = dbli(n)*(dpm2 - dzero*dp)*dtmp
              dp2pri = (2.0d0*dzero*dppr - dbli(nnp1)*dp)*dtmp
              drat = dp/dppr
c
c           cubically-convergent iterative improvement of root
c
              dzeri = dzero - drat*(1.0d0 + drat*dp2pri/(2.0d0*dppr))
              if ( abs(dzeri - dzero).le.tol )  go to 30
              dzero = dzeri
   20         continue
              ierr = 9
   30     continue
          ddifx = ddif*dzero
          ab(k) = dsum-ddifx
          wt(k) = 2.0d0*(1.d0 - dzero*dzero)/(dbli(n)*dpm2)**2*ddif
          i = n - k + 1
          ab(i) = dsum + ddifx
          wt(i) = wt(k)
   40 continue
      if ( mod(n,2).eq.0 )  return
      ab(n/2+1) = dsum
      dprod = n
      do 50 k=1,n-1,2
   50 dprod = dbli(n-k-1)*dprod/dbli(n-k)
      wt(n/2+1) = 2.0d0/dprod**2*ddif
      return
      end
