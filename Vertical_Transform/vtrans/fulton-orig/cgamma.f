      subroutine cgeval( rgamma, rtb, pb, pt, tol, n, c, psihat, ierr )
      integer            n, ierr
      double precision   rgamma, rtb, pb, pt, tol, c(0:n), psihat(0:n)
c
c   Purpose
c
c       Computes the eigenvalues (phase speeds) of the vertical
c       structure problem for constant  gamma = dT/dlogp + kappa*T.
c
c   Arguments
c
c       input
c
c           rgamma  gas constant  R  times the static stability  gamma
c
c           rtb     gas constant  R  times the surface temperature  Tb
c
c           pb, pt  bottom and top pressures
c
c           tol     tolerance for computing the phase speeds  c(k):
c                       tol.gt.0    use absolute tolerance  tol
c                       tol.lt.0    use relative tolerance -tol
c                       tol.eq.0    use relative tolerance  1.0e-22
c
c           n       number of vertical modes desired
c
c       output
c
c           c       phase speeds of the vertical modes
c
c           psihat  normalization constants for the vertical modes
c
c           ierr    error flag (zero on normal return)
c
c   Method
c
c       Newton's method is used to solve the trancendental equations
c       which determine the phase speeds.  The normalization constants
c       for the vertical structure functions are computed by a
c       construction involving the green's function.
c       Double precision is used for all internal calculuations.
c
c   Error conditions
c
c       ierr = 1    rgamma.le.0.0
c       ierr = 2       rtb.le.0.0
c       ierr = 3        pt.le.0.0
c       ierr = 4        pb.le.pt
c       ierr = 5         n.le.0
c       ierr.lt.0   Newton's method failed to converge for vertical
c                   mode  k = -(ierr+1)  (the results for modes
c                   k = 0  to  k = -(ierr+2)  should be correct)
c
c   Language    Fortran 77
c
c   History     Written by Scott R. Fulton in March 1984
c
c
      integer           i, k, nitmax, nstart
      double precision  alpha, beta, cold, ctol, cval
      double precision  eps, mu, nu, pi, z, zc, zt
      double precision  f0, fn, f0p, fnp, g0, gn
      parameter         ( nitmax = 10 )
      double precision  zero, one, two, three, four
      parameter         ( zero = 0.0d0,   one = 1.0d0 )
      parameter         ( two  = 2.0d0, three = 3.0d0, four = 4.0d0 )
c
c   eigenvalue relations for external and internal modes
c
      f0( mu )=-mu*(one-alpha/(two*mu*mu))*sinh(mu*zt)-beta*cosh(mu*zt)
      fn( nu )= nu*(one+alpha/(two*nu*nu))*sin( nu*zt)-beta*cos( nu*zt)
c
c   derivatives of eigenvalue relations for external and internal modes
c
      f0p( mu ) = -(one + alpha/(two*mu*mu) + beta*zt)*sinh( mu*zt )
     2               - mu*zt*(one - alpha/(two*mu*mu))*cosh( mu*zt )
      fnp( nu ) =  (one - alpha/(two*nu*nu) + beta*zt)*sin(  nu*zt )
     2               + nu*zt*(one + alpha/(two*nu*nu))*cos(  nu*zt )
c
c   Green's functions for external and internal modes
c
      g0( z, mu ) = cosh( mu*z ) - (alpha/mu)*sinh( mu*z )
      gn( z, nu ) = cos(  nu*z ) - (alpha/nu)*sin(  nu*z )
c
c   argument checks
c
      ierr = 0
      if ( rgamma.le.zero )  ierr = 1
      if (    rtb.le.zero )  ierr = 2
      if (     pt.le.zero )  ierr = 3
      if (     pb.le.pt   )  ierr = 4
      if (      n.le.0    )  ierr = 5
      if ( ierr.ne.0 )  return
c
c   set up the necessary constants
c
      beta  = rgamma/rtb
      alpha = one/two - beta
      if ( tol.gt.zero )  ctol = tol
      if ( tol.eq.zero )  eps  = 1.0d-22
      if ( tol.lt.zero )  eps  = abs( tol )
      zt = log( pb/pt )
      zc = two*beta/alpha
c
c   calculate the eigenvalue and normalization constant (external mode)
c
      if ( zt.gt.zc )  then
          cval = sqrt( rtb*(one - exp( -zt )) )
          mu = sqrt( one/four - rgamma/cval**2 )
          do 10 i=1,nitmax
              cold = cval
              mu = mu - f0( mu )/f0p( mu )
              cval = sqrt( rgamma/(one/four - mu*mu) )
              if ( tol.le.zero )  ctol = eps*cval
              if ( abs( cval - cold ).le.ctol )  go to 20
   10     continue
          ierr = -1
          return
   20     c(0) = cval
          psihat(0) = sqrt( -two*mu*g0(zero,mu)/(f0p(mu)*g0(zt,mu)) )
          nstart = 1
      else if ( zt.eq.zc )  then
          nu = zero
          c(0) = sqrt( four*rgamma )
          psihat(0) = sqrt(three*alpha/
     &                ((one-two*alpha)*(one+two*alpha+four*alpha**2)))
          nstart = 1
      else
          nstart = 0
      end if
c
c   calculate the eigenvalue and normalization constant (internal modes)
c
      pi = acos( -one )
      do 50 k=nstart,n
          nu = k*pi/zt
          cval = sqrt( rgamma/(one/four + nu*nu) )
          do 30 i=1,nitmax
              cold = cval
              nu = nu - fn( nu )/fnp( nu )
              cval = sqrt( rgamma/(one/four + nu*nu) )
              if ( tol.le.zero )  ctol = eps*cval
              if ( abs( cval - cold ).le.ctol )  go to 40
   30     continue
          ierr = -(k + 1)
          return
   40     c(k) = cval
          psihat(k) = sqrt( two*nu*gn(zero,nu)/(fnp(nu)*gn(zt,nu)) )
   50 continue
      return
      end
      subroutine cgefun( k, np, p, rgamma, rtb, pb, pt, c, psihat, psi )
      integer            k, np
      double precision   p(np), rgamma, rtb, pb, pt, c(0:0), psihat(0:0)
      double precision   psi(np)
c
c   Purpose
c
c       computes the eigenfunctions of the vertical structure
c       problem for constant  gamma = dT/dlogp + kappa*T.
c
c   Arguments
c
c       dimension   p(np), c(0:n), psihat(0:n), psi(np)
c
c       input
c
c           k       index of the vertical structure function desired
c
c           np      number of pressure levels
c
c           p       array of pressure levels
c
c           rgamma  gas constant  R  times the static stability  gamma
c
c           rtb     gas constant  R  times the surface temperature  Tb
c
c           pb, pt  bottom and top pressures
c
c           c       phase speeds of the vertical modes
c
c           psihat  normalization constants for the vertical modes
c
c       output
c
c           psi     psi(i)  is the value of vertical structure
c                   function  k  at  p = p(i)
c
c   Use
c
c       First call  cgeval  to compute  c  and  psihat.
c       Then  call  cgefun  to compute  psi  as desired.
c
c   Language    Fortran 77
c
c   History     Written by Scott R. Fulton in March 1984
c
c
      integer  i
      double precision  alpha, beta, dp, mu, nu, z, zc, zt
      double precision  dpb, drgam, g0, gn
      double precision  one, two, four
      parameter         ( one = 1.0d0, two = 2.0d0, four = 4.0d0 )
c
c   green's functions for external and internal modes
c
      g0( z, mu ) = cosh( mu*z ) - (alpha/mu)*sinh( mu*z )
      gn( z, nu ) = cos(  nu*z ) - (alpha/nu)*sin(  nu*z )
c
c   set up the necessary constants
c
      if ( k.lt.0 )  return
      drgam = rgamma
      beta  = drgam/rtb
      alpha = one/two - beta
      dpb = pb
      zt = log( dpb/pt )
      zc = two*beta/alpha
      dp = dpb - pt
c
c   calculate the eigenfunction (vertical structure function) values
c
      if ( k.eq.0 .and. zt.gt.zc )  then
          mu = sqrt( one/four - drgam/c(0)**2 )
          do 10 i=1,np
   10     psi(i) = sqrt( dp/p(i) )*psihat(0)*g0( log(dpb/p(i)), mu )
      else if ( k.eq.0 .and. zt.eq.zc )  then
          do 20 i=1,np
   20     psi(i) = sqrt( dp/p(i) )*psihat(0)*(one - alpha*log(dpb/p(i)))
      else
          nu = sqrt( drgam/c(k)**2 - one/four )
          do 30 i=1,np
   30     psi(i) = sqrt( dp/p(i) )*psihat(k)*gn( log(dpb/p(i)), nu )
      end if
      return
      end
