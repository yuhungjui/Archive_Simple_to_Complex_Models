      block data vtdoc
c
c   Revision History
c
c       05/01/84:  original version
c       02/22/99:  updated for use under UNIX (double precision)
c
c *********************** Introduction to  VTRANS **********************
c
c      VTRANS  is a collection of routines for solving the vertical
c          structure problem and performing vertical transforms.
c
c   Definitions
c
c       The vertical structure problem is
c
c                    d  (  1   dpsi)   psi
c           l(psi) = -- (-----*----) + --- = 0     (0.lt.pt.le.p.le.pb),
c                    dp (sigma  dp )   c*c
c
c           dpsi                     dpsi   p*sigma
c           ---- = 0  at  p = pt,    ---- + -------*psi = 0  at  p = pb.
c            dp                       dp      R*T
c
c       Here  p  is pressure,  psi  is the vertical structure function,
c       c  is the phase speed,  sigma is the static stability,  R  is
c       the gas constant,  T  is the basic state temperature, and  pt
c       and  pb  are the top and bottom pressures, respectively.
c
c       The  routine  vtset  solves this problem (approximately)
c       using a basic state temperature profile supplied by the user.
c       The resulting eigenvalues  c(k)  (k=0, ..., n)  are real and
c       the corresponding eigenfunctions  psi(k,p)  (k=0, ..., n)
c       are orthonormal  in the vertical inner product
c
c                        1
c           ( u, v ) = ----- * integral (p=pt to p=pb) of u(p)*v(p),
c                      pb-pt
c
c       where  u  and  v  denote any functions of  p.  The forward
c       transform of a function function  f(p)  is defined by
c
c           fhat(k) = ( f(p), psi(k,p) )
c
c       and the corresponding inverse transform is
c
c           f(p) = sum (k=0 to n) of fhat(k)*psi(k,p).
c
c       This transform pair has the property that for any two functions
c       f  and  g  with transforms  fhat  and  ghat,  respectively,
c       their vertical inner product  is given by
c
c           ( f, g ) = sum (k=0 to n) of fhat(k)*ghat(k)
c
c       and their energy inner product is given by
c
c           ( f, l(g) ) = sum (k=0 to n) of fhat(k)*ghat(k)/c(k)**2
c
c       Thus if  u(p)  is the wind speed, the kinetic energy associated
c       with vertical mode  k  is  uhat(k)**2,  and if  phi  is the
c       perturbation geopotential, the available potential energy
c       associated with vertical mode  k  is  (phihat(k)/c(k))**2.
c
c       Routines from this package compute numerical approximations
c       to these transforms for functions supplied at discrete pressure
c       levels utilizing a continuous function fit to the specified
c       data. This continuous fit takes the form of a linear combination
c       of vertical basis functions (Chebyshev polynomials modified to
c       satisfy the boundary conditions), with the coefficients in this
c       linear combination referred to as the 'vbf' coefficients of the
c       function.  Thus a forward transform of a function is computed
c       from discrete data by first obtaining the vbf coefficients
c       (routine  vproj)  and then computing the spectral coefficients
c       (routine  vtran  with the matrix  tfor).  Similarly, an inverse
c       transform is computed from these spectral coefficients by first
c       obtaining the vbf coefficients (routine  vtran  with the matrix
c       tinv)  and then evaluating that representation at any desired
c       pressure levels (routine  veval).  This may be represented as:
c
c                     vproj                  vtran  (tfor)
c             *------------------->     *--------------------->
c      physical                     vbf                      spectral
c   representation             representation             representation
c   **************             **************             **************
c    p(i), f(p(i))                 fvbf(j)                    fhat(k)
c    (i=1, ..., m)              (j=0, ..., n)              (k=0, ..., n)
c                     veval                  vtran  (tinv)
c             <-------------------*     <---------------------*
c
c   Contents
c
c       User-level routines
c
c           vtset   solves the vertical structure problem for a
c                   temperature profile given at discrete pressures
c                   and computes the transform matrices  tfor  and  tinv
c
c           vproj   computes the vbf representation of a function
c
c           vtran   computes the forward (vbf to spectral) and
c                   inverse (spectral to vbf) vertical transforms
c
c           veval   evaluates the vbf representation of a function
c
c       Internal routines
c
c           rrgset  sets up the temperature and stability approximations
c
c           rrgval  evaluates the stability approximation
c
c           funset  sets up a continuous approximation to function data
c
c           funval  evaluates the continuous approximation to a function
c
c           vbfval  evaluates the vertical basis functions
c
c           csset   sets up a cubic spline interpolate
c
c           csval   evaluates a cubic spline interpolate
c
c           gaussl  computes gauss-legendre weights and abscissas
c
c   Required Library Routine
c
c       rsg     real symmetric generalized eigenproblem solver (eispack)
c
c   Language    Fortran 77
c
c   Reference   Fulton, S. R. and W. H. Schubert, 1985:  
c               Vertical normalmode transforms:  theory and application.  
c               Monthly Weather Review, vol. 113, 647-658.
c
c   History     Written in April 1983 by Scott R. Fulton
c
c
      end
      subroutine vtset( m, p, rt, n, w, work, con, c, tfor, tinv, ierr )
      integer           m, n, ierr
      double precision  p(m), rt(m,5), w(0:n,0:1), work(0:n,0:n)
      double precision  con(4), c(0:n), tfor(0:n,0:n), tinv(0:n,0:n)
c
c
c   Purpose
c
c       Given a basic state temperature profile this routine solves the
c       vertical structure problem and generates the transform matrices
c
c   Arguments
c
c       input
c
c           m       number of pressure levels at which the basic
c                   state temperature profile is supplied (m.ge.3)
c
c           p       array of pressure levels (in increasing order).
c                   the minimum and maximum pressure levels input
c                   become the top and bottom pressures, respectively.
c                   note:  the top pressure must be positive.
c
c           rt      array containing the corresponding (absolute)
c                   temperatures  T  multiplied by the gas constant  R:
c                       rt(i,1) = R*T( p(i) )   (i = 1, ..., m)
c
c           n       desired spectral truncation (index of the last
c                   vertical mode)--must have  n.ge.sqrt( m+4 )
c
c       work space
c
c           w       array of length at least  2*(n+1)
c
c           work    array of length at least  (n+1)**2
c
c       output
c
c           rt      the storage locations  rt(i,j)  (i=1, ..., m)  for
c                   j=2, ..., 5  contain parameters defining the cubic
c                   spline fit to  R*T  and the Chebyshev polynomial
c                   fit to the reciprocal of the static stability
c                   (the input values  rt(i,1)  are unchanged).
c
c           work    contains the vertical inner products
c                   of the vertical basis functions
c
c           con     contains constants to be used by vtrans routines
c
c           c       contains the phase speeds--in the same
c                   units as sqrt( rt )--in decreasing order
c
c           tfor    forward transform matrix (to be used by  vtran)
c
c           tinv    inverse transform matrix (to be used by  vtran)
c                   note:  the kth column of  tinv  consists of the vbf
c                   coefficients of vertical structure function  psi(k),
c                   so that function may be evaluated by calling
c                       veval( tinv(0,k), md, m, p, con, n, w, psival )
c                   (see the comments in the routine  veval)
c
c           ierr    integer error flag  (zero on normal return)
c
c           note:   upon return, the cubic spline fit to  r*t  may be
c                   evaluated at any pressure level  pl  as
c                       value = csval( pl, md, m, p, rt, 0 )
c                   (where  md  is the derivative desired, e.g.  md=0
c                   for the value, md=1 for the derivative with respect
c                   to  p, etc.), and the reciprocal of the static
c                   stability  sigma  may be evaluated at  pl  as
c                       value = pl**2*rrgval( pl, con, m, rt, ierr )
c                   (where  ierr=0  is returned if  sigma  is positive)
c
c   Method
c
c       The vertical structure problem is solved by the Rayleigh-Ritz-
c       Galerkin method using polynomial basis functions (Chebyshev
c       polynomials modified to satisfy the boundary conditions).  The
c       inner products of the basis functions are computed by numerical
c       quadrature, using a least-squares Chebyshev polynomial fit to
c       the reciprocal of the static stability  gamma.  The resulting
c       generalized eigenproblem is solved by the eispack routine  rsg.
c
c   Timing
c
c       With a typical (tropical) basic state temperature profile
c       (supplied at 23 pressure levels), execution time was about
c       0.00016*n**2  seconds on the ncar cray1 for  8 .le. n .le. 64.
c
c   Error Conditions
c
c       ierr =  1   m  too small  (m  must be at least  3)
c
c       ierr =  2   pressure levels not in increasing order
c
c       ierr =  3   top pressure not positive
c
c       ierr =  4   surface temperature not positive
c
c       ierr =  5   stability not strictly positive on the domain
c
c       ierr =  6   n  too small  (n  must be at least  sqrt( m+4 ))
c
c       ierr =  7   vertical inner product matrix not positive definite
c
c       ierr =  8   energy   inner product matrix not positive definite
c
c       ierr =  9   problems in computing weights and abscissas
c
c       ierr = 10+i convergence failed for phase speed  i
c                   (phase speeds  0, 1, ..., i-1  are okay)
c
c
      double precision rrgval
      integer          i, ig, j, k, ndeg, ng, ngauss
      double precision beta, dp, fact, pb, pg, pt, rrg, rrgb, rtb
      double precision wg, wt
      common  /rrgcom/ ndeg, ngauss
      save    /rrgcom/
c
c ************************** initializations ***************************
c
c   set up the temperature and stablity approximations
c
      ierr = 0
      call rrgset( m, p, rt, ierr )
      if ( ierr.ne.0 )  return
c
c   set up the constants array
c
      pt = p(1)
      if ( pt.le.0.0d0 )  ierr = 3
      pb = p(m)
      dp = pb - pt
      rtb = rt(m,1)
      if ( rtb.le.0.0d0 )  ierr = 4
      con(1) = pt
      con(2) = pb
      con(3) = dp
      rrgb = rrgval( pb, con, m, rt, ierr )
      if ( ierr.ne.0 )  return
      beta = 1.0d0/(rtb*rrgb)
      con(4) = beta
c
c ************* inner products of vertical basis functions *************
c
c   compute the gauss-legendre abscissas and weights
c
      ng = n + 2 + (ndeg+1)/2
      if ( mod( ndeg, 2 ).eq.0 )  ng = ng + 1
      if ( 2*ng.gt.(n+1)**2 )  ierr = 6
      if ( ierr.ne.0 )  return
      call gaussl( ng, pt, pb, tinv(0,0), tinv(ng,0), ierr )
      if ( ierr.ne.0 )  return
      do 10 i=0,n
      do 10 j=i,n
      tfor(i,j) = 0.0d0
   10 work(i,j) = 0.0d0
c
c   compute the integrals using gauss-legendre quadrature (ng  points)
c
      do 30 ig=1,ng
          pg = tinv(   ig-1,0)
          wg = tinv(ng+ig-1,0)
          rrg = rrgval( pg, con, m, rt, ierr )
          if ( ierr.ne.0 )  return
          wt = wg*rrg*pg**2
          call vbfval( pg, 1, con, n, w )
          do 20 i=0,n
          do 20 j=i,n
          tfor(i,j) = tfor(i,j) + wt*w(i,1)*w(j,1)
   20     work(i,j) = work(i,j) + wg*w(i,0)*w(j,0)
   30 continue
c
c   add in the boundary term and scale the matrices
c
      call vbfval( pb, 0, con, n, w )
      fact = beta*rrgb*pb
      do 40 i=0,n
      do 40 j=i,n
      tfor(i,j) = (tfor(i,j) + fact*w(i,0)*w(j,0))/dp
   40 work(i,j) =  work(i,j)/dp
c
c ************** solution of the generalized eigenproblem **************
c
      call rsg( n+1, n+1, tfor, work, c, 1, tinv, w(0,0), w(0,1), ierr )
      if ( ierr.ne.0 )  then
          if ( ierr.eq.7*n+8 )  then
              ierr = 7
          else
              ierr = 9 + ierr
          end if
          return
      end if
c
c   compute the phase speeds from the eigenvalues
c
      do 50 i=0,n
          if ( c(i).le.0.0d0 )  ierr = 8
          if ( ierr.ne.0 )  return
          c(i) = 1.0d0/sqrt( c(i) )
   50 continue
c
c   compute the forward transformation matrix
c
      do 60 i=0,n
      do 60 j=i+1,n
   60 work(j,i) = work(i,j)
      do 80 i=0,n
      do 80 j=0,n
          tfor(i,j) = 0.0d0
          do 70 k=0,n
   70     tfor(i,j) = tfor(i,j) + tinv(k,i)*work(k,j)
   80 continue
      return
      end
      subroutine vproj( m,p,f,in,ip,ng,pg,wg,prom,n,con,w,fvbf,ierr )
      integer           m, in, ip, ng, n, ierr
      double precision  p(m), f(m,4), pg(ng), wg(ng), prom(0:n,0:n)
      double precision  con(4), w(0:n), fvbf(0:n)
c
c
c   Purpose
c
c       computes the vbf representation of a function
c
c   Arguments
c
c       input
c
c                   ************** definition of function **************
c
c           m       number of pressure levels at which the function is
c                   specified (m  must be at least 3)
c
c           p       array of pressure levels at which the function is
c                   specified (must be in increasing order, extending
c                   from the top pressure to the bottom pressure as
c                   defined in the call to  vtset)
c
c           f       array containing the corresponding function values
c                   (fval(i,1)  is the value at  p(i), i = 1, ..., m)
c
c           in      integer specifying whether or not to integrate:
c                       in.gt.0    integrate the function with respect
c                                  to  -log(p)  before projecting it
c                                  (the integral is zero at p = pb)
c                       in.le.0    do not integrate the function
c
c                   ************* definition of projection *************
c                   the function is projected using a least-squares fit
c                   in the norm generated by the discrete inner product
c                     (u,v) = sum (i=1 to ng) of wg(i)*u(pg(i))*v(pg(i))
c                   where  pg  and  wg  may be supplied or generated.
c
c           ip      integer specifying the inner product desired:
c                       ip = 1    Chebyshev inner product
c                       ip = 2    vertical inner product
c                       ip = 3    user-supplied inner product
c                   any other value of  ip  implies that the projection
c                   was defined in a previous call to  vproj  and that
c                   ng, pg, wg  and  prom  are unchanged from that call.
c
c           ng      number of gaussian pressure levels
c
c           pg      array for or containing the gaussian pressure levels
c                   (generated by  vproj  when  ip = 1  or  ip = 2)
c
c           wg      array for or containing the corresponding weights
c                   (generated by  vproj  when  ip = 1  or  ip = 2)
c
c                   ********** examples of typical projections *********
c                   ip=1, ng.ge.2*n:  approximate uniform fit
c                   ip=2, ng.ge.2*n:  continuous least-squares fit
c
c           prom    array for or containing the projection matrix
c                   (generated by  vproj  when  ip = 1, 2  or  3)
c
c           n       spectral truncation (index of last vertical mode)
c
c           con     vector of constants produced by  vtset
c
c       work space
c
c           w       vector of length at least  n+1
c
c       output
c
c           f       the storage locations  f(i,j)  (i=1, ..., m)  for
c                   j = 1, 2, 3, 4  contain parameters defining the
c                   cubic spline fit to the function that was projected,
c                   i.e., the function as input or its integral.
c
c           pg, wg  contain levels and weights generated (ip = 1  or  2)
c
c           prom    contains the matrix generated (ip = 1, 2  or  3)
c
c           fvbf    vbf coefficients of the function:  fvbf(j)  is
c                   the coefficient of vertical basis function  j
c                   in the vbf expansion of the function
c
c           ierr    integer error flag  (zero on normal return)
c
c           note:   upon return, the cubic spline fit to the function
c                   that was projected may be evaluated at any pressure
c                   level  pl  as    value = funval( pl, md, m, p, f ),
c                   where  md  is the derivative desired, e.g.  md=0
c                   for the value, md=1 for the derivative with respect
c                   to  p, etc.).
c
c   Method
c
c       A cubic spline interpolate of the function is set up (this is
c       integrated if requested) and projected onto the vertical basis
c       functions using a least-squares (Galerkin) projection in the
c       specified inner product.
c
c   Error Conditions
c
c       ierr = 1    m  too small  (m must be at least  3)
c
c       ierr = 2    pressure levels not in increasing order
c
c       ierr = 3    minimum and maximum pressures do not match the
c                   top and surface pressures defined in  vtset
c
c
      double precision funval
      integer          i, ig, j, k
      double precision dp, fact, pb, pibyng, pt, wf
c
c   initializations and argument checks
c
      ierr = 0
      pt = con(1)
      pb = con(2)
      dp = con(3)
      if ( p(1).ne.pt .or. p(m).ne.pb )  ierr = 3
      if ( ng.le.0 )  ierr = 4
      if ( ierr.ne.0 )  return
      if ( ip.lt.1 .or. ip.gt.3 )  go to 100
c
c ***** set up the inner product and associated projection matrix ******
c
c   compute the abscissas and weights defining the inner product
c
      if ( ip.eq.1 )  then
          pibyng = acos( -1.0d0 )/ng
          do 10 ig=1,ng
              pg(ig) = pb - dp*(1.0d0 + cos( (ig - 0.5d0)*pibyng ))/2
              wg(ig) = 1.0d0
   10     continue
      else if ( ip.eq.2 )  then
          call gaussl( ng, pt, pb, pg, wg, ierr )
          if ( ierr.ne.0 )  return
          do 20 ig=1,ng
   20     wg(ig) = wg(ig)/dp
      end if
c
c   generate the projection matrix
c
      do 30 i=0,n
      do 30 j=i,n
   30 prom(i,j) = 0.0d0
      do 50 ig=1,ng
          call vbfval( pg(ig), 0, con, n, w )
          do 40 i=0,n
          do 40 j=i,n
   40     prom(i,j) = prom(i,j) + wg(ig)*w(i)*w(j)
   50 continue
c
c   compute its cholesky decomposition
c
      do 70 i=0,n
      do 70 j=i,n
          fact = prom(i,j)
          do 60 k=0,i-1
   60     fact = fact - prom(i,k)*prom(j,k)
          if ( i.ne.j )  then
              prom(j,i) = fact/prom(i,i)
          else if ( fact.gt.0.0d0 )  then
              prom(i,i) = sqrt( fact )
          else
              ierr = 8
              return
          end if
   70 continue
c
c ******************* calculation of the projection ********************
c
c   compute the right-hand side of the linear system
c
  100 call funset( m, p, f, in, con, ierr )
      if ( ierr.ne.0 )  return
      do 110 i=0,n
  110 fvbf(i) = 0.0d0
      do 130 ig=1,ng
          call vbfval( pg(ig), 0, con, n, w )
          wf = wg(ig)*funval( pg(ig), 0, m, p, f )
          do 120 i=0,n
  120     fvbf(i) = fvbf(i) + wf*w(i)
  130     continue
c
c   solve the linear system for the vbf coefficients
c
      fvbf(0) = fvbf(0)/prom(0,0)
      do 150 i=1,n
          do 140 j=0,i-1
  140     fvbf(i) = fvbf(i) - prom(i,j)*fvbf(j)
          fvbf(i) = fvbf(i)/prom(i,i)
  150 continue
      fvbf(n) = fvbf(n)/prom(n,n)
      do 170 i=n-1,0,-1
          do 160 j=i+1,n
  160     fvbf(i) = fvbf(i) - prom(j,i)*fvbf(j)
          fvbf(i) = fvbf(i)/prom(i,i)
  170 continue
      return
      end
      subroutine vtran( f, n, tmat, g )
      integer           n
      double precision  f(0:n), tmat(0:n,0:n), g(0:n)
c
c
c   Purpose
c
c       computes a vertical transform of a function
c       between its vbf and spectral represnetations
c
c   Arguments
c
c       input
c
c           f       vector containing input coefficients of function--
c                   for forward transform:  vbf      coefficients  fvbf
c                   for inverse transform:  spectral coefficients  fhat
c
c           n       spectral truncation (index of last vertical mode)
c
c           tmat    transform matrix produced by  vtset--
c                   for forward transform:  tfor
c                   for inverse transform:  tinv
c
c       output
c
c           g       vector containing output coefficients of function--
c                   for forward transform:  spectral coefficients  fhat
c                   for inverse transform:  vbf      coefficients  fvbf
c
c
      integer          i, j
c
      do 20 i=0,n
          g(i) = 0.0d0
          do 10 j=n,0,-1
   10     g(i) = g(i) + tmat(i,j)*f(j)
   20 continue
      return
      end
      subroutine veval( fvbf, md, m, p, con, n, w, fval )
      integer           md, m, n
      double precision  fvbf(0:n), p(m), con(4), w(0:n,2), fval(m,2)
c
c   Purpose
c
c       evaluates the vbf representation of a function
c
c   Arguments
c
c       dimension   fvbf(0:n), p(m), con(4), w(0:n,md+1), fval(m,md+1)
c
c       input
c
c           fvbf    vbf coefficients of the function
c
c           md      order of the highest derivative desired:
c                       md = 0    evaluate the function itself
c                       md = 1    evaluate the function and derivative
c
c           m       number of output pressure levels
c
c           p       array of output pressure levels (in any order)
c
c           con     vector of constants produced by  vtset
c
c           n       spectral truncation (index of last vertical mode)
c
c       work space
c
c           w       vector of length at least  (md+1)*(n+1)
c
c       output
c
c           fval    vector containing the values requested:
c                       fval(i,1)  is the value      at  p(i)
c                       fval(i,2)  is the derivative at  p(i)  (if md=1)
c                   for  i = 1, ..., m
c   Note:  the kth column of  tinv  consists of the vbf coefficients of
c          vertical structure function  psi(k),  so that function may be
c          evaluated by calling  veval  with  fvbf = tinv(0,k).
c
c
      integer  i, j
c
c   evaluate the vbf representation at the desired  p  values
c
      do 30 j=1,m
          call vbfval( p(j), md, con, n, w )
          fval(j,1) = 0.0d0
          do 10 i=n,0,-1
   10     fval(j,1) = fval(j,1) + fvbf(i)*w(i,1)
          if ( md.eq.0 )  go to 30
          fval(j,2) = 0.0d0
          do 20 i=n,0,-1
   20     fval(j,2) = fval(j,2) + fvbf(i)*w(i,2)
   30 continue
      return
      end
      subroutine rrgset( m, p, rt, ierr )
      integer            m, ierr
      double precision   p(m), rt(m,5)
c
c   sets up the temperature and stability approximations
c
      integer          j, k, ndeg, ngauss
      double precision dp, fg, kappa, pb, pg, pibyng, rgamma
      double precision sg, tk, tkm1, tkm2, wt
      double precision csval
      parameter        ( kappa = 0.2859d0 )
      common           /rrgcom/ ndeg, ngauss
      save             /rrgcom/
      data  ndeg, ngauss / 0, 0 /
c
c   set up the cubic spline interpolate of  rt
c
      ierr = 0
      call csset( m, p, rt, 3, 3, ierr )
      if ( ierr.ne.0 )  return
c
c   set up the Chebyshev polynomial fit to  1/(r*gamma)
c
      if ( ndeg.lt.0 )  then
          ndeg = abs( ndeg )
      else
          ndeg = m - 1
      end if
      if ( ngauss.lt.0 )  then
          ngauss = abs( ngauss )
      else
          ngauss = 2*ndeg
      end if
      pb = p(m)
      dp = pb - p(1)
      do 10 k=0,ndeg
   10 rt(k+1,5) = 0.0d0
      if ( ngauss.le.0 )  stop
      pibyng = acos( -1.0d0 )/ngauss
      do 30 j=1,ngauss
          sg = cos( (j - 0.5d0)*pibyng )
          pg = pb - dp*(sg + 1.0d0)/2.0d0
          rgamma = kappa*csval( pg, 0, m, p, rt, j )
     2             -  pg*csval( pg, 1, m, p, rt, 2 )
          if ( rgamma.le.0.0d0 )  then
              ierr = 5
              print 9000, pg, rgamma
 9000         format('0error in  rrgset--stability not positive ',
     2               'at  p =',1pe9.2,':  r*gamma =',e9.2/)
              rgamma = 1.0d0
          end if
          fg = 1.0d0/rgamma
          tkm2 = 1.0d0
          tkm1 = sg
          rt(1,5) = rt(1,5) + fg
          rt(2,5) = rt(2,5) + fg*tkm1
          do 20 k=2,ndeg
              tk = 2*sg*tkm1 - tkm2
              rt(k+1,5) = rt(k+1,5) + fg*tk
              tkm2 = tkm1
              tkm1 = tk
   20     continue
   30 continue
      if ( ierr.ne.0 )  return
      rt(1,5) = rt(1,5)/ngauss
      wt = 2.0d0/ngauss
      do 40 k=1,ndeg
   40 rt(k+1,5) = wt*rt(k+1,5)
      return
      end
      function rrgval( pval, con, m, rt, ierr )
      integer          m, ierr
      double precision rrgval, pval, con(3), rt(m,5)
c
c   evaluates the Chebyshev polynomial fit to  1/(r*gamma)  at  pval
c
      integer          k, ndeg, ngauss
      double precision dp, pb, s, tk, tkm1, tkm2
      common  /rrgcom/ ndeg, ngauss
      save    /rrgcom/
c
      ierr = 0
      pb = con(2)
      dp = con(3)
      s = 2*(pb - pval)/dp - 1.0d0
      rrgval = rt(1,5) + rt(2,5)*s
      tkm2 = 1.0d0
      tkm1 = s
      do 10 k=2,ndeg
          tk = 2*s*tkm1 - tkm2
          rrgval = rrgval + rt(k+1,5)*tk
          tkm2 = tkm1
          tkm1 = tk
   10 continue
      if ( rrgval.gt.0.0d0 )  return
      ierr = 5
      print 9000, pval, rrgval
 9000 format('0error in  rrgval--stability not positive ',
     2       'at  p =',1pe9.2,':  1/(r*gamma) =',e9.2/)
      return
      end
      subroutine funset( m, p, f, in, con, ierr )
      integer            m, in, ierr
      double precision   p(m), f(m,4), con(4)
c
c   sets up a continuous approximation to a function given discretely
c
      integer          i
      double precision a0, a1, a2, a3, beta, d0, d1, d2, d3
      double precision pb, pt, slope1, slopen
      double precision csval
c
c   preliminaries
c
      ierr = 0
      pt   = con(1)
      pb   = con(2)
      beta = con(4)
      if ( p(1).ne.pt .or. p(m).ne.pb )  ierr = 3
      if ( ierr.ne.0 )  return
      if ( in.le.0 )  go to 30
c
c   integrate the function with respect to  -log(p)
c
      call csset( m, p, f, 3, 3, ierr )
      if ( ierr.ne.0 )  return
      slope1 = -f(1,1)/pt
      slopen = -f(m,1)/pb
      f(m,1) = 0.0d0
      do 10 i=m-1,1,-1
          a0 = f(i,1) 
     &       - p(i)*(f(i,2)-p(i)*(f(i,3)/2.0d0-p(i)*f(i,4)/6.0d0))
          a1 = f(i,2) - p(i)*(f(i,3) - p(i)*f(i,4)/2.0d0)
          a2 = (f(i,3) - p(i)*f(i,4))/4.0d0
          a3 = f(i,4)/18.0d0
          d0 = log( p(i+1)/p(i) )
          d1 = p(i+1)    - p(i)
          d2 = p(i+1)**2 - p(i)**2
          d3 = p(i+1)**3 - p(i)**3
          f(i,1) = f(i+1,1) + a0*d0 + a1*d1 + a2*d2 + a3*d3
   10 continue
c
c   set up the interpolate and add in the constant of integration
c
      f(1,2) = slope1
      f(m,2) = slopen
      call csset( m, p, f, 1, 1, ierr )
      a0 = -csval( pb, 0, m, p, f, 0 )
      do 20 i=1,m
   20 f(i,1) = f(i,1) + a0
      return
c
c   set up the cubic spline interpolate without integrating
c
   30 f(1,2) = 0.0d0
      f(m,2) = -beta*f(m,1)/pb
      call csset( m, p, f, 1, 1, ierr )
      return
      end
      function funval( pval, md, m, p, f )
      integer          md, m
      double precision funval, pval, p(m), f(m,4), csval
c
      funval = csval( pval, md, m, p, f, 0 )
      return
      end
      subroutine vbfval( pval, md, con, n, chi )
      integer            md, n
      double precision  pval, con(4), chi(0:n,0:1)
c
c   evaluates vertical basis functions:  Chebyshev polynomials w/b.c.'s
c
      integer          j
      double precision a, b, beta, dp, fact, pb, pt, rm1j
      double precision s, tj, tjp1, tjp2
c
      pt   = con(1)
      pb   = con(2)
      dp   = con(3)
      beta = con(4)
      s  = 2*(pb - pval)/dp - 1.0d0
      if ( beta.eq.0.0d0 )  stop
      fact = 2*pb/(beta*dp)
      tj   = 1.0d0
      tjp1 = s
      rm1j = 1.0d0
      do 30 j=0,n
          tjp2 = 2*s*tjp1 - tj
          b = (j+2)**2
          a = rm1j + b*(1.0d0 + fact*(1.0d0 + rm1j))
          chi(j,0) = tjp2 - a - b*s
          if ( md.le.0 )  go to 20
          if ( abs(s).lt.1.0d0 .and. pt.lt.pval .and. pval.lt.pb )  then
              chi(j,1) = -(2/dp)*((j+2)*(tjp1-s*tjp2)/(1.0d0 - s*s) - b)
          else if ( (s.le.-1.0d0.or.pval.ge.pb).and.mod(j,2).eq.0 ) then
              chi(j,1) = 4*b/dp
          else
              chi(j,1) = 0.0d0
          end if
   20     tj   = tjp1
          tjp1 = tjp2
          rm1j = -rm1j
   30 continue
      return
      end
      subroutine csset( n, x, f, ibc1, ibcn, ierr )
      integer           n, ibc1, ibcn, ierr
      double precision x(n), f(n,4)
c
c   Purpose
c
c       sets up a cubic spline to interpolate specified data
c
c   Arguments
c
c       input
c
c           n       number of data points  (n.ge.2)
c
c           x       abscissas (in increasing order)
c
c           f       on input,  f(i,1)  contains the ordinate
c                   corresponding to  x(i)  (i=1, ..., n).
c                   with some boundary conditions, additional
c                   information must be supplied in  f(1,2)
c                   and/or  f(n,2)  as described below.
c
c           ibc1,   integers specifying the boundary conditions
c           ibcn    to be applied at  x(1)  and  x(n), respectively.
c                   ibc1 (ibcn)  is interpreted as follows:
c                       value    condition
c                         1    first derivative at endpoint supplied
c                              in  f(1,2)  (f(n,2))
c                         2    second derivative at endpoint supplied
c                              in  f(1,2)  (f(n,2))
c                         3    'not-a-knot' condition:  continuity
c                              of the third derivative imposed at
c                              x(2)  (x(n-1)) --requires  n.ge.3
c                         4    slope at endpoint estimated by fitting
c                              a cubic polynomial to the first (last)
c                              four data points--requires  n.ge.4
c
c       output
c
c           f       contains the polynomial coefficients of the cubic
c                   spline.  the input values  f(i,1)  (i=1, ..., n)
c                   remain unchanged.
c
c           ierr    error flag  (zero on normal return)
c
c   Error Conditions
c
c       ierr = 1    n  too small for choice of boundary conditions
c       ierr = 2    abscissas not in increasing order
c
c   Use
c
c       First call  csset  to set up the cubic spline.  Then use
c       the function  csval  to evaluate the spline or any of its
c       derivatives as desired.
c
c   Method
c
c       A tridiagonal system is set up from the requirement of
c       continuity of the second derivative at the interior
c       breakpoints plus the boundary conditions.  The system is
c       solved by Gaussian elimination for the first derivatives
c       at the breakpoints, from which the coefficients of the
c       cubic spline as a piecewise polynomial are obtained.
c
c   History
c
c       Written by Scott R. Fulton (February, 1982)
c       Based on the routine  cubspl  by Carl DeBoor
c
c   Reference
c
c       DeBoor, Carl (1978):  A Practical Guide to Splines.
c               Springer-Verlag, New York, 392 pp.
c
c
      integer          i, jbc1, jbcn
      double precision dd(3), fact
c
c   argument checks
c
      ierr = 0
      jbc1 = ibc1
      jbcn = ibcn
      if ( ibc1.lt.1 .or. ibc1.gt.4 )  jbc1 = 3
      if ( ibcn.lt.1 .or. ibcn.gt.4 )  jbcn = 3
      if ( n.lt.max( 2 , jbc1 , jbcn ) )  then
          ierr = 1
          return
      end if
      do 10 i=2,n
          if ( x(i-1).ge.x(i) )  then
              ierr = 2
              return
          end if
   10 continue
c
c   compute  x  differences and first divided differences of  f
c
      do 20 i=2,n
          f(i,3) = x(i) - x(i-1)
          f(i,4) = (f(i,1) - f(i-1,1))/f(i,3)
   20 continue
c
c   set up first equation from the boundary condition at  x(1)
c
      if ( jbc1.eq.1 )  then
          f(1,4) = 1.0d0
          f(1,3) = 0.0d0
      else if ( jbc1.eq.2 )  then
          f(1,4) = 2.0d0
          f(1,3) = 1.0d0
          f(1,2) = 3.0d0*f(2,4) - f(2,3)*f(1,2)/2.0d0
      else if ( jbc1.eq.3 )  then
          f(1,4) = f(3,3)
          f(1,3) = f(2,3) + f(3,3)
          f(1,2) = ( f(3,3)*(f(2,3) + 2.0d0*f(1,3))*f(2,4)
     2             + f(2,3)*f(2,3)*f(3,4) )/f(1,3)
      else
          f(1,4) = 1.0d0
          f(1,3) = 0.0d0
          dd(2) = (f(3,4) - f(2,4))/(f(2,3) + f(3,3))
          dd(3) = (f(4,4) - f(3,4))/(f(3,3) + f(4,3))
          dd(3) = ( dd(3) - dd(2) )/(f(2,3) + f(3,3) + f(4,3))
          f(1,2) = f(2,4) - f(2,3)*(dd(2) - (f(2,3)+f(3,3))*dd(3))
      end if
c
c   set up equations for continuity of second derivative at interior
c   breakpoints  x  and do forward pass of gaussian elimination
c
      do 30 i=2,n-1
          fact = -f(i+1,3)/f(i-1,4)
          f(i,2) = fact*f(i-1,2)+3.0d0*(f(i+1,3)*f(i,4)+f(i,3)*f(i+1,4))
          f(i,4) = fact*f(i-1,3)+2.0d0*(f(i,3) + f(i+1,3))
   30 continue
c
c   set up last equation from the boundary condition at  x(n)
c   and complete forward pass of gaussian elimination
c
      if ( jbcn.eq.2 )  then
          f(n,2) = 3.0d0*f(n,4) + f(n,3)*f(n,2)/2.0d0
          fact   = -1.0d0/f(n-1,4)
          f(n,4) =  fact*f(n-1,3) + 2.0d0
          f(n,2) = (fact*f(n-1,2) + f(n,2))/f(n,4)
      else if ( jbcn.eq.3 )  then
          fact   =  f(n-1,3) + f(n,3)
          dd(1)  = (f(n-1,1) - f(n-2,1))/f(n-1,3)
          f(n,2) = ( f(n-1,3)*(f(n,3) + 2.0d0*fact)*f(n,4)
     2             + f(n,3)*f(n,3)*dd(1) )/fact
          fact   = -fact/f(n-1,4)
          f(n,4) = (fact + 1.0d0)*f(n-1,3)
          f(n,2) = (fact*f(n-1,2) + f(n,2))/f(n,4)
      else if ( jbcn.eq.4 )  then
          dd(3) = (f(n-2,1) - f(n-3,1))/f(n-2,3)
          dd(2) = (f(n-1,1) - f(n-2,1))/f(n-1,3)
          dd(3) = ( dd(2) - dd(3))/(f(n-2,3) + f(n-1,3))
          dd(2) = (f(n,4) - dd(2))/(f(n-1,3) + f(n  ,3))
          dd(3) = ( dd(2) - dd(3))/(f(n-2,3) + f(n-1,3) + f(n,3))
          f(n,2) = f(n,4) + f(n,3)*(dd(2) + (f(n-1,3)+f(n,3))*dd(3))
      end if
c
c   do backward pass of gaussian elimination
c
      do 40 i=n-1,1,-1
   40 f(i,2) = (f(i,2) - f(i,3)*f(i+1,2))/f(i,4)
c
c   compute polynomial coefficients of the cubic spline
c
      do 50 i=2,n
          fact  =  f(i,3)
          dd(1) = (f(i,1) - f(i-1,1))/fact
          dd(3) =  f(i-1,2) + f(i,2) - 2.0d0*dd(1)
          f(i-1,3) = 2.0d0*(dd(1) - f(i-1,2) - dd(3))/fact
          f(i-1,4) = 6.0d0*dd(3)/(fact*fact)
   50 continue
      return
      end
      function csval( xval, md, n, x, f, isw )
      integer          md, n, isw
      double precision csval, xval, x(n), f(n,4)
c
c
c   Purpose
c
c       evaluates the cubic spline set up by the subroutine  csset
c
c   Arguments
c
c       dimension   x(n), f(n,4)
c
c       input
c
c           xval    abscissa at which the spline is to be evaluated
c
c           md      integer specifying the derivative desired:
c                       md = 0    spline value
c                       md = 1    first  derivative
c                       md = 2    second derivative
c                       md = 3    third  derivative
c
c           n,x,f   unchanged from call to  csset
c
c           isw     integer switch to help locate the interval
c                   containing  xval:
c                       isw.gt.0    ascending search
c                       isw.eq.0    bisection search
c                       isw.lt.0    descending search
c                   if  isw = +1 (-1), the search starts with the
c                   first (last) interval.  thus when interpolating
c                   at an increasing (decreasing) sequence of
c                   points  xval = y(j)  (j=1, ..., k), set
c                   isw = +j (-j)  for maximum efficiency.
c
c       output
c
c           csval   returns the requested spline value or
c                   derivative.  extrapolation is used when  xval
c                   is outside the interval  (x(1),x(n)),  and the
c                   third derivative is taken to be continuous
c                   from the right at the breakpoints.
c
c   Use
c
c       First call  csset  to set up the cubic spline.  Then use
c       the function  csval  to evaluate the spline or any of its
c       derivatives as desired.
c
c   Method
c
c       An index  i  is determined as described above such that
c       x(i).le.xval.lt.x(i+1).  The desired value is then obtained
c       using the polynomial representation of the cubic spline
c       on that interval.
c
c   History
c
c       Written by Scott R. Fulton (February, 1982)
c       Based on the routine  ppvalu  by Carl DeBoor
c
c   Reference
c
c       DeBoor, Carl (1978):  A Practical Guide to Splines.
c               Springer-Verlag, New York, 392 pp.
c
c
      integer          i, ih, im, j, k
      double precision dx, fact
      save             i
      data             i / 1 /
c
c   return zero for fourth- and higher-order derivatives
c
      j = max( md+1 , 1 )
      if ( j.gt.4 )  then
          csval = 0.0d0
          return
      end if
c
c   find the  x  interval containing  xval
c
      if ( isw )  10, 30, 80
c
c   1. descending search
c
   10 i = max( i, 1 )
      i = min( i, n-1 )
      if ( isw.eq.-1 .or. xval.ge.x(i+1) )  i = n-1
      i = i+1
   20 i = i-1
      if ( x(i).gt.xval .and. i.gt.1 )  go to 20
      go to 100
c
c   2. bisection search
c
   30 i  = 1
      ih = n
   40 if ( ih-i.le.1 )  go to 100
      im = (i + ih)/2
      if ( xval - x(im) )  50, 60, 70
   50 ih = im
      go to 40
   60 i  = im
      ih = im
      go to 40
   70 i = im
      go to 40
c
c   3. ascending search
c
   80 i = max( i, 1 )
      i = min( i, n-1 )
      if ( isw.eq.+1 .or. x(i).gt.xval )  i = 1
      i = i-1
   90 i = i+1
      if ( xval.ge.x(i+1) .and. i+1.lt.n )  go to 90
c
c   evaluate the cubic spline
c
  100 csval = f(i,4)
      if ( j.eq.4 )  return
      fact = 4-j
      dx = xval - x(i)
      do 200 k=3,j,-1
          csval = f(i,k) + dx*csval/fact
          fact  = fact - 1.0d0
  200 continue
      return
      end
