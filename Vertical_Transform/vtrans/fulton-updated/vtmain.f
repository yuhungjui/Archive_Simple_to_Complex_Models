      program vtrans
c
c   Revision History
c
c       05/21/84:  original version (used for MWR paper 1985)
c       02/22/99:  updated for use under UNIX
c
c   Disclaimer
c
c       This driver program for  vtrans  is not claimed to be very
c       understandable, well tested, or even useful--use at your risk.
c
c   Purpose
c
c       This program uses  vtrans  routines to solve the vertical
c       problem and compute vertical transforms for temperature
c       profiles and functions specified at discrete pressure levels.
c
c   Input
c
c       The input to this program consists of lines which specify the
c       data to be used and what to do with it.  Directive lines are
c       recognized by their first three characters (in columns 1-3)
c       and may contain various other parameters.  Data lines begin
c       with at least one blank (or '?' if the line is to be ignored).
c       Each line read is echoed on the output, and unrecognized lines
c       are treated as comments and ignored.
c
c       Directive   Parameters
c       ---------   ----------
c
c       temperature Specifies a basic state temperature profile from
c       ---         data.  Columns 21-80 become a label for the basic
c                   state.  This line must be followed by data lines
c                   giving the temperature data  (p(i), t(i)), one data
c                   point per line in 2f10 format, with pressure in any
c                   units (be consistent throughout) and temperature in
c                   either celsius or kelvin.  The minimum and maximum
c                   pressures input become the top and bottom pressures.
c
c       gamma       Specifies an analytical basic state of constant
c       ---         static stability  gamma = kappa*t - p*dt/dp.
c                   parameter   meaning     columns  format   default
c                     gamma   (in kelvin)    21-30     f       23.79
c                     tbot    surface temp.  31-40     f       29.38
c                     ptop    top    press.  41-50     f       100.0
c                     pbot    bottom press.  51-60     f      1010.0
c                     dp      p increment    61-70     f        ptop
c                   Here  dp  is used in generating discrete temperature
c                   data from this analytically specified basic state.
c
c   **** note:  only one basic state may be specified at a time ****
c
c       top         Redefines the top pressure for temp. and function
c       ---         parameter   meaning     columns  format   default
c                      ptop   new top press  21-30     f      old ptop
c                   Note:  can only increase top pressure (not decrease)
c
c       nsolve      Solves the vertical structure problem numerically
c       ---         parameter   meaning     columns  format   default
c                       n     truncation     16-20     i         16
c                      ndeg   poly. degree   21-30     i        nt-1
c                     ngauss  gauss points   31-40     i       2*ndeg
c                   where  nt  is the number of temperature data points.
c                   (ndeg  and  ngauss  refer to the Chebyshev fit to
c                   the stability profile--see also  stability  below)
c
c       asolve      Solves the vertical structure problem analytically
c       ---         parameter   meaning     columns  format   default
c                      n      truncation     16-20     i         32
c                      tol    c tolerance    21-30     f       10e-10
c
c   **** Note:  Only one solution method may be used at a time ****
c
c       rnsolve     Solves reference problem (for errors) numerically
c       ---         parameter   meaning     columns  format   default
c                      nref   ref. trunc.    16-20     i         96
c                      ndeg   poly. degree   21-30     i        nt-1
c                     ngauss  gauss points   31-40     i       2*ndeg
c                   where  nt  is the number of temperature data points.
c                   (ndeg  and  ngauss  refer to the cyebyshev fit to
c                   the stability profile--see also  stability  below)
c
c       rasolve     Solves reference problem (for errors) analytically
c       ---         parameter   meaning     columns  format   default
c                      nref   ref. trunc.    16-20     i         96
c                      tol    c tolerance    21-30     f    (rel)10e-15
c
c   **** Note:  Only one solution method may be used at a time ****
c
c       output      Specifies output pressure levels as follows:
c       ---         parameter   meaning     columns  format   default
c                     pmin    min pressure   21-30     f        ptop
c                     pmax    max pressure   31-40     f        pbot
c                     dp      increment      41-50     f        ptop
c                   if pmin=-1.0 the temperature data pressures are used
c                   if pmin=-2.0 the function    data pressures are used
c                   if pmin=-3.0 this line must be followed by data
c                   lines giving the desired output pressure levels,
c                   one value per line, in f10 format.
c
c       solution    Prints solution of the vertical structure problem
c       ---         (phase speeds and vertical structure functions)
c                   parameter   meaning     columns  format   default
c                      md      derivative       20     i         0
c                      k1      start mode    21-30     i         0
c                      k2      stop  mode    31-40     i         n
c                      k3      increment     41-50     i         1
c                   where  md = 0  means print  values  of function
c                          md = 1  means print  d/dp    of function
c                          md = 2  means print -d/dlogp of function
c
c       stability   Tests the fit of a polynomial stability profile
c       ---         parameter   meaning     columns  format   default
c                      ndeg   poly. degree   21-30     i        nt-1
c                     ngauss  gauss points   31-40     i       2*ndeg
c                   where  nt  is the number of temperature data points.
c
c       errors      Prints errors in the (numerical) solution of the
c       ---         vertical structure problem, compared to the refernce
c                   solution (see  vrn  and  vra  above).  parameters:
c                   k1, k2, k3  exactly as for  solution  (see above).
c
c       reconstruct Reconstructs the vertical structure functions from
c       ---         their values at discrete levels and prints the mode-
c                   by-mode energetics of reconstruction.  parameters:
c                   k1, k2, k3  exactly as for  solution  (see above).
c                   must be followed by data lines giving the desired
c                   levels to be used (one per line, f format).  These
c                   lines must be followed by a line giving a character
c                   label for the levels (beginning in column 1).
c
c       function    Specifies a function to be transformed from data.
c       ---         Columns 21-80 become a label for the function.
c                   This line must be followed by data lines giving the
c                   function data  (p(i), f(i)), one data point per line
c                   in 2f10 format (use same pressure units as before).
c
c       evaluate    Specifies the index of a vertical structure function
c       ---         to be evaluated at the output pressure levels for
c                   for subsequent projection and reconstruction.
c                   parameter   meaning     columns  format   default
c                       k    v.s.f. index    21-30     i         0
c
c   **** Note:  Only one function may be defined at a time ****
c
c       project     Projects the function onto the vertical basis.
c       ---         parameter   meaning     columns  format   default
c                       ip   inner product   21-30     i     1 (cheby.)
c                       ng   num. of pnts.   31-40     i      2*(nf+n)
c                   where  nf  is the number of function data points
c                   (see vtrans documentation for meanings of ip and ng)
c                   use  ip = 3  for collocation projection of  f  data.
c
c       integrate   Projects the integral of the function onto the
c       ---         vertical basis (integral set equal to zero at pbot).
c                   Parameters:  same as for  project  (see above)
c                   To multiply the function by  kappa  or  rgas  before
c                   integrating, put 'kappa' or 'rgas ' in col. 16-20.
c
c       physical    Prints physical space profiles of temperature,
c       ---         stability (sigma and gamma) and the function.
c                   Parameters:  md, k1, k2, k3  (as for  solution
c                   above).  Here  k1, k2, k3  specify the vertical
c                   modes used in reconstructing the function--the
c                   default is to use all modes.
c
c       spectral    Prints phase speeds, spectral coefficients of the
c       ---         function and mode-by-mode energetics.  Parameters:
c                   k1, k2, k3  exactly as for  solution  (see above).
c
c       end         Ends program execution
c       ---
c
c   Language    Fortran 77
c
c   History     Written by Scott R. Fulton in April 1984
c
c
      integer             nmax        , ntmax      , nfmax      
      parameter         ( nmax   =  32, ntmax  = 47, nfmax = 47)
      integer             npmax       , nrmax      , ngmax
      parameter         ( npmax  =  92, nrmax  = 96, ngmax = 99 )
      integer             nipmax      , nppage     , npp
      parameter         ( nipmax = 158, nppage = 17, npp   = 21 )
      integer           i, ierr, in, ip, ipage, j, k
      integer           m(npp), me(0:nmax,npp), md, mode, m1, m2
      integer           n, ncol, ndeg, ndl, nf, ng, ngauss
      integer           ngl, nip, np, npages, nr, nt, n1, n2, n3
      double precision  cerr, cerrl, con(4), dp, efsum, ensum
      double precision  err2, erra, erri, errr, factor, fp, fs
      double precision  gamma, gc, gs, kappa, pbot, perrl
      double precision  pmax, pmin, psib, psirb, ptop, relerr
      double precision  rgamma, rgas, rtb
      double precision  sigma, tbot, tc, time, tinf, tol, ts, tzero
      double precision  c (0: nmax), psihat(0: nmax), fvbf(0:nmax)
      double precision  cr(0:nrmax), psihar(0:nrmax), fhat(0:nmax)
      double precision  ef(0:nmax), en(0:nmax)
      double precision  tfor(0: nmax,0: nmax), tinv(0: nmax,0: nmax)
      double precision  rfor(0:nrmax,0:nrmax), rinv(0:nrmax,0:nrmax)
      double precision  work(0:nrmax,0:nrmax), prom(0: nmax,0: nmax)
      double precision  pt(ntmax), rt(ntmax,9), pf(nfmax)
      double precision  fdata(nfmax), f(nfmax,4)
      double precision  w(0:nrmax,2), psierr(0:nrmax)
      double precision  pip(nipmax), wip(nipmax), vsf(npmax,nppage)
      double precision  p(npmax), fval(npmax,2), finteg(npmax)
      double precision  pg(ngmax), wg(ngmax), psi(ngmax), psir(ngmax)
      double precision  pnew, ttop, ftop
      double precision  csval, adquad, rrgval, funval
      real              second
      character         line*80, first, first3*3, label(0:2)*8, labelf*60
      character         labeln*18, labelp(3)*11, labelr*10, labelt*60
      logical           anavsp, anavsr
      external          fgrand
      common            /fgradc/ con, rt, nt
      common   /rrgcom/ ndeg, ngauss
      common   /gamcom/ rgamma
      parameter         ( rgas = 287.0d0, kappa = 0.2859d0 )
      double precision  zero, one
      parameter         ( zero = 0.0d0, one = 1.0d0 )
      data       label  / '  values', '   d/dp ', '-d/dlogp' /
      data       labelp / 'chebyshev', 'vertical', 'collocation' /
c
c   initializations
c
      n  = 0
      nt = 0
      nf = 0
      ip = 0
      np = 0
      ptop = zero
      pbot = zero
c
c   read the next data line and decide what to do
c
   10 read '(a)', line
      print '('' input--'',a)', line
   20 first3 = line(1:3)
      if ( first3.eq.'end' )  stop
      if ( first3.eq.'tem' )  go to 100
      if ( first3.eq.'gam' )  go to 150
      if ( first3.eq.'top' )  go to 170
      if ( first3.eq.'nso' )  go to 200
      if ( first3.eq.'aso' )  go to 210
      if ( first3.eq.'rns' )  go to 220
      if ( first3.eq.'ras' )  go to 230
      if ( first3.eq.'out' )  go to 250
      if ( first3.eq.'sol' )  go to 300
      if ( first3.eq.'sta' )  go to 400
      if ( first3.eq.'err' )  go to 500
      if ( first3.eq.'rec' )  go to 540
      if ( first3.eq.'fun' )  go to 600
      if ( first3.eq.'eva' )  go to 650
      if ( first3.eq.'pro' )  go to 700
      if ( first3.eq.'int' )  go to 710
      if ( first3.eq.'phy' )  go to 800
      if ( first3.eq.'spe' )  go to 900
      go to 10
c
c   input the basic state temperature profile
c
  100 labelt = line(21:80)
      nt = 0
  110 read '(a)', line
      print '('' input--'',a)', line
c         first = line(1:1)
          read ( line, '(a1)' )  first
          if ( first.eq.'?' )  then
              print *,'(this data point ignored)'
              go to 110
          end if
          if ( first.ne.' ' )  go to 120
          if ( nt.ge.ntmax )  then
              print *, 'too many data points'
              stop
          end if
          nt = nt + 1
          read ( line, 1100 )  pt(nt), rt(nt,1)
 1100     format(2f10.0)
      go to 110
  120 call sorter( nt, pt, rt )
      ptop = pt( 1)
      pbot = pt(nt)
      if ( rt(nt,1).lt.100.0d0 )  then
          tzero = 273.15d0
          print *, 'temperatures are assumed to be in celsius'
      else
          tzero = zero
          print *, 'temperatures are assumed to be in kelvin'
      end if
      do 130 i=1,nt
  130 rt(i,1) = rgas*(rt(i,1) + tzero)
      rgamma = zero
      rtb    = rt(nt,1)
      go to 20
c
c   input the constant gamma basic state
c
  150 read ( line, '(20x,5f10.0)' )  gamma, tbot, ptop, pbot, dp
      if ( gamma.le.zero )  gamma =  23.79d0
      if (  tbot.le.zero )   tbot =  29.38d0
      if (  ptop.le.zero )   ptop =  100.0d0
      if (  pbot.le.zero )   pbot = 1010.0d0
      if (    dp.le.zero )     dp =   ptop
      call plevel( ptop, pbot, dp, pbot, nt, pt )
      if ( nt.gt.ntmax )  then
          print *, 'too many pressure levels'
          stop
      end if
      ptop = pt( 1)
      pbot = pt(nt)
      write ( labelt, 5150 )  gamma, tbot, ptop, pbot
 5150 format('gamma =',f6.2,'  tb =',f7.2,
     2        '  pt =',f8.1,'  pb =',f9.1)
      if ( tbot.lt.100.0d0 )  then
          tzero = 273.15
      else
          tzero = zero
      end if
      tbot = tbot + tzero
      tinf = gamma/kappa
      do 160 i=1,nt
  160 rt(i,1) = rgas*(tinf + (tbot - tinf)*(pt(i)/pbot)**kappa)
      call csset( nt, pt, rt, 3, 3, ierr )
      if ( ierr.ne.0 )  print 9000, 'csset', ierr
 9000 format('0********** error in ',a,':  ierr =',i4,' **********')
      rgamma = rgas*gamma
      rtb    = rt(nt,1)
      go to 10
c
c   redefine the top pressure
c
  170 read ( line, '(20x,f10.0)' )  pnew
      if ( rgamma.gt.zero )  then
          print *, 'use  gamma  line instead'
          go to 10
      else if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      end if
      if ( pnew.le.zero )  pnew = ptop
      if ( pnew.lt.ptop )  then
          print *, 'cannot decrease top pressure'
          go to 10
      end if
      ptop = pnew
      j = 0
  172 j = j + 1
      if ( pt(j+1).le.ptop )  go to 172
      ttop = csval( ptop, 0, nt, pt, rt, 0 )
      print 2170, 'new temperature data'
 2170 format('0********** ',a,' **********'/)
      pt(1)   = ptop
      rt(1,1) = ttop
      print '(1x,f10.1,f10.2)', pt(1), rt(1,1)/rgas - tzero
      nt = nt - j + 1
      do 174 i=2,nt
          pt(i)   = pt(i+j-1)
          rt(i,1) = rt(i+j-1,1)
          print '(1x,f10.1,f10.2)', pt(i), rt(i,1)/rgas - tzero
  174 continue
      if ( nf.le.0 )  go to 10
      if ( ip.le.0 )  then
          print *, 'must project or integrate the function first'
          go to 10
      end if
      j = 0
  176 j = j + 1
      if ( pf(j+1).le.ptop )  go to 176
      ftop = funval( ptop, in, nf, pf, f )
      if ( in.gt.0 )  ftop = -ptop*ftop
      print 2170, 'new function data'
      pf(1)    = ptop
      fdata(1) = ftop
      print '(1x,f10.1,1pe17.5)', pf(1), fdata(1)
      nf = nf - j + 1
      do 178 i=2,nf
          pf(i)    = pf(i+j-1)
          fdata(i) = fdata(i+j-1)
          print '(1x,f10.1,1pe17.5)', pf(i), fdata(i)
  178 continue
      go to 10
c
c   solve the vertical structure problem (numerically)
c
  200 read ( line, '(15x,i5,2i10)' )  n, ndeg, ngauss
      if ( 1.gt.n .or. n.gt.nmax )  n = 16
      if ( ndeg+1.gt.5*ntmax )  then
          print *, 'the polynomial degree requested is too large'
          stop
      end if
      if (   ndeg.ne.0 )    ndeg = -ndeg
      if ( ngauss.ne.0 )  ngauss = -ngauss
c     time = second( work )
      time = second( )
      call vtset( nt, pt, rt, n, w, work, con, c, tfor, tinv, ierr )
c     time = second( work ) - time
      time = second( ) - time
      if ( ierr.ne.0 )  print 9000, 'vtset', ierr
      print 2200, time, 'numerical', n, 'vertical structure'
 2200 format('0cpu time:',f8.3,' seconds for ',a,
     2       ' solution of  n =',i3,2x,a,' problem')
      print 2210, ndeg, ngauss
 2210 format(' (stability fit:  polynomial of degree',i4,' using',i4,
     2       ' gauss points)')
      anavsp = .false.
      go to 10
c
c   solve the vertical structure problem (analytically)
c
  210 read ( line, '(15x,i5,f10.0)' )  n, tol
      if ( rgamma.le.zero )  then
          print *, 'cannot do it--basic state specified by data'
          go to 10
      end if
      if ( 1.gt.n .or. n.gt.nmax )  n = nmax
      if ( tol.eq.zero )  tol = -1.0d-10
c     time = second( work )
      time = second( )
      call cgeval( rgamma, rtb, pbot, ptop, tol, n, c, psihat, ierr )
c     time = second( work ) - time
      time = second( ) - time
      if ( ierr.ne.0 )  print 9000, 'cgeval', ierr
      print 2200, time, 'analytical', n, 'vertical structure'
      anavsp = .true.
      go to 10
c
c   solve the reference vertical structure problem (numerically)
c
  220 read ( line, '(15x,i5,2i10)' )  nr, ndeg, ngauss
      if ( 1.gt.nr .or. nr.gt.nrmax )  nr = nrmax
      if ( ndeg+1.gt.5*ntmax )  then
          print *, 'the polynomial degree requested is too large'
          stop
      end if
      if (   ndeg.ne.0 )    ndeg = -ndeg
      if ( ngauss.ne.0 )  ngauss = -ngauss
c     time = second( work )
      time = second( )
      call vtset( nt, pt, rt, nr, w, work, con, cr, rfor, rinv, ierr )
c     time = second( work ) - time
      time = second( ) - time
      if ( ierr.ne.0 )  print 9000, 'vtset', ierr
      ng = nr + 3
      call gaussl( ng, ptop, pbot, pg, wg, ierr )
      if ( ierr.ne.0 )  print 9000, 'gaussl', ierr
      print 2200, time, 'numerical', nr, 'reference'
      print 2210, ndeg, ngauss
      anavsr = .false.
      go to 10
c
c   solve the reference vertical structure problem (analytically)
c
  230 read ( line, '(15x,i5,f10.0)' )  nr, tol
      if ( rgamma.le.zero )  then
          print *, 'cannot do it--basic state specified by data'
          go to 10
      end if
      if ( 1.gt.nr .or. nr.gt.nrmax )  nr = nrmax
c     time = second( work )
      time = second( )
      call cgeval( rgamma, rtb, pbot, ptop, tol, nr, cr, psihar, ierr )
c     time = second( work ) - time
      time = second( ) - time
      if ( ierr.ne.0 )  print 9000, 'cgeval', ierr
      print 2200, time, 'analytical', nr, 'reference'
      ng = ngmax
      call gaussl( ng, ptop, pbot, pg, wg, ierr )
      if ( ierr.ne.0 )  print 9000, 'gaussl', ierr
      anavsr = .true.
      go to 10
c
c   generate the output pressure levels
c
  250 read ( line, '(20x,3f10.0)' )  pmin, pmax, dp
      if ( ptop.le.zero .or. pbot.le.zero )  then
          print *, 'must specify basic state first'
          go to 10
      end if
      if ( pmin.eq.-one )  then
          np = nt
          do 251 i=1,np
  251     p(i) = pt(i)
      else if ( pmin.eq.-2.0d0 )  then
          if ( nf.le.0 )  then
              print *, 'must specify basic state first'
              go to 10
          end if
          np = nf
          do 252 i=1,np
  252     p(i) = pf(i)
      else if ( pmin.eq.-3.0d0 )  then
          np = 0
  253     read '(a)', line
          print '('' input--'',a)', line
c             first = line(1:1)
              read ( line, '(a1)' )  first
              if ( first.ne.' ' )  go to 254
              np = np + 1
              read ( line, 1100 )  p(np)
          go to 253
  254     if ( np.gt.npmax )  print *, 'too many output pressure levels'
          go to 20
      else
          if ( pmin.le.zero )  pmin = ptop
          if ( pmax.le.zero )  pmax = pbot
          if (   dp.le.zero )    dp = ptop
          call plevel( pmin, pmax, dp, pbot, np, p )
          if ( np.gt.npmax )  print *, 'too many output pressure levels'
      end if
      go to 10
c
c   print the solution of the vertical structure problem
c
  300 read ( line, '(19x,i1,3i10)' )  k, n1, n2, n3
      if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      else if ( np.le.0 )  then
          print *, 'must specify the output pressure levels first'
          go to 10
      end if
      k  = max( k, 0 )
      k  = min( k, 2 )
      md = min( k, 1 )
      if ( 0.gt.n1 .or. n1.gt.n )  n1 = 0
      if ( 1.gt.n2 .or. n2.gt.n )  n2 = n
      if ( n3.eq.0 )  n3 = 1
      ncol = (n2 - n1)/n3 + 1
      npages = (ncol - 1)/nppage + 1
      labeln = '                  '
      do 350 ipage=1,npages
          m1 = n1 + (ipage - 1)*nppage*n3
          m2 = min( n2, m1 + (nppage - 1)*n3 )
          ncol = (m2 - m1)/n3 + 1
          do 330 j=1,ncol
              mode = m1 + (j-1)*n3
              if ( anavsp )  then
                  call  cgefun( mode, np, p, rgamma, rtb, pbot, ptop,
     2                          c, psihat, vsf(1,j) )
              else
                  i = mode*(n + 1)
                  call veval( tinv(i,0), md, np, p, con, n, w, fval )
                  do 310 i=1,np
  310             vsf(i,j) = fval(md*np+i,1)
              end if
              if ( k.eq.2 )  then
                  do 320 i=1,np
  320             vsf(i,j) = -p(i)*vsf(i,j)
              end if
              m(j) = mode
  330     continue
          if ( npages.gt.1 )  write ( labeln, 5300 ) ipage, npages
          print 2300, labeln, labelt
          if ( anavsp )  print 2310
          print 2320, label(k), (m(j), j=1,ncol )
          print '(10x)'
          do 340 i=1,np
  340     print '(1x,f8.1,1x,17f7.3)', p(i), (vsf(i,j), j=1,ncol)
          print '(''0  c(k) = '',17f7.3)', (c(m(j)), j=1,ncol)
  350 continue
 5300 format('     (page',i2,' of',i2,')')
 2300 format('1 '//'0solution of the vertical structure problem',a,
     3           //'0basic state:  ',a60)
 2310 format(15x,'(the vertical structure problem ',
     2           'was solved analytically)')
 2320 format(//21x,a,' of the vertical structure functions',
     4       //'0     p',3x,17(2x,'k=',i2:1x))
      print 2330, n
 2330 format('0(phase speeds  c  in meters per second, ',
     2       'spectral truncation  n =',i4,')')
      go to 10
c
c   test the stability fit
c
  400 ndl = ndeg
      ngl = ngauss
      read ( line, '(20x,2i10)' )  ndeg, ngauss
      if ( np.le.0 )  then
          print *, 'must specify the output pressure levels first'
          go to 10
      else if ( ndeg+1.gt.5*ntmax )  then
          print *, 'the polynomial degree requested is too large'
          go to 10
      end if
      if (   ndeg.ne.0 )    ndeg = -ndeg
      if ( ngauss.ne.0 )  ngauss = -ngauss
      call rrgset( nt, pt, rt, ierr )
      if ( ierr.ne.0 )  print 9000, 'rrgset', ierr
      con(2) = pbot
      print 2400, labelt, ndeg, ngauss
 2400 format('1 '//'0test of stability profile fit',
     2           //'0basic state:  ',a60,
     3           //'0description:  chebyshev least-squares fit to  ',
     4                            '1/gamma  using',
     5            /'               polynomial of degree',i4,' with',i4,
     6                            ' gauss points',
     7           //15x,'temperature',13x,'error',16x,'gamma',
     8            /6x,'p',5x,'spline',2x,'chebyshev',3x,'absolute',
     9             3x,'relative',3x,'spline',2x,'chebyshev'/)
      err2 = zero
      erri = zero
      relerr = 1.0d-06
      finteg(np) = adquad( p(np), pbot, fgrand, relerr, j, w, ierr )
      if ( ierr.ne.0 )  print *, 'integral from ', p(i), 'to', pbot
      if ( ierr.eq.0 )  then
      do 410 i=np-1,1,-1
          finteg(i) = finteg(i+1)
     2              + adquad( p(i), p(i+1), fgrand, relerr, j, w, ierr )
          if ( ierr.ne.0 ) print *, 'integral from ', p(i), 'to', p(i+1)
  410 continue
      end if
      do 420 i=1,np
          ts = csval( p(i), 0, nt, pt, rt, 0 )
          gs = kappa*ts - p(i)*csval( p(i), 1, nt, pt, rt, 2 )
          ts = ts/rgas - tzero
          gs = gs/rgas
          gc = rrgval( p(i), con, nt, rt, ierr )
          if ( ierr.ne.0 )  print 9000, 'rrgval', ierr
          gc = one/(rgas*gc)
          tc = (p(i)/con(2))**kappa*(rtb + finteg(i))
          tc = tc/rgas - tzero
          erra = tc - ts
          errr = erra/(ts + tzero)
          err2 = err2 + erra**2
          erri = max( erri, abs( erra ) )
          print 2410, p(i), ts, tc, erra, errr, gs, gc
 2410     format(1x,f8.1,f10.2,f9.2,f12.5,1pe11.2,0p,f10.3,f9.3)
  420 continue
      err2 = sqrt( err2/np )
      print 2420, err2, erri
 2420 format('0',10x,'rms absolute error:',f10.5,
     2          /11x,'max absolute error:',f10.5)
      if ( (ndeg.ne.ndl .or. ngauss.ne.ngl) .and. ndl.ne.0 )  then
          ndeg   = -ndl
          ngauss = -ngl
          call rrgset( nt, pt, rt, ierr )
          if ( ierr.ne.0 )  print 9000, 'rrgset', ierr
          con(2) = pbot
      end if
      go to 10
c
c   compute and print errors in the vertical structure problem solution
c
  500 read ( line, '(20x,3i10)' )  n1, n2, n3
      if ( anavsp )  then
          print *, 'no errors--the problem was solved analytically'
          go to 10
      else if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      else if ( nr.le.0 )  then
          print *, 'must solve the reference problem first'
          go to 10
      end if
      if ( 0.gt.n1 .or. n1.gt.n )  n1 = 0
      if ( 1.gt.n2 .or. n2.gt.n )  n2 = n
      if ( n3.eq.0 )  n3 = 1
      do 520 k=n1,n2,n3
          i = k*(n + 1)
          call veval( tinv(i,0), 0, ng,  pg, con, n, w, psi )
          call veval( tinv(i,0), 0, 1, pbot, con, n, w, psib )
          if ( anavsr )  then
              call cgefun(k,ng, pg,rgamma,rtb,pbot,ptop,cr,psihar,psir )
              call cgefun(k,1,pbot,rgamma,rtb,pbot,ptop,cr,psihar,psirb)
          else
              i = k*(nr + 1)
              call veval( rinv(i,0), 0, ng,  pg, con, nr, w, psir  )
              call veval( rinv(i,0), 0, 1, pbot, con, nr, w, psirb )
          end if
          j = sign( one, psib/psirb )
          psierr(k) = zero
          do 510 i=1,ng
  510     psierr(k) = psierr(k) + wg(i)*(psir(i) - j*psi(i))**2
          psierr(k) = sqrt( psierr(k)/(pbot - ptop) )
  520 continue
      labelr = 'analytical'
      if ( .not.anavsr )  write ( labelr, '(''   n='',i3,2x)' )  nr
      print 2500, labelt, n, labelr
 2500 format('1 '//'0errors in the numerical solution ',
     2             'of the vertical structure problem',
     3           //'0basic state:  ',a60,
     4           //9x,' ----------- eigenvalue  c(k) ------------ ',
     5             1x,'eigenfunction psi(k)'/9x,' spectral truncation',
     6             2x,' relative difference',2x,' norm of difference ',
     7            /4x,'k',6x,'n=',i3,3x,a10,6x,'value',4x,'log10',
     8             8x,'value',4x,'log10'/)
      do 530 k=n1,n2,n3
          cerr = (c(k) - cr(k))/cr(k)
          if ( cerr.ne.zero )  cerrl = log10( abs( cerr ) )
          if ( cerr.eq.zero )  cerrl = zero
          if ( psierr(k).ne.zero )  perrl = log10( abs( psierr(k) ) )
          if ( psierr(k).eq.zero )  perrl = zero
          print '(1x,i4,2x,2f11.5,1pe13.2, 0pf7.2, 1pe15.2,  0pf7.2)',
     2               k, c(k), cr(k), cerr, cerrl, psierr(k), perrl
  530 continue
      go to 10
c
c   reconstruct vertical modes from their values at discrete levels
c
  540 read ( line, '(20x,3i10)' )  n1, n2, n3
      if ( anavsp )  then
          print *, 'this option available for numerical solution only'
          go to 10
      else if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      end if
      if ( 0.gt.n1 .or. n1.gt.n )  n1 = 0
      if ( 1.gt.n2 .or. n2.gt.n )  n2 = n
      if ( n3.eq.0 )  n3 = 1
      nf = 0
  550 read '(a)', line
      print '('' input--'',a)', line
c         first = line(1:1)
          read ( line, '(a1)' )  first
          if ( first.eq.'?' )  then
              print *,'(this level ignored)'
              go to 550
          end if
          if ( first.ne.' ' )  go to 560
          if ( nf.ge.nfmax )  then
              print *, 'too many data points'
              stop
          end if
          nf = nf + 1
          read ( line, '(f80.0)' )  pf(nf)
          if ( nf.gt.1 .and. pf(nf).le.pf(nf-1) )  then
              print *, 'levels must be in increasing order'
              stop
          end if
      go to 550
  560 labelf = line(1:60)
      nip = 2*(nf + n)
      ncol = (n2 - n1)/n3 + 1
      npages = (ncol - 1)/npp + 1
      first  = '1'
      labeln = '                  '
      do 590 ipage=1,npages
          m1 = n1 + (ipage - 1)*npp*n3
          m2 = min( n2, m1 + (npp - 1)*n3 )
          ncol = (m2 - m1)/n3 + 1
          do 580 j=1,ncol
              mode = m1 + (j-1)*n3
              i = mode*(n + 1)
              call veval( tinv(i,0), 0, nf, pf, con, n, w, f )
              call vproj(nf,pf,f,0,1,nip,pip,wip,prom,n,con,w,fvbf,ierr)
              if ( ierr.ne.0 )  print 9000, 'vproj', ierr
              call vtran( fvbf, n, tfor, fhat )
              do 570 k=0,n
  570         me(k,j) = nint( 100.0d0*fhat(k)**2 )
              m(j) = mode
  580     continue
          if ( npages.gt.1 )  write ( labeln, 5300 ) ipage, npages
          print 2550, first, labeln, nf, labelf, (m(j), j=1,ncol )
          print '(10x)'
          do 585 k=0,n
  585     print '( 1x, i3, f7.2, 5x, 21i5 )',
     2                 k , c(k), (me(k,j), j=1,ncol)
          first = ' '
  590 continue
 2550 format(a,9x///'0percent energy for vertical structure function ',
     2              'reconstruction',5x,a,
     3             /'0from values at',i4,'  discrete levels:  ',a,
     4            //'0  k  c(k)   mode:',i3,20i5)
      print 2330, n
      nf = 0
      go to 10
c
c   input the function to be transformed
c
  600 labelf = line(21:80)
      nf = 0
  610 read '(a)', line
      print '('' input--'',a)', line
c         first = line(1:1)
          read ( line, '(a1)' )  first
          if ( first.eq.'?' )  then
              print *,'(this data point ignored)'
              go to 610
          end if
          if ( first.ne.' ' )  go to 620
          if ( nf.ge.nfmax )  then
              print *, 'too many data points'
              stop
          end if
          nf = nf + 1
          read ( line, 1100 )  pf(nf), fdata(nf)
      go to 610
  620 call sorter( nf, pf, fdata )
      go to 20
c
c   evaluate a vertical structure function at discrete levels
c
  650 read ( line, '(20x,i10)' )  mode
      if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      else if ( np.le.0 )  then
          print *, 'must specify the output pressure levels first'
          go to 10
      else if ( np.gt.nfmax )  then
          print *, 'too many pressure levels for input data'
          go to 10
      end if
      write ( labelf, 5650 )  mode
 5650 format('vertical structure function',i4,'  from discrete data')
      nf = np
      do 660 i=1,nf
  660 pf(i) = p(i)
      if ( anavsp )  then
          call  cgefun( mode, nf, pf, rgamma, rtb, pbot, ptop,
     2                  c, psihat, fdata )
      else
          i = mode*(n + 1)
          call veval( tinv(i,0), 0, nf, pf, con, n, w, fdata )
      end if
      print 2650, mode
 2650 format('0internally generated function data:',
     2      /'0       p        psi(',i2,')'/)
      print '(1x,f10.1,f13.5)', ( pf(i), fdata(i), i=1,nf )
      go to 10
c
c   project the function (or its integral) onto the vertical basis
c
  700 in = 0
      factor = one
      go to 720
  710 in = 1
      if ( line(16:20).eq.'kappa' )  then
          factor = kappa
      else if ( line(16:20).eq.'rgas ' )  then
          factor = rgas
      else
          factor = one
      end if
  720 read ( line, '(20x,2i10)' )  ip, nip
      if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      end if
      if ( ip.le.0 )  ip = 1
      if ( ip.eq.1 .and. nip.le.0 )  nip = 2*(nf + n)
      if ( ip.eq.2 .and. nip.le.0 )  nip = 2*(nf + n)
      if ( ip.ge.3 )  then
          nip = nf
          if ( nip.gt.nipmax )  stop
          do 730 i=1,nip
              pip(i) = pf(i)
              wip(i) = one
  730     continue
      end if
      do 740 i=1,nf
  740 f(i,1) = factor*fdata(i)
      if ( factor.ne.one )  print 2700, factor
 2700 format('0********** function values multiplied by',1pe12.5,
     2       ' **********')
c     time = second( work )
      time = second( )
      call vproj( nf,pf,f,in,ip,nip,pip,wip,prom,n,con,w,fvbf,ierr )
c     time = second( work ) - time
      time = second( ) - time
      if ( ierr.ne.0 )  print 9000, 'vproj', ierr
      if ( in.eq.0 )  print 2710, time, 'function', n, labelp(ip), nip
      if ( in.eq.1 )  print 2710, time, 'integral', n, labelp(ip), nip
 2710 format('0cpu time:',f8.3,' seconds for projection of ',a,
     2       '  (n =',i3')',
     3       /' (projection inner product:  ',a,' with',i4,' points)')
      go to 10
c
c   compute and print physical space profiles
c
  800 read ( line, '(19x,i1,3i10)' )  k, n1, n2, n3
      if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      else if ( ip.le.0 )  then
          print *, 'must project or integrate the function first'
          go to 10
      else if ( np.le.0 )  then
          print *, 'must specify the output pressure levels first'
          go to 10
      end if
      k  = max( k, 0 )
      k  = min( k, 2 )
      md = min( k, 1 )
      if ( 0.gt.n1 .or. n1.gt.n )  n1 = 0
      if ( 1.gt.n2 .or. n2.gt.n )  n2 = n
      if ( n3.eq.0 )  n3 = 1
      if ( n1.ne.0 .or. n2.ne.n .or. n3.ne.1 )  then
          call vtran( fvbf, n, tfor, fhat )
          do 810 i=0,n
  810     en(i) = zero
          do 820 i=n1,n2,n3
  820     en(i) = fhat(i)
          call vtran( en, n, tinv, ef )
      else
          do 830 i=0,n
  830     ef(i) = fvbf(i)
      end if
      call veval( ef, md, np, p, con, n, w, fval )
      if ( k.eq.2 )  then
          do 840 i=1,np
  840     fval(md*np+i,1) = -p(i)*fval(md*np+i,1)/factor
      end if
      print 2800, labelt, labelf
 2800 format('1 '//'0physical space profiles',
     2           //'0basic state:  ',a60,
     3           //'0function:     ',a60)
      if ( in.gt.0 )  print 8000
 8000 format(15x,'(integrated with respect to  -log(p) )')
      print 2810, label(k)
 2810 format('0 '//21x,'basic  state',15x,a,' of function',
     2           //6x,'p',9x,'t',6x,'gamma',6x,'sigma',
     3             7x,'cubic spline',4x,'projection'/)
      err2 = zero
      erri = zero
      do 850 i=1,np
          ts    = csval( p(i), 0, nt, pt, rt, 0 )/rgas - tzero
          gamma = rrgval( p(i), con, nt, rt, ierr )
          if ( ierr.ne.0 )  print 9000, 'rrgval', ierr
          gamma = one/gamma
          sigma = gamma/p(i)**2
          gamma = gamma/rgas
          fs = funval( p(i), md, nf, pf, f )
          if ( k.eq.2 )  fs = -p(i)*fs/factor
          fp = fval(md*np+i,1)
          err2 = err2 + (fp - fs)**2
          erri = max( erri, abs( fp - fs) )
          print '( 1x, 0pf8.1, f10.2, f9.2, 1pe13.3, e16.4, e15.4 )',
     2                  p(i) ,   ts , gamma,  sigma,   fs ,   fp
  850 continue
      err2 = sqrt( err2/np )
      print 2820, err2, erri
 2820 format('0',38x,'rms projection error:',1pe12.4,
     2          /39x,'max projection error:',e12.4)
      if ( n1.ne.0 .or. n2.ne.n .or. n3.ne.1 )  print 2830, n1, n2, n3
 2830 format('0(vertical modes',i4,'  to',i4,'  step',i4,
     2       '  used in reconstructing the projection)')
      if ( tzero.eq.zero )  print 2840,  'kelvin', n
      if ( tzero.ne.zero )  print 2840, 'celsius', n
 2840 format('0(t  in ',a,',  gamma  in kelvin, ',
     2       'spectral truncation  n =',i4,')')
      go to 10
c
c   compute and print spectral space profiles
c
  900 read ( line, '(20x,3i10)' )  n1, n2, n3
      if ( n.le.0 )  then
          print *, 'must solve the vertical structure problem first'
          go to 10
      else if ( ip.le.0 )  then
          print *, 'must project or integrate the function first'
          go to 10
      end if
      if ( 0.gt.n1 .or. n1.gt.n )  n1 = 0
      if ( 1.gt.n2 .or. n2.gt.n )  n2 = n
      if ( n3.eq.0 )  n3 = 1
      call vtran( fvbf, n, tfor, fhat )
      ensum = zero
      efsum = zero
      do 910 k=0,n
          if ( in.eq.0 )  en(k) =  fhat(k)**2
          if ( in.gt.0 )  en(k) = (fhat(k)/c(k))**2
          ensum = ensum + en(k)
          if ( in.eq.0 )  ef(k) =  fhat(k)
          if ( in.gt.0 )  ef(k) =  fhat(k)/c(k)**2
          efsum = efsum + abs( ef(k) )
  910 continue
      if ( ensum.eq.zero )  ensum = one
      if ( efsum.eq.zero )  efsum = one
      print 2900, labelt, labelf
 2900 format('1 '//'0spectral space profiles',
     2           //'0basic state:  ',a60,
     3           //'0function:     ',a60)
      if ( in.gt.0 )  print 8000
      print 2910
 2910 format(2x/'0',33x,'perturbation energy',2x,'forcing  efficiency',
     2      /3x,'k',6x,'c(k)',8x,'fhat(k)',9x,'value',3x,'percent',
     3                                     6x,'value',3x,'percent'/)
      print 2920, ( k, c(k), fhat(k), en(k), 100*en(k)/ensum,
     2                ef(k), 100*abs( ef(k) )/efsum, k=n1,n2,n3 )
 2920 format(1x,i3,f11.3,1pe16.5,e14.3,0pf7.2,1pe14.3,0pf7.2)
      print 2930, ensum, efsum
 2930 format('0',25x,'totals:',1p,e12.3,e21.3)
      if ( in.eq.0 )  print 2940, ' ', ' '
      if ( in.gt.0 )  print 2940, '/c(k)**2', '/c(k)**2'
 2940 format('0perturbation energy = fhat(k)**2',a,
     2      /' forcing  efficiency = fhat(k)   ',a)
      print 2950, n
 2950 format('0(phase speeds  c  in meters per second, ',
     2       'spectral truncation  n =',i4,')')
      go to 10
      end
      subroutine sorter( n, x, y )
      integer            n
      double precision   x(n), y(n)
c
c   sorts pairs (x(i),y(i)) into order of increasing  x  (bubble sort)
c
      integer           i, j, k
      double precision  temp
c
      i = 1
   10 if ( i.ge.n )  return
      i = i + 1
      j = i
   20 if ( x(j-1).eq.x(j) )  then
          do 30 k=j+1,n
              x(k-1) = x(k)
              y(k-1) = y(k)
   30     continue
          n = n - 1
          i = i - 1
          go to 20
      end if
      if ( x(j-1).le.x(j) )  go to 10
      temp   = x(j-1)
      x(j-1) = x(j)
      x(j)   = temp
      temp   = y(j-1)
      y(j-1) = y(j)
      y(j)   = temp
      j = j-1
      if ( j.gt.1 )  go to 20
      go to 10
      end
      subroutine plevel( pmin, pmax, dp, pbot, np, p )
c
c   routine to generate the output pressure levels
c
      integer           np, i
      double precision  pmin, pmax, dp, pbot, p(1), p1, p2, p3
      double precision  zero, one
      parameter         ( zero = 0.0d0, one = 1.0d0 )
c
      p1 = min( pmin, pmax )
      p2 = max( pmin, pmax )
      p3 = abs( dp )
      if ( p3.eq.zero )  p3 = p2 - p1
      if ( p3.eq.zero )  p3 = one
      np = (p2 - p1)/p3 + one + 1.0d-6
      do 10 i=1,np
   10 p(i) = p1 + (i-1)*p3
      if ( p(np).ge.pbot )  return
      np = np+1
      p(np) = pbot
      return
      end
      function fgrand( pval )
      double precision  fgrand, pval
c
c   integrand to give  r*t  from  1/(r*gamma)
c
      double precision  kappa, rgamma, rrgval
      parameter         ( kappa = 0.2859d0 )
      integer           ierr
      integer           ntmax
      parameter         ( ntmax =  47 )
      integer           nt
      double precision  con(4), rt(ntmax,9)
      common   /fgradc/ con, rt, nt
c
      rgamma = rrgval( pval, con, nt, rt, ierr )
      if ( ierr.ne.0 )  stop
      fgrand = (con(2)/pval)**kappa/(pval*rgamma)
      return
      end
