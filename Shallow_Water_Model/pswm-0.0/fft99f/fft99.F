      subroutine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
      integer          ifax(13), inc, jump, n, lot, isign
      double precision a(n), work(n), trigs(3*n/2+1)
c
c Purpose      Performs multiple Fast Fourier Transforms.  This package
c              will perform a number of simultaneous real/half-complex
c              periodic Fourier transforms or corresponding inverse
c              transforms, i.e.  given a set of real data vectors, the
c              package returns a set of 'half-complex' Fourier
c              coefficient vectors, or vice versa.  The length of the
c              transforms must be an even number greater than 4 that has
c              no other factors except possibly powers of 2, 3, and 5.
c              This is an all fortran version of the craylib package
c              that is mostly written in cal.
c
c              The package fft99f contains several user-level routines:
c
c            subroutine fftfax
c                an initialization routine that must be called once
c                before a sequence of calls to the fft routines
c                (provided that n is not changed).
c
c            subroutines fft99 and fft991
c                two fft routines that return slightly different
c                arrangements of the data in gridpoint space.
c
c Usage        Let n be of the form 2**p * 3**q * 5**r, where p .ge. 1,
c              q .ge. 0, and r .ge. 0.  Then a typical sequence of
c              calls to transform a given set of real vectors of length
c              n to a set of 'half-complex' Fourier coefficient vectors
c              of length n is
c
c                   integer          ifax(13), inc, isign, jump, m, n
c                   double precision a(m*(n+2)), trigs(3*n/2+1),
c                  +                 work(m*(n+1))
c
c                   call fftfax (n, ifax, trigs)
c                   call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
c
c              See the individual write-ups for fftfax, fft99, and
c              fft991 below, for a detailed description of the
c              arguments.
c
c History      The package was written by Clive Temperton at ecmwf in
c              November, 1978.  It was modified, documented, and tested
c              for NCAR by Russ Rew in September, 1980.
c              Converted to single/double precision by Scott Fulton in
c              November 1990 (also added save statements)
c
c-----------------------------------------------------------------------
c
c subroutine fftfax (n,ifax,trigs)
c
c Purpose      A set-up routine for fft99 and fft991.  It need only be
c              called once before a sequence of calls to the fft
c              routines (provided that n is not changed).
c
c Argument     ifax(13), trigs(3*n/2+1)
c Dimensions
c
c Arguments
c
c On Input     n
c               an even number greater than 4 that has no prime factor
c               greater than 5.  N is the length of the transforms (see
c               the documentation for fft99 and fft991 for the
c               definitions of the transforms).
c
c              ifax
c               an integer array.  The number of elements actually used
c               will depend on the factorization of n.  Dimensioning
c               ifax for 13 suffices for all n less than a million.
c
c              trigs
c               a floating point array of dimension 3*n/2 if n/2 is
c               even, or 3*n/2+1 if n/2 is odd.
c
c On Output    ifax
c               contains the factorization of n/2.  ifax(1) is the
c               number of factors, and the factors themselves are stored
c               in ifax(2),ifax(3),...  If fftfax is called with n odd,
c               or if n has any prime factors greater than 5, ifax(1)
c               is set to -99.
c
c              trigs
c               an array of trignometric function values subsequently
c               used by the fft routines.
c
c-----------------------------------------------------------------------
c
c subroutine fft991 (a,work,trigs,ifax,inc,jump,n,m,isign)
c                       and
c subroutine fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
c
c Purpose      Perform a number of simultaneous real/half-complex
c              periodic Fourier transforms or corresponding inverse
c              transforms, using ordinary spatial order of gridpoint
c              values (fft991) or explicit cyclic continuity in the
c              gridpoint values (fft99).  Given a set
c              of real data vectors, the package returns a set of
c              'half-complex' Fourier coefficient vectors, or vice
c              versa.  The length of the transforms must be an even
c              number that has no other factors except possibly powers
c              of 2, 3, and 5.  These version of fft991 and fft99 are
c              optimized for use on the cray-1.
c
c Argument     a(m*(n+2)), work(m*(n+1)), trigs(3*n/2+1), ifax(13)
c Dimensions
c
c Arguments
c
c On Input     a
c               an array of length m*(n+2) containing the input data
c               or coefficient vectors.  This array is overwritten by
c               the results.
c
c              work
c               a work array of dimension m*(n+1)
c
c              trigs
c               an array set up by fftfax, which must be called first.
c
c              ifax
c               an array set up by fftfax, which must be called first.
c
c              inc
c               the increment (in words) between successive elements of
c               each data or coefficient vector (e.g.  inc=1 for
c               consecutively stored data).
c
c              jump
c               the increment (in words) between the first elements of
c               successive data or coefficient vectors.  On the cray-1,
c               try to arrange data so that jump is not a multiple of 8
c               (to avoid memory bank conflicts).  For clarification of
c               inc and jump, see the examples below.
c
c              n
c               the length of each transform (see definition of
c               transforms, below).
c
c              m
c               the number of transforms to be done simultaneously.
c
c              isign
c               = +1 for a transform from Fourier coefficients to
c                    gridpoint values.
c               = -1 for a transform from gridpoint values to Fourier
c                    coefficients.
c
c On Output    a
c               If isign = +1, and m coefficient vectors are supplied
c               each containing the sequence:
c
c               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
c
c               then the result consists of m data vectors each
c               containing the corresponding n+2 gridpoint values:
c
c               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
c               for fft99, x(n-1),x(0),x(1),x(2),...,x(n-1),x(0).
c                   (explicit cyclic continuity)
c
c               When isign = +1, the transform is defined by:
c                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c                 and i=sqrt (-1)
c
c               If isign = -1, and m data vectors are supplied each
c               containing a sequence of gridpoint values x(j) as
c               defined above, then the result consists of m vectors
c               each containing the corresponding Fourier cofficients
c               a(k), b(k), 0 .le. k .le n/2.
c
c               When isign = -1, the inverse transform is defined by:
c                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
c                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
c
c               A call with isign=+1 followed by a call with isign=-1
c               (or vice versa) returns the original data.
c
c               Note: the fact that the gridpoint values x(j) are real
c               implies that b(0)=b(n/2)=0.  For a call with isign=+1,
c               it is not actually necessary to supply these zeros.
c
c Examples      Given 19 data vectors each of length 64 (+2 for explicit
c               cyclic continuity), compute the corresponding vectors of
c               Fourier coefficients.  The data may, for example, be
c               arranged like this:
c
c first data   a(1)=    . . .                a(66)=             a(70)
c vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
c
c second data  a(71)=   . . .                                  a(140)
c vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
c
c               and so on.  Here inc=1, jump=70, n=64, m=19, isign=-1,
c               and fft99 should be used (because of the explicit cyclic
c               continuity).
c
c               Alternatively the data may be arranged like this:
c
c                first         second                          last
c                data          data                            data
c                vector        vector                          vector
c
c                 a(1)=         a(2)=                           a(19)=
c
c                 x(63)         x(63)       . . .               x(63)
c        a(20)=   x(0)          x(0)        . . .               x(0)
c        a(39)=   x(1)          x(1)        . . .               x(1)
c                  .             .                               .
c                  .             .                               .
c                  .             .                               .
c
c               in which case we have inc=19, jump=1, and the remaining
c               parameters are the same as before.  In either case, each
c               coefficient vector overwrites the corresponding input
c               data vector.
c
c-----------------------------------------------------------------------
c
c     subroutine "fft99" - multiple fast real periodic transform
c     corresponding to old scalar routine fft9
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis and welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input and output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data 'vector'
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
c         i.e. explicit cyclic continuity; (n+2) locations required
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
      integer i, ia, ib, ibase, igo, ink, j, jbase, k
      integer l, la, m, nfax, nh, nx
c
      nfax=ifax(1)
      if ( nfax.le.0 ) return
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
c  50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=ia
      do 100 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in cyclic boundary points
  110 continue
      ia=1
      ib=n*inc+1
cdir$ ivdep
      do 120 l=1,lot
      a(ia)=a(ib)
      a(ib+inc)=a(ia+inc)
      ia=ia+jump
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
      end
