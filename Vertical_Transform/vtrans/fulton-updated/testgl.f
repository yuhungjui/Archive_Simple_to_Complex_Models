      program testgl
c
c   Main program to test the Gauss-Legendre subroutine  gaussl
c
      integer            nmax
      parameter          ( nmax = 100 )
      integer            j, n, ierr
      double precision   a, b, x(nmax), w(nmax)
c
c   Loop on the number of Gauss points
c
      a = -1.0
      b = +1.0
c     do 10 n=1,16
      do 10 n=24,96,24
          call  gaussl( n, a, b, x, w, ierr )
          print 2000, n, a, b
 2000     format(/'Gauss-Legendre abscissas and weights from gaussl:',
     &           /'n =',i4,'  a =',f5.2,'  b =',f5.2,
     &          //'  j   x(j)   w(j)')
          if ( ierr.eq.0 )  then
              print '(i3,2f20.15)', ( j, x(j), w(j), j=1,n )
          else
              print *, 'gaussl error:  ierr =', ierr
          end if
   10 continue
      end
