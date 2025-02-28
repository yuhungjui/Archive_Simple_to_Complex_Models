      subroutine fft99a(a,work,trigs,inc,jump,n,lot)
      integer           inc, jump, n, lot
      SINGLE_OR_DOUBLE  a(n), work(n), trigs(3*n/2+1)
c
c     subroutine fft99a - preprocessing step for fft99, isign=+1
c     (spectral to gridpoint transform)
c
      integer           ia, iabase, ib, ibbase, ink
      integer           ja, jabase, jb, jbbase, k, l, nh, nx
      SINGLE_OR_DOUBLE  c, s, two
      parameter         ( two = 2.0EorD0 )
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
cdir$ ivdep
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
c
c     remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cdir$ ivdep
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+
     *    (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-
     *    (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
cdir$ ivdep
      do 40 l=1,lot
      work(ja)=two*a(ia)
      work(ja+1)=-two*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
c
   50 continue
      return
      end
