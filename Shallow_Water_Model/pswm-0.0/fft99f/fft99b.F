      subroutine fft99b(work,a,trigs,inc,jump,n,lot)
      integer           inc, jump, n, lot
      double precision  work(n), a(n), trigs(3*n/2+1)
c
c     subroutine fft99b - postprocessing step for fft99, isign=-1
c     (gridpoint to spectral transform)
c
      integer           ia, iabase, ib, ibbase, ink 
      integer           ja, jabase, jb, jbbase, k, l, nh, nx
      double precision  c, s, scale, zero, one, half, two
      parameter         ( zero = 0.0d0, one = 1.0d0 )
      parameter         ( half = 0.5d0, two = 2.0d0 )
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) and a(n/2)
      scale=one/dble(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
cdir$ ivdep
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=zero
      a(jb+inc)=zero
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
c
c     remaining wavenumbers
      scale=half*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
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
      a(ja)=scale*((work(ia)+work(ib))
     *   +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))
     *   -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
      scale=two*scale
cdir$ ivdep
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
c
   50 continue
      return
      end
