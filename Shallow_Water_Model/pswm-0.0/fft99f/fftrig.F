      subroutine fftrig(trigs,n,mode)
      integer           n, mode
      double precision  trigs(3*n/2+1)
      integer           i, imode, l, la, nh, nn
      double precision  angle, del, pi, one, two, half
      parameter         ( one = 1.0d0 )
      parameter         ( two = 2.0d0, half = 0.5d0 )
c
      pi=two*asin(one)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/dble(nn)
      l=nn+nn
      do 10 i=1,l,2
      angle=half*dble(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=half*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=half*dble(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=half*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=dble(i-1)*del
      trigs(la+i)=two*sin(angle)
   30 continue
      return
   40 continue
      del=half*del
      do 50 i=2,n
      angle=dble(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end
