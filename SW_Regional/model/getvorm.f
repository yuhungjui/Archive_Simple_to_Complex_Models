      subroutine getvorm(nx,my,lev,wmamp,vorm,divm)
      dimension vorm(nx,lev,my),divm(nx,lev,my)
c
      data aa,bb,dd,ee/0.,2.,0.,1./
c      wmamp=1.e-6
      pi=atan(1.0)*4.
c
      im=nx/2
      jm=my/2
c
c      clat=my-(my/2.)+0.0001
c      clon=nx-(nx/2.)+0.0001
      clat=120.00001
      clon=240.00001
      a=50.
      b=30.
      rf=5.0
      thda=0.
      pi=4.*atan(1.0)
      rot=thda*pi/180.
      amp=-7.5e-6
c      amp=0.
c
      do 24 j=1,my
      do 24 i=1,nx
        vorm(i,1,j)=-(aa-bb*sin(2.*pi*j/my))
     &              *(dd-ee*sin(2.*pi*i/nx))*wmamp
c        vorm(i,1,j)=-tanh(pi*(j-jm)/my)*tanh(pi*(i-im)/nx)*wmamp
        x=(i-clon)*1.0
        y=(j-clat)*1.0
        arg=atan2(y,x)
        ra=(x**2+y**2)**0.5
        ax=ra*cos(arg+rot)/a
        by=ra*sin(arg+rot)/b
        r=(ax**2+by**2)**0.5
c        r=((x/a)**2+(y/b)**2)**0.5
        if(r.eq.0)then
          divm(i,1,j)=1.*amp
        else if(r.gt.0. .and. r.lt.1.)then
          divm(i,1,j)=(1.-exp(-rf/r*exp(1./(r-1.))))*amp
        else
          divm(i,1,j)=0.
        endif
 24   continue
c
      return
      end
