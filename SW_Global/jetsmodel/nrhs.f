      subroutine nrhs(iter)
c
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
c  working array
c

c--------------------------------------------llming

      integer temp(55),iter

      real ratio
      real force(0:int(its))
c      real force(0:4001001)

c------------------------------------------------------

      dimension gg(nx,lev,my),hh(nx,lev,my),ws1(mlmax,2,lev)
     &,         ws2(mlmax,2,lev),wk(nx,my),wks(mlmax,2)
c
c  compute nonlinear term for VOR and DIV eqn in physical space
c
      do 10 k=1,lev
      do 10 j=1,my
      do 10 i=1,nx
        gg(i,k,j)=-1.0*(vor(i,k,j)+cor(j))*uu(i,k,j)
        hh(i,k,j)=-1.0*(vor(i,k,j)+cor(j))*vv(i,k,j)
 10   continue
C
C  Calculate spectral tendency for VORTICITY eqn,  -alpha(GG, HH)
C
      call alpha2s(jtrun,mlmax,nx,my,lev,gg,hh,weight,cim
     *, onocos,poly,dpoly,vorten)

c---------------------------------------add viscosity     ! llming

      do j=1,my,1
         do i=1,nx,1
            gg(i,1,j)=2.4848e75*vor(i,1,j)    ! viscosity coefficient  efolding=0.5 days
         enddo
      enddo

      call  tranrs(jtrun,mlmax,nx,my,lev,poly,weight,gg,ws1)

      do i=1,mlmax,1
         vorten(i,1,1)=vorten(i,1,1)-ws1(i,1,1)*epsd16(i)   ! -  correct
         vorten(i,2,1)=vorten(i,2,1)-ws1(i,2,1)*epsd16(i)   ! -
      enddo

c----------------------------------------


c----------------------------------------add Ekman drag    ! llming

      do j=1,my,1
         do i=1,nx,1
c            gg(i,1,j)=2.314815e-8*vor(i,1,j)    ! viscosity coefficient  efolding=500 days
            gg(i,1,j)=1.15741e-8*vor(i,1,j)    ! viscosity coefficient  efolding=1000 days
c            gg(i,1,j)=7.71605e-9*vor(i,1,j)    ! viscosity coefficient  efolding=1500 days
c            gg(i,1,j)=3.858025e-9*vor(i,1,j)    ! viscosity coefficient  efolding=3000 days
         enddo
      enddo

      call  tranrs(jtrun,mlmax,nx,my,lev,poly,weight,gg,ws1)

      do i=1,mlmax,1
         vorten(i,1,1)=vorten(i,1,1)-ws1(i,1,1)   ! -
         vorten(i,2,1)=vorten(i,2,1)-ws1(i,2,1)   ! -
      enddo

c--------------------------------------------------


c--------------------------------------------------   add forcing  ! llming

c      ratio=(1-dt/172800.)/(1+dt/172800.)
c
c      force(0)=0.
c      force(iter)=force(iter-1)*ratio+
c     +            sqrt(1.-ratio**2)*ran(iter)*(5.e-13)    !! forcing amplitude
c
c      do i=2,56,1
c         temp(i)=mlsort(i,56)
c         vorten(temp(i),1,1)=vorten(temp(i),1,1)+force(iter)
c         vorten(temp(i),2,1)=vorten(temp(i),2,1)+force(iter)
c      enddo

c--------------------------------------------------

c
C IF BAROTROPIC CALCULATION, NO TENDENCY FOR D AND PHI
C
      if(baro)then
        do 20 i=1,mlmax*2*lev
          divten(i,1,1)=0.0
          phiten(i,1,1)=0.0
 20     continue
        return
      endif
C
C  Calculate spectral tendency for DIVERGENCE eqn, alpha(HH, -GG)
C                                                  - L2(II + PPHY)
      do 30 i=1,nx*my*lev
        hh(i,1,1)=-1.0*hh(i,1,1)
 30   continue
C
      call alpha2s(jtrun,mlmax,nx,my,lev,hh,gg,weight,cim
     *, onocos,poly,dpoly,divten)
C
      do 40 k=1,lev
      do 40 j=1,my
      do 40 i=1,nx
        gg(i,k,j)=(uu(i,k,j)**2+vv(i,k,j)**2)*0.5*radsq*onocos(j)
     *           +(hmean+phi(i,k,j)-topo(i,j))
c     *           +phi(i,k,j)
 40   continue
c
      call  tranrs(jtrun,mlmax,nx,my,lev,poly,weight,gg,ws1)
C
      do 50 k=1,lev
      do 50 i=1,mlmax
        divten(i,1,k)=divten(i,1,k)+ws1(i,1,k)*eps4(i)
        divten(i,2,k)=divten(i,2,k)+ws1(i,2,k)*eps4(i)
 50   continue
c
c add friction
c
      do j=1,my
      do i=1,nx
        gg(i,1,j)=fricv*vor(i,1,j)
        hh(i,1,j)=fricd*div(i,1,j)
      enddo
      enddo
      call  tranrs(jtrun,mlmax,nx,my,lev,poly,weight,gg,ws1)
      call  tranrs(jtrun,mlmax,nx,my,lev,poly,weight,hh,ws2)
      do i=1,mlmax*2
        vorten(i,1,1)=vorten(i,1,1)-ws1(i,1,1)
        divten(i,1,1)=divten(i,1,1)-ws2(i,1,1)
      enddo
C
C    Calculate spectral tendecny for MASS eqn
C
      do 60 k=1,lev
      do 60 j=1,my
      do 60 i=1,nx
        gg(i,k,j)=(hmean+phi(i,k,j)-topo(i,j))*uu(i,k,j)
        hh(i,k,j)=(hmean+phi(i,k,j)-topo(i,j))*vv(i,k,j)
 60   continue
      call alpha2s(jtrun,mlmax,nx,my,lev,gg,hh,weight,cim
     *, onocos,poly,dpoly,phiten)
c
      do i=1,mlmax*2
        phiten(i,1,1)=-1.0*phiten(i,1,1)
      enddo
c
c-add mass forcing
c for mass
c      amp=-9.58333e-4
c      amp=-20.e-4
c for vorticity
      voramp=1.e-4
c      amp=0.
      t0=0.5*24.
      t1=0.*24.
      t2=t0-t1
      w0=1./(0.5*86400.)
      amp=voramp*w0
      xtau=tau-t1
      if(xtau.lt.0.)then
        gm=-1.
      else
        gm=abs(2.*xtau/t2-1.)
      endif
      if(gm.eq.0.)then
        ft=1.
      else if(gm.gt.0. .and. gm.lt.1.)then
        ft=1.-exp(-320./gm*exp(1./(gm-1.)))
      else
        ft=0.
      endif
      print *,'ft=',ft
      do 323 i=1,mlmax*2
c        phiten(i,1,1)=phiten(i,1,1)+fspc(i,1)*amp*ft
        vorten(i,1,1)=vorten(i,1,1)+fspc(i,1)*amp*ft
 323  continue
c
      return
      end
