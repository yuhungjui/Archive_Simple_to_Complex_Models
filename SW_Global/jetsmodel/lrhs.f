      subroutine lrhs
C-----------------------------------------------------------------------
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
c  working array
c
      dimension gg(nx,lev,my),hh(nx,lev,my)
      dimension ws1(mlmax,2,lev)
c
c  compute linear term for VOR and DIV eqn in physical space
c
      do 10 k=1,lev
      do 10 j=1,my
      do 10 i=1,nx
        gg(i,k,j)=-1.0*cor(j)*uu(i,k,j)
        hh(i,k,j)=-1.0*cor(j)*vv(i,k,j)
 10   continue
C
C  Calculate spectral tendency for VORTICITY eqn,  -alpha(GG, HH)
C
      call alpha2s(jtrun,mlmax,nx,my,lev,gg,hh,weight,cim
     *, onocos,poly,dpoly,vorten)

c---------------------------------------add viscosity     ! llming
      do j=1,my,1
         do i=1,nx,1
            gg(i,1,j)=6.5*vor(i,1,j)    ! viscosity coefficient
         enddo
      enddo

      call  tranrs(jtrun,mlmax,nx,my,lev,poly,weight,gg,ws1)

      do i=1,mlmax,1
         vorten(i,1,1)=vorten(i,1,1)+ws1(i,1,1)*eps4(i)
         vorten(i,2,1)=vorten(i,2,1)+ws1(i,2,1)*eps4(i)
      enddo

c      do i=1,mlmax*2
c         vorten(i,1,1)=vorten(i,1,1)+ws1(i,1,1)*eps4(i)
c      enddo
c----------------------------------------

      if(baro)then
        do 20 i=1,mlmax*2*lev
          divten(i,1,1)=0.0
          phiten(i,1,1)=0.0
 20     continue
        return
      endif
c
C  Calculate spectral tendency for DIVERGENCE eqn, alpha(HH, -GG)
C                                                  - L2(II + PPHY)
C
      do 30 i=1,nx*my*lev
        hh(i,1,1)=-1.0*hh(i,1,1)
 30   continue
C
      call alpha2s(jtrun,mlmax,nx,my,lev,wd2,wd1,weight,cim
     *, onocos,poly,dpoly,divten)
c
      do 40 k=1,lev
      do 40 i=1,mlmax
        divten(i,1,k)=divten(i,1,k)+pspcnow(i,1,k)*eps4(i)
        divten(i,2,k)=divten(i,2,k)+pspcnow(i,2,k)*eps4(i)
 40   continue
C
C    Calculate spectral tendecny for MASS eqn
C
      do 50 i=1,mlmax*2*lev
        phiten(i,1,1)=-hmean*dspcnow(i,1,1)
 50   continue
C
      return
      end
