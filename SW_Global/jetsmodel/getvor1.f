      subroutine getvor1(nx,my,lev,vorg1,div,sinl)

******************************************************************************

      dimension vorg1(nx,lev,my),div(nx,lev,my),sinl(my)

c      local array

      integer nvx
      parameter(nvx=1)
      dimension xlon(nx),xlat(my)
      dimension rlon(nvx),rlat(nvx),rmax(nvx),voramp(nvx)

c      data rlon/115./, rlat/32./, rmax/500./, voramp/-5.0e-5/
c      data rlon/120./, rlat/25./, rmax/500./, voramp/-5.0e-5/
c      data rlon/120./, rlat/25./, rmax/500./, voramp/-5.0e-5/ ! High Pressure

      data rlon/90./,rlat/20./
      data rmax/800./,voramp/3.e-4/    ! 1 Low Pressure
c      data rmax/300./,voramp/1.67e-4/    ! 1 Low Pressure

      radx=6.37e3
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r

      xlon(1)=0.
      do 10 i=2,nx
        xlon(i)=xlon(1)+float(i-1)*360./nx
 10   continue
      do 15 j=1,my
        xlat(j)=asin(sinl(j))*r2d
 15   continue



      do n=1,nvx,1 

      ix=rlon(n)/(360./nx)+1.001
      jy=(rlat(n)-xlat(1))/(180./my)+1.001
      print *,'h  ix=',ix,' jy=',jy,' rlat=',rlat,' rlon=',rlon

      jb=jy-30
      je=jy+30
      ib=ix-30
      ie=ix+30

      do 200 j=1,my
      do 200 i=1,nx
        vorg1(i,1,j)=0.
        div(i,1,j)=0.
 200  continue

      a=rmax(nvx)
      b=a
      do 101 j=jb,je
      do 101 i=ib,ie
c      do 101 j=1,my
c      do 101 i=1,nx
        dx=radx*cos(rlat(n)*d2r)*(rlon(n)-xlon(i))*d2r
        dy=radx*(rlat(n)-xlat(j))*d2r
        vorg1(i,1,j)=voramp(n)*exp(-(dx/a)**2-(dy/b)**2)
 101  continue
      print *,'vorg1(ix,1,jy)=',vorg1(ix,1,jy)

      enddo

c
c      rf=15.
c      do j=1,my
c        do i=1,nx
c          dx=radx*cos(rlat*d2r)*(rlon-xlon(i))*d2r
c          dy=radx*(rlat-xlat(j))*d2r
c          r=(dx*dx+dy*dy)**0.5
c          rr=r/rmax
c          if(rr.eq.0.)then
c            vorg(i,1,j)=1.*voramp
c          else if(rr.gt.0. .and. rr.lt.1.)then
c            vorg(i,1,j)=(1.-exp(-rf/rr*exp(1./(rr-1.))))*voramp
c          else
c            vorg(i,1,j)=0.
c          endif
c        enddo
c      enddo
c      print *,'vorg(ix,1,jy)=',vorg(ix,1,jy)
c
      return
      end
