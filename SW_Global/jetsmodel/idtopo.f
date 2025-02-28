      subroutine idtopo(nx,my,sinl,topo)

******************************************************************************
      dimension sinl(my),topo(nx,my)
c
      dimension xlon(nx),xlat(my)
      dimension rlat(2),rlon(2),rotang(2),topamp(2),rxx(2),ryy(2)
      dimension rx(2),ry(2)
c
      data rlon/101.,115./,rlat/0.,1.75/,rotang/-30.,50./
      data rx/300.,400./,ry/800.,700./,topamp/1200.,1500./
      data rxx/40.,80./,ryy/300.,260./,topamp/1200.,1500./
c
******************************************************************************
      radx=6.37e3           ! earth radius
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
c
      xlon(1)=0.
      do 10 i=2,nx
         xlon(i)=xlon(1)+float(i-1)*360./nx
 10   continue
      do 15 j=1,my
         xlat(j)=asin(sinl(j))*r2d
 15   continue                           
c
      do n=1,2      !  two topography
         rang=rotang(n)*d2r
         print *,' topo rlat=',rlat(n),' rlon=',rlon(n),rotang(n)
c        if(n.eq.1)then
c
        rf=5.
        do j=1,my
        do i=1,nx
           dx=radx*cos(rlat(n)*d2r)*(rlon(n)-xlon(i))*d2r
           dy=radx*(rlat(n)-xlat(j))*d2r
           thda=atan2(dy,dx)
           r=(dx*dx+dy*dy)**0.5
c          rmax=((rx(n)*cos(thda+rang))**2+(ry(n)*sin(thda+rang))**2)**0.5
c          rr=r/rmax
           aa=r*cos(thda+rang)/rx(n)
           bb=r*sin(thda+rang)/ry(n)
           rr=(aa**2+bb**2)**0.5

           if(rr.eq.0.)then
              topo(i,j)=1.*topamp(n)
           else if(rr.gt.0. .and. rr.lt.1.)then
              topo(i,j)=(1.-exp(-rf/rr*exp(1./(rr-1.))))*topamp(n)
           else
              if(n.eq.1)topo(i,j)=0.
           endif

        enddo
        enddo
c      else
c        do j=1,my
c        do i=1,nx
c          dx=radx*cos(rlat(n)*d2r)*(rlon(n)-xlon(i))*d2r
c          dy=radx*(rlat(n)-xlat(j))*d2r
c          r=(dx*dx+dy*dy)**0.5
c          thda=atan2(dy,dx)
cc          rmax=((rxx(n)*cos(thda+rang))**2+(ryy(n)*sin(thda+rang))**2)**0.5
cc          rr=r/rmax
c          aa=r*cos(thda+rang)/rxx(n)
c          bb=r*sin(thda+rang)/ryy(n)
c          rr=(aa**2+bb**2)**0.5
c          if(rr.ge.0. .and. rr.lt.2.5)then
c          topo(i,j)=topamp(n)*exp(-((dx/rxx(n))*cos(thda+rang))**2
c     &                            -((dy/ryy(n))*sin(thda+rang))**2)
cc          topo(i,j)=topamp(n)*exp(-rr)
c          else
c           if(n.eq.1)topo(i,j)=0.
c          endif
c        enddo
c        enddo
c
c      endif
      enddo
c
****************************************************************************

      nxmy8=nx*my*8
      open(99,file='../DATAOUT/TOPO3.OUT',
     &     form='unformatted',status='unknown',access='direct',
     &     recl=nxmy8)
      write(99,rec=1)topo
      close(99)

****************************************************************************
c
c      if(1.eq.1)stop
      return
      end
