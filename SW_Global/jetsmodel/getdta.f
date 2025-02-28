      subroutine getdta(nx,my,lev,ut,vt,phi,sinl,rad,cosl)
      dimension ut(nx,lev,my),vt(nx,lev,my),phi(nx,lev,my)
     &,         sinl(my)
c local array
      dimension xlon(nx),xlat(my)
      data rlon/115./, rlat/32./, rmax/500./, hamp/300./
c
      radx=6.37e3
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
      ix=rlon/(360./nx)+1.001
      jy=(rlat-xlat(1))/(180./my)+1.001
      print *,'h  ix=',ix,' jy=',jy,' rlat=',rlat,' rlon=',rlon
c
      jb=jy-30
      je=jy+30
      ib=ix-30
      ie=ix+30
c
      do 200 j=1,my
      do 200 i=1,nx
        ut(i,1,j)=0.
        vt(i,1,j)=0.
        phi(i,1,j)=0.
 200  continue
c
      a=rmax
      b=a
      do 101 j=jb,je
      do 101 i=ib,ie
c      do 101 j=1,my
c      do 101 i=1,nx
        dx=radx*cos(rlat*d2r)*(rlon-xlon(i))*d2r
        dy=radx*(rlat-xlat(j))*d2r
        phi(i,1,j)=hamp*exp(-(dx/a)**2-(dy/b)**2)*9.806
        cc=exp(-(dx/a)**2-(dy/b)**2)
 101  continue
      print *,'phi(ix,1,jy)=',phi(ix,1,jy),hamp
c
      return
      end
