      subroutine bogusuv(nx,my,lev,np,ntyp,rad,rlat,rlon,sinl
     *,                  vm,rvm,rm2,b,ut,vt,center,ix,jy,cosl)
      dimension ut(nx,lev,my),vt(nx,lev,my),center(5000,2,np)
     *, rlat(np),rlon(np),vm(np),rvm(np),ix(np),jy(np)
     *, sinl(my),cosl(my)
c
c  working array
c
      dimension xlon(nx),xlat(my)
c
      radx=rad/1000.
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
      do 20 n=1,ntyp
c
      print *,'vm=',vm(n),' rvm=',rvm(n),' rm2=',rm2,' b=',b
      print *,'rlon=',rlon(n),' rlat=',rlat(n)
c
      ix(n)=rlon(n)/(360./nx)+1.001
      jy(n)=(rlat(n)-xlat(1))/(180./my)+1.001
      print *,'ix=',ix(n),' jy=',jy(n)
      center(1,1,n)=ix(n)
      center(1,2,n)=jy(n)
c      center(1,1,n)=rlon(n)
c      center(1,2,n)=rlat(n)
      jb=jy(n)-20
      je=jy(n)+20
      ib=ix(n)-20
      ie=ix(n)+20
      do 100 j=jb,je
      xx=cosl(j)/rad
      do 100 i=ib,ie
        dx=radx*cos(rlat(n)*d2r)*(rlon(n)-xlon(i))*d2r
        dy=radx*(rlat(n)-xlat(j))*d2r
        thda=atan2(dy,dx)
        r=(dx*dx+dy*dy)**0.5
        if(r.le.rm2)then
          fac=r/rm2
          ubog= (vm(n)*(r/rvm(n))*exp((1.-(r/rvm(n))**b)/b))
     *         *sin(thda)
          vbog=-(vm(n)*(r/rvm(n))*exp((1.-(r/rvm(n))**b)/b))
     *         *cos(thda)

          ut(i,lev,j)=ut(i,lev,j)*fac+ubog*(1.-fac)*xx
          vt(i,lev,j)=vt(i,lev,j)*fac+vbog*(1.-fac)*xx
        endif
 100  continue
c
  20  continue
c
      return
      end
