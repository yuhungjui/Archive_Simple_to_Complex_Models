      subroutine gettopo(nxr,myr,topo)
c
      dimension topo(nxr,myr)
c
      data clon/120./,clat/120./,tradx/12./,trady/14./,tamp/3500./
c
c  ktopo=0: elips or circle
c       =2: strip
      data ktopo/2/
c
      do 110 i=1,nxr*myr
      topo(i,1)=0.
 110  continue
c
      if(ktopo.eq.0)then
c- elips
      ixy=50
      ix=clon
      jy=clat
      print *,'clat,clon=',clat,clon
      a=tradx
      b=trady
      do 101 j=jy-ixy,jy+ixy
      do 101 i=ix-ixy,ix+ixy
        x=1.0*(i-clon)
        y=1.0*(j-clat)
        topo(i,j)=tamp-((tamp/a)*sqr(x*x+y*y))*9.806
c        topo(i,j)=tamp*exp(-(x/a)**2-(y/b)**2)*9.806
 101  continue
c
      else if(ktopo.eq.2)then
c- strip
      wide=trady/2.
      pir2=4.*atan(1.0)/2.
      do 8 j=1,myr
      if(abs(j-clat).le.wide)then
      do 6 i=1,nxr
          topo(i,j)=tamp*cos(pir2*(j-clat)/wide)*9.806
 6    continue
      endif
 8    continue
c
      endif
      call qmaxn3(topo,'regional','output sgeo',1,1,1,nxr,myr,1)
      write(41,rec=1)topo
      close(41)
c
      return
      end
