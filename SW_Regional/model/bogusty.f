      subroutine bogusty(nx,my,dx,ut,vt)
      dimension ut(nx,1,my),vt(nx,1,my)
c
      data yj/120.001/,xi/120.001/
c
      data v1/0.000518/,v2/0.002718/,v3/0.000218/
      data r1/52.5/,r2/62.5/,r3/120/
      dy=dx
c
      do 93 j=1,my
      do 93 i=1,nx
         ut(i,1,j)=0.
         vt(i,1,j)=0.
 93   continue
                                   
      do 1789 j=1,my
      do 1789 i=1,nx
         ubog=0.
         vbog=0.
         x=dx*(xi-i)/1000.
         y=dy*(yj-j)/1000.
         thda=atan2(y,x)
         r=(x*x+y*y)**0.5
         if(r.ge.r3)then
c            ubog=0.5*(-1*(v2-v1)*(r1/r)**2-(v3-v2)*(r2/r)**2
c     &           +v3*(r3/r)**2)*r*1000.*sin(thda)
c            vbog=0.5*(-1*(v2-v1)*(r1/r)**2-(v3-v2)*(r2/r)**2
c     &           +v3*(r3/r)**2)*(-1)*r*1000.*cos(thda)
         else if(r.ge.r2)then
              ubog=0.5*(v3-(v2-v1)*(r1/r)**2-(v3-v2)*(r2/r)**2)
     &             *r*1000.*sin(thda)
              vbog=0.5*(v3-(v2-v1)*(r1/r)**2-(v3-v2)*(r2/r)**2)
     &             *(-1)*r*1000.*cos(thda)
         else if(r.ge.r1)then
              ubog=0.5*(v2-(v2-v1)*(r1/r)**2)*r*1000.*sin(thda)
              vbog=0.5*(v2-(v2-v1)*(r1/r)**2)*(-1)*r*1000.*cos(thda)
         else
              ubog=0.5*v1*r*1000.*sin(thda)
              vbog=0.5*v1*(-1)*r*1000.*cos(thda)
         endif
c
      ut(i,1,j)=ubog
      vt(i,1,j)=vbog
 1789 continue
c
      return
      end
