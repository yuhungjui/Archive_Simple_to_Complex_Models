      subroutine getvorm2(nx,my,lev,dx,vorm,divm)
      dimension vorm(nx,lev,my),divm(nx,lev,my)
c
c      data yj/512.001/,xi/512.001/
      data yj/300.001/,xi/300.001/
c
      data r1/9.5/,r2/52.5/,r3/62.5/,r4/120/
      data d1/2.5/,d2/2.5/,d3/2.5/,d4/15.0/
      data v1/0.015918/,v2/0.000518/,v3/0.002718/,v4/0.000218/
      data v5/-0.000082/
      dy=dx
c
      do 421 j10=1,my
      do 421 i10=1,nx
      vorm(i10,1,j10)=0.
      divm(i10,1,j10)=0.
 421  continue
c                                   
      do 1789 j=1,my
      do 1789 i=1,nx
         vorbog=0.
         s1=0.
         s2=0.
         x=dx*(xi-i)/1000.
         y=dy*(yj-j)/1000.
         thda=atan2(y,x)
         r=(x*x+y*y)**0.5
         if(r.ge.(r4+d4))then
           vorbog=v5
         else if(r.ge.(r4-d4))then
              s1=(r-r4+d4)/(2*d4)
              s2=(r4+d4-r)/(2*d4)
              vorbog=v4*(1-3*s1**2+2*s1**3)+v5*(1-3*s2**2+2*s2**3)
         else if(r.ge.(r3+d3))then
              vorbog=v4
         else if(r.ge.(r3-d3))then
              s1=(r-r3+d3)/(2*d3)
              s2=(r3+d3-r)/(2*d3)
              vorbog=v3*(1-3*s1**2+2*s1**3)+v4*(1-3*s2**2+2*s2**3)
         else if(r.ge.(r2+d2))then
              vorbog=v3
         else if(r.ge.(r2-d2))then
              s1=(r-r2+d2)/(2*d2)
              s2=(r2+d2-r)/(2*d2)
              vorbog=v2*(1-3*s1**2+2*s1**3)+v3*(1-3*s2**2+2*s2**3)
         else if(r.ge.(r1+d1))then
              vorbog=v2
         else if(r.ge.(r1-d1))then
              s1=(r-r1+d1)/(2*d1)
              s2=(r1+d1-r)/(2*d1)
              vorbog=v1*(1-3*s1**2+2*s1**3)+v2*(1-3*s2**2+2*s2**3)
         else
              vorbog=v1
         endif
c
         vorm(i,1,j)=vorbog
 1789 continue
c
      return
      end
