      subroutine bogusuv2(nx,my,dx,ut,vt)
      dimension ut(nx,1,my),vt(nx,1,my)
c
      data b/0.7/,rm2/2000./,yj/112.501/,xi/97.501/,vm/5./,rvm/100./
c1024      data b/1.0/,rm2/2000./,yj/112.501/,xi/127.501/,vm/10./,rvm/100./
c      data b/1.0/,rm2/2000./,yj/112.501/,xi/120.001/,vm/30./,rvm/100./
c
      dy=dx
c
c      do 90 j=1,my
c      do 90 i=1,nx
c        ut(i,1,j)=0.
c        vt(i,1,j)=0.
c 90   continue
c
      do 100 j=1,my
      do 100 i=1,nx
        x=dx*(xi-i)/1000.
        y=dy*(yj-j)/1000.
        thda=atan2(y,x)
        r=(x*x+y*y)**0.5
        if(r.le.rm2)then
          fac=r/rm2
          ubog= (vm*(r/rvm)*exp((1.-(r/rvm)**b)/b))*sin(thda)
          vbog=-(vm*(r/rvm)*exp((1.-(r/rvm)**b)/b))*cos(thda)
c
          ut(i,1,j)=ut(i,1,j)+ubog
          vt(i,1,j)=vt(i,1,j)+vbog
c          ut(i,1,j)=ut(i,1,j)*fac+ubog*(1.-fac)
c          vt(i,1,j)=vt(i,1,j)*fac+vbog*(1.-fac)
        endif
 100  continue
c
      return
      end
