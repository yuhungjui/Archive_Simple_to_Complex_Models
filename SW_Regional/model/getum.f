c
      subroutine getum(nx,my,lev,um,vm)
      dimension um(nx,lev,my),vm(nx,lev,my)
c      data u0,v0/-2.0,4.0/
      data u0,v0/-1.0,2.0/
      pi=4.*atan(1.0)
      pi2=2.*atan(1.0)
      do 20 j=1,my
      do 20 i=1,nx
        um(i,1,j)=u0
        vm(i,1,j)=v0
 20   continue
      return
      end
