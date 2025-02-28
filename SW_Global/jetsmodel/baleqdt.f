      subroutine baleqdt
c
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
c  working array
c
      dimension gg(nx,lev,my),hh(nx,lev,my),phispc(mlmax,2,lev)
c
c  solve balence eqation to get geopotential
c  first compute vorticity and divergence at grid point
c
c      call trandv(jtrun,mlmax,nx,my,lev,uu,vv,weight,cim
c     *, onocos,poly,dpoly,vspcnow,dspcnow)
c
c      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
c      call transr(jtrun,mlmax,nx,my,lev,poly,dspcnow,div)
c
      do 10 k=1,lev
      do 10 j=1,my
      do 10 i=1,nx
        gg(i,k,j)=-(cor(j)+vor(i,k,j))*vv(i,k,j)
        hh(i,k,j)= (cor(j)+vor(i,k,j))*uu(i,k,j)
 10   continue
c
c  get right hand side of balence equation
c
      call alpha2s(jtrun,mlmax,nx,my,lev,gg,hh,weight,cim
     *, onocos,poly,dpoly,phispc)
c
      do 25 k=1,lev
      do 20 i=2,mlmax
        phispc(i,1,k)=phispc(i,1,k)/eps4(i)
        phispc(i,2,k)=phispc(i,2,k)/eps4(i)
 20   continue
      phispc(1,1,k)=0.
      phispc(1,2,k)=0.
 25   continue
c
      call transr(jtrun,mlmax,nx,my,lev,poly,phispc,phi)
c
c  compute geopotential
c
      do 30 k=1,lev
      do 30 j=1,my
      do 30 i=1,nx
        phi(i,k,j)=phi(i,k,j)-(uu(i,k,j)**2+vv(i,k,j)**2)*0.5
     &             *radsq*onocos(j)
 30   continue
c
      return
      end
