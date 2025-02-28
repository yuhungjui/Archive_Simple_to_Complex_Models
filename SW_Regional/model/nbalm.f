      subroutine nbalm
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      include '../include/rspec.h'
      include '../include/rfftcom.h'
c
      dimension w1(nxr,levr,myr),w2(nxr,levr,myr)
     1,         ws(jtrunxy2,2,levr)
c
      do 100 j=1,myr
      do 100 i=1,nxr
        w1(i,1,j)= (corr(i,j)+vorm(i,1,j))*vm(i,1,j)
        w2(i,1,j)=-(corr(i,j)+vorm(i,1,j))*um(i,1,j)
 100  continue
c
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,w1,w2,ws)
c
      do 53 k=1,levr
      do 50 mn= 2, jtrunxy
      ms=mn+jtrunxy
      ws(mn,1,k)=-ws(mn,1,k)/cmn(mn)
      ws(mn,2,k)=-ws(mn,2,k)/cmn(mn)
      ws(ms,1,k)=-ws(ms,1,k)/cmn(mn)
      ws(ms,2,k)=-ws(ms,2,k)/cmn(mn)
 50   continue
c      ws(1,1,k)=0.
c      ws(1,2,k)=0.
c      ws(1+jtrunxy,1,k)=0.
c      ws(1+jtrunxy,2,k)=0.
 53   continue
c
      call rtransr(jtrunx,jtruny,nxr,myr,levr,ws,ttm)
c
      do 120 j=1,myr
      do 120 i=1,nxr
        ttm(i,1,j)=ttm(i,1,j)-(um(i,1,j)**2+vm(i,1,j)**2)*0.5
c     &             +sgeor(i,j)  
 120  continue
c
      return
      end
