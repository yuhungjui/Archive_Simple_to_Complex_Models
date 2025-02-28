c
      subroutine regrhs
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      include '../include/rspec.h'
c
c  working array
c
      dimension uu(nxr,levr,myr),vv(nxr,levr,myr)
     1,         w1(nxr,myr),ws1(jtrunxy2,2,levr)
     *,         tmpspc(jtrunxy2,2,levr)
c
      max2lev=jtrunxy2*2*levr
c
      if(.not. linear)then
c vorticity equation
      do 82 j=1,myr
      do 82 i=1,nxr
        uu(i,1,j)=-um(i,1,j)*rvorr(i,1,j)
     &            -utr(i,1,j)*(corr(i,j)+rvorr(i,1,j))
     &            -utr(i,1,j)*rvorr(i,1,j)
        vv(i,1,j)=-vm(i,1,j)*rvorr(i,1,j)
     &            -vtr(i,1,j)*(corr(i,j)+rvorr(i,1,j))
     &            -vtr(i,1,j)*rvorr(i,1,j)
 82   continue
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,uu,vv,vortenr)
c divergence equation
      do 81 j=1,myr
      do 81 i=1,nxr
        w1(i,j)=(utr(i,1,j)**2+vtr(i,1,j)**2)*0.5+ttr(i,1,j)
        vv(i,1,j)=-1.0*vv(i,1,j)
 81    continue
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,vv,uu,divtenr)
      call rtranrs(jtrunx,jtruny,nxr,myr,1,w1,ws1(1,1,1))
      do 83 mn= 1, jtrunxy
        ms=mn+jtrunxy
        ws1(mn,1,1)=-ws1(mn,1,1)*cmn(mn)
        ws1(mn,2,1)=-ws1(mn,2,1)*cmn(mn)
        ws1(ms,1,1)=-ws1(ms,1,1)*cmn(mn)
        ws1(ms,2,1)=-ws1(ms,2,1)*cmn(mn)
 83   continue
      do 84 i=1,jtrunxy2*2*levr
        divtenr(i,1,1)=divtenr(i,1,1)-ws1(i,1,1)
 84   continue
C    Calculate spectral tendecny for MASS eqn
      do 85 j=1,myr
      do 85 i=1,nxr
        uu(i,1,j)=-ttr(i,1,j)*um(i,1,j)
     &            -(ttm(i,1,j)-topo(i,j))*utr(i,1,j)
     &            -ttr(i,1,j)*utr(i,1,j)
        vv(i,1,j)=-ttr(i,1,j)*vm(i,1,j)
     &            -(ttm(i,1,j)-topo(i,j))*vtr(i,1,j)
     &            -ttr(i,1,j)*vtr(i,1,j)
 85   continue
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,uu,vv,temtenr)
c
      else
c vorticity equation
      do 20 j=1,myr
      do 20 i=1,nxr
        uu(i,1,j)=-um(i,1,j)*rvorr(i,1,j)
     &            -utr(i,1,j)*(corr(i,j)+rvorr(i,1,j))
        vv(i,1,j)=-vm(i,1,j)*rvorr(i,1,j)
     &            -vtr(i,1,j)*(corr(i,j)+rvorr(i,1,j))
 20   continue
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,uu,vv,vortenr)
c divergence equation
      do 21 j=1,myr
      do 21 i=1,nxr
        w1(i,j)=(utr(i,1,j)**2+vtr(i,1,j)**2)*0.5+ttr(i,1,j)
        vv(i,1,j)=-1.0*vv(i,1,j)
 21    continue
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,vv,uu,divtenr)
      call rtranrs(jtrunx,jtruny,nxr,myr,1,w1,ws1(1,1,1))
      do 22 mn= 1, jtrunxy
        ms=mn+jtrunxy
        ws1(mn,1,1)=-ws1(mn,1,1)*cmn(mn)
        ws1(mn,2,1)=-ws1(mn,2,1)*cmn(mn)
        ws1(ms,1,1)=-ws1(ms,1,1)*cmn(mn)
        ws1(ms,2,1)=-ws1(ms,2,1)*cmn(mn)
 22   continue
      do 23 i=1,jtrunxy2*2*levr
        divtenr(i,1,1)=divtenr(i,1,1)-ws1(i,1,1)
 23   continue
C    Calculate spectral tendecny for MASS eqn
      do 24 j=1,myr
      do 24 i=1,nxr
        uu(i,1,j)=-ttr(i,1,j)*um(i,1,j)
     &            -(ttm(i,1,j)-topo(i,j))*utr(i,1,j)
        vv(i,1,j)=-ttr(i,1,j)*vm(i,1,j)
     &            -(ttm(i,1,j)-topo(i,j))*vtr(i,1,j)
 24   continue
      call ralpha2s( nxr,myr,levr,jtrunx,jtruny,jtrunxy2
     1            , pir,dox,doy,uu,vv,temtenr)
c
      endif
c
c add friction
      do 31 j=1,myr
      do 31 i=1,nxr
        uu(i,1,j)=ck(i,j)*rvorr(i,1,j)
        vv(i,1,j)=ck(i,j)*rdivr(i,1,j)
 31   continue
      call rtranrs(jtrunx,jtruny,nxr,myr,levr,uu,tmpspc)
      call rtranrs(jtrunx,jtruny,nxr,myr,levr,vv,ws1)
      do 32 m=1,max2lev
          vortenr(m,1,1)=vortenr(m,1,1)-tmpspc(m,1,1)
          divtenr(m,1,1)=divtenr(m,1,1)-ws1(m,1,1)
 32   continue
      do 41 j=1,myr
      do 41 i=1,nxr
        uu(i,1,j)=ck(i,j)*ttr(i,1,j)
 41   continue
      call rtranrs(jtrunx,jtruny,nxr,myr,levr,uu,tmpspc)
      do 42 m=1,max2lev
          temtenr(m,1,1)=temtenr(m,1,1)-tmpspc(m,1,1)
 42   continue
c
      return
      end
