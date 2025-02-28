      subroutine rgetrdy
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      include '../include/rspec.h'
      include '../include/rfftcom.h'
c
      call gettopo(nxr,myr,topo)
c
      call getum(nxr,myr,levr,um,vm)
c      call rtrandv(jtrunx,jtruny,nxr,myr,levr,pir,dox,doy
c     1,            um,vm,vormspc,divmspc)
c      call rtransr(jtrunx,jtruny,nxr,myr,levr,vormspc,vorm)
c      call rtransr(jtrunx,jtruny,nxr,myr,levr,divmspc,divm)
c      call rtranuv(jtrunx,jtruny,nxr,myr,levr,pir,dox,doy,cmn
c     1            ,vormspc,divmspc,um,vm)
      call nbalm
      do 90 j=1,myr
      do 90 i=1,nxr
        ttm(i,1,j)=ttm(i,1,j)+hmean*gravr
cc        utr(i,1,j)=um(i,1,j)
c        vtr(i,1,j)=vm(i,1,j)
 90   continue
c
c      call getvorm1(nxr,myr,levr,dx,vorm,divm)
c      call rtranrs(jtrunx,jtruny,nxr,myr,levr,vorm,vormspc)
c      call rtranrs(jtrunx,jtruny,nxr,myr,levr,divm,divmspc)
c      call rtranuv(jtrunx,jtruny,nxr,myr,levr,pir,dox,doy,cmn
c     1            ,vormspc,divmspc,utr,vtr)
c      call bogusty(nxr,myr,dx,utr,vtr)
      call bogusuv(nxr,myr,dx,utr,vtr)
c      call bogusuv2(nxr,myr,dx,utr,vtr)
c      call bogusbanduv(nxr,myr,dx,utr,vtr)
      call rtrandv(jtrunx,jtruny,nxr,myr,levr,pir,dox,doy
     1,            utr,vtr,vornowr,divnowr)
      call rtransr(jtrunx,jtruny,nxr,myr,levr,vornowr,rvorr)
      call rtransr(jtrunx,jtruny,nxr,myr,levr,divnowr,rdivr)
      call nbal
c
      if(.not. linear)then
        do 10 j=1,myr
        do 10 i=1,nxr
        utr(i,1,j)=utr(i,1,j)+um(i,1,j)
        vtr(i,1,j)=vtr(i,1,j)+vm(i,1,j)
        ttr(i,1,j)=ttr(i,1,j)+ttm(i,1,j)
 10     continue
      endif
c
      call rtranrs(jtrunx,jtruny,nxr,myr,levr,ttr,temnowr)
      call rtrandv(jtrunx,jtruny,nxr,myr,levr,pir,dox,doy
     1,            utr,vtr,vornowr,divnowr)
      call rtransr(jtrunx,jtruny,nxr,myr,levr,vornowr,rvorr)
      call rtransr(jtrunx,jtruny,nxr,myr,levr,divnowr,rdivr)
c compute pv
c      if(.not. linear)then
c      do 212 j=1,myr
c      do 212 i=1,nxr
c        pvr(i,1,j)=(rvorr(i,1,j)+corr(i,j))/ttr(i,1,j)
c 212  continue
c      else
      do 213 j=1,myr
      do 213 i=1,nxr
        pvr(i,1,j)=(rvorr(i,1,j)+corr(i,j))/(ttr(i,1,j)+hmean*gravr)
 213  continue
c      endif
c
      do 210 ml=1,jtrunxy2*2*levr
        voroldr(ml,1,1)=vornowr(ml,1,1)
        divoldr(ml,1,1)=divnowr(ml,1,1)
        temoldr(ml,1,1)=temnowr(ml,1,1)
 210  continue
c
      call outflds(0.)
c
      return
      end
