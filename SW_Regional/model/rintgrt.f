      subroutine rintgrt
c
c  time integration use leapfog scheme
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      include '../include/rspec.h'
c
      dimension strspc(jtrunxy2,2,levr)
      logical forward
      data forward/.true./
c
      max2lev=jtrunxy2*2*levr
c
      iend=tauer*3600./dtr+0.001
      iout=tauor*3600./dtr+0.001
c
      print *,' begin leapfog time integration'
c
      do 100 iter=1,iend
c
      taux=iter*dtr/3600.
      tau=taux
c
      print *,'forcast tau=',taux
c
      call regrhs
c
      k=1
      temtenr(1,1,k)=0.
      temtenr(1,2,k)=0.
      temtenr(1+jtrunxy,1,k)=0.
      temtenr(1+jtrunxy,2,k)=0.
      vortenr(1,1,k)=0.
      vortenr(1,2,k)=0.
      vortenr(1+jtrunxy,1,k)=0.
      vortenr(1+jtrunxy,2,k)=0.
      divtenr(1,1,k)=0.
      divtenr(1,2,k)=0.
      divtenr(1+jtrunxy,1,k)=0.
      divtenr(1+jtrunxy,2,k)=0.
c
      if(forward)then
c
c  leapfog time integration
c
        do 10 i=1,max2lev
          vornowr(i,1,1)=dtr*vortenr(i,1,1)+voroldr(i,1,1)
          divnowr(i,1,1)=dtr*divtenr(i,1,1)+divoldr(i,1,1)
          temnowr(i,1,1)=dtr*temtenr(i,1,1)+temoldr(i,1,1)
 10    continue
c
        dta=dtr*2.0
        forward=.false.
c
      else
c
c  leapfog time integration
c
        do 20 i=1,max2lev
          vortenr(i,1,1)=dta*vortenr(i,1,1)+voroldr(i,1,1)
          divtenr(i,1,1)=dta*divtenr(i,1,1)+divoldr(i,1,1)
          temtenr(i,1,1)=dta*temtenr(i,1,1)+temoldr(i,1,1)
 20     continue
c
c  do robert time filter
c
        do 30 i=1,max2lev
          voroldr(i,1,1)= vornowr(i,1,1)+tfiltr*(voroldr(i,1,1)
     *                   -2.0*vornowr(i,1,1)+vortenr(i,1,1))
          divoldr(i,1,1)= divnowr(i,1,1)+tfiltr*(divoldr(i,1,1)
     *                   -2.0*divnowr(i,1,1)+divtenr(i,1,1))
          temoldr(i,1,1)= temnowr(i,1,1)+tfiltr*(temoldr(i,1,1)
     *                   -2.0*temnowr(i,1,1)+temtenr(i,1,1))
          vornowr(i,1,1)=vortenr(i,1,1)
          divnowr(i,1,1)=divtenr(i,1,1)
          temnowr(i,1,1)=temtenr(i,1,1)
30    continue
c
      endif
      if(dohd)call rhdiffu( dta,myr,nxr,jtrunxy,levr,hfiltr
     1                    , vornowr,divnowr,temnowr
     2                    , cmn)
c
c  transform spectrum coefficient to grid point
c
      call rtransr(jtrunx,jtruny,nxr,myr,levr,vornowr,rvorr)
      call rtransr(jtrunx,jtruny,nxr,myr,levr,divnowr,rdivr)
      call rtransr(jtrunx,jtruny,nxr,myr,levr,temnowr,ttr)
      call rtranuv(jtrunx,jtruny,nxr,myr,levr,pir,dox,doy,cmn
     1            ,vornowr,divnowr,utr,vtr)
c
c compute stream function
c      strspc(1,1,1)=0.
c      strspc(1,2,1)=0.
c      strspc(1+jtrunxy,1,1)=0.
c      strspc(1+jtrunxy,2,1)=0.
c      do 32 mn= 2, jtrunxy
c      ms=mn+jtrunxy
c      strspc(mn,1,1)=-vornowr(mn,1,1)/cmn(mn)
c      strspc(mn,2,1)=-vornowr(mn,2,1)/cmn(mn)
c      strspc(ms,1,1)=-vornowr(ms,1,1)/cmn(mn)
c      strspc(ms,2,1)=-vornowr(ms,2,1)/cmn(mn)
c 32   continue
c      call rtransr(jtrunx,jtruny,nxr,myr,levr,strspc,pvr)
c compute pv
c      if(.not. linear)then
c      do 32 j=1,myr
c      do 32 i=1,nxr
c        pvr(i,1,j)=(rvorr(i,1,j)+corr(i,j))/ttr(i,1,j)
c 32   continue
c      else
      do 213 j=1,myr
      do 213 i=1,nxr
        pvr(i,1,j)=(rvorr(i,1,j)+corr(i,j))/(ttr(i,1,j)+hmean*gravr)
 213  continue
c      endif
c
c      call boundary(nxr,myr,levr,ttr,rvorr,rdivr)
c
      if(mod(iter,iout).eq.0)then
        call outflds(taux)
      endif
c
 100  continue
c
      print *,' finish leapfor time integration'
C
      return
      end
