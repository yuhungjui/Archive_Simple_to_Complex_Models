      subroutine rcons
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rfftcom.h'
c
      namelist /modlist/tauer,tauor,dtr,dx,dy,tfiltr,hmean
     *,                bplane,blat,hfiltr,ckd,dohd,linear
c
      pir  = 4.0*atan(1.0)
      d2r=pir/180.
c
      open(1,file='../rnamlsts',status='unknown')
      read(1,modlist,end=5)
 5    continue
      print modlist
c
c  initialize ifax,ifay and trigsx,trigsy for rfftmlt routine
c
c for cray
c
      call fftfax (nxr,ifax,trigsx)
      call fftfax (myr,ifay,trigsy)
c
      dox=nxr*dx
      doy=myr*dy
      print *,dox,doy
c
      pidx=2.*pir/dox
      pidy=2.*pir/doy
      do 50 n=1,jtruny
      do 50 m=1,jtrunx
	mn=m+(n-1)*jtrunx
	cmn(mn)=((m-1)*pidx)**2+((n-1)*pidy)**2
 50   continue
c
      if(bplane)then
      jh=myr/2
      f0  =2.0*omegar*sin(blat*d2r)
      beta=2.0*omegar*cos(blat*d2r)/radr
      do 100 j=1,myr
      do 100 i=1,nxr
        corr(i,j)=f0+beta*(j-jh)*dy
 100  continue
      else
      f0  =2.0*omegar*sin(blat*d2r)
      do 102 j=1,myr
      do 102 i=1,nxr
        corr(i,j)=f0
 102  continue
      endif
c
      do 105 j=1,myr
      do 105 i=1,nxr
        ck(i,j)=ckd
 105  continue
c
c      hmean=hmean*gravr
c
      nxmy=nxr*myr*8
      open(21,file='../dat/H.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(22,file='../dat/U.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(23,file='../dat/V.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(24,file='../dat/VOR.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(25,file='../dat/PV.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(26,file='../dat/DIV.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(35,file='../dat/HM.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      open(36,file='../dat/UM.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      open(37,file='../dat/VM.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      open(38,file='../dat/VORM.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      open(40,file='../dat/DIVM.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      open(41,file='../dat/TOPO.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
c
      return
      end
