      program plot1
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
cleon      data lon,lat/120,120/
cleon      data lon,lat/512,512/
      data lon,lat/275,325/
c
      namelist /modlist/tauer,tauor,dtr,dx,dy,tfiltr,hmean
     *,                bplane,blat,hfiltr,ckd,dohd,linear
c
      dimension p(nxr,myr),u(nxr,myr),v(nxr,myr)
     1,         vo(nxr,myr),di(nxr,myr),s(nxr,myr)
      dimension xw(im,jm),xwu(im,jm),xwv(im,jm)
      real kx(ll),ky(ll)
      character*4 head
      character*11 lab1
      character*2 xlab,ylab
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
      open(41,file='../dat/TOPO.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      open(16,file='../dat/track.dat',form='unformatted'
     *,status='unknown')
c
      open(1,file='../rnamlsts',status='unknown')
      read(1,modlist,end=5)
 5    continue
      print modlist
c
      call opngks
c
      read(41,rec=1)topo
      close(41)
c
      itauo=0
      nn=1
  51  nrec=itauo+1
      if (nrec.gt.((tauer/tauor)+1.001))go to 155
c      itau=(nrec-1)*tauor
      titau=(nrec-1)*tauor
      print *,' '
      print *,'#################################'
      print *,' '
      print *,'Reading          tau=',titau
      read(21,rec=nrec)p
      read(22,rec=nrec)u
      read(23,rec=nrec)v
      read(24,rec=nrec)vo
      read(25,rec=nrec)s
      read(26,rec=nrec)di
      itauo=nrec
c
      xlab='KM'
      ylab='KM'
c
      head='TOPO'
      indexplot=1
      tinc=500.0
      if (indexplot.eq.1)then
      print *,'Plotting Topography'
      call datacut(nxr,myr,topo,lon,lat,im,jm,xw)
      do i=1,im
         do j=1,jm
            xw(i,j)=xw(i,j)/9.80616
         enddo
      enddo
      call set(0.15,0.9,0.15,0.9,-200.,200.,-200.,200.,1)
      call gslwsc(2.0)
      CALL CONREC(xw,im,im,jm,0.,0.,tinc,-1,0,0)
      call labmod('(i4)','(i4)',4,4,14,14,0,0,0)
      call periml(15,2,15,2)
      endif
c
      head='GEOP'
      indexplot=0
      finc=500.
      if (indexplot.eq.1)then
      print *,'Plotting GEOP at tau=',titau
      call datacut(nxr,myr,p,lon,lat,im,jm,xw)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call label(head,titau,lab1)
      call plchhq(0.52,0.75,lab1,0.012,0.,0.)
      call set(0.15,0.9,0.15,0.9,0.,3000.,0.,3000.,1)
      call gslwsc(1.0)
      CALL CONREC(xw,im,im,jm,0.,0.,finc,-1,0,0)
      endif
c
      head=' PV '
      indexplot=1
      finc=20.0
      if (indexplot.eq.1)then
      print *,'Plotting  PV  at tau=',titau
      call datacut(nxr,myr,s,lon,lat,im,jm,xw)
      pmax=0
      do i=1,im
         do j=1,jm
            xw(i,j)=xw(i,j)*10e8
            if (xw(i,j).gt.pmax)then
               pmax=xw(i,j)
c               print *, pmax
               kx(nn)=i
               ky(nn)=j
            endif
         enddo
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call label(head,titau,lab1)
      call plchhq(0.52,0.95,lab1,0.012,0.,0.)
      call plchhq(0.52,0.05,xlab,0.012,0.,0.)
      call plchhq(0.05,0.51,ylab,0.012,0.,0.)
      call set(0.15,0.9,0.15,0.9,-80.,80.,-80.,80.,1)
      call gslwsc(1.0)
c      CALL CONREC(xw,im,im,jm,0.,0.,finc,-1,0,0)
      CALL CPCNRC(xw,im,im,jm,0.,400.,finc,-1,0,0)
      print *,' after conrec'
      endif
c
      nn=nn+1
c
      head='WIND'
      indexplot=0
      if (indexplot.eq.1)then
      print *,'Plotting WIND at tau=',itau
      call datacut(nxr,myr,u,lon,lat,im,jm,xwu)
      call datacut(nxr,myr,v,lon,lat,im,jm,xwv)
      call set(0.15,0.9,0.15,0.9,0.,1.,0.,1.,1)
      call gslwsc(1.0)
      CALL VELVCT(xwu,im,xwv,im,im,jm,5.,30.,-1,40,0,SPV)
      endif
c
      call frame
      go to 51
c
 155  continue
c
      do i=1,ll
         write(16)kx(i),ky(i)
      enddo
c
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
c
      call clsgks
c
      end
c
      subroutine datacut(nx,my,x,ix,jy,im,jm,xw)
      dimension x(nx,my),xw(im,jm)
      ib=ix-20
      ie=ix+20
      jb=jy-20
      je=jy+20
      ib=ix-(im-1)/2
      ie=ix+(im-1)/2
      jb=jy-(jm-1)/2
      je=jy+(jm-1)/2
      do 10 i=ib,ie
      do 10 j=jb,je
      xw(i-ib+1,j-jb+1)=x(i,j)
 10   continue
      return
      end
c
      subroutine label(head,titau,lab1)
c
      character*4 head
      character*11 lab1
c
      if (titau.eq.0)then
         write(lab1,1001)head
 1001    format(a4,2x,'0.0HR')
      else
          write(lab1,1000)head,titau
c 1000     format(a4,1x,i3.3,'HR')
 1000     format(a4,1x,f4.1,'HR')
      endif
c
      return
      end
