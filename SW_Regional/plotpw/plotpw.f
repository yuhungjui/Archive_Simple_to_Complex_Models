      program plotpvwind
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      parameter (ncl=6)
      parameter (ipror=50)
      data lon,lat/450,450/
      data umeanflow,vmeanflow/-0.0,0.0/
c      data iprox,jproy/275,350/
      data iprox,jproy/250,350/
c
      namelist /modlist/tauer,tauor,dtr,dx,dy,tfiltr,hmean
     *,                bplane,blat,hfiltr,ckd,dohd,linear
c
      dimension u(nxr,myr),v(nxr,myr),s(nxr,myr)
      dimension xw(im,jm),xwu(im,jm),xwv(im,jm)
      real kx(ll),ky(ll),rd(ncl),rdtopo(ncl)
      integer ishade(ncl+1),ishadetopo(ncl+1)
      character*4 head
      character*11 lab1
      character*2 xlab,ylab
      character*4 llb(ncl+2),llbtopo(ncl+2)
      data rd/10.,20.,30.,40.,50.,60./
      data rdtopo/100.,1000.,1500.,2000.,2500.,3000./
      data ishade/0,9,10,11,12,13,17/
      data ishadetopo/1,3,4,5,6,7,8/
c
      nxmy=nxr*myr*8
      open(22,file='../dat/U.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(23,file='../dat/V.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(25,file='../dat/PV.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
c      open(41,file='../dat/taiwan.out',form='unformatted'
c     *     ,access='direct',recl=nxmy,status='unknown')
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
      call dfclrs
      do i=1,ncl
         write(llb(i+1),'(i4)')int(rd(i))
      enddo
      llb(1)='    '
      llb(ncl+2)='    '
      do i=1,ncl
         write(llbtopo(i+1),'(i4)')int(rdtopo(i))
      enddo
      llbtopo(1)='    '
      llbtopo(ncl+2)='    '
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
      read(22,rec=nrec)u
      read(23,rec=nrec)v
      read(25,rec=nrec)s
      itauo=nrec
c
      xlab='KM'
      ylab='KM'
c
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call bckclr(0,0,1)
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
      call plchhq(0.87,0.82,'UNIT M',0.012,0.,0.)
      call set(0.15,0.8,0.15,0.8,-300.,300.,-300.,300.,1)
      call gslwsc(2.0)
c      CALL CONREC(xw,im,im,jm,0.,0.,tinc,-1,0,0)
      call pl2d(0.15,0.8,0.15,0.8,xw,im,jm,rdtopo,ncl,0,
     +          ishadetopo,0,-9999)
      call set(0.15,0.8,0.5,0.8,-300.,300.,-300.,300.,1)
      call plotbar(llbtopo,3.3,.3,1.0)
      endif
c--------------------------------
      head='TOPO'
      indexplot=0
      tinc=1000.0
      if (indexplot.eq.1)then
      print *,'Plotting Topography'
      call datacut(nxr,myr,topo,lon,lat,im,jm,xw)
      do i=1,im
         do j=1,jm
            xw(i,j)=xw(i,j)/9.80616
            if (xw(i,j).gt.pmax)pmax=xw(i,j)
         enddo
      enddo
      call set(0.15,0.8,0.15,0.8,-300.,300.,-300.,300.,1)
      call cpsetc('ILT',' ')
      call cpcnrc(xw,im,im,jm,100.,4000.,1000.,-1,0,0)
      endif
c--------------------------------
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
      call gslwsc(1.0)
      call label(head,titau,lab1)
      call plchhq(0.48,0.75,lab1,0.012,0.,0.)
      call plchhq(0.48,0.05,xlab,0.012,0.,0.)
      call plchhq(0.05,0.48,ylab,0.012,0.,0.)
      call plchhq(0.90,0.48,'UNIT 10e-8 1/s',0.012,0.,0.)
      call set(0.15,0.8,0.15,0.8,-300.,300.,-300.,300.,1)
      call gslwsc(1.0)
c      CALL CONREC(xw,im,im,jm,0.,0.,finc,-1,0,0)
c      CALL CPCNRC(xw,im,im,jm,0.,400.,finc,-1,0,0)
      call pl2d(0.15,0.8,0.15,0.8,xw,im,jm,rd,ncl,0,
     +          ishade,0,-9999)
c      if (nrec.eq.1)then
         xx=float(iprox-lon+(im-1)/2+1)
         yyb=float(jproy-ipror-lat+(jm-1)/2+1)
         yyt=float(jproy+ipror-lat+(jm-1)/2+1)
         call gslwsc(3.0)
c         call line(xx,yyb,xx,yyt)
c      endif
      call set(0.15,0.8,0.15,0.45,-300.,300.,-300.,300.,1)
      call plotbar(llb,3.3,.3,1.0)
      call set(0.15,0.8,0.15,0.8,-300.,300.,-300.,300.,1)
      call gslwsc(1.0)
      call labmod('(i4)','(i4)',4,4,14,14,0,0,0)
      call periml(6,2,6,2)
      print *,' after pl2d'
      endif
c
      head='WIND'
      indexplot=1
      if (indexplot.eq.1)then
      print *,'Plotting WIND at tau=',titau
      call datacut(nxr,myr,u,lon,lat,im,jm,xwu)
      call datacut(nxr,myr,v,lon,lat,im,jm,xwv)
      call set(0.15,0.8,0.15,0.8,1.,float(im),1.,float(jm),1)
      call gslwsc(1.0)
c      CALL VELVCT(xwu,im,xwv,im,im,jm,5.,30.,-1,40,0,SPV)
c
      i=0
      j=0
      do i=16,im-15,15
         do j=16,jm-15,15
            xwu(i,j)=xwu(i,j)+umeanflow
            xwv(i,j)=xwv(i,j)+vmeanflow
            call windbar(float(i),float(j),
     +                   xwu(i,j),xwv(i,j),1,1.5,199)
         enddo
      enddo
c
      endif
c
      nn=nn+1
c
      call frame
      go to 51
c
 155  continue
c
      print *,'TOPO Max. = ',pmax
      do i=1,ll
         write(16)kx(i),ky(i)
      enddo
      close(16)
c
      close(22)
      close(23)
      close(25)
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
