      program plottrack
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      data lon,lat/120,120/
c
      namelist /modlist/tauer,tauor,dtr,dx,dy,tfiltr,hmean
     *,                bplane,blat,hfiltr,ckd,dohd,linear
c
      dimension p(nxr,myr),u(nxr,myr),v(nxr,myr)
     1,         vo(nxr,myr),di(nxr,myr),s(nxr,myr)
      dimension xw(im,jm),xwu(im,jm),xwv(im,jm)
      real kx(ll),ky(ll),pmax
      character*4 head
      character*10 lab1
c
c***  idatasource=  0   ,  read from PV.OUT
c***               else ,  read from track.dat
c
      idatasource=1
c
      if (idatasource.eq.0)then
          go to 1915
      else
          go to 1916
      endif
c
 1915 nxmy=nxr*myr*8
      open(25,file='../dat/PV.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(41,file='../dat/TOPO.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      go to 1917
c
 1916 nxmy=nxr*myr*8
      open(16,file='../dat/track.dat',form='unformatted'
     *,status='unknown')
      open(41,file='../dat/TOPO.OUT',form='unformatted'
     *     ,access='direct',recl=nxmy,status='unknown')
      do i=1,ll
      read(16)kx(i),ky(i)
      enddo
 1917 continue
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
      head='TOPO'
      indexplot=1
      tinc=1000.0
      if (indexplot.eq.1)then
      print *,'Plotting Topography'
      call datacut(nxr,myr,topo,lon,lat,im,jm,xw)
      do i=1,im
         do j=1,jm
            xw(i,j)=xw(i,j)/9.80616
         enddo
      enddo
      call set(0.2,0.8,0.2,0.8,0.,1.,0.,1.,1)
c-0330 call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      call label(head,itau,lab1)
c      call plchhq(0.5,0.85,lab1,0.012,0.,0.)
c-0330 call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call set(0.2,0.8,0.2,0.8,0.,1.,0.,1.,1)
c      call plchhq(0.32,0.85,'F  LINEAR',0.015,0.,0.)
      call plchhq(0.55,0.9,'F  NONLINEAR',0.015,0.,0.)
      call gslwsc(3.0)
      call plchhq(0.06,0.55,'KM',0.012,0.,0.)
      call plchhq(0.55,0.12,'KM',0.012,0.,0.)
c
      call set(0.35,0.7,0.35,0.7,0.,3000.,0.,3000.,1)
      call gslwsc(1.0)
      CALL CONREC(xw,im,im,jm,0.,0.,tinc,-1,0,0)
      call gslwsc(3.0)
      call labmod('(i4)','(i4)',4,4,14,14,0,0,0)
      call periml(3,5,3,5)
      endif
c
      itauo=0
      nn=1
c
      if (idatasource.eq.0)then
  51  nrec=itauo+1
      if (nrec.gt.((tauer/tauor)+1.001))go to 155
      itau=(nrec-1)*tauor
      print *,'Reading tau=',itau
      read(25,rec=nrec)s
c
      itauo=nrec
      pmax=0
c
      print *,'Judging PV maximum position at tau=',itau
      call datacut(nxr,myr,s,lon,lat,im,jm,xw)
      do i=1,im
         do j=1,jm
            xw(i,j)=xw(i,j)*1.e8
            if (xw(i,j).gt.pmax) then
               pmax=xw(i,j)
               kx(nn)=i
               ky(nn)=j
            endif
         enddo
      enddo
c
      print *,nn
      nn=nn+1
c
      go to 51
 155  continue
c
      endif
c
      do i=1,ll
      print *,i,kx(i),ky(i)
      enddo
c
      call set(0.35,0.7,0.35,0.7,1.,61.,1.,61.,1)
      call gslwsc(4.0)
      call curve(kx,ky,ll)
c      call frame
c
c
      if (idatasource.eq.0)then
         close(25)
      else
         close(16)
      endif
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
      subroutine label(head,itau,lab1)
c
      character*4 head
      character*10 lab1
c
      if (itau.eq.0)then
         write(lab1,1001)head
 1001    format(a4,1x,'000HR')
      else
          write(lab1,1000)head,itau
 1000     format(a4,1x,i3.3,'HR')
      endif
c
      return
      end
