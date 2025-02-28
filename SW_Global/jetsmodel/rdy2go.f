      subroutine rdy2go

      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'

      real sum(my)   ! llming

      dimension xlon(nx),xlat(my),rforce(nx,my)

      integer x,temp(mlmax)    !  llming
      real ran1(mlmax),ran2(mlmax)       !  llming
      real vorg(nx,lev,my),w2(nx,my)
      real vorg1(nx,lev,my),vorg2(nx,lev,my)

      nxmy=nx*my

c---------------------------------------llming

      open (43,file='1random.dat',status='unknown',
     +      form='formatted')  ! llming
      open (44,file='2random.dat',status='unknown',
     +      form='formatted')  ! llming

      do i=1,mlmax,1
         vspcnow(i,1,1)=0.
         vspcnow(i,2,1)=0.
      enddo

c------------------------------------------------------------

      do i=1,451,1
         read(43,*)ran1(i)
         read(44,*)ran2(i)
      enddo

      x=0

      do i=36,46,1
         do j=1,i,1
            temp(j)=mlsort(j,i)
c            vspcnow(temp(j),1,1)=ran1(x+j)*(2.72e-6)-(1.36e-6)  !
c            vspcnow(temp(j),2,1)=ran2(x+j)*(2.72e-6)-(1.36e-6)  !
            vspcnow(temp(j),1,1)=ran1(x+j)*(5.e-6)-(2.5e-6)  !
            vspcnow(temp(j),2,1)=ran2(x+j)*(5.e-6)-(2.5e-6)  !

c            vspcnow(temp(j),1,1)=ran1(x+j)*(6.e-6)-(3.e-6)   ! del16_2
c            vspcnow(temp(j),2,1)=ran2(x+j)*(6.e-6)-(3.e-6)
         enddo
            x=x+i
            print*,'x=',x
      enddo

      do j=1,my,1              ! llming
         do i=1,nx,1
            div(i,1,j)=0.
         enddo
      enddo

      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,div,dspcnow)

c------------------------------------------------------------

      if(bogus) then
        call bogusuv(nx,my,lev,np,ntyp,rad,rlat,rlon,sinl,vm,rvm,rm2
     &              ,bb,uu,vv,center,ix,jy,cosl)
        print *,'aft bogusuv'
      endif

c------------------------------------------------------------

      call getvor1(nx,my,lev,vorg1,div,sinl)
      call getvor2(nx,my,lev,vorg2,div,sinl)

      do i=1,nx,1
         do j=1,my,1
            vorg(i,1,j)=vorg1(i,1,j)+vorg2(i,1,j)
         enddo
      enddo

      open(21,file='../DATAOUT/VORG.OUT',access='direct'
     &, form='unformatted',recl=nx*my*8,status='unknown')

      do 100 j=1,my,1
      do 100 i=1,nx,1
         w2(i,j)=vorg(i,1,j)
 100  continue
      write(21,rec=nrec)w2

c---------------------------------------------------------------

      do i=1,nx,1
         do j=1,my,1
            vor(i,1,j)=vor(i,1,j)+vorg(i,1,j)
         enddo
      enddo

      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,vor,vspcnow)
      call tranuv(jtrun,mlmax,nx,my,lev,radsq,onocos,eps4,cim
     *                 ,poly,dpoly,vspcnow,dspcnow,uu,vv)

      call baleqdt
c      print *,'baledt phi=',phi(307,1,325)/grav

c-------------------------------------------------------------








c--------------------------------------------------llming    go

      do j=1,my,1
         xxc=rad/cosl(j)
         do i=1,nx,1
            uu(i,1,j)=uu(i,1,j)*xxc
            vv(i,1,j)=vv(i,1,j)*xxc
         enddo
      enddo

c----------------------------------------llming

      call uven0

c----------------------------------------------------------llming

      do j=1,my,1
         sum(j)=0.
      enddo

      do j=1,my,1
         do i=1,nx,1
            sum(j)=sum(j)+uu(i,1,j)
         enddo
         zaw(j,0)=sum(j)/nx
      enddo

c----------------------------------------------------------llming   

c--------------------------------------------------llming   return

      do j=1,my,1
         xxc=rad/cosl(j)
         do i=1,nx,1
            uu(i,1,j)=uu(i,1,j)/xxc
            vv(i,1,j)=vv(i,1,j)/xxc
         enddo
      enddo

c----------------------------------------llming

c----------------------------------------------------------llming
 
      open(31,file='3random.dat',status='unknown')

      do i=1,int(its),1
         read(31,*)ran(i)
      enddo

c----------------------------------------------------------llming










c---------------
      if(ktop.eq.0)then
        do 225 i=1,nxmy
          topo(i,1)=0.
 225    continue

      else if(ktop.eq.1)then
c      call idtopo(nx,my,sinl,topo)
      nxmy8=nx*my*8
      open(99,file='../DATAOUT/TOPO3.OUT',
     &     form='unformatted',status='unknown',access='direct',
     &     recl=nxmy8)
        read(99,rec=1)topo
        close(99)

c*********************************************************************

      else if(ktop.eq.2)then
        lenc=nx*my*8
        open(61,file='../topollming8.dat'
     &,   access='direct',form='unformatted',recl=lenc
     &,   status='unknown')
        read(61,rec=1)topo
        close(61)

        open(62,file='../DATAOUT/TOPO3.OUT'
     &,   access='direct',form='unformatted',recl=lenc
     &,   status='unknown')
        write(62,rec=1)topo
        close(62)
      endif

c********************************************************************

c      call qmaxn3 (topo,'topo',' ',1,1,1,nx,my,1)
      do i=1,nxmy
        topo(i,1)=topo(i,1)*grav
      enddo
      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,topo,topospc)
      do 226 i=1,mlmax
        topospc(i,1)=topospc(i,1)*eps4(i)
        topospc(i,2)=topospc(i,2)*eps4(i)
 226  continue
c
c  transform grid point to spectrum
c
      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,phi,pspcnow)
c
      call trandv(jtrun,mlmax,nx,my,lev,uu,vv,weight,cim
     *, onocos,poly,dpoly,vspcnow,dspcnow)
c
c  transform spectrun vorticity to grid point
c  transform spectrum divergence to grid point
c
      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
      call transr(jtrun,mlmax,nx,my,lev,poly,dspcnow,div)
c
c  set initial spectrum coefficient to old
c
      do 20 i=1,mlmax*2*lev
        vspcold(i,1,1)=vspcnow(i,1,1)
        dspcold(i,1,1)=dspcnow(i,1,1)
        pspcold(i,1,1)=pspcnow(i,1,1)
 20   continue
c
      call outflds(0.)
c
      pi=4.*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
      xlon(1)=0.
      do 10 i=2,nx
        xlon(i)=xlon(1)+float(i-1)*360./nx
 10   continue
      do 15 j=1,my
        xlat(j)=asin(sinl(j))*r2d
 15   continue
c
      do 99 j=1,my
      do 99 i=1,nx
        rforce(i,j)=0.
 99   continue
      clat=3.
      clon=108.
c      print *,'rforce clat,clon=',clat,clon
      a=100.
      b=100.
      rf=30.
      thda=0.
      rot=thda*pi/180.
      radx=rad/1000.
      do 101 j=1,my
      do 101 i=1,nx
        dx=radx*cos(clat*d2r)*(clon-xlon(i))*d2r
        dy=radx*(clat-xlat(j))*d2r
        arg=atan2(dy,dx)
        ra=(dx**2+dy**2)**0.5
        aa=ra*cos(arg+rot)/a
        bb=ra*sin(arg+rot)/b
        r=(aa**2+bb**2)**0.5
        if(r.eq.0)then
          rforce(i,j)=1.
        else if(r.gt.0. .and. r.lt.1.)then
          rforce(i,j)=1.-exp(-rf/r*exp(1./(r-1.)))
        else
          rforce(i,j)=0.
        endif
 101  continue
c      call qmaxn3 (rforce,'forc',' ',1,1,1,nx,my,1)
      call  tranrs(jtrun,mlmax,nx,my,1,poly,weight,rforce,fspc)
c      
c      nxmy8=nx*my*8
c      open(98,file='../DATAOUT/FORCE.OUT',access='direct'
c     &, form='unformatted',recl=nxmy8,status='unknown')
c      write(98,rec=1)rforce
c      close(98)
c      if(1.eq.1)stop

      return
      end
