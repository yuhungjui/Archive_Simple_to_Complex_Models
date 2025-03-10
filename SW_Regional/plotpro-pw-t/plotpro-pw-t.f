      program plotpropwt
c
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=1, IWKID=1)
      parameter (ncl=4)
      parameter (ipror=50)
      data lon,lat/275,325/
c      data iprox,jproy/275,350/
      data iprox,jproy/250,350/
      data indexwind/1/
c
      namelist /modlist/tauer,tauor,dtr,dx,dy,tfiltr,hmean
     *,                bplane,blat,hfiltr,ckd,dohd,linear
c
      dimension u(nxr,myr),v(nxr,myr),s(nxr,myr)
      dimension xw(im,jm),xwu(im,jm),xwv(im,jm)
      real kx(ll),ky(ll),rd(ncl),rdtopo(ncl)
      real side(ll,(2*ipror+1)),xwork(ll)
      real uwind(ll,(2*ipror+1)),uwork(ll)
      real vwind(ll,(2*ipror+1)),vwork(ll)
      real yt,yb
      real xcor(ll)
      integer ishade(ncl+1),ishadetopo(ncl+1)
      character*4 head
      character*11 lab1
      character*2 xlab,ylab
      character*4 llb(ncl+2),llbtopo(ncl+2)
      CHARACTER*16 LDASH(1)
      DATA  LDASH(1) /'$$''$$''$$''$$''$$''$'/
      data rd/20.,30.,40.,50./
      data rdtopo/500.,1500.,2500.,3500./
      data ishade/0,9,11,12,15/
      data ishadetopo/1,3,4,5,6/
c
      nxmy=nxr*myr*8
      open(22,file='../dat/U.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(23,file='../dat/V.OUT',access='direct'
     &, form='unformatted',recl=nxmy,status='unknown')
      open(25,file='../dat/PV.OUT',access='direct'
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
c      call opngks
      CALL GOPKS (IERRF, ISZDM)
      CALL GOPWK (IWKID, LUNIT, IWTYPE)
      CALL GACWK (IWKID)
c
      read(41,rec=1)topo
      close(41)
c
      itauo=0
      nn=1
  51  nrec=itauo+1
      if (nrec.gt.((tauer/tauor)+1.001))go to 155
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
c Put PV data to side(i,j) array
c
      do j=1,myr
            if (abs(j-jproy).le.ipror)then
            side(nrec,(j-jproy+ipror+1))=s(iprox,j)
            uwind(nrec,(j-jproy+ipror+1))=u(iprox,j)
            vwind(nrec,(j-jproy+ipror+1))=v(iprox,j)
            endif
      enddo
c      print *,'Time Step = ',nrec
c
c      call frame
      go to 51
c
 155  continue
c
c Prepare and draw two work arrays
c
      do j=1,(2*ipror+1)
         imod=mod(j,10)
         if (imod.eq.1)then
c
c         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
         do i=1,ll
            xcor(i)=(i-1)*tauor
            xwork(i)=side(i,j)*10.e8
         enddo
c
c Plot time-series
c
         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c            CALL AGSETI('DASH/LENGTH.',16)
c            CALL DISPLA(2,0,1)
            jj=int(j/10)
c            print *,jj
c            yt=1.-.08*float(jj)-.02
c            yb=1.-.08*float(jj+1)
            yt=.12+.08*float(jj+1)-.02
            yb=.12+.08*float(jj)
c            print *,yb,yt
            if (jj.ne.0)then
            CALL SET(.2,.8,YB,YT,0.,12.,-10.,110.,1)
c            CALL ANOTAT('Time','Time',3,0,1,LDASH)
c            call ezy(side(1,j),61,'$')
c            call ezxy(xcor,xwork,ll,' $')
            call curve(xcor,xwork,ll)
c            CALL LABMOD ('(F5.0)','(F5.0)',9,9,9,9,2,2,0)
            CALL LABMOD ('(I5)','(F5.0)',9,9,9,9,3,2,0)
            call periml(12,2,4,1)
                if (indexwind.eq.1)then
                print *,'#################################'
                print *,' '
                print *,'Plotting Windbar at point = ('
     +              ,iprox,',',(j+jproy-ipror-1),' ) '
                do ix=1,ll
                   if (mod(ix,5).eq.1)then
                   call windbar(xcor(ix),30.,
     +             uwind(ix,j),vwind(ix,j),1,1.5,100)
                   endif
                enddo
                endif
            endif
            if (jj.eq.0)then
            CALL SET(.2,.8,YB,YT,0.,12.,-10.,110.,1)
c            CALL ANOTAT('Time','Time',3,0,1,LDASH)
c            call ezy(side(1,j),61,'Time Evolution$')
c            call ezxy(xcor,xwork,ll,'Time Evolution$')
            call curve(xcor,xwork,ll)
c            CALL LABMOD ('(F5.0)','(F5.0)',9,9,9,9,2,2,0)
            CALL LABMOD ('(I5)','(F5.0)',9,9,9,9,3,2,0)
            call periml(12,2,4,1)
                if (indexwind.eq.1)then
                print *,'#################################'
                print *,' '
                print *,'Plotting Windbar at point = ('
     +              ,iprox,',',(j+jproy-ipror-1),' ) '
                do ix=1,ll
                   if (mod(ix,5).eq.1)then
                   call windbar(xcor(ix),30.,
     +             uwind(ix,j),vwind(ix,j),1,1.5,100)
                   endif
                enddo
                endif
            call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
            call gslwsc(1.0)
            call plchhq(0.5,0.08,'Time',0.012,0.,0.)
            call plchhq(0.8,0.08,'(Hour)',0.010,0.,0.)
            call plchhq(0.10,0.55,'PV',0.012,0.,0.)
            call plchhq(0.08,0.98,'(*10e8 1/s)',0.010,0.,0.)
            call plchhq(0.10,0.94,'North',0.010,0.,0.)
            call plchhq(0.10,0.16,'South',0.010,0.,0.)
c            call plchhq(0.5,0.04,'Northwest Side',0.022,0.,0.)
            call plchhq(0.5,0.04,'Southwest Side',0.022,0.,0.)
            endif
         endif
      enddo
c      call frame
c
      close(22)
      close(23)
      close(25)
c
c      call clsgks
      CALL GDAWK (IWKID)
      CALL GCLWK (IWKID)
      CALL GCLKS
c
      stop
      end
