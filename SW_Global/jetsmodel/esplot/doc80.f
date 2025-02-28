      program vormap 
c=====================================================================
c	open file, read data, and plot data
c=====================================================================
cllming      include '../include/paramp.h'

      integer nx,my,lev,itime,nxmy4
      parameter(nx=64,my=64,lev=1,itime=1)

      real pv(1:nx,1:my,1:lev,1:itime),um(1:nx,1:my,1:lev),
     +     up(1:nx,1:my,1:lev,1:itime),vor(1:nx,1:my,1:lev,1:itime),
     +     phi(1:nx,1:my,1:lev,1:itime),vp(1:nx,1:my,1:lev,1:itime),
     +     tt(1:nx,1:my,1:itime),h(1:nx,1:my,1:lev,1:itime),
     +     vm(1:nx,1:my,1:lev),hm(1:nx,1:my,1:lev),
     +     pvm(1:nx,1:my,1:lev),vorm(1:nx,1:my,1:lev),
     +     phim(1:nx,1:my,1:lev),strm(1:nx,1:my,1:lev),
     +     str(1:nx,1:my,1:lev,1:itime),
     +     ttm(1:nx,1:my),track_x(lev,itime),track_y(lev,itime),
     +     tau(nx,my,lev,itime),tilt(itime),vel(nx,my,lev,itime)

      real vortemp(1:nx,1:my,1:lev,1:itime)
      real enplot(1:nx,1:my,1:lev,1:itime)

      character title*13

c--------------------------------------------------
      title(1:13)= 'Vorticity'  !_!
c-----------------------------------------------------------

c-----read data-------------------------------------------------------llming

c      nxmy4=nx*my*4
c      open(22,file='VOR.dat',form='unformatted',
c     +     status='unknown',access='direct',recl=nxmy4)

      open(22,file='enplot80.dat',form='formatted',status='unknown')

      read(22,*),((enplot(i,j,1,1),i=1,nx),j=1,my)

      do i=1,nx,1
         do j=1,my,1
               enplot(i,j,1,1)=enplot(i,j,1,1)*1.0e4
c               enplot(i,j,1,1)=enplot(i,j,1,1)*1.0e1
         enddo
      enddo

c---------------------------------------------------------------llming

c--------------------------------------------------------
cllming      call ana_field(pv,10000.0,nx,my,lev,itime)
cllming      call ana_field(vor,10000.0,nx,my,lev,itime)
c------------------------------------------------------------
c     NCARPLOT !!!
c--------------------------------------------------------------
      CALL OPNGKS

c      do 20 m= iplot,ipend,ito

      do 20 m=1,itime,1
c      do 20 m=1,1,1

      call plot(enplot(1,1,1,m),track_x(1,m),track_y(1,m)
     +,nx,my,m,title,0)  !_! title , "1 center", "0 no center"

      CALL FRAME
 20   continue
      CALL CLSGKS
c-----------------------------------------------------
      print *, 'plot is OK !!'
c---------------------------------------------------------------
      STOP
      END

c=======================================================   
      subroutine plot(f,track_x,track_y,nx,my,istep,title,idot)
c-----------------------------------------------------------
c     f : value need to plot (2D)
c     track_x, track_y : center (1 point)
c     istep : time step
c     idot : 1: need to plot center
c            0: no center
c------------------------------------------------                                                                                
c******** parameter for pl2d ***********
ca      parameter(ncl=16) ! for str !
c      parameter(ncl=3) ! for PV & VOR!
c      parameter(ncl=12) ! for PV & VOR!
      parameter(ncl=9) ! for PV & VOR!
ca      parameter(ncl=9) ! for tt !

      character llb(ncl+2)*4
      real rd(ncl)
      integer ishade(ncl+1)

c***************************************
      integer nx,my,istep,idot
      real f(1:nx,1:my)
      real track_x,track_y
      character label*24,num*5,title*13
c---------------------------------------------------

c-------  set for title  ----------------------------
       label(1:13)=title  !_!
       label(14:16)=' t='  !_!
       label(22:24)=' hr'  !_!
                            
c*******    data for pl2d **************
c
c      data rd     /10.,20.,50.,80.,100.,200./ !for pv!
c      data rd     /5.,10.,20.,40.,80.,160./ !for pv!
c      data rd     /1.,10.,100.,1000.,2000.,5000./ !for pv!
c      data rd     /1.,10.,100.,1000.,2000./ !for pv!
c      data rd     /1000.,2000./ !for pv!
c      data rd     /1000.,2000.,5000./ !for pv!
c      data rd     /1.,2.,9./ !for pv!
c      data rd     /1.,10.,100.,1000.,2000.,3000.,4000.,
c     +                   5000.,6000.,7000.,8000.,9000./ !for pv!
      data rd     /1.,10.,20.,50.,100.,200.,
     +                   500.,1000.,2000./ !for pv!
c      data ishade /1,37,42,44,47,50,12/  ! for pv!
c      data ishade /1,37,42,44,47,50,55/  ! for pv!
c      data ishade /1,35,39,43,47,51,55/  ! for pv!   llming
c      data ishade /1,1,1,1,39,47,55/  ! for pv!
c      data ishade /1,1,1,1,39,47/  ! for pv!
c      data ishade /1,39,47/  ! for pv!
c      data ishade /1,39,47,55/  ! for pv!
c      data ishade /1,1,1,1,39,47,47,47,47,47,47,47,47/  ! for pv!
      data ishade /1,1,1,1,39,39,39,47,47,47/  ! for pv!
c      data ishade /1,1,1,1,1,1,39,47,47,47,47,47,47/  ! for pv!

ca      data rd /-50.,-40.,-30.,-20.,-15.,-10.,-5.,
ca     +         0.,5.,10.,15.,20.,25.,30.,35.,40./ !for str!
ca      data ishade /28,190,188,116,113,84,87,90,1,131,
ca     +             133,58,55,14,150,98,96/  ! for str!

ca      data rd /-4.,-3.,-2.,-1.,0.,1.,2.,3.,4./ !for tt!
ca      data ishade /28,90,86,80,1,1,129,102,14,149/  ! for tt!
c-----------------------------------------------------------                                                                                
      do i=1,ncl
        write(llb(i+1),'(i4)')int(rd(i))
      enddo
      llb(1)='    '
      llb(ncl+2)='    '
      spval=-9999.
                                                                                
c******************************************
c     FOR PLOT !!!!
c******************************************
c--- for title
         write(num,45) 1920*(istep-1)   !_!               !!!!!!!!!!!!!!!!!!!!
cllming   45    format(f5.2)
   45    format(i5)
         label(17:21)=num
         print *,label(1:24)
c---                                                                                
         call dfclrs
         call bckclr(0,0,1)
         call pcseti('FN',22)
         call pcseti('CC',55)
         call cpsetc('ilt-informational label text',' ')
         call set(0.15,0.85,0.15,0.85,1.,nx,1.,my,1)
         call pl2d(0.15,0.85,0.15,0.85,f(1,1),nx,my,rd,ncl,1,
     +             ishade,0,spval)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++040420
c---- for center

         if (idot .eq. 1)then
         call NGDOTS(track_x,track_y,1,2.,98)
         endif
c----
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         call plotbar(llb,3.1,.2,1.0)  !_! plot color bar
cllming         call cpcnrc(temp,nnn,nnn,nnn,0.,100.,10.,-1,-1,-341)
         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c         call set(0.1,0.9,0.1,0.9,1.,dom1,1.,dom1,1)
c         call velvct(temp3,nnn1,temp4,nnn1,nnn1,nnn1,
c         call set(0.1,0.9,0.1,0.9,1.,dom1,1.,dom1,1)
c         call velvct(temp3,nnn1,temp4,nnn1,nnn1,nnn1, 
c     +        -50.,50.,-1,0,0,999.)
c      do j=1,nnn1
c       do i=1,nnn1
c         xxx=real(i)
c         yyy=real(j)
c         if (vel2(i,j).GE.25.) then
c      call  windbar(xxx,yyy,temp3(i,j),temp4(i,j),1,2.0,199)
c         endif
c      enddo
c      enddo

         call set(0.15,0.85,0.15,0.85,1.,64.,1.,64.,1)
         call gslwsc(2.)
         call tick4(24,16,24,16)
c         call perim(1,64,1,64)
         call perim(1,1,1,1)
         call tick4(12,8,12,8)
         call line(1.,1.,64.,64.)

         call gslwsc(1.)
         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
         call plchhq(.5,.975,'Energy Spectrum',24.,0.,0.)
         call plchhq(.5,.915,'m',24.,0.,0.)
cc        call plchhq(.86,.08,'vor',18.,0.,0.)

         call plchhq(.15,.87,'0',18.,0.,0.)
         call plchhq(.261,.87,'10',18.,0.,0.)
         call plchhq(.372,.87,'20',18.,0.,0.)
         call plchhq(.483,.87,'30',18.,0.,0.)
         call plchhq(.594,.87,'40',18.,0.,0.)
         call plchhq(.705,.87,'50',18.,0.,0.)
         call plchhq(.816,.87,'60',18.,0.,0.)
c         call plchhq(.849,.87,'63',18.,0.,0.)
c         call plchhq(.5,.27,'180E',18.,0.,0.)
c         call plchhq(.7,.27,' 90W',18.,0.,0.)
c         call plchhq(.9,.27,'  0 ',18.,0.,0.)

         call plchhq(.07,.5,'n',24.,0.,0.)
cc         call plchhq(.95,.088,'*10e-4',18.,0.,0.)

         call plchhq(.12,.15,'0',18.,0.,0.)
         call plchhq(.12,.261,'10',18.,0.,0.)
         call plchhq(.12,.372,'20',18.,0.,0.)
         call plchhq(.12,.483,'30',18.,0.,0.)
         call plchhq(.12,.594,'40',18.,0.,0.)
         call plchhq(.12,.705,'50',18.,0.,0.)
         call plchhq(.12,.816,'60',18.,0.,0.)
c         call plchhq(.06,'90N',20.,0.,0.)

      RETURN
      END
                                                                                

