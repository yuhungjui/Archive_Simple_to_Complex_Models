      program reconv


c----------------------------------------------------------c
c  Program           : Regional Spectral Convective Model  c
c  Original Author   : Prof. Hung-Chi Kuo                  c
c  Reprogram         : Yu-Ming Tsai                        c
c  Date              : 2008 Dec. 17                        c
c  Institute         : CSUAS, Fort Collins, CO             c
c----------------------------------------------------------c


c---------------------------------------------------------------------------c
c  This is a spectral convection model with Ooyama's thermodynamic scheme.  c
c  Galerikin (x dir.) and Tau (z dir.) methods are used.                    c
c  The accompany include file 'c2o.f' could be found at                     c
c  the end of this program file.                                            c
c---------------------------------------------------------------------------c


           include 'c2o.f'


      call psetup                  !!! problem setup

      call consetup                !!! constant setup

      call set_lanczos_filter      !!! settings for using Lanczos filter




c----------------------------------------------------------
c  Select an example to run.
c----------------------------------------------------------


c----------------------------------------------------------
c  ex.1) Dry hydrostatic adjustment experiment.
c     (Hydrostatic adjustment by acoustic wave.)
c     This experiment will generate lots of acoustic waves,
c  so remember to integrate with the small enough time step (i.e. 0.075 sec).
c  Because of the rapid propagating acoustic waves will hit the boundary
c  and then reflect back, do not integrate this experiment over "150 seconds".
c----------------------------------------------------------

      call initial_cond1



c----------------------------------------------------------
c  ex.2) Rising dry bubble in a hydrostatic atmosphere.
c----------------------------------------------------------

c      call initial_cond2



c----------------------------------------------------------
c  ex.3) Rising moist bubble in a hydrostatic atmosphere.
c----------------------------------------------------------

c      call initial_cond3



c----------------------------------------------------------
c  ex.4) Condensation with a specified tilt updraft.
c----------------------------------------------------------

c      call initial_cond4







      call basic_field             !!! set the basic state


           time=0.0
      call drive (time)            !!! output the initial fields

      call run                     !!! begin RK4 integration

           stop 'normal end'
           end





c************************************************
      subroutine psetup

      include 'c2o.f'

c----------------------------------------------------------
c  Open binary output files.
c      mxmz8=(mnx+1)*(mnz+1)*8
c      open(21,file='./dataout/entropy.dat',access='direct',
c     +     form='unformatted',recl=mxmz8,status='unknown')   !!! 1
c      open(22,file='./dataout/drydensity.dat',access='direct',
c     +     form='unformatted',recl=mxmz8,status='unknown')   !!! 2
c      open(23,file='./dataout/moistdensity.dat',access='direct',
c     +     form='unformatted',recl=mxmz8,status='unknown')   !!! 3
c      open(24,file='./dataout/uu.dat',access='direct',
c     +     form='unformatted',recl=mxmz8,status='unknown')   !!! 4
c      open(25,file='./dataout/ww.dat',access='direct',
c     +     form='unformatted',recl=mxmz8,status='unknown')   !!! 5
c----------------------------------------------------------

c----------------------------------------------------------
c  Open ASC II output files.
      open(31,file='./dataout/entropy_asc.dat',status='unknown')      !!! 1
      open(32,file='./dataout/drydensity_asc.dat',status='unknown')   !!! 2
      open(33,file='./dataout/moistdensity_asc.dat',status='unknown') !!! 3
      open(34,file='./dataout/uu_asc.dat',status='unknown')           !!! 4
      open(35,file='./dataout/ww_asc.dat',status='unknown')           !!! 5
      open(36,file='./dataout/temp_asc.dat',status='unknown')         !!! 6
      open(37,file='./dataout/pressure_asc.dat',status='unknown')     !!! 7
      open(38,file='./dataout/etac_asc.dat',status='unknown')         !!! 8
c----------------------------------------------------------


c----------------------------------------------------------
c  Set spatial and time frame of the model.
c----------------------------------------------------------

      dox   = 2500.0       ! X domain size (meters)
      doz   = 2500.0       ! Y domain size (meters)

      timax =    3.0       ! Total integration time (seconds)
      delt  =    0.075     ! Time integration step for resolution 33*33 (seconds)
c      delt  =   0.01875   ! Time integration step for resolution 65*65 (seconds)
      ito   =   10         ! Output time step interval (step)

      tpmax =    2.5       ! max temperature anomaly (K)


      return
      end





c************************************************
      subroutine consetup

      include 'c2o.f'

c------------------------------------------------
c  Reference:
c  Wallace and Hobbs, 2006:
c  Atmospheric Science, An Introductory Survey
c------------------------------------------------

      cw    =   4218.    !! Specific heat of liquid water at 0 C              (J K^-1 kg^-1)
      ci    =   2106.    !! Specific heat of ice at 0 C                       (J K^-1 kg^-1)
      xlw0  =2501000.    !! Latent heat of condensation at 0 C                (J kg^-1)
      xli0  =2380000.    !! Doesn't be used in the model 
      rv    =    461.51  !! Gas constant for water vapor                      (J K^-1 kg^-1)
      ra    =    287.05  !! Gas constant for dry air                          (J K^-1 kg^-1)
      cvv   =   1463.    !! Specific heat of water vapor at constant volume   (J K^-1 kg^-1)
      cpv   =   1952.    !! Specific heat of water vapor at constant pressure (J K^-1 kg^-1)
      cva   =    717.    !! Specific heat of dry air at constant volume       (J K^-1 kg^-1)
      cpa   =   1004.    !! Specific heat of dry air at constant pressure     (J K^-1 kg^-1)
      t0    =    273.15  !! Reference temperature 0 C                         (K)
      xi0   = 100000./(ra*t0) !! Reference dry air density = 1.27538          (kg m^-3)
      e0    =    610.7        !! Used in function ew(t)                       
      ggg   =      9.81       !! Gravity at sea level                         (m s^-2)
      xlamd0=xlamd(t0)        !! Specific entropy gain by vaporizing condensate at 0 C
      seta0 = seta(t0)        !! Mass density of saturated vapor at 0 C       (kg m^-3)
      dddt  =      0.05       !! Adjustable scaling constant                  (K)

      return
      end





c************************************************
      subroutine basic_field

      include 'c2o.f'

      real term(-1:mnx,0:mnz,9),imnx(13),fmnx(3*mnx/2+1),
     &     phy(0:mnz),spec(0:nz),work((mnz+2)*(mnx+2))

      call fftfax(mnx,imnx,fmnx)

      do k=1,9,1
         do j=0,mnz,1
            do i=-1,mnx,1
               term(i,j,k)=0.0   !!! initialize matrix to 0.
            enddo
         enddo
      enddo

      do k=1,8,1
         do j=0,mnz,1
            bar(j,k)=0.0   !!! initialize matrix to 0.
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            term(i,j,1)=sigs(i,j)   !!! these vlaues are from initial condition (spectral).
            term(i,j,2)=xis (i,j)
            term(i,j,3)=etas(i,j)
            term(i,j,4)=uds (i,j)
            term(i,j,5)=wds (i,j)
         enddo
      enddo

c      print*,'term(12,12,3)s=',term(12,12,3)

      call delxz

      do k=1,5,1
         do i=-1,nx,1

            do j=0,nz,1
               spec(j)=term(i,j,k)
            enddo
            
            call cctstp(mnz+1,z,nz,0.0,doz,spec,phy)   !!! spc --> phy
           
            do j=0,mnz,1
               term(i,j,k)=phy(j)
            enddo

         enddo
      enddo

      do k=1,5,1
         call fft99(term(-1,0,k),work,fmnx,imnx,1,mnx+2,mnx,mnz+1,+1)   !!! spc --> phy
      enddo

c      print*,'term(12,12,3)_2p=',term(12,12,3)

      larry=(mnx+2)*(mnz+1)

c  Description: subroutine thermodiagnost(sigg(1),xi(2),eta(3),
c                               temp(6),p_sigg(7),p_xi(8),p_eta(9),mx)

      call thermodiagnost(term(-1,0,1),term(-1,0,2),term(-1,0,3),
     &     term(-1,0,6),term(-1,0,7),term(-1,0,8),term(-1,0,9),larry)

c  Description: subroutine compute_pressure(xi(2),eta(3),temp(6),p_sigg(7),etac,mx)
      call compute_pressure(term(-1,0,2),term(-1,0,3),term(-1,0,6),
     &     term(-1,0,7),term(-1,0,8),larry)

      do j=0,mnz,1

c  Finally used in drive before plot
        
         bar(j,1)=term(-1,j,1)            !!! entropy
         bar(j,2)=term(-1,j,2)            !!! xi
         bar(j,3)=term(-1,j,3)            !!! eta
         bar(j,6)=term(-1,j,6)            !!! temp
         bar(j,7)=term(-1,j,7)/100.       !!! pressure
         bar(j,8)=bar(j,3)/seta(bar(j,6)) !!!
         
 
         do i=1,mnx,1
            ttbar(i,j)=term(i,j,6)
            xibar(i,j)=term(i,j,2)
         enddo

      enddo
 
c      print*,'bar(16,3)=',bar(16,3)

      return
      end





c************************************************
      subroutine initial_cond1

c----------------------------------------------------------
c  ex.1) Dry hydrostatic adjustment experiment.
c     (Hydrostatic adjustment by acoustic wave.)
c     This experiment will generate lots of acoustic waves,
c  so remember to integrate with the small enough time step (i.e. 0.075 sec).
c  Because of the rapid propagating acoustic waves will hit the boundary
c  and then reflect back, do not integrate this experiment over "150 seconds".
c----------------------------------------------------------

      include 'c2o.f'

      real tini(-1:mx,0:mz+1,5),work((mz+2)*(mx+2)),wk(-1:mx)
      real temp(1:mx,1:mz)

      call fftfax(mx,imx,fmx)
      call mctset(mz,imz,fmz)

      pi=acos(-1.0)

      call colpts(dox,doz,mx,mz,zj(0),xxi(-1))

      do k=1,5,1
         do i=-1,mx,1
            do j=0,mz+1,1
               tini(i,j,k)=0.0
            enddo
         enddo
      enddo

      z01  =  doz/2.
      h1   =  50.
      x01  =  dox/2.
      hx1  =  200.
      z02  =  doz/2.
      hz2  =  200.

      do i=-1,mx,1
         do j=0,mz,1
            tini(i,j,2)=xi0*0.9                      !!! = 1.147846 (kg m^-3)
            tini(i,j,3)=seta(293.15)*0.000000005     !!! eta
            tini(i,j,4)=0.0                          !!! u
            tini(i,j,5)=0.0                          !!! v

            tt(i,j)    =293.15-ggg/ra*zj(j)          !!! basic state temp. profile.
            pp(i,j)    =tini(i,j,2)*ra*tt(i,j)
         enddo
      enddo

c      print*,'t0           =', t0
c      print*,'efun(293.15) =', efun(293.15 )
c      print*,'seta(293.15) =', seta(293.15 )
c      print*,'tini(16,16,3)_1=', tini(16,16,3)   !!! 1

      do i=-1,mx,1            !!! Gaussian shape heating
         do j=0,mz,1
            tt(i,j)=tt(i,j)+tpmax*exp(-((xxi(i)-x01)/hx1)**2)
     &      *exp(-(( zj(j)-z02)/hz2)**2)
         enddo
      enddo


      larry=(mx+2)*(mz+1)
 
c      print*,'tini(16,16,3)_2=',tini(16,16,3)   !!! 2

      call compute_entropy(tini(-1,0,2),tini(-1,0,3),
     &                     tini(-1,0,1),tt(-1,0),larry)

c      print*,'tini(18,24,3)_3=',tini(18,24,3)   !!! 3

      do k=1,5,1
         call mctranf(mz,tini(-1,0,k),+1,mx+2,mx+2,1,imz,fmz,work)  !!! phy --> spc
         call fft99(tini(-1,0,k),work,fmx,imx,1,mx+2,mx,nz+1,-1)    !!! phy --> spc
      enddo

c      print*,'tini(18,24,3)_4s=',tini(18,24,3)   !!! 4

      do i=-1,nx,1
         do j=0,nz,1
            sigs(i,j)=tini(i,j,1)
            xis (i,j)=tini(i,j,2)
            etas(i,j)=tini(i,j,3)
            uds (i,j)=tini(i,j,4)
            wds (i,j)=tini(i,j,5)
         enddo
      enddo

c------------------------------------------------
c  Fianlly we get the initial spectral coefficients of
c  sigs, xis, etas, uds, wds.
c------------------------------------------------

      return
      end





c************************************************
      subroutine initial_cond2

c----------------------------------------------------------
c  ex.2) Rising dry bubble in a hydrostatic atmosphere.
c----------------------------------------------------------

      include 'c2o.f'

      real tini(-1:mx,0:mz+1,5),work((mz+2)*(mx+2)),wk(-1:mx)
      real temp(1:mx,1:mz)

      call fftfax(mx,imx,fmx)
      call mctset(mz,imz,fmz)

      pi=acos(-1.0)

      call colpts(dox,doz,mx,mz,zj(0),xxi(-1))

      do k=1,5,1
         do i=-1,mx,1
            do j=0,mz+1,1
               tini(i,j,k)=0.0
            enddo
         enddo
      enddo

      z01  =  doz/2.
      h1   =  50.
      x01  =  dox/2.
      hx1  =  200.
      z02  =  doz/4.
      hz2  =  200.

      do i=-1,mx,1
         do j=0,mz,1
            tt(i,j)    =293.15-ggg/cpa*zj(j)            !!! basic state temp. profile.
            tini(i,j,2)=xi0*(tt(i,j)/293.15)**(cva/ra)  !!! (kg m^-3)
            tini(i,j,3)=seta(293.15)*0.0000001        !!! eta
            pp(i,j)    =tini(i,j,2)*ra*tt(i,j)
            tini(i,j,4)=0.0                             !!! u
            tini(i,j,5)=0.0                             !!! v
         enddo
      enddo

c      print*,'t0           =', t0
c      print*,'efun(293.15) =', efun(293.15 )
c      print*,'seta(293.15) =', seta(293.15 )
c      print*,'tini(16,16,3)_1=', tini(16,16,3)   !!! 1

      do i=-1,mx,1            !!! Gaussian shape heating
         do j=0,mz,1
            tt(i,j)=tt(i,j)+tpmax*exp(-((xxi(i)-x01)/hx1)**2)
     &      *exp(-(( zj(j)-z02)/hz2)**2)

            tini(i,j,2)=pp(i,j)/(ra*tt(i,j))
         enddo
      enddo

c      print*,'tt(15,15)=',tt(15,15)

      larry=(mx+2)*(mz+1)
 
c      print*,'tini(16,16,3)_2=',tini(16,16,3)   !!! 2

      call compute_entropy(tini(-1,0,2),tini(-1,0,3),
     &                     tini(-1,0,1),tt(-1,0),larry)

c      print*,'tini(18,24,1)_3=',tini(18,24,1)   !!! 3

      do k=1,5,1
         call mctranf(mz,tini(-1,0,k),+1,mx+2,mx+2,1,imz,fmz,work)  !!! phy --> spc
         call fft99(tini(-1,0,k),work,fmx,imx,1,mx+2,mx,nz+1,-1)    !!! phy --> spc
      enddo

c      print*,'tini(18,24,3)_4s=',tini(18,24,3)   !!! 4

      do i=-1,nx,1
         do j=0,nz,1
            sigs(i,j)=tini(i,j,1)
            xis (i,j)=tini(i,j,2)
            etas(i,j)=tini(i,j,3)
            uds (i,j)=tini(i,j,4)
            wds (i,j)=tini(i,j,5)
         enddo
      enddo

c------------------------------------------------
c  Fianlly we get the initial spectral coefficients of
c  sigs, xis, etas, uds, wds.
c------------------------------------------------

      return
      end





c************************************************
      subroutine initial_cond3

c----------------------------------------------------------
c  ex.3) Rising moist bubble in a hydrostatic atmosphere.
c----------------------------------------------------------

      include 'c2o.f'

      real tini(-1:mx,0:mz+1,5),work((mz+2)*(mx+2)),wk(-1:mx),
     &     rtb(0:mz)
      real temp(1:mx,1:mz)

      call fftfax(mx,imx,fmx)
      call mctset(mz,imz,fmz)

      pi=acos(-1.0)

      call colpts(dox,doz,mx,mz,zj(0),xxi(-1))

      do k=1,5,1
         do i=-1,mx,1
            do j=0,mz+1,1
               tini(i,j,k)=0.0
            enddo
         enddo
      enddo

      z01  =  doz/2.
      h1   =  50.
      x01  =  dox/2.
      hx1  =  200.
      z02  =  doz/4.
      hz2  =  200.
      ha   = 1500.

      do i=-1,mx,1
         do j=0,mz,1

c  rtb is the mixing ratio of total water density with respect to the dry air density.

            rtb(j)     =0.0135*exp(-zj(j)/ha)               !!! 
            tt(i,j)    =293.15-ggg/cpa*zj(j)                !!! basic state temp. profile.
            tini(i,j,2)=xi0*(tt(i,j)/293.15)**(cva/ra)
     &                 *exp(rv/ra*(0.0135-rtb(j)))          !!! (kg m^-3)
            tini(i,j,3)=tini(i,j,2)*rtb(j)                  !!! eta
            pp(i,j)    =tini(i,j,2)*(ra+rtb(j)*rv)*tt(i,j)  !!! pressure
            tini(i,j,4)=0.0                                 !!! u
            tini(i,j,5)=0.0                                 !!! v

         enddo
      enddo

c      print*,'t0           =', t0
c      print*,'efun(293.15) =', efun(293.15 )
c      print*,'seta(293.15) =', seta(293.15 )
c      print*,'tini(16,16,3)_1=', tini(16,16,3)   !!! 1

      do i=-1,mx,1            !!! Gaussian shape heating
         do j=0,mz,1
            tt(i,j)=tt(i,j)+tpmax*exp(-((xxi(i)-x01)/hx1)**2)
     &      *exp(-(( zj(j)-z02)/hz2)**2)

            tini(i,j,2)=pp(i,j)/((ra+rtb(j)*rv)*tt(i,j))
            tini(i,j,3)=tini(i,j,2)*rtb(j)
         enddo
      enddo

c      print*,'tt(15,15)=',tt(15,15)

      larry=(mx+2)*(mz+1)
 
c      print*,'tini(16,16,3)_2=',tini(16,16,3)   !!! 2

      call compute_entropy(tini(-1,0,2),tini(-1,0,3),
     &                     tini(-1,0,1),tt(-1,0),larry)

c      print*,'tini(18,24,1)_3=',tini(18,24,1)   !!! 3

      do k=1,5,1
         call mctranf(mz,tini(-1,0,k),+1,mx+2,mx+2,1,imz,fmz,work)  !!! phy --> spc
         call fft99(tini(-1,0,k),work,fmx,imx,1,mx+2,mx,nz+1,-1)    !!! phy --> spc
      enddo

c      print*,'tini(18,24,3)_4s=',tini(18,24,3)   !!! 4

      do i=-1,nx,1
         do j=0,nz,1
            sigs(i,j)=tini(i,j,1)
            xis (i,j)=tini(i,j,2)
            etas(i,j)=tini(i,j,3)
            uds (i,j)=tini(i,j,4)
            wds (i,j)=tini(i,j,5)
         enddo
      enddo

c------------------------------------------------
c  Fianlly we get the initial spectral coefficients of
c  sigs, xis, etas, uds, wds.
c------------------------------------------------

      return
      end





c************************************************
      subroutine initial_cond4

c----------------------------------------------------------
c  ex.4) Condensation with a specified tilt updraft.
c----------------------------------------------------------

      include 'c2o.f'

      real tini(-1:mx,0:mz+1,5),work((mz+2)*(mx+2)),wk(-1:mx),
     &     rtb(0:mz)
      real temp(1:mx,1:mz)

      call fftfax(mx,imx,fmx)
      call mctset(mz,imz,fmz)

      pi=acos(-1.0)

      call colpts(dox,doz,mx,mz,zj(0),xxi(-1))

      do k=1,5,1
         do i=-1,mx,1
            do j=0,mz+1,1
               tini(i,j,k)=0.0
            enddo
         enddo
      enddo

      x01  =  dox/2.
      hx1  =  200.
      z02  =  doz/2.
      hz2  =  200.
      ha   = 1500.
      umax =    5.0
      phy  = pi*(-30./180.)
      x01  = x01-doz/2.*tan(phy)

      do i=-1,mx,1
         do j=0,mz,1
            xx         =x01+tan(phy)*zj(j)
c  rtb is the mixing ratio of total water density with respect to the dry air density.
            rtb(j)     =0.01275*exp(-zj(j)/ha)              !!! 
            tt(i,j)    =293.15-ggg/cpa*zj(j)                !!! basic state temp. profile.
            tini(i,j,2)=xi0*(tt(i,j)/293.15)**(cva/ra)
     &                 *exp(rv/ra*(0.01275-rtb(j)))         !!! (kg m^-3)
            tini(i,j,3)=tini(i,j,2)*rtb(j)                  !!! eta
            pp(i,j)    =tini(i,j,2)*(ra+rtb(j)*rv)*tt(i,j)  !!! pressure
            tini(i,j,4)=umax*sin(phy)*sin(pi*zj(j)/doz)
     &                 *exp(-((xxi(i)-xx)/hx1)**2)          !!! u
            tini(i,j,5)=umax*cos(phy)*sin(pi*zj(j)/doz)                                 
     &                 *exp(-((xxi(i)-xx)/hx1)**2)          !!! w
         enddo
      enddo

c      print*,'t0           =', t0
c      print*,'efun(293.15) =', efun(293.15 )
c      print*,'seta(293.15) =', seta(293.15 )
c      print*,'tini(16,16,3)_1=', tini(16,16,3)   !!! 1

c      do i=-1,mx,1            !!! Gaussian shape heating
c         do j=0,mz,1
c            tt(i,j)=tt(i,j)+tpmax*exp(-((xxi(i)-x01)/hx1)**2)
c     &      *exp(-(( zj(j)-z02)/hz2)**2)
c
c            tini(i,j,2)=pp(i,j)/((ra+rtb(j)*rv)*tt(i,j))
c            tini(i,j,3)=tini(i,j,2)*rtb(j)
c         enddo
c      enddo

c      print*,'tt(15,15)=',tt(15,15)

      larry=(mx+2)*(mz+1)
 
c      print*,'tini(16,16,3)_2=',tini(16,16,3)   !!! 2

      call compute_entropy(tini(-1,0,2),tini(-1,0,3),
     &                     tini(-1,0,1),tt(-1,0),larry)

c      print*,'tini(18,24,1)_3=',tini(18,24,1)   !!! 3

      do k=1,5,1
         call mctranf(mz,tini(-1,0,k),+1,mx+2,mx+2,1,imz,fmz,work)  !!! phy --> spc
         call fft99(tini(-1,0,k),work,fmx,imx,1,mx+2,mx,nz+1,-1)    !!! phy --> spc
c   LOT=nz+1 because truncates at wavenumber nz.
      enddo

c      print*,'tini(18,24,3)_4s=',tini(18,24,3)   !!! 4

      do i=-1,nx,1
         do j=0,nz,1
            sigs(i,j)=tini(i,j,1)
            xis (i,j)=tini(i,j,2)
            etas(i,j)=tini(i,j,3)
            uds (i,j)=tini(i,j,4)
            wds (i,j)=tini(i,j,5)
         enddo
      enddo

c------------------------------------------------
c  Fianlly we get the initial spectral coefficients of
c  sigs, xis, etas, uds, wds.
c------------------------------------------------

      return
      end





c************************************************
      subroutine colpts(dox,doz,mx,mz,zj,xi)

c------------------------------------------------
c  Set the regular grid points in the X direction,
c  and irregular collocation grid points in the Z 
c  direction in the model
c------------------------------------------------

      real dox,doz,zj(0:mz),xi(-1:mx)

      pi=acos(-1.0)

      do i=0,mx,1
         xi(i)=float(i)*dox/float(mx)
      enddo

      do j=0,mz,1
         zj(j)=0.5*doz*(cos(j*pi/mz)+1.0)
      enddo


      xi(-1)=xi(mx-1)

      return
      end





c************************************************
      subroutine run

      include 'c2o.f'

      real rk(-1:nx)
      save rk

      pi=acos(-1.0)
 
      call waveno(rk(-1))

      itmax=int(timax/delt)


      do istep=1,itmax,1

         time=delt*float(istep)

         call rk4(rk(-1),time)   !!! 4th order Runge-Kutta method

         kk=mod(istep,ito)
         if (kk .eq. 0) then
            call drive(time)

c            print*,'time       =',time,'  seconds'
c            print*,'spectral coeffients.'
c            print*,'sigs(16,16)=',sigs(16,16)
c            print*,'xis(16,16) =',xis(16,16)
c            print*,'etas(16,16)=',etas(16,16)
c            print*,'uds(16,16) =',uds(16,16)
c            print*,'wds(16,16) =',wds(16,16)
c            print*,'                       '

         endif

      enddo
 
         print*,'                       '
         print*,'timax=',timax
         print*,'delt =',delt
         print*,'itmax=',itmax
         print*,'                       '
      return
      end





c************************************************
      subroutine rk4 (rk,time)

c------------------------------------------------
c  time intergration scheme in the spectral space.
c  Galerikin in X direction
c  Tau in Z direction
c------------------------------------------------

      include 'c2o.f'

      real rk(-1:nx)
      real tv(-1:nx,0:nz,1:5),tf(-1:nx,0:nz,1:5)
      save tv,tf

c------------------1-----------------------------

      call rhs(time)
      call lanc_filter

      do j=0,nz,1
         do i=-1,nx,1
            tv(i,j,1)=sigs(i,j)   !!! spectral coefficients
            tv(i,j,2)= xis(i,j)
            tv(i,j,3)=etas(i,j)
            tv(i,j,4)= uds(i,j)
            tv(i,j,5)= wds(i,j)

            tf(i,j,1)= f(i,j,1)   !!! weighting 1
            tf(i,j,2)= f(i,j,2)
            tf(i,j,3)= f(i,j,3)
            tf(i,j,4)= f(i,j,4)
            tf(i,j,5)= f(i,j,5)
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            sigs(i,j)=tv(i,j,1)+0.5*delt* f(i,j,1)
            xis (i,j)=tv(i,j,2)+0.5*delt* f(i,j,2)
            etas(i,j)=tv(i,j,3)+0.5*delt* f(i,j,3)
            uds (i,j)=tv(i,j,4)+0.5*delt* f(i,j,4)
            wds (i,j)=tv(i,j,5)+0.5*delt* f(i,j,5)
         enddo
      enddo

      call boundary(wds(-1,0))

c------------------2-----------------------------

      call rhs(time)
      call lanc_filter

      do k=1,5,1
         do j=0,nz,1
            do i=-1,nx,1
               tf(i,j,k)=tf(i,j,k)+f(i,j,k)*2.0   !!! weighting 2
            enddo
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            sigs(i,j)=tv(i,j,1)+0.5*delt* f(i,j,1)
            xis (i,j)=tv(i,j,2)+0.5*delt* f(i,j,2)
            etas(i,j)=tv(i,j,3)+0.5*delt* f(i,j,3)
            uds (i,j)=tv(i,j,4)+0.5*delt* f(i,j,4)
            wds (i,j)=tv(i,j,5)+0.5*delt* f(i,j,5)
         enddo
      enddo
            
      call boundary (wds(-1,0))

c------------------3-----------------------------

      call rhs(time)
      call lanc_filter

      do k=1,5,1
         do j=0,nz,1
            do i=-1,nx,1
               tf(i,j,k)=tf(i,j,k)+f(i,j,k)*2.0   !!! weighting 2
            enddo
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            sigs(i,j)=tv(i,j,1)+delt* f(i,j,1)
            xis (i,j)=tv(i,j,2)+delt* f(i,j,2)
            etas(i,j)=tv(i,j,3)+delt* f(i,j,3)
            uds (i,j)=tv(i,j,4)+delt* f(i,j,4)
            wds (i,j)=tv(i,j,5)+delt* f(i,j,5)
         enddo
      enddo
            
      call boundary (wds(-1,0))

c------------------4-----------------------------

      call rhs(time)
      call lanc_filter

      do k=1,5,1
         do j=0,nz,1
            do i=-1,nx,1
               tf(i,j,k)=tf(i,j,k)+f(i,j,k)   !!! weighting 1
            enddo
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            sigs(i,j)=tv(i,j,1)+delt* tf(i,j,1)/6.
            xis (i,j)=tv(i,j,2)+delt* tf(i,j,2)/6.
            etas(i,j)=tv(i,j,3)+delt* tf(i,j,3)/6.
            uds (i,j)=tv(i,j,4)+delt* tf(i,j,4)/6.
            wds (i,j)=tv(i,j,5)+delt* tf(i,j,5)/6.
         enddo
      enddo
            
      call boundary (wds(-1,0))

c------------------end---------------------------


      return
      end





c************************************************
      subroutine rhs (time)

      include 'c2o.f'

      real term(-1:mx,0:mz+1,15),work((mz+2)*(mx+2))

      pi=acos(-1.0)

      do k=1,15,1
         do j=0,mz+1,1
            do i=-1,mx,1
               term(i,j,k)=0.0
            enddo
         enddo
      enddo

      do j=0,mz+1,1
         do i=-1,mx,1
            bf(i,j) =0.0
            pxp(i,j)=0.0
            pzp(i,j)=0.0
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            term(i,j,1) = sigs(i,j)   !!! spectral coefficients
            term(i,j,2) = xis (i,j)
            term(i,j,3) = etas(i,j)
            term(i,j,4) = uds (i,j)
            term(i,j,5) = wds (i,j)

            term(i,j,6) = sigs(i,j)
            term(i,j,7) = xis(i,j)
            term(i,j,8) = etas(i,j)

            term(i,j,9) = sigs(i,j)
            term(i,j,10)= xis(i,j)
            term(i,j,11)= etas(i,j)
         enddo
      enddo

      do k=6,8,1
         call fder99(term(-1,0,k),work,1,mx+2,nx,nz+1,dox)   !!! derivative in spectral space
      enddo

      do k=9,11,1
         call mcsder(nz,0.0,doz,term(-1,0,k),nx+2,mx+2,1,work)  !!! derivative in spectral space
      enddo

      do k=1,11,1
         call fft99(term(-1,0,k),work,fmx,imx,1,mx+2,mx,nz+2,+1)  !! spc --> phy
      enddo

      do k=1,11,1
         call mctranf(mz,term(-1,0,k),-1,mx+2,mx+2,1,imz,fmz,work)  !! spc --> phy
      enddo

      do j=0,mz,1
       do i=-1,mx,1
        term(i,j,2)=(term(i,j,2)+abs(term(i,j,2))+0.00001)/2.
        term(i,j,3)=(term(i,j,3)+abs(term(i,j,3))+0.00000001)/2.
       enddo
      enddo

      larry=(mx+2)*(mz+1)

      call thermodiagnost(term(-1,0,1),term(-1,0,2),term(-1,0,3),
     & tt(-1,0),p1(-1,0),p2(-1,0),p3(-1,0),larry )  !! thermodynamic diagnose in phy space
 
      do j=0,mz,1
         do i=-1,mx,1
            pxp(i,j)=p1(i,j)*term(i,j,6)+p2(i,j)*term(i,j,7)+
     &               p3(i,j)*term(i,j,8)
            bf(i,j) =p1(i,j)*term(i,j,9)+p2(i,j)*term(i,j,10)+
     &               p3(i,j)*term(i,j,11)+
     &               (term(i,j,2)+term(i,j,3))*ggg
         enddo
      enddo

      do k=1,5,1
       m=k+5
       l=k+10
       do j=0,mz,1
        do i=-1,mx,1
         term(i,j,m)=term(i,j,k)*term(i,j,4)/
     &               (term(i,j,2)+term(i,j,3))
         term(i,j,l)=term(i,j,k)*term(i,j,5)/
     &               (term(i,j,2)+term(i,j,3))
        enddo
       enddo
      enddo

      do i=-1,mx,1
         bf(i,0)   =0.0
         bf(i,mz)  =0.0
         bf(i,1)   =0.0
         bf(i,mz-1)=0.0
      enddo

      do k=6,15,1
         call mctranf(mz,term(-1,0,k),+1,mx+2,mx+2,1,imz,fmz,work)  !! phy --> spc
      enddo

      do k=6,15,1
         call fft99(term(-1,0,k),work,fmx,imx,1,mx+2,mx,nz+1,-1)  !! phy --> spc
      enddo

      call mctranf(mz,pxp(-1,0),+1,mx+2,mx+2,1,imz,fmz,work)  !! phy --> spc
      call fft99(pxp(-1,0),work,fmx,imx,1,mx+2,mx,nz+1,-1)    !! phy --> spc
      call mctranf(mz,bf(-1,0),+1,mx+2,mx+2,1,imz,fmz,work)   !! phy --> spc
      call fft99(bf(-1,0),work,fmx,imx,1,mx+2,mx,nz+1,-1)     !! phy --> spc

      do k=5,10,1
         call fder99(term(-1,0,k),work,1,mx+2,nx,nz+1,dox)
      enddo
      do k=11,15,1
         call mcsder(nz,0.0,doz,term(-1,0,k),nx+2,mx+2,1,work)
      enddo
     
      do k=1,5,1
       m=k+5
       l=k+10
        do j=0,nz,1
         do i=-1,nx,1
          f(i,j,k)=(term(i,j,m)+term(i,j,l)) *(-1.0)
         enddo
        enddo
       enddo

       do j=0,nz,1
        do i=-1,nx,1
         f(i,j,4)=f(i,j,4)-pxp(i,j)
         f(i,j,5)=f(i,j,5)-bf(i,j)
        enddo
       enddo
 
       return
       end





c************************************************
      subroutine drive(time)

c------------------------------------------------
c  Deal with all the output processes of the model.
c------------------------------------------------

c------------------------------------------------
c   term(i,j,1):  sigg, entropy
c   term(i,j,2):  xis , dry air density
c   term(i,j,3):  eta , total airborne moisture
c   term(i,j,4):  uu  , horizontal velocity
c   term(i,j,5):  ww  , vertical velocity
c   term(i,j,6):  t   , temperature
c   term(i,j,7):  p   , pressure 
c   term(i,j,8):  etac, saturated total airborn moisture density
c   term(i,j,9):  div , divergence
c------------------------------------------------

      include 'c2o.f'

      real term(-1:mnx,0:mnz,9),imnx(13),fmnx(3*mnx/2+1),
     &   phy(0:mnz),spec(0:nz),work((mnz+2)*(mnx+2)),
     &   termm(0:mnx,0:mnz,9)

      real tp(0:mnx,0:mnz)
      real temp(1:mnx+1,1:mnz+1,16)


      call fftfax(mnx,imnx,fmnx)

      do k=1,9,1
         do j=0,mnz,1
            do i=-1,mnx,1
               term(i,j,k)=0.0   !!!   Initialize matrix
            enddo
         enddo
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            term(i,j,1) = sigs(i,j)   !!! spectral coefficients
            term(i,j,2) = xis (i,j)
            term(i,j,3) = etas(i,j)
            term(i,j,4) = uds (i,j)
            term(i,j,5) = wds (i,j)
         enddo
      enddo

      call delxz

      do k=1,5,1
         do i=-1,nx,1

            do j=0,nz,1
               spec(j)=term(i,j,k)
            enddo
            
            call cctstp(mnz+1,z,nz,0.0,doz,spec,phy)   !!! spc --> phy
           
            do j=0,mnz,1
               term(i,j,k)=phy(j)
            enddo

         enddo
      enddo

      do k=1,5,1
         call fft99(term(-1,0,k),work,fmnx,imnx,1,mnx+2,mnx,mnz+1,+1)  !! spc --> phy
      enddo

      larry=(mnx+2)*(mnz+1)

c  Thermodynamic diagnose in physical space
c  Description: subroutine thermodiagnost(sigg(1),xi(2),eta(3),
c                               temp(6),p_sigg(7),p_xi(8),p_eta(9),mx)
      call thermodiagnost(term(-1,0,1),term(-1,0,2),term(-1,0,3),
     &     term(-1,0,6),term(-1,0,7),term(-1,0,8),term(-1,0,9),larry)

c  Description: subroutine compute_pressure(xi(2),eta(3),temp(6),p_sigg(7),etac,mx)
      call compute_pressure(term(-1,0,2),term(-1,0,3),term(-1,0,6),
     &     term(-1,0,7),term(-1,0,8),larry)

c      call diverg(term(-1,0,4),term(-1,0,5),term(-1,0,9),
c     &     mnx,mnz,dox,doz)


      z01  =  doz/2.
      h1   =  50.
      x01  =  dox/2.
      hx1  =  200.
      z02  =  doz/2.
      hz2  =  200.

      do i=0,mnx,1
         x(i)=0.+(dox/float(mnx))*float(i)
      enddo

      do j=0,mnz,1
         z(j)=0.+(doz/float(mnz))*float(j)
      enddo

      do j=0,mnz,1            !!! Gaussian shape heating
         do i=0,mnx,1
            tp(i,j)=tpmax*exp(-((x(i)-x01)/hx1)**2)
     &      *exp(-(( z(j)-z02)/hz2)**2)
         enddo
      enddo


c  Minus background profile
      do j=0,mnz,1
         do i=0,mnx,1
            termm(i,j,1)=term(i,j,1)-bar(j,1)              !!! entropy
            termm(i,j,2)=term(i,j,2)-bar(j,2)              !!! xi
            termm(i,j,3)=term(i,j,3)-bar(j,3)              !!! eta
            termm(i,j,4)=term(i,j,4)                       !!! uu
            termm(i,j,5)=term(i,j,5)                       !!! ww
            termm(i,j,6)=term(i,j,6)-bar(j,6)              !!! temperature(t')
c            termm(i,j,6)=term(i,j,6)-bar(j,6)-tp(i,j)      !!! temperature(t'')
            termm(i,j,7)=term(i,j,7)/100.-bar(j,7)         !!! pressure
            termm(i,j,8)=term(i,j,8)                       !!! etac
            termm(i,j,9)=term(i,j,9)                       !!! divergence
         enddo
      enddo

c  The reason to do this is because
c  when use -r8 to compile, we cannot wite a file from index 0.
      do j=1,mnz+1,1
         do i=1,mnx+1,1
            temp(i,j,1)=termm(i-1,j-1,1)
            temp(i,j,2)=termm(i-1,j-1,2)
            temp(i,j,3)=termm(i-1,j-1,3)
            temp(i,j,4)=termm(i-1,j-1,4)
            temp(i,j,5)=termm(i-1,j-1,5)
            temp(i,j,6)=termm(i-1,j-1,6)
            temp(i,j,7)=termm(i-1,j-1,7)
            temp(i,j,8)=termm(i-1,j-1,8)
            temp(i,j,9)=termm(i-1,j-1,9)
         enddo
      enddo

            print*,'                       '
            print*,'time       =',time,'  seconds'
            print*,'physical space values.'
            print*,'sigs(16,16)=',temp(16,16,1)
            print*,'sigs(22,22)=',temp(22,22,1)
            print*,'xis (16,16)=',temp(16,16,2)
            print*,'xis (22,22)=',temp(22,22,2)
            print*,'etas(16,16)=',temp(16,16,3)
            print*,'etas(22,22)=',temp(22,22,3)
            print*,'uds (16,16)=',temp(16,16,4)
            print*,'uds (22,22)=',temp(22,22,4)
            print*,'wds (16,16)=',temp(16,16,5)
            print*,'wds (22,22)=',temp(22,22,5)
            print*,'temp(16,16)=',temp(16,16,6)
            print*,'temp(22,22)=',temp(22,22,6)
            print*,'pres(16,16)=',temp(16,16,7)
            print*,'pres(22,22)=',temp(22,22,7)
            print*,'etac(16,16)=',temp(16,16,8)
            print*,'etac(22,22)=',temp(22,22,8)

      nrec=int(time/delt)/ito+1
      do j=1,mnz+1,1
         do i=1,mnx+1,1

c  Binary write.
c            write(21,rec=nrec)temp(i,j,1)   !!! entropy
c            write(22,rec=nrec)temp(i,j,2)   !!! dry air density
c            write(23,rec=nrec)temp(i,j,3)   !!! moist density
c            write(24,rec=nrec)temp(i,j,4)   !!! u
c            write(25,rec=nrec)temp(i,j,5)   !!! w

c  Asc II write.
             write(31,*)temp(i,j,1)   !!! entropy
             write(32,*)temp(i,j,2)   !!! dry air density
             write(33,*)temp(i,j,3)   !!! moist density
             write(34,*)temp(i,j,4)   !!! u
             write(35,*)temp(i,j,5)   !!! w
             write(36,*)temp(i,j,6)   !!! temp
             write(37,*)temp(i,j,7)   !!! pressure
             write(38,*)temp(i,j,8)   !!! condensed liquid water

         enddo
      enddo


      return
      end
            















c************************************************
      subroutine delxz

c------------------------------------------------
c  Calculate all the regular grid points for the output.
c------------------------------------------------

      include 'c2o.f'

      dx=dox/float(mnx)
      dz=doz/float(mnz)

      z(0)=0.0
      do i=1,mnz,1
         z(i)=z(i-1)+dz
      enddo

      x(0)=0.0
      do j=1,mnx,1
         x(j)=x(j-1)+dx
      enddo

      x(-1)=x(mnx)

      return
      end





c************************************************
      subroutine set_lanczos_filter

      include 'c2o.f'

      pi=acos(-1.0)

      fitc(0)=1.0

      do i=1,nz,1
         aa=i*pi/nz
         fitc(i)=sin(aa)/aa
      enddo

      fitc2(0)=1.0

      do i=1,mz,1
         aa=i*pi/mz
         fitc2(i)=sin(aa)/aa
      enddo

      fitf(-1)=1.0
      fitf(0) =1.0

      do i=1,nx/2
         aa=i*2*pi/nx
         fitf(i*2-1) = sin(aa)/aa
         fitf(i*2)   = fitf(i*2-1)
      enddo

      do j=0,nz,1
         do i=-1,nx,1
            filt(i,j)=fitf(i)*fitc(j)
         enddo
      enddo

      return
      end





c************************************************
      subroutine lanc_filter

      include 'c2o.f'

      do k=1,5,1
         do j=0,nz,1
            do i=-1,nx,1
               f(i,j,k)=f(i,j,k)*filt(i,j)
            enddo
         enddo
      enddo

      return
      end





c************************************************
      subroutine boundary(wk)

      include 'c2o.f'

      common sum1(-1:nx),sum2(-1:nx)
      real wk(-1:nx,0:nz)

      do i=-1,nx,1
         sum1(i)=0.0
         sum2(i)=0.0
      enddo

      wku=0.0
      wkb=0.0

      do i=-1,nx,1
         sign=1.0
         do j=0,nz-2,1
            sum1(i)=sum1(i)+wk(i,j)
            sum2(i)=sum2(i)+wk(i,j)*sign
            sign=sign*(-1.0)
         enddo
      enddo

      wk(-1,nz-1)=(sum2(-1)-sum1(-1))*0.5
c     &             +(wku-wkb)*0.5
      wk(-1,nz)  =(sum2(-1)+sum1(-1))*(-0.5)
c     &             -(wku-wkb)*0.5
     
      do i=0,nx,1
         wk(i,nz-1)=(sum2(i)-sum1(i))*0.5
         wk(i,nz)  =(sum2(i)+sum1(i))*(-0.5)
      enddo

      return
      end





c************************************************
      subroutine waveno(rk)

      include 'c2o.f'

      real rk(-1:nx)

      pi=acos(-1.0)
      ax=2.0*pi/dox

      do i=-1,nx,1
         ii=(i+1)/2
         rk(i)=(float((i+1)/2)*ax)**2
      enddo

      return
      end





c************************************************
      subroutine compute_entropy(xi,eta,sigg,t,mx)
      real xi(1:mx),eta(1:mx),sigg(1:mx),t(1:mx)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      do i=1,mx,1
         ssig1  =  sigma1(xi(i),eta(i),t(i))
         ssig2  =  sigma2(xi(i),eta(i),t(i))
         if ((eta(i)-seta(t(i))) .le. 0.) then
            sigg(i)=ssig1
         else
            sigg(i)=ssig2
         endif
      enddo
      return
      end





c************************************************
      subroutine thermodiagnost(sigg,xi,eta,t,p1,p2,p3,mx)

      real xi(1:mx),eta(1:mx),sigg(1:mx),p1(1:mx),p2(1:mx),
     &     p3(1:mx),t(1:mx),xi0
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/ xlamd0,seta0
      common/constemp/ dddt
 
      do i=1,mx,1
         call compute_t1_t2(xi(i),eta(i),sigg(i),t1,t2)
         call compute_p_coeff(xi(i),eta(i),t1,t2,
     &     p11,p12,p13,p21,p22,p23)
         t(i)=max(t1,t2)
         xomeg1=(1+tanh((t1-t2)/dddt))/2.
         xomeg2=1-xomeg1
         p1(i)=p11*xomeg1 + p21*xomeg2
         p2(i)=p12*xomeg1 + p22*xomeg2
         p3(i)=p13*xomeg1 + p23*xomeg2
      enddo
      return
      end





c************************************************
      subroutine compute_t1_t2(xi,eta,sigg,t1,t2)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/ xlamd0,seta0
      real xi,eta,sigg,t1,t2,dt,tn,tm

      t1=t0*(xi/xi0)**( xi*ra/(xi*cva+eta*cvv))*
     &   (eta/seta0)**(eta*rv/(xi*cva+eta*cvv))*
     &   exp((sigg-eta*xlamd0)/(xi*cva+eta*cvv))

      tm=t0+10.0
      icount=4

      do i=1,icount,1
         tn=tm
         g1=xi*cva*alog(tn/t0)-xi*ra*alog(xi/xi0)+
     &      eta*cfun(tn)+dfun(tn)-sigg
         g2=xi*cva/tn+eta*dcfun(tn)+ddfun(tn)
         g3=eta*ddcfun(tn)+dddfun(tn)-xi*cva/tn/tn
         dt= -1.0 * g1 / g2 * (1.0-0.5*g1*g3/g2/g2)
         tm = tn+dt
         if (tm .le. 100.) then
            tm=tm+222.
c            i=0     !!!!!!!!! problem
         endif
      enddo
      t2=tm
      return
      end





c************************************************
      subroutine compute_p_coeff(xi,eta,t1,t2,p11,p12,p13,p21,p22,p23)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/ xlamd0,seta0
 
      p11=p1sig(xi,eta,t1)
      p12=p1xi (xi,eta,t1)
      p13=p1eta(xi,eta,t1)
      p21=p2sig(xi,eta,t2)
      p22=p2xi (xi,eta,t2)
      p23=p2eta(xi,eta,t2)

      return
      end





c************************************************
      subroutine compute_pressure(xi,eta,tt,p,etac,mx)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/ xlamd0,seta0
      real xi(1:mx),eta(1:mx),tt(1:mx),p(1:mx),etac(1:mx)
      real psuetac
 
      do i=1,mx,1
         etav=min(eta(i),seta(tt(i)))
         psuetac=eta(i)-etav
         etac(i)=(psuetac+abs(psuetac))*0.5
         p(i)=(xi(i)*ra+etav*rv)*tt(i)
      enddo
      return
      end





c************************************************
      function efun(t)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real efun,t
      efun=ew(t)
      return
      end
c************************************************
      function ew(t)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real ew,t,aa
      aa=(t/t0)**((cpv-cw)/rv)*exp((cpv-cw)/rv*(t0-t)/t)
      ew=e0*exp(xlw0/rv/t/t0*(t-t0))*aa
      return
      end
c************************************************
      function xl(t)
c------------------------------------------------
c  L(T) Latent heat of condensation as a function of temprature.
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real xl,t
      xl=xlw(t)
      return
      end
c************************************************
      function dxl(t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.6)
c  Kirchhoff equation,
c  Theoretically, cw is also a function of temprature,
c  but right now we ignore this relation because
c  it is nearly constant.
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real dxl,t
      dxl=dxlw(t)
      return
      end
c************************************************
      function xlw(t)
c------------------------------------------------
c  L(T) Latent heat of condensation as a function of temprature.
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real xlw,t
      xlw=(cpv-cw)*(t-t0)+xlw0
      return
      end
c************************************************
      function dxlw(t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.6)
c  Kirchhoff equation,
c  Theoretically, cw is also a function of temprature,
c  but right now we ignore this relation because
c  it is nearly constant.
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real dxlw,t
      dxlw=cpv-cw
      return
      end
c************************************************
      function seta(t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.8)
c  Mass density of saturated vapor.
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real seta,t
      seta=efun(t)/rv/t
      return
      end
c************************************************
      function xlamd(t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.9)
c  Specific entropy gain by vaporizing condensate.
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real xlamd,t
      xlamd=xl(t)/t
      return
      end
c************************************************
      function sigma1(xi,eta,t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.16)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real sigma1,xi,eta,t
      sigma1=(xi*cva+eta*cvv)*alog(t/t0)+eta*xlamd0
     &  -xi*ra*alog(xi/xi0)-eta*rv*alog(eta/seta0)
      return
      end
c************************************************
      function sigma2(xi,eta,t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.17)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real sigma2,xi,eta,t
      sigma2=xi*cva*alog(t/t0)
     &  -xi*ra*alog(xi/xi0)+eta*cfun(t)+dfun(t)
      return
      end
c************************************************
      function cfun(t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.11b)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real cfun,t
      cfun=cvv*alog(t/t0)-rv*alog(seta(t)/seta0)-xl(t)/t+xlamd0
      return
      end
c************************************************
      function dfun(t)
c------------------------------------------------
c  Ooyama (1990) Equation (3.10)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real dfun,t
      dfun=efun(t)*xl(t)/rv/t/t
      return
      end
c************************************************
      function ddfun(t)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real ddfun,t
      ddfun=dxl(t)*efun(t)/rv/t/t+
     &      xl(t)*efun(t)/rv/t/t/t*(xl(t)/rv/t-2.0)
      return
      end
c************************************************
      function dcfun(t)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real dcfun,t
      dcfun=(cpv-dxl(t))/t
      return
      end
c************************************************
      function dddfun(t)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real dddfun,t
      dddfun=3.0*efun(t)/rv/t/t/t*(xl(t)/rv/t-1.0)*
     &       (dxl(t)-2.0*xl(t)/t)+
     &       efun(t)*(xl(t)/rv/t/t)**3-
     &       efun(t)*dxl(t)/rv/t/t/t
      return
      end
c************************************************
      function ddcfun(t)
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real ddcfun,t
      ddcfun=(cpv-dxl(t))/t/t*(-1.0)
      return
      end
c************************************************
      function p1sig(xi,eta,t)
c------------------------------------------------
c  Kuo and Cheng (1999) Equation (B.4a)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real p1sig,xi,eta,t
      p1sig=(xi*ra+eta*rv)/(xi*cva+eta*cvv)*t
      return
      end
c************************************************
      function p1xi(xi,eta,t)
c------------------------------------------------
c  Kuo and Cheng (1999) Equation (B.4b)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real p1xi,xi,eta,t
      p1xi=ra*t+
     &  (ra*(1.+alog(xi/xi0))-cva*alog(t/t0))*p1sig(xi,eta,t)
      return
      end
c************************************************
      function p1eta(xi,eta,t)
c------------------------------------------------
c  Kuo and Cheng (1999) Equation (B.4c)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real p1eta,xi,eta,t
      p1eta=rv*t+(rv*(1.0+alog(eta/seta0))-cvv*alog(t/t0)-xlamd0)*
     &  p1sig(xi,eta,t)
      return
      end
c************************************************
      function p2sig(xi,eta,t)
c------------------------------------------------
c  Kuo and Cheng (1999) Equation (B.8a)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real p2sig,xi,eta,t
      p2sig=(xi*ra+dfun(t))*t/(xi*cva+(eta*dcfun(t)+ddfun(t))*t)
      return
      end
c************************************************
      function p2xi(xi,eta,t)
c------------------------------------------------
c  Kuo and Cheng (1999) Equation (B.8b)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real p2xi,xi,eta,t
      p2xi=ra*t+
     &  (ra*(1.+alog(xi/xi0))-cva*alog(t/t0))*p2sig(xi,eta,t)
      return
      end
c************************************************
      function p2eta(xi,eta,t)
c------------------------------------------------
c  Kuo and Cheng (1999) Equation (B.8c)
c------------------------------------------------
      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa
      common/const1/xlamd0,seta0
      real p2eta,xi,eta,t
      p2eta=cfun(t)*p2sig(xi,eta,t)*(-1.0)
      return
      end






c************************************************ program end.








c------------------------------------------------
c  The accompany include file "c2o.f" (for reference)
c  The following scripts should be saved as another file: c2o.f
c  to be pratically used.





c      parameter ( mx=48, mz=48, nx=mx/3*2, nz=mz/3*2 )
c      parameter ( mnx=nx, mnz=nz )

c      common/model1/ sigs(-1:nx,0:nz), xis(-1:nx,0:nz),
c     & etas(-1:nx,0:nz), uds(-1:nx,0:nz), wds(-1:nx,0:nz),
c     & imx(13), imz(13), fmx(3*mx/2+1), fmz(2*mz),
c     & f(-1:nx,0:nz,1:5)

c      common/pressu/ pp(-1:mx,0:mz), pxp(-1:mx,0:mz+1),
c     & pzp(-1:mx,0:mz+1), bf(-1:mx,0:mz+1),p1(-1:mx,0:mz),
c     & p2(-1:mx,0:mz),p3(-1:mx,0:mz)

c      common/temperature/ tt(-1:mx,0:mz)

c      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa,
c     & ggg

c      common/const1/ xlamd0,seta0

c      common/constemp/ dddt

c      common/basic/ bar(0:mnz,8), xibar(-1:mnx,0:mnz),
c     &                            ttbar(-1:mnx,0:mnz)

c      common/ex/ zi,h
c      common/job/ dox, doz, dzz
c      common/mmm/ zj(0:mz), xxi(-1:mx)
c      common/diffusion/ edx, edz
c      common/ou/ x(-1:mnx), z(0:mnz)
c      common/filter/ fitf(-1:nx), fitc(0:nz), filt(-1:nx,0:nz),
c     & fitc2(0:mz)

c      common/timeset/ timax,delt
c      common/output/ ito
c      common/temppert/ tpmax





c************************************************ The End, Be Happy !!















