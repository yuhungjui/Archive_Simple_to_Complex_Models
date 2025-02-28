      program wind_curve

      integer my,days
      parameter(my=96,days=3001)
      real zaw(my,days),zawsum(my),tazaw(my),lati(my)
      real xcoord,ycoord
      real temp,zawmax,zawmin,tazawmax,tazawmin
      character*5 xlabel(9),ylabel(7)

      open(21,file='zaw.dat')
      open(22,file='tazaw.dat')

      do i=1,my,1
         do j=1,days,1
            read(21,*)zaw(i,j)
         enddo
      enddo

      temp=zaw(1,1)
      do i=1,my,1
         do j=1,days,1
            if (zaw(i,j) .gt. temp) then
               temp=zaw(i,j)
            endif
         enddo
      enddo
      zawmax=temp
      print*,'zawmax  =',zawmax

      temp=zaw(1,1)
      do i=1,my,1
         do j=1,days,1
            if (zaw(i,j) .lt. temp) then
               temp=zaw(i,j)
            endif
         enddo
      enddo
      zawmin=temp
      print*,'zawmin  =',zawmin
      

      do i=1,my,1
         zawsum(i)=0.
         do j=1,days,1
            zawsum(i)=zawsum(i)+zaw(i,j)
         enddo
         tazaw(i)=zawsum(i)/real(days)
         write(22,*)tazaw(i)
      enddo

      temp=tazaw(1)
      do i=1,my,1
            if (tazaw(i) .gt. temp) then
               temp=tazaw(i)
            endif
      enddo
      tazawmax=temp
      print*,'tazawmax=',tazawmax

      temp=tazaw(1)
      do i=1,my,1
            if (tazaw(i) .lt. temp) then
               temp=tazaw(i)
            endif
      enddo
      tazawmin=temp
      print*,'tazawmin=',tazawmin

      do i=1,my,1
         lati(i)= -90.+ 180./95.*i
      enddo

c*******************************************************************

      call opngks

      call gscr(1,0,1.,1.,1.)      !   white background
      call gslwsc(1.)              !   thickness
      call gscr(1,1,0.,0.,0.)      !   black words
      call set(0.15,0.85,0.15,0.85,-20.,20.,-90.,90.,1)
      call perim(4,2,2,3)

c      call gscr(1,1,1.,0.,0.)
      call curve(tazaw(1),lati(1),96)
      call gscr(1,1,0.35,0.35,0.35)
      call line(0.,-90.,0.,90.)
      call line(-20.,0.,20.,0.)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call gscr(1,1,0.,0.,0.)
      call pcseti('FN',22)

      data xlabel /"-20","-15","-10"," -5","  0","  5"," 10",
     +             " 15"," 20"/ 
      data ylabel /"90S","60S","30S"," EQ","30N","60N","90N"/

      do i=1,9,1
      xcoord=0.152+(i-1)*0.0875
      ycoord=0.10
      call plchhq(xcoord,ycoord,xlabel(i),24.,0,0)
      enddo

      do j=1,7,1
      xcoord=0.108
      ycoord=0.15+(j-1)*0.1167
      call plchhq(xcoord,ycoord,ylabel(j),24.,0,0)
      enddo

      call pcseti('FN',22)
      call plchhq(0.5,0.04,"Velocity (m/s)",30.,0.,0)
      call plchhq(0.025,0.5,"Latitude",30.,90.,0)

c      call plchhq(0.5,0.94,"Time-mean zonal-mean zonal wind",26.,0.,0)
c      call plchhq(0.40,0.88,"Time Step=",22.,0.,0)
c      call plchhq(0.60,0.88,time,24.,0.,0)

c      call plchhq(0.575,0.27,"N=64",24.,0.,0)
c      call plchhq(0.615,0.23,"dt=0.0001",24.,0.,0)
c      call plchhq(0.65,0.19,"domain=[-1,1]",24.,0.,0)

      call frame
      call clsgks

      stop
      end
