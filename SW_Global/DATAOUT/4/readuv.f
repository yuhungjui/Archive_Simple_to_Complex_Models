      program readwind
 
      integer nx,my,time
      parameter(nx=192,my=96,time=4)
      real u(nx,time,my)
      real v(nx,time,my)
      real sqv(nx,time,my),sqsum(time),vrms(time)

      open (21,file='U.dat',access='direct',
     +      form='unformatted',recl=nx*my*4,status='unknown')
      open (22,file='V.dat',access='direct',
     +      form='unformatted',recl=nx*my*4,status='unknown')

c      open (24,file='uasc.dat',status='unknown')
c      open (25,file='vasc.dat',status='unknown')

      do k=1,time
         read(21,rec=k),((u(i,k,j),i=1,nx),j=1,my)
         read(22,rec=k),((v(i,k,j),i=1,nx),j=1,my)
      enddo

c      do i=1,nx,1
c         do j=1,my,1
c            write(24,*)  u(i,1,j)
c            write(25,*)  v(i,1,j)
c         enddo
c      enddo

      do i=1,nx,1
         do j=1,my,1
            do k=1,time
               sqv(i,k,j)=u(i,k,j)**2+v(i,k,j)**2
            enddo
         enddo
      enddo

      do k=1,2,1
         sqsum(k)=0.
      enddo
      do i=1,nx,1
         do j=1,my,1
            do k=1,time
               sqsum(k)=sqsum(k)+sqv(i,k,j)
            enddo
         enddo
      enddo

      do k=1,time
         vrms(k)=sqrt(sqsum(k)/real(nx*my))
      enddo

      print*,'vrms( 0)=',vrms(1)
      print*,'vrms(20)=',vrms(2)
      print*,'vrms(40)=',vrms(3)
      print*,'vrms(60)=',vrms(4)

      stop
      end
