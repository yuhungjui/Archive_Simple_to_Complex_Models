      program transformation_8bytes_to_4bytes 

      integer n,nx,my,nxmy8,nxmy4
      parameter(n=4,nx=192,my=96)    ! change

      double precision phi(nx,my,n)
      double precision vor(nx,my,n)
      double precision   u(nx,my,n)
      double precision   v(nx,my,n)
      double precision   pv(nx,my,n)
      double precision topo(nx,my)
 
      real phi2(nx,my,n)   
      real vor2(nx,my,n)   
      real   u2(nx,my,n)   
      real   v2(nx,my,n)   
      real  pv2(nx,my,n)   
      real topo2(nx,my)   
 
      nxmy4=nx*my*4
      nxmy8=nx*my*8

c      open(21,file='./8/PHI.OUT',access='direct',
c     +  form='unformatted',recl=nxmy8,status='unknown')

      open(22,file='./8/VOR.OUT',access='direct',
     +  form='unformatted',recl=nxmy8,status='unknown')
      open(23,file='./8/U.OUT',access='direct',
     +  form='unformatted',recl=nxmy8,status='unknown')
      open(24,file='./8/V.OUT',access='direct',
     +  form='unformatted',recl=nxmy8,status='unknown')

c      open(25,file='./8/PV.OUT',access='direct',
c     +  form='unformatted',recl=nxmy8,status='unknown')
c      open(26,file='./8/TOPO3.OUT',access='direct',
c     +  form='unformatted',recl=nxmy8,status='unknown')

      do k=1,n,1
c         read(21,rec=k) ,((phi(i,j,k),i=1,nx),j=1,my)
         read(22,rec=k) ,((vor(i,j,k),i=1,nx),j=1,my)
         read(23,rec=k) ,((u(i,j,k),i=1,nx),j=1,my)
         read(24,rec=k) ,((v(i,j,k),i=1,nx),j=1,my)
c         read(25,rec=k) ,((pv(i,j,k),i=1,nx),j=1,my)
      enddo

c      read(26,rec=1) ,((topo(i,j),i=1,nx),j=1,my)

      print *,'vor',vor(96,48,1)
      print *,'u',u(96,48,1)
      print *,'v',v(96,48,1)
   
      do i=1,nx,1
         do j=1,my,1
            do k=1,n,1
c               phi2(i,j,k)=phi(i,j,k)
               vor2(i,j,k)=vor(i,j,k)
               u2(i,j,k)=u(i,j,k)
               v2(i,j,k)=v(i,j,k)
c               pv2(i,j,k)=pv(i,j,k)
c               topo2(i,j)=topo(i,j)
            enddo
         enddo
      enddo

c      open(31,file='./4/PHI.dat',access='direct',
c     +  form='unformatted',recl=nxmy4,status='unknown')

      open(32,file='./4/VOR.dat',access='direct',
     +  form='unformatted',recl=nxmy4,status='unknown')
      open(33,file='./4/U.dat',access='direct',
     +  form='unformatted',recl=nxmy4,status='unknown')
      open(34,file='./4/V.dat',access='direct',
     +  form='unformatted',recl=nxmy4,status='unknown')

c      open(35,file='./4/PV.dat',access='direct',
c     +  form='unformatted',recl=nxmy4,status='unknown')
c      open(36,file='./4/TOPO3.dat',access='direct',
c     +  form='unformatted',recl=nxmy4,status='unknown')

      do k=1,n,1
c         write(31,rec=k),((phi2(i,j,k),i=1,nx),j=1,my) 
         write(32,rec=k),((vor2(i,j,k),i=1,nx),j=1,my) 
         write(33,rec=k),((u2(i,j,k),i=1,nx),j=1,my) 
         write(34,rec=k),((v2(i,j,k),i=1,nx),j=1,my) 
c         write(35,rec=k),((pv2(i,j,k),i=1,nx),j=1,my) 
      enddo

c      write(36,rec=1),((topo2(i,j),i=1,nx),j=1,my) 

      stop
      end

