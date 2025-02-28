      subroutine uven0

c---------------------------------------------------------------
c  select wanted wavenumber to calculate kinetic energy.  
c
c  written by LLMing 
c  2005 Aug. 21
c
c---------------------------------------------------------------

      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'

      real u0spc(mlmax,2,lev),v0spc(mlmax,2,lev)    !!!! llming
      real uch(nx,lev,my),vch(nx,lev,my)              !!!! llming
      real en(nx,my)
      real ensum,enave(mlmax)

      open(45,file='enave0.dat',status='unknown',form='formatted')


c----------------- mlmax=1 ---------------------------------------

      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,uu,u0spc)
      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,vv,v0spc)

      do i=2,mlmax,1
         do j=1,2,1
            u0spc(i,j,1)=0.
            v0spc(i,j,1)=0.
         enddo
      enddo

      call transr(jtrun,mlmax,nx,my,lev,poly,u0spc,uch)
      call transr(jtrun,mlmax,nx,my,lev,poly,v0spc,vch)

      do i=1,nx,1
         do j=1,my,1
            en(i,j)=0.5*(uch(i,1,j)**2+vch(i,1,j)**2)
         enddo
      enddo

      ensum=0.
      do i=1,nx,1
         do j=1,my,1
            ensum=ensum+en(i,j)
         enddo
      enddo

      enave(1)=ensum/real(nx*my)

c----------------- mlmax=1 ---------------------------------------




c----------------- mlmax=2 - 2079 --------------------------------

      do k=2,2079,1   

      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,uu,u0spc)
      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,vv,v0spc)

         do i=1,k-1,1
            do j=1,2,1
               u0spc(i,j,1)=0.
               v0spc(i,j,1)=0.
            enddo
         enddo

         do i=k+1,mlmax,1
            do j=1,2,1
               u0spc(i,j,1)=0.
               v0spc(i,j,1)=0.
            enddo
         enddo

         call transr(jtrun,mlmax,nx,my,lev,poly,u0spc,uch)
         call transr(jtrun,mlmax,nx,my,lev,poly,v0spc,vch)

         do i=1,nx,1
            do j=1,my,1
               en(i,j)=0.5*(uch(i,1,j)**2+vch(i,1,j)**2)
            enddo
         enddo

         ensum=0.
         do i=1,nx,1
            do j=1,my,1
               ensum=ensum+en(i,j)
            enddo
         enddo

         enave(k)=ensum/real(nx*my)

      enddo    

c----------------- mlmax=2 - 2079 --------------------------------



c----------------- mlmax= 2080 -----------------------------------

      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,uu,u0spc)
      call tranrs(jtrun,mlmax,nx,my,lev,poly,weight,vv,v0spc)

      do i=1,2079,1
         do j=1,2,1
            u0spc(i,j,1)=0.
            v0spc(i,j,1)=0.
         enddo
      enddo

      call transr(jtrun,mlmax,nx,my,lev,poly,u0spc,uch)
      call transr(jtrun,mlmax,nx,my,lev,poly,v0spc,vch)

      do i=1,nx,1
         do j=1,my,1
            en(i,j)=0.5*(uch(i,1,j)**2+vch(i,1,j)**2)
         enddo
      enddo

      ensum=0.
      do i=1,nx,1
         do j=1,my,1
            ensum=ensum+en(i,j)
         enddo
      enddo

      enave(2080)=ensum/real(nx*my)

c----------------- mlmax=2080  -------------------------------


      do k=1,2080,1
         write(45,*) enave(k)
      enddo

      return
      end
