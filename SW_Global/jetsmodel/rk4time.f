      subroutine rk4time
c
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
c  working array
c
      dimension ws2(mlmax,2,lev,3)

c--------------------------------------------------llming

      integer tmark
      real sum(my)
      open(65,file='zaw.dat',status='unknown')

c------------------------------------------------------
c
      print *,' begin runkut4 time integration'
c
      dta=dt*0.5
      max2lev=mlmax*2*lev
c
      iend=taue*3600./dt+0.001
      iout=tauo*3600./dt+0.001
c
      do 100 iter=1,iend
c
      taux=iter*dt/3600.
c
      print *,'forcast tau=',taux
c
      if(linear)then
        call lrhs
      else
        call nrhs(iter)           ! llming
      endif
      do 10 i=1,max2lev
        vspcold(i,1,1)=vspcnow(i,1,1)
        dspcold(i,1,1)=dspcnow(i,1,1)
        pspcold(i,1,1)=pspcnow(i,1,1)
        ws2(i,1,1,1)=vorten(i,1,1)
        ws2(i,1,1,2)=divten(i,1,1)
        ws2(i,1,1,3)=phiten(i,1,1)
 10   continue
C
      do 20 i=1,max2lev
      vspcnow(i,1,1)=vspcold(i,1,1)+dta*vorten(i,1,1)
      dspcnow(i,1,1)=dspcold(i,1,1)+dta*divten(i,1,1)
      pspcnow(i,1,1)=pspcold(i,1,1)+dta*phiten(i,1,1)
 20   continue
c
      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
      call transr(jtrun,mlmax,nx,my,lev,poly,dspcnow,div)
      call transr(jtrun,mlmax,nx,my,lev,poly,pspcnow,phi)
      call tranuv(jtrun,mlmax,nx,my,lev,radsq,onocos,eps4,cim
     *, poly,dpoly,vspcnow,dspcnow,uu,vv)
      if(linear)then
        call lrhs
      else
        call nrhs(iter)       ! llming
      endif
C
      do 30 i=1,max2lev
        ws2(i,1,1,1)=ws2(i,1,1,1)+vorten(i,1,1)*2.0
        ws2(i,1,1,2)=ws2(i,1,1,2)+divten(i,1,1)*2.0
        ws2(i,1,1,3)=ws2(i,1,1,3)+phiten(i,1,1)*2.0
 30   continue
C
      do 40 i=1,max2lev
      vspcnow(i,1,1)=vspcold(i,1,1)+dta*vorten(i,1,1)
      dspcnow(i,1,1)=dspcold(i,1,1)+dta*divten(i,1,1)
      pspcnow(i,1,1)=pspcold(i,1,1)+dta*phiten(i,1,1)
 40   continue
c
      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
      call transr(jtrun,mlmax,nx,my,lev,poly,dspcnow,div)
      call transr(jtrun,mlmax,nx,my,lev,poly,pspcnow,phi)
      call tranuv(jtrun,mlmax,nx,my,lev,radsq,onocos,eps4,cim
     *, poly,dpoly,vspcnow,dspcnow,uu,vv)
C
      if(linear)then
        call lrhs
      else
        call nrhs(iter)       ! llming
      endif
C
      do 50 i=1,max2lev
        ws2(i,1,1,1)=ws2(i,1,1,1)+vorten(i,1,1)*2.0
        ws2(i,1,1,2)=ws2(i,1,1,2)+divten(i,1,1)*2.0
        ws2(i,1,1,3)=ws2(i,1,1,3)+phiten(i,1,1)*2.0
 50   continue
C
      do 60 i=1,max2lev
      vspcnow(i,1,1)=vspcold(i,1,1)+dt*vorten(i,1,1)
      dspcnow(i,1,1)=dspcold(i,1,1)+dt*divten(i,1,1)
      pspcnow(i,1,1)=pspcold(i,1,1)+dt*phiten(i,1,1)
 60   continue
c
      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
      call transr(jtrun,mlmax,nx,my,lev,poly,dspcnow,div)
      call transr(jtrun,mlmax,nx,my,lev,poly,pspcnow,phi)
      call tranuv(jtrun,mlmax,nx,my,lev,radsq,onocos,eps4,cim
     *, poly,dpoly,vspcnow,dspcnow,uu,vv)
C
      if(linear)then
        call lrhs
      else
        call nrhs(iter)            !  llming
      endif
c
      do 70 i=1,max2lev
        ws2(i,1,1,1)=ws2(i,1,1,1)+vorten(i,1,1)
        ws2(i,1,1,2)=ws2(i,1,1,2)+divten(i,1,1)
        ws2(i,1,1,3)=ws2(i,1,1,3)+phiten(i,1,1)
 70   continue
C
      do 80 i=1,max2lev
        vspcnow(i,1,1)=vspcold(i,1,1)+dt*ws2(i,1,1,1)/6.0
        dspcnow(i,1,1)=dspcold(i,1,1)+dt*ws2(i,1,1,2)/6.0
        pspcnow(i,1,1)=pspcold(i,1,1)+dt*ws2(i,1,1,3)/6.0
 80   continue
c
c  do diffusion
c
      if(dohdi)call hdiffu(dt)
c
c  transform spectrum coefficient to grid point
c
      call transr(jtrun,mlmax,nx,my,lev,poly,vspcnow,vor)
      call transr(jtrun,mlmax,nx,my,lev,poly,dspcnow,div)
      call transr(jtrun,mlmax,nx,my,lev,poly,pspcnow,phi)
      call tranuv(jtrun,mlmax,nx,my,lev,radsq,onocos,eps4,cim
     *, poly,dpoly,vspcnow,dspcnow,uu,vv)

c------------------------------------------------! llming
c   calculate the zonal average of zonal wind.


      if (mod(iter,2880).eq.0) then

         tmark=iter/2880*10

         do j=1,my,1
            sum(j)=0.
         enddo

         do j=1,my,1
            do i=1,nx,1
               sum(j)=sum(j)+uu(i,1,j)
            enddo
            zaw(j,tmark)=sum(j)/nx
         enddo

      endif

c------------------------------------------------------

c
      print *,'outflds vor =',vor(96,1,48)
c
c  compute stream function
c
      do 81 i=2,mlmax
        strspc(i,1)=-vspcnow(i,1,1)/eps4(i)
        strspc(i,2)=-vspcnow(i,2,1)/eps4(i)
 81   continue
      strspc(1,1)=0.
      strspc(1,2)=0.
      call transr(jtrun,mlmax,nx,my,1,poly,strspc,str)
c
      if(bogus)call track(nx,my,lev,np,ntyp,vor,ix,jy,iter
     *,  center,sinl)
c
c------------------------------------------------------llming

      if(iter .eq. 864000) then
      do j=1,my,1
         do k=0,3000,10
            write(65,*)zaw(j,k)
         enddo
      enddo
      endif

c--------------------------------------------------------llming
   

      if(mod(iter,iout) .eq. 0) then
      call outflds(taux)
      endif
c     +  call uven
      

c
 100  continue
c
c    write()center
c
      print *,' fisnsh runkut4 time integration'
C
      return
      end
