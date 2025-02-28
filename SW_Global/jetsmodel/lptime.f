      subroutine lptime
c
c  time integration use leapfog scheme
c
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'

      logical forward
      data forward/.true./

c--------------------------------------------------llming

      integer tmark
      real sum(my)
      open(65,file='zaw.dat',status='unknown')

c--------------------------------------------------------
c
c      real dta    ! llming    
c
c
      print *,'                                                   '
      print *,'c--------- begin leapfog time integration --------c'
      print *,'                                                   '
c
      max2lev=mlmax*2*lev
c
      iend=taue*3600./dt+0.001
      iout=tauo*3600./dt+0.001
c
      do 100 iter=1,iend
c
      taux=iter*dt/3600.
      tau=taux
c
cllming      print *,'forcast tau=',taux
c
      if(linear)then
        call lrhs
      else
        call nrhs(iter)
      endif
      vorten(1,1,1)=0.
      divten(1,1,1)=0.
      phiten(1,1,1)=0.
      vorten(1,2,1)=0.
      divten(1,2,1)=0.
      phiten(1,2,1)=0.
c
      if(forward)then
c
c  leapfog time integration
c
        do 10 i=1,max2lev
          vspcnow(i,1,1)=dt*vorten(i,1,1)+vspcold(i,1,1)
          dspcnow(i,1,1)=dt*divten(i,1,1)+dspcold(i,1,1)
          pspcnow(i,1,1)=dt*phiten(i,1,1)+pspcold(i,1,1)
 10    continue
c
        dta=dt*2.0
        forward=.false.
c
      else
c
c  leapfog time integration
c
        do 20 i=1,max2lev
          vorten(i,1,1)=dta*vorten(i,1,1)+vspcold(i,1,1)
          divten(i,1,1)=dta*divten(i,1,1)+dspcold(i,1,1)
          phiten(i,1,1)=dta*phiten(i,1,1)+pspcold(i,1,1)
 20     continue
c
c  do robert time filter
c
        do 30 i=1,max2lev
          vspcold(i,1,1)= vspcnow(i,1,1)+tfilter*(vspcold(i,1,1)
     *                   -2.0*vspcnow(i,1,1)+vorten(i,1,1))
          dspcold(i,1,1)= dspcnow(i,1,1)+tfilter*(dspcold(i,1,1)
     *                   -2.0*dspcnow(i,1,1)+divten(i,1,1))
          pspcold(i,1,1)= pspcnow(i,1,1)+tfilter*(pspcold(i,1,1)
     *                   -2.0*pspcnow(i,1,1)+phiten(i,1,1))
          vspcnow(i,1,1)=vorten(i,1,1)
          dspcnow(i,1,1)=divten(i,1,1)
          pspcnow(i,1,1)=phiten(i,1,1)
30    continue
c
      endif
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

c--------------------------------------------------llming    go

      do j=1,my,1
         xxc=rad/cosl(j)
         do i=1,nx,1
            uu(i,1,j)=uu(i,1,j)*xxc
            vv(i,1,j)=vv(i,1,j)*xxc
         enddo
      enddo

c------------------------------------------------! llming
c
c   calculate the zonal average of zonal wind.

      if (mod(iter,int(1.*86400./dt)).eq.0) then

         tmark=iter/int(1.*86400./dt)*1

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

c-------------------------------------------------------

      if (mod(iter,int(its)) .eq. 0)  then
         call uven80
      endif

c--------------------------------------------------llming    return

      do j=1,my,1
         xxc=rad/cosl(j)
         do i=1,nx,1
            uu(i,1,j)=uu(i,1,j)/xxc
            vv(i,1,j)=vv(i,1,j)/xxc
         enddo
      enddo

c----------------------------------------llming


cllming      print *,'lptime vor',vor(128,1,64)    ! llming

cllming      print *,'lptime uu=',uu(128,1,64)*rad/cosl(64)    ! llming
cllming      print *,'lptime vv=',vv(128,1,64)*rad/cosl(64)    ! llming

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

cllming      if(bogus)call track(nx,my,lev,np,ntyp,vor,ix,jy,iter
cllming     *,  center,sinl)

c------------------------------------------------------llming
     
      if(iter .eq. int(its)) then
      do j=1,my,1
         do k=0,int(taue/24.),1
            write(65,*)zaw(j,k)
         enddo
      enddo
      endif

c--------------------------------------------------------llming
c
      if (mod(iter,iout) .eq. 0)  then
         call outflds(taux)
      endif

c
 100  continue
c
      print *,' finish leapfor time integration'
c
      return
      end
