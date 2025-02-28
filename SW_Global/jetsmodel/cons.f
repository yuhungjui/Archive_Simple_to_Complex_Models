      subroutine cons
      include '../include/param.h'
      include '../include/const.h'
      include '../include/fftcom.h'
c
      rad=6.37e6
      radsq=rad*rad
      grav=9.81
      omega=7.292e-5
      cp=1004.5
      rgas=287.
      radsq = rad**2
      pi=4.0*atan(1.0)
      capa=1.0/3.5
      its=taue*3600./dt
c
c  build pointer arrays for locating zonal and total wavenumber
c  values in the one-dimensional spherical harmonic arrays.
c
      call sortml(jtrun,mlmax,msort,lsort,mlsort)
c
      do 150 ml=1,mlmax
      rl= lsort(ml)
      rm= msort(ml)-1
      rlm= rl-1.0
      if(msort(ml).eq.1) rm= 0.0
      if(lsort(ml).eq.1) rlm= 0.0
      eps4(ml)= rl*rlm/radsq
      epsd16(ml)= (eps4(ml))**8
      cim(ml)= rm
  150 continue
c
c  gaussian quadrature weights and latitudes
c
      one= 1.0
      onem= -one
      call gausl3(my,onem,one,weight,sinl)
c
      jm2= my/2
cdir$ ivdep
      do 155 j=1,jm2
      sinl(my+1-j)= -sinl(j)
      weight(my+1-j)= weight(j)
      onocos(j)= 1.0/(1.0-sinl(j)*sinl(j))
      onocos(my+1-j)= onocos(j)
      cosl(j)= 1.0/sqrt(onocos(j))
      cosl(my+1-j)= cosl(j)
  155 continue
c
c  coriolis parameter for each latitude
c
      do 30 j=1,my
      cor(j)= 2.0*omega*sinl(j)
  30  continue
c
      call fftfax(nx,ifax,trigs)
c
c  associated legendre polynomials and their derivatives
c
      call lgndr(jm2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
c
      return
      end
