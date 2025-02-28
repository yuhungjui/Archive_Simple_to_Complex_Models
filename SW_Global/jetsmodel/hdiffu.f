      subroutine hdiffu(dta)
c
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
c  compute diffusion coefficients
c
      hfilt2=hfilt
      if(ktop.eq.2)then       ! about topography
        hfilt2=hfilt*100.
      else if(ktop.eq.1)then
        hfilt2=hfilt*20.
      endif
c
c  difuse vorticity and divergence fields
c
      k=1
      do 30 ml=1,mlmax
        c1=1.+dta*hfilt*eps4(ml)**2
        c2=1.+dta*hfilt2*eps4(ml)**2
        vspcnow(ml,1,k)=vspcnow(ml,1,k)/c1
        vspcnow(ml,2,k)=vspcnow(ml,2,k)/c1
        dspcnow(ml,1,k)=dspcnow(ml,1,k)/c2
        dspcnow(ml,2,k)=dspcnow(ml,2,k)/c2
        pspcnow(ml,1,k)=pspcnow(ml,1,k)/c1
        pspcnow(ml,2,k)=pspcnow(ml,2,k)/c1
 30   continue
c
      return
      end
