      subroutine rhdiffu( dta,my,nx,mlmax,lev,hfilt
     1                  , vornow,divnow,temnow
     2                  , eps4)
c

      dimension vornow(mlmax*2,2,lev),divnow(mlmax*2,2,lev)
     2        , temnow(mlmax*2,2,lev)
     3        , eps4(mlmax)
c
      do 100 k=1,lev
c
c  compute diffusion coefficients
c
      vdiffu =hfilt
c
c  difuse vorticity and divergence fields
c  diffuse moisture and temperature fields
c
      do 30 ml=1,mlmax
        mn=ml+mlmax
        c1=1.+dta*vdiffu*eps4(ml)**2
        c2=1.+dta*vdiffu*eps4(ml)**2
        c3=1.+dta*vdiffu*eps4(ml)**2
        vornow(ml,1,k)=vornow(ml,1,k)/c1
        vornow(ml,2,k)=vornow(ml,2,k)/c1
        vornow(mn,1,k)=vornow(mn,1,k)/c1
        vornow(mn,2,k)=vornow(mn,2,k)/c1
        divnow(ml,1,k)=divnow(ml,1,k)/c2
        divnow(ml,2,k)=divnow(ml,2,k)/c2
        divnow(mn,1,k)=divnow(mn,1,k)/c2
        divnow(mn,2,k)=divnow(mn,2,k)/c2
        temnow(ml,1,k)=temnow(ml,1,k)/c2
        temnow(ml,2,k)=temnow(ml,2,k)/c2
        temnow(mn,1,k)=temnow(mn,1,k)/c2
        temnow(mn,2,k)=temnow(mn,2,k)/c2
 30   continue
c
 100  continue
c
      return
      end
