      subroutine boundary(nxr,myr,lev,ttr,utr,vtr)
      dimension ttr(nxr,lev,myr),utr(nxr,lev,myr),vtr(nxr,lev,myr)
c
      k=1
c
      do 10 j=1,myr
       uavg=(utr(1,k,j)+utr(nxr,k,j))*0.5
       vavg=(vtr(1,k,j)+vtr(nxr,k,j))*0.5
       tavg=(ttr(1,k,j)+ttr(nxr,k,j))*0.5
       utr(1,k,j)=uavg
       utr(nxr,k,j)=uavg
       vtr(1,k,j)=vavg
       vtr(nxr,k,j)=vavg
       ttr(1,k,j)=tavg
       ttr(nxr,k,j)=tavg
 10   continue
      do 20 i=2,nxr-1
       uavg=(utr(i,k,1)+utr(i,k,myr))*0.5
       vavg=(vtr(i,k,1)+vtr(i,k,myr))*0.5
       tavg=(ttr(i,k,1)+ttr(i,k,myr))*0.5
       utr(i,k,1)=uavg
       utr(i,k,myr)=uavg
       vtr(i,k,1)=vavg
       vtr(i,k,myr)=vavg
       ttr(i,k,1)=tavg
       ttr(i,k,myr)=tavg
 20   continue
c
      return
      end
