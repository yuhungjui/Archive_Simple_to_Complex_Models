c
      subroutine outflds(taux)
      include '../include/rparam.h'
      include '../include/rconst.h'
      include '../include/rgrid.h'
      dimension p(nxr,myr),u(nxr,myr),v(nxr,myr)
     1,         vo(nxr,myr),di(nxr,myr),s(nxr,myr)
c
      if(taux.eq.0.)then
        do 90 j=1,myr
        do 90 i=1,nxr
          p(i,j)=ttm(i,1,j)/gravr
          u(i,j)=um(i,1,j)
          v(i,j)=vm(i,1,j)
          vo(i,j)=vorm(i,1,j)
          di(i,j)=divm(i,1,j)
 90    continue
        write(35,rec=1)p
        write(36,rec=1)u
        write(37,rec=1)v
        write(38,rec=1)vo
        write(40,rec=1)di
        close(35)
        close(36)
        close(37)
        close(38)
        close(40)
      endif
c
      do 100 j=1,myr
      do 100 i=1,nxr
        p(i,j)=ttr(i,1,j)/gravr
        u(i,j)=utr(i,1,j)
        v(i,j)=vtr(i,1,j)
        vo(i,j)=rvorr(i,1,j)
        s(i,j)=pvr(i,1,j)
        di(i,j)=rdivr(i,1,j)
 100  continue
c
      call qmaxn3(p,'regional','output p',1,1,1,nxr,myr,1)
      call qmaxn3(u,'regional','output u',1,1,1,nxr,myr,1)
      call qmaxn3(v,'regional','output v',1,1,1,nxr,myr,1)
c
      itauo=taux+0.01
      nxmyr=nxr*myr
      nrec=taux/tauor+1.001
      write(21,rec=nrec)p
      write(22,rec=nrec)u
      write(23,rec=nrec)v
      write(24,rec=nrec)vo
      write(25,rec=nrec)s
      write(26,rec=nrec)di
c
      return
      end
