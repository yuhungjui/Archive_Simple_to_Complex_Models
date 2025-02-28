      subroutine gettopo(nxr,myr,topo)
c
      parameter (nxr1=200,myr1=400)
c
      dimension topo(nxr,myr),topoori(nxr1,myr1)
      dimension xlon(nxr1,myr1),xlat(nxr1,myr1)
      dimension topot(nxr1/2,myr1/2)
c
      data lon,lat/450,450/
      indextopo=1
c
      nxmy=nxr*myr*8
      open(61,file='../dat/terll.dat'
     &, form='formatted',status='old')
c
 111  format(3x,f8.0,7x,f8.5,7x,f8.4,4x)
      do i=1,200
      do j=1,400
      read(61,111)topoori(i,j),xlat(i,j),xlon(i,j)
c      print *,topoori(i,j),xlon(i,j),xlat(i,j)
      enddo
      enddo
c
      do i=1,nxr
      do j=1,myr
      topo(i,j)=0.
      enddo
      enddo
c
  99  continue
      close(61)
c
      if (indextopo.eq.1)then
      nnx=nxr1/2
      mmy=myr1/2
      do i=1,nnx
      do j=1,mmy
      topot(i,j)=topoori(i*2,j*2)
      enddo
      enddo
c
      ib=lon-(nnx-1)/2
      ie=lon+(nnx-1)/2+1
      jb=lat-(mmy-1)/2
      je=lat+(mmy-1)/2+1
c      print *,ib,ie,jb,je
      do i=ib,ie
      do j=jb,je
      topo(i,j)=topot(i-ib+1,j-jb+1)*9.80616
c      print *,topo(i,j)
      enddo
      enddo
c
      endif
c
      write(41,rec=1)topo
      close(41)
c
      return
c
      end
