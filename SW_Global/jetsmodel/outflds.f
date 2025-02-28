      subroutine outflds(taux)
c
c  output grid point value to a file
c
      include '../include/param.h'
      include '../include/const.h'
      include '../include/grid.h'
c
c  working array
c
      dimension w2(nx,my),w3(nx,my)
c
      itau=taux
      nrec=itau/tauo+1
      print *,'outfild tau=',itau,' nrec=',nrec
c

      do 90 j=1,my
      xxc=rad/cosl(j)
      do 90 i=1,nx
         w2(i,j)=uu(i,1,j)*xxc     !!!
         w3(i,j)=vv(i,1,j)*xxc     !!!
 90   continue 
      write(11,rec=nrec)w2
      write(12,rec=nrec)w3
      print *,'outflds uu,vv=',w2(49,58),w3(145,37)    ! llming

c
      do 100 j=1,my
      do 100 i=1,nx
        w2(i,j)=vor(i,1,j)
        w3(i,j)=phi(i,1,j)/grav
 100  continue
      write(13,rec=nrec)w2
      write(14,rec=nrec)w3
      print *,'outflds vor=',w2(49,58),w2(145,37)    ! llming


c
c  compute potential vorticity
c
      do 101 j=1,my
      do 101 i=1,nx
        w3(i,j)=(vor(i,1,j)+cor(j))/((hmean+phi(i,1,j)-topo(i,j))/grav)
 101  continue
      write(15,rec=nrec)w3

      return
      end
