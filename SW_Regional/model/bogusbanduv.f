      subroutine bogusbanduv(nx,my,dx,ut,vt)
c
      dimension ut(nx,1,my),vt(nx,1,my)
c
      data wild/50000./,yj/85.001/,xi/127.501/,vormax/0.00003/
      wind=0.
      dy=dx
c
      do 1829 j=1,my
      do 1829 i=1,nx
         rrr=(j-yj)*dy
         rr=(rrr*rrr)**0.5
         if (rr.le.wild)then
         wind=0.5*vormax*rrr
         ut(i,1,j)=ut(i,1,j)+wind
c       ut(i,1,j)=wind
         endif
 1829 continue
c
      return
      end
