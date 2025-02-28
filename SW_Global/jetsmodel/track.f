      subroutine track(nx,my,lev,np,ntyp,vor,ix,jy,ic,center,sinl)
      dimension center(5000,2,np),vor(nx,lev,my),sinl(my)
     *, ix(np),jy(np)
      data ixyrange/12/
c
      pi=4.0*atan(1.0)
      r2d=180./pi
      ip=ic+1
c
      do 20 n=1,ntyp
c
      ib=ix(n)-ixyrange
      ie=ix(n)+ixyrange
      jb=jy(n)-ixyrange
      je=jy(n)+ixyrange
      vormax=vor(ix(n),lev,jy(n))
c
      do 10 j=jb,je
      do 10 i=ib,ie
        if(vor(i,lev,j).gt.vormax)then
          ix(n)=i
          jy(n)=j
          vormax=vor(i,lev,j)
        endif
 10   continue
c
      center(ip,1,n)=ix(n)
      center(ip,2,n)=jy(n)
c      center(ip,1,n)=(ix(n)-1)*(360./nx)
c      center(ip,2,n)=(jy(n)-1)*(180./my)+asin(sinl(1))*r2d
c
 20   continue
c
      return
      end
