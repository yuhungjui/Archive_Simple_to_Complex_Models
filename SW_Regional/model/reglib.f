      subroutine rtransr(jtrunx,jtruny,nx,my,lev,s,r)
      include '../include/rparamt.h'
      include '../include/rfftcom.h'
      dimension r(nx,lev,my),s(jtrunx*jtruny*2,2,lev)
      dimension cc(nxt+3,myt),dd(myt+3,nxt),work(nxt*myt,2)
c
      jtrunxy=jtrunx*jtruny
      jtrunx2=jtrunx*2
c
      do 68 k=1,lev
c
      do 60 i=1,(nx+3)
      do 60 j=1,my
	cc(i,j)=0.
   60 continue
      do 61 i=1,(my+3)
      do 61 j=1,nx
	dd(i,j)=0.
   61 continue
c
      do 62 n=1,jtruny
	nn=n*2-1
	np=nn+1
      do 62 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)= (s(mn,1,k)+s(ms,1,k))*0.5
	dd(np,mm)= (s(mn,2,k)-s(ms,2,k))*0.5
	dd(nn,mp)= (s(mn,2,k)+s(ms,2,k))*0.5
	dd(np,mp)=-(s(mn,1,k)-s(ms,1,k))*0.5
   62 continue
c
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,jtrunx2,+1)
c
      do 64 j=1,my
      do 64 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	cc(mm,j)=dd(j,mm)
	cc(mp,j)=dd(j,mp)
   64 continue
c
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,my,+1)
c
      do 66 j=1,my
      do 66 i=1,nx
	r(i,k,j)=cc(i,j)
   66 continue
c
   68 continue
c
      return
      end
      subroutine rtranrs(jtrunx,jtruny,nx,my,lev,r,s)
      include '../include/rparamt.h'
      include '../include/rfftcom.h'
      dimension r(nx,lev,my),s(jtrunx*jtruny*2,2,lev)
      dimension cc(nxt+3,myt),dd(myt+3,nxt),work(nxt*myt,2)
c
      jtrunxy=jtrunx*jtruny
      jtrunx2=jtrunx*2
c
      do 68 k=1,lev
c
      do 58 i=1,(nx+3)
      do 58 j=1,my
	cc(i,j)=0.
   58 continue
      do 59 i=1,(my+3)
      do 59 j=1,nx
	dd(i,j)=0.
   59 continue
c
      do 60 j=1,my
      do 60 i=1,nx
	cc(i,j)=r(i,k,j)
   60 continue
c
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,my,-1)
c
      do 64 j=1,my
      do 62 m=1,jtrunx
         mm= 2*m-1
	 mp= mm+1
	 dd(j,mm) = cc(mm,j)
	 dd(j,mp) = cc(mp,j)
   62 continue
   64 continue
c
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,jtrunx2,-1)
c
      do 66 n=1,jtruny
	nn= 2*n-1
	np= nn+1
      do 66 m=1,jtrunx
	mm= 2*m-1
	mp= mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	s(mn,1,k) = dd(nn,mm)-dd(np,mp)
	s(mn,2,k) = dd(np,mm)+dd(nn,mp)
	s(ms,1,k) = dd(nn,mm)+dd(np,mp)
	s(ms,2,k) =-dd(np,mm)+dd(nn,mp)
   66 continue
c
   68 continue
c
      return
      end
      subroutine rtrngra(jtrunx,jtruny,nx,my,pi,dox,doy,s,dxpl,dypl)
c
      include '../include/rparamt.h'
      include '../include/rfftcom.h'
c
      dimension s(jtrunx*jtruny*2,2),dxpl(nx,my),dypl(nx,my)
      dimension cc(nxt+3,myt),dd(myt+3,nxt),work(nxt*myt,2)
c
      odx=2.*pi/dox
      ody=2.*pi/doy
      jtrunxy=jtrunx*jtruny
      jtrunx2=jtrunx*2
      jtruny2=jtruny*2
c da/dx
      do 70 i=1,(nx+3)
      do 70 j=1,my
	cc(i,j)=0.
   70 continue
      do 71 i=1,(my+3)
      do 71 j=1,nx
	dd(i,j)=0.
   71 continue
c y-direction
      do 72 n=1,jtruny
	nn=n*2-1
	np=nn+1
      do 72 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
        dd(nn,mm)= (s(mn,1)+s(ms,1))*0.5
        dd(np,mm)= (s(mn,2)-s(ms,2))*0.5
        dd(nn,mp)= (s(mn,2)+s(ms,2))*0.5
        dd(np,mp)=-(s(mn,1)-s(ms,1))*0.5
   72 continue
c
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,jtrunx2,+1)
c x-direction
      do 74 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
        wm1=(m-1)*odx
      do 74 j=1,my
        cc(mm,j)=-dd(j,mp)*wm1
        cc(mp,j)= dd(j,mm)*wm1
   74 continue
c
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,my,+1)
c
      do 76 j=1,my
      do 76 i=1,nx
        dxpl(i,j)=cc(i,j)
   76 continue
c da/dy
      do 60 i=1,(nx+3)
      do 60 j=1,my
	cc(i,j)=0.
   60 continue
      do 61 i=1,(my+3)
      do 61 j=1,nx
	dd(i,j)=0.
   61 continue
c
      do 62 n=1,jtruny
	nn=n*2-1
	np=nn+1
      do 62 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
        cc(mm,nn)= (s(mn,1)+s(ms,1))*0.5
        cc(mm,np)= (s(mn,2)-s(ms,2))*0.5
        cc(mp,nn)= (s(mn,2)+s(ms,2))*0.5
        cc(mp,np)=-(s(mn,1)-s(ms,1))*0.5
   62 continue
c
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,jtruny2,+1)
c
      do 64 n=1,jtruny
	nn=n*2-1
	np=nn+1
        wn1=(n-1)*ody
      do 64 i=1,nx
	dd(nn,i)=-cc(i,np)*wn1
	dd(np,i)= cc(i,nn)*wn1
   64 continue
c
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,nx,+1)
c
      do 66 j=1,my
      do 66 i=1,nx
        dypl(i,j)=dd(j,i)
   66 continue
c
      return
      end
      subroutine rtranuv1(jtrunx,jtruny,nx,my,lev,pi,dox,doy,cmn
     1,                  sv,sd,uu,vv)
c
      dimension sd(jtrunx*jtruny*2,2,lev)
     1,         sv(jtrunx*jtruny*2,2,lev)
      dimension uu(nx,lev,my),vv(nx,lev,my)
      dimension cmn(jtrunx*jtruny)
c
      include '../include/rparamt.h'
      include '../include/rfftcom.h'
      dimension cc(nxt+3,myt),dd(myt+3,nxt)
      dimension work(nxt*myt,2)
c
      odx=2.*pi/dox
      ody=2.*pi/doy
      jtrunxy=jtrunx*jtruny
      jtrunx2=jtrunx*2
      jtruny2=jtruny*2
c
      do 666 k=1,lev
c
      do 70 i=1,(nx+3)
      do 70 j=1,my
	cc(i,j)=0.
   70 continue
      do 71 i=1,(my+3)
      do 71 j=1,nx
	dd(i,j)=0.
   71 continue
c
      n=1
      m=1
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)=(sd(mn,1,k)+sd(ms,1,k))*0.5
	dd(np,mm)=0.
	dd(nn,mp)=0.
	dd(np,mp)=0.
c
        n=1
	nn=n*2-1
	np=nn+1
      do 72 m=2,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)=-(sd(mn,1,k)+sd(ms,1,k))*0.5/cmn(mn)
	dd(np,mm)=-(sd(mn,2,k)-sd(ms,2,k))*0.5/cmn(mn)
	dd(nn,mp)=-(sd(mn,2,k)+sd(ms,2,k))*0.5/cmn(mn)
	dd(np,mp)= (sd(mn,1,k)-sd(ms,1,k))*0.5/cmn(mn)
   72 continue
      do 73 n=2,jtruny
	nn=n*2-1
	np=nn+1
      do 73 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)=-(sd(mn,1,k)+sd(ms,1,k))*0.5/cmn(mn)
	dd(np,mm)=-(sd(mn,2,k)-sd(ms,2,k))*0.5/cmn(mn)
	dd(nn,mp)=-(sd(mn,2,k)+sd(ms,2,k))*0.5/cmn(mn)
	dd(np,mp)= (sd(mn,1,k)-sd(ms,1,k))*0.5/cmn(mn)
   73 continue
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,jtrunx2,+1)
c      do 74 m=1,jtrunx
      do 74 m=2,jtrunx
	mm=m*2-1
	mp=mm+1
        wm1=(m-1)*odx
      do 74 j=1,my
	cc(mm,j)=-dd(j,mp)*wm1
	cc(mp,j)= dd(j,mm)*wm1
   74 continue
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,my,+1)
      do 76 j=1,my
      do 76 i=1,nx
	uu(i,k,j)=cc(i,j)
   76 continue
c
      do 78 i=1,(nx+3)
      do 78 j=1,my
	cc(i,j)=0.
   78 continue
      do 79 i=1,(my+3)
      do 79 j=1,nx
	dd(i,j)=0.
   79 continue
c
      n=1
      m=1
	nn=n*2-1
	np=nn+1
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	cc(mm,nn)=(sv(mn,1,k)+sv(ms,1,k))*0.5
	cc(mm,np)=0.
	cc(mp,nn)=0.
	cc(mp,np)=0.
c
        n=1
	nn=n*2-1
	np=nn+1
      do 80 m=2,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	cc(mm,nn)=-(sv(mn,1,k)+sv(ms,1,k))*0.5/cmn(mn)
	cc(mm,np)=-(sv(mn,2,k)-sv(ms,2,k))*0.5/cmn(mn)
	cc(mp,nn)=-(sv(mn,2,k)+sv(ms,2,k))*0.5/cmn(mn)
	cc(mp,np)= (sv(mn,1,k)-sv(ms,1,k))*0.5/cmn(mn)
 80   continue
      do 82 n=2,jtruny
	nn=n*2-1
	np=nn+1
      do 82 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	cc(mm,nn)=-(sv(mn,1,k)+sv(ms,1,k))*0.5/cmn(mn)
	cc(mm,np)=-(sv(mn,2,k)-sv(ms,2,k))*0.5/cmn(mn)
	cc(mp,nn)=-(sv(mn,2,k)+sv(ms,2,k))*0.5/cmn(mn)
	cc(mp,np)= (sv(mn,1,k)-sv(ms,1,k))*0.5/cmn(mn)
 82   continue
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,jtruny2,+1)
c      do 84 n=1,jtruny
      do 84 n=2,jtruny
	nn=n*2-1
	np=nn+1
        wn1=(n-1)*ody
      do 84 i=1,nx
	dd(nn,i)=-cc(i,np)*wn1
	dd(np,i)= cc(i,nn)*wn1
   84 continue
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,nx,+1)
      do 86 j=1,my
      do 86 i=1,nx
	uu(i,k,j)=uu(i,k,j)-dd(j,i)
   86 continue
c--------------
      do 90 i=1,(nx+3)
      do 90 j=1,my
	cc(i,j)=0.
   90 continue
      do 91 i=1,(my+3)
      do 91 j=1,nx
	dd(i,j)=0.
   91 continue
c
      n=1
      m=1
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)=(sv(mn,1,k)+sv(ms,1,k))*0.5
	dd(np,mm)=0.
	dd(nn,mp)=0.
	dd(np,mp)=0.
c
      n=1
      nn=n*2-1
      np=nn+1
      do 92 m=2,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)=-(sv(mn,1,k)+sv(ms,1,k))*0.5/cmn(mn)
	dd(np,mm)=-(sv(mn,2,k)-sv(ms,2,k))*0.5/cmn(mn)
	dd(nn,mp)=-(sv(mn,2,k)+sv(ms,2,k))*0.5/cmn(mn)
	dd(np,mp)= (sv(mn,1,k)-sv(ms,1,k))*0.5/cmn(mn)
   92 continue
      do 93 n=2,jtruny
	nn=n*2-1
	np=nn+1
      do 93 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n -1)*jtrunx
	ms=mn+jtrunxy
	dd(nn,mm)=-(sv(mn,1,k)+sv(ms,1,k))*0.5/cmn(mn)
	dd(np,mm)=-(sv(mn,2,k)-sv(ms,2,k))*0.5/cmn(mn)
	dd(nn,mp)=-(sv(mn,2,k)+sv(ms,2,k))*0.5/cmn(mn)
	dd(np,mp)= (sv(mn,1,k)-sv(ms,1,k))*0.5/cmn(mn)
   93 continue
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,jtrunx2,+1)
c      do 94 m=1,jtrunx
      do 94 m=2,jtrunx
	mm=m*2-1
	mp=mm+1
        wm1=(m-1)*odx
      do 94 j=1,my
	cc(mm,j)=-dd(j,mp)*wm1
	cc(mp,j)= dd(j,mm)*wm1
   94 continue
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,my,+1)
      do 96 j=1,my
      do 96 i=1,nx
	vv(i,k,j)=cc(i,j)
   96 continue
c
      do 98 i=1,(nx+3)
      do 98 j=1,my
	cc(i,j)=0.
   98 continue
      do 99 i=1,(my+3)
      do 99 j=1,nx
	dd(i,j)=0.
   99 continue
c
      n=1
      m=1
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	cc(mm,nn)=(sd(mn,1,k)+sd(ms,1,k))*0.5
	cc(mm,np)=0.
	cc(mp,nn)=0.
	cc(mp,np)=0.
c
        n=1
	nn=n*2-1
	np=nn+1
      do 100 m=2,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n-1)*jtrunx
	ms=mn+jtrunxy
	cc(mm,nn)=-(sd(mn,1,k)+sd(ms,1,k))*0.5/cmn(mn)
	cc(mm,np)=-(sd(mn,2,k)-sd(ms,2,k))*0.5/cmn(mn)
	cc(mp,nn)=-(sd(mn,2,k)+sd(ms,2,k))*0.5/cmn(mn)
	cc(mp,np)= (sd(mn,1,k)-sd(ms,1,k))*0.5/cmn(mn)
 100  continue
      do 102 n=2,jtruny
	nn=n*2-1
	np=nn+1
      do 102 m=1,jtrunx
	mm=m*2-1
	mp=mm+1
	mn=m+(n -1)*jtrunx
	ms=mn+jtrunxy
	cc(mm,nn)=-(sd(mn,1,k)+sd(ms,1,k))*0.5/cmn(mn)
	cc(mm,np)=-(sd(mn,2,k)-sd(ms,2,k))*0.5/cmn(mn)
	cc(mp,nn)=-(sd(mn,2,k)+sd(ms,2,k))*0.5/cmn(mn)
	cc(mp,np)= (sd(mn,1,k)-sd(ms,1,k))*0.5/cmn(mn)
 102  continue
      call rfftmlt(cc,work,trigsx,ifax,1,nx+3,nx,jtruny2,+1)
c      do 104 n=1,jtruny
      do 104 n=2,jtruny
	nn=n*2-1
	np=nn+1
        wn1=(n-1)*ody
      do 104 i=1,nx
	dd(nn,i)=-cc(i,np)*wn1
	dd(np,i)= cc(i,nn)*wn1
  104 continue
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,nx,+1)
      do 106 j=1,my
      do 106 i=1,nx
	vv(i,k,j)=vv(i,k,j)+dd(j,i)
  106 continue
c
 666  continue
c
      return
      end
      subroutine rtrandv(jtrunx,jtruny,nx,my,lev,pi,dox,doy,ru,rv
     *,                  sv,sd)
      include '../include/rparamt.h'
      include '../include/rfftcom.h'
      dimension ru(nx,lev,my),sd(jtrunx*jtruny*2,2,lev)
      dimension rv(nx,lev,my),sv(jtrunx*jtruny*2,2,lev)
      dimension dd(myt+3,nxt),work(nxt*myt,2)
      dimension aa(nxt+3,myt),bb(nxt+3,myt),cc(myt+3,nxt)
c
      odx=2.*pi/dox
      ody=2.*pi/doy
      jtrunxy=jtrunx*jtruny
      jtrunx2=jtrunx*2
      jtruny2=jtruny*2
c
      do 168 k=1,lev
c
      do 58 i=1,nx+3
      do 58 j=1,my
	aa(i,j)=0.
	bb(i,j)=0.
   58 continue
      do 59 i=1,my+3
      do 59 j=1,nx
	dd(i,j)=0.
	cc(i,j)=0.
   59 continue
c
      do 60 j=1,my
      do 60 i=1,nx
	aa(i,j)=rv(i,k,j)
	bb(i,j)=ru(i,k,j)
   60 continue
c
      call rfftmlt(aa,work,trigsx,ifax,1,nx+3,nx,my,-1)
      call rfftmlt(bb,work,trigsx,ifax,1,nx+3,nx,my,-1)
c
      do 62 j=1,my
      do 62 m=1,jtrunx
         mm= 2*m-1
	 mp= mm+1
         wm1=m-1
	 dd(j,mm) =-aa(mp,j)*wm1*odx
	 dd(j,mp) = aa(mm,j)*wm1*odx
	 cc(j,mm) =-bb(mp,j)*wm1*odx
	 cc(j,mp) = bb(mm,j)*wm1*odx
   62 continue
c
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,jtrunx2,-1)
      call rfftmlt(cc,work,trigsy,ifay,1,my+3,my,jtrunx2,-1)
c
      do 66 n=1,jtruny
	nn= 2*n-1
	np= nn+1
      do 66 m=1,jtrunx
	mm= 2*m-1
	mp= mm+1
	mn= m+(n -1)*jtrunx
	ms= mn+jtrunxy
	sv(mn,1,k) = dd(nn,mm)-dd(np,mp)
	sv(mn,2,k) = dd(np,mm)+dd(nn,mp)
	sv(ms,1,k)=  dd(nn,mm)+dd(np,mp)
	sv(ms,2,k)= -dd(np,mm)+dd(nn,mp)
	sd(mn,1,k) = cc(nn,mm)-cc(np,mp)
	sd(mn,2,k) = cc(np,mm)+cc(nn,mp)
	sd(ms,1,k)=  cc(nn,mm)+cc(np,mp)
	sd(ms,2,k)= -cc(np,mm)+cc(nn,mp)
   66 continue
c
      do 8 i=1,(nx+3)
      do 8 j=1,my
        aa(i,j)=0.
        bb(i,j)=0.
   8  continue
      do 9 i=1,(my+3)
      do 9 j=1,nx
        dd(i,j)=0.
        cc(i,j)=0.
   9  continue
      do 67 j=1,my
      do 67 i=1,nx
	dd(j,i)=-ru(i,k,j)
	cc(j,i)= rv(i,k,j)
   67 continue
c
      call rfftmlt(dd,work,trigsy,ifay,1,my+3,my,nx,-1)
      call rfftmlt(cc,work,trigsy,ifay,1,my+3,my,nx,-1)
c
      do 68 i=1,nx
      do 68 n=1,jtruny
         nn= 2*n-1
	 np= nn+1
         wn1=n-1
	 aa(i,nn) =-dd(np,i)*wn1*ody
	 aa(i,np) = dd(nn,i)*wn1*ody
	 bb(i,nn) =-cc(np,i)*wn1*ody
	 bb(i,np) = cc(nn,i)*wn1*ody
   68 continue
c
      call rfftmlt(aa,work,trigsx,ifax,1,nx+3,nx,jtruny2,-1)
      call rfftmlt(bb,work,trigsx,ifax,1,nx+3,nx,jtruny2,-1)
c
      do 69 n=1,jtruny
	nn= 2*n-1
	np= nn+1
      do 69 m=1,jtrunx
	mm= 2*m-1
	mp= mm+1
	mn= m+(n -1)*jtrunx
	ms= mn+jtrunxy
	sv(mn,1,k) = sv(mn,1,k)+(aa(mm,nn)-aa(mp,np))
	sv(mn,2,k) = sv(mn,2,k)+(aa(mm,np)+aa(mp,nn))
	sv(ms,1,k) = sv(ms,1,k)+(aa(mm,nn)+aa(mp,np))
	sv(ms,2,k) = sv(ms,2,k)-(aa(mm,np)-aa(mp,nn))
	sd(mn,1,k) = sd(mn,1,k)+(bb(mm,nn)-bb(mp,np))
	sd(mn,2,k) = sd(mn,2,k)+(bb(mm,np)+bb(mp,nn))
	sd(ms,1,k) = sd(ms,1,k)+(bb(mm,nn)+bb(mp,np))
	sd(ms,2,k) = sd(ms,2,k)-(bb(mm,np)-bb(mp,nn))
 69   continue
c
  168 continue
c
      return
      end
      subroutine rtranuv(jtrunx,jtruny,nx,my,lev,pi,dox,doy,cmn
     1,                  sv,sd,uu,vv)
c
      dimension sd(jtrunx*jtruny*2,2,lev)
     1,         sv(jtrunx*jtruny*2,2,lev)
      dimension uu(nx,lev,my),vv(nx,lev,my)
      dimension cmn(jtrunx*jtruny)
c
      include '../include/rparamt.h'
      dimension wsd(jtrunxt*jtrunyt*2,2,levt)
     1,         wsv(jtrunxt*jtrunyt*2,2,levt)
      dimension dddx(nxt,myt),dddy(nxt,myt)
     1,         dvdx(nxt,myt),dvdy(nxt,myt)
c
      jtrunxy=jtrunx*jtruny
c
      do 666 k=1,lev
c
      wsd(1,1,k)=sd(1,1,k)
      wsd(1,2,k)=sd(1,2,k)
      wsd(jtrunxy+1,1,k)=sd(jtrunxy+1,1,k)
      wsd(jtrunxy+1,2,k)=sd(jtrunxy+1,2,k)
      do 50 mn= 2, jtrunxy
      ms=mn+jtrunxy
      wsd(mn,1,k)=-sd(mn,1,k)/cmn(mn)
      wsd(mn,2,k)=-sd(mn,2,k)/cmn(mn)
      wsd(ms,1,k)=-sd(ms,1,k)/cmn(mn)
      wsd(ms,2,k)=-sd(ms,2,k)/cmn(mn)
 50   continue
      wsv(1,1,k)=sv(1,1,k)
      wsv(1,2,k)=sv(1,2,k)
      wsv(jtrunxy+1,1,k)=sv(jtrunxy+1,1,k)
      wsv(jtrunxy+1,2,k)=sv(jtrunxy+1,2,k)
      do 51 mn= 2, jtrunxy
      ms=mn+jtrunxy
      wsv(mn,1,k)=-sv(mn,1,k)/cmn(mn)
      wsv(mn,2,k)=-sv(mn,2,k)/cmn(mn)
      wsv(ms,1,k)=-sv(ms,1,k)/cmn(mn)
      wsv(ms,2,k)=-sv(ms,2,k)/cmn(mn)
 51   continue
c
      call rtrngra(jtrunx,jtruny,nx,my,pi,dox,doy,wsd(1,1,k)
     *,            dddx,dddy)
      call rtrngra(jtrunx,jtruny,nx,my,pi,dox,doy,wsv(1,1,k)
     *,            dvdx,dvdy)
c
      do 60 j=1,my
      do 60 i=1,nx
        uu(i,k,j)=dddx(i,j)-dvdy(i,j)
        vv(i,k,j)=dddy(i,j)+dvdx(i,j)
 60   continue
c
 666  continue
c
      return
      end
      subroutine ralpha2s( nx,my,lev,jtrunx,jtruny,mlmax
     1                  , pi,dox,doy,vdmer,vdzon,divten)
c
      dimension vdmer(nx,lev,my),vdzon(nx,lev,my)
      dimension divten(mlmax,2,lev)
c
      include '../include/rparamt.h'
      include '../include/rfftcom.h'
      dimension work(nxt*myt,2)
      dimension gg(nxt+3,myt),hh(myt+3,nxt)
c
      odx=2.*pi/dox
      ody=2.*pi/doy
      jtrunx2=jtrunx*2
      jtruny2=jtruny*2
      jtrunxy=jtrunx*jtruny
c
      do 200 k=1,lev
c
      do 125 i=1,(nx+3)
      do 125 j=1,my
        gg(i,j)=0.
 125  continue
      do 126 i=1,(my+3)
      do 126 j=1,nx
        hh(i,j)=0.
 126  continue
c
      do 130 j=1,my
      do 130 i=1,nx
        gg(i,j)= vdmer(i,k,j)
 130  continue
c
      call rfftmlt(gg,work,trigsx,ifax,1,nx+3,nx,my,-1)
      do 140 j=1,my
      do 140 m=1,jtrunx
        mm=2*m-1
        mp=mm+1
        hh(j,mm)=-gg(mp,j)*(m-1)*odx
        hh(j,mp)= gg(mm,j)*(m-1)*odx
 140  continue
c
      call rfftmlt(hh,work,trigsy,ifay,1,my+3,my,jtrunx2,-1)
      do 150 n=1,jtruny
        nn=2*n-1
        np=nn+1
      do 150 m=1,jtrunx
        mm=2*m-1
        mp=mm+1
        mn= m+(n -1)*jtrunx
        ms= mn+jtrunxy
        divten(mn,1,k) = hh(nn,mm)-hh(np,mp)
        divten(mn,2,k) = hh(np,mm)+hh(nn,mp)
        divten(ms,1,k) = hh(nn,mm)+hh(np,mp)
        divten(ms,2,k) =-hh(np,mm)+hh(nn,mp)
 150  continue
c
      do 155 i=1,(nx+3)
      do 155 j=1,my
        gg(i,j)=0.
 155  continue
      do 156 i=1,(my+3)
      do 156 j=1,nx
        hh(i,j)=0.
 156  continue
c
      do 160 j=1,my
      do 160 i=1,nx
        hh(j,i)= vdzon(i,k,j)
 160  continue
c
      call rfftmlt(hh,work,trigsy,ifay,1,my+3,my,nx,-1)
      do 170 i=1,nx
      do 170 n=1,jtruny
        nn= 2*n-1
        np= nn+1
        gg(i,nn) =-hh(np,i)*(n-1)*ody
        gg(i,np) = hh(nn,i)*(n-1)*ody
 170  continue
c
      call rfftmlt(gg,work,trigsx,ifax,1,nx+3,nx,jtruny2,-1)
c
      do 180 n=1,jtruny
        nn= 2*n-1
        np= nn+1
      do 180 m=1,jtrunx
        mm= 2*m-1
        mp= mm+1
        mn= m+(n -1)*jtrunx
        ms= mn+jtrunxy
        divten(mn,1,k) = divten(mn,1,k)+(gg(mm,nn)-gg(mp,np))
        divten(mn,2,k) = divten(mn,2,k)+(gg(mm,np)+gg(mp,nn))
        divten(ms,1,k) = divten(ms,1,k)+(gg(mm,nn)+gg(mp,np))
        divten(ms,2,k) = divten(ms,2,k)-(gg(mm,np)-gg(mp,nn))
 180  continue
c
 200  continue
c
      return
      end
