      subroutine pl2d (vl,vr,vb,vt,rp,m,n,rd,ncl,nd,
     *                   ishad,ity,dummy)
c
c  input:
c        vl,vr,vb,vt : frame domain
c        rp(m,n) : input data array
c	 rd(ncl) : contour line distribution
c        nd      : contour level use index 
c                  ex. 0  - no contour line to be drawn
c                      +n - contour line to be drawn, and
c                           a labeled line for "n" contour lines
c                      1  - each contour line labeled
c                      2  - a unlabeled line next to a labeled line
c                      -n - as "+n", but dash line to be drawn for
c                           "rd(i)" negative 
c        ishad(ncl+1)  : shading density (color) index array
c                        ishad = 0 -- no shading
c                        ishad > 0 -- solid line or area fill 
c                        ishad < 0 -- dot fill 
c                       * abs(ishad) increasing ,then density increasing
c                       * abs(ishad) =< 15 for ity not equal zero
c        ity     : shading type index (-2 to +2)
c                  index= 0 -- color full fill
c                         1 -- color parallel line fill or dot
c                         2 -- color cross line fill or dot
c                        -1 -- parallel line fill or dot
c                        -2 -- cross line fill or dot
c                      
c   provider : Chin-Tzu Fong        
c
      parameter (nz=40)   ! ncl must be not greater than nz
      dimension rp(m,n)
      parameter (nama=50000000,ncra=100000)
      dimension iama(nama),rwrk(95000),iwrk(95000),xcra(ncra),
     +          ycra(ncra),iara(nz),igra(nz)
      dimension rd(ncl),ishad(ncl+1),iasf(13)
      data iasf/13*1/
      external shader
      external drawcl
      external filland
c
      call gslwsc(1.0)
c
      if(nz .lt. ncl) then
      print*,' fatal error in "plot2d", ncl too large (must be not'
      print*,' greater than nz(=40)'
      end if
      call gsasf(iasf)
c
      if(ity.lt.0)then
      call gsplci(1)
      call gstxci(1)
      call gsfais(0)
      call cpseti('LBC',0)
      end if
c
      call param (ishad,ncl,ity)
c
      call cpsetr('VPL',vl)
      call cpsetr('VPR',vr)
      call cpsetr('VPB',vb)
      call cpsetr('VPT',vt)
      call cpsetr('VPS',0.)
      call cpseti('CLS',0)
      call cpseti('CFB',1)
      call cpseti('NCL',ncl)
      call arinam(iama,nama)
      call cprect(rp,m,m,n,rwrk,95000,iwrk,95000)
      call cppkcl(rp,rwrk,iwrk)
      call cpsetc('ILT',' ')
      call cpsetr('HLL',2.)
      call cpsetc('HLT','$ZDV$')
c     call cpsetc('HLT','H''L')
      if(ity.eq.0)call cpsetc('HLT',' ')
      call cpseti('HLB',3)
c     call cpseti('HLB',0)
      call CPSETR('SPV - SPECIAL VALUE',dummy)   
      call cpseti('HLO',7)
      call cpseti('GIL',5)
c     call cpseti('NSD',-2)
      call cpseti('NLZ',1)
c
      do 101 i=1,ncl
c      call cpseti('LLP',3)
      call cpseti('LLP',0)
      call cpseti('LLB',0)
      call cpseti('RWC',200)
      call cpseti('LLO',0)
      call cpsetr('CWM',1.0)
      call cpsetr('LLS',0.015)
      call cpseti('PAI',i)
      call cpsetr('CLV',rd(i))
      call cpsetr('CLL',2.)
      if(rd(i).eq.0.)call cpsetr('CLL',4.)
      call cpseti('CLU',1)
c
      if(nd.eq.0)then
      call cpseti('CLU',0)
      go to 50
      end if
      if(mod(i-1,nd).eq.0)call cpseti('CLU',3)
      if(nd.lt.0.and.rd(i).lt.0.)
     * call cpsetc('CLD','$$''''$$''''')
c
 50   call cpseti('AIA',i+1)
      call cpseti('AIB',i)
c
  101 continue
c
      ifilland= 0           ! swich of the land area colorred
c
      call cplbam(rp,rwrk,iwrk,iama)
      call cplbdr(rp,rwrk,iwrk)
      call cpclam(rp,rwrk,iwrk,iama)
      call arscam(iama,xcra,ycra,ncra,iara,igra,ncl,shader)
      call cpcldm(rp,rwrk,iwrk,iama,drawcl)
c
      if(ifilland.eq.1)then
         call mapstc('OU','PO')
         call mapint
         call mapbla(iama)
         call cplbam(rp,rwrk,iwrk,iama)
         call arscam(iama,xcra,ycra,ncra,iara,igra,ncl,filland)
      endif
 999  return
      end



      subroutine drawcl(xcs,ycs,ncs,iai,iag,nai)
      dimension xcs(*),ycs(*),iai(*),iag(*)
      idr=1
      do 101 i=1,nai
      if(iai(i).lt.0)idr=0
  101 continue
      if(idr.ne.0)call curved(xcs,ycs,ncs)
      return
      end 

c
      subroutine shader(xcs,ycs,ncs,iai,iag,nai)
      parameter (nz=40)
      dimension xcs(*),ycs(*),iai(*),iag(*)
      dimension dst(100000),ind(100000)
      dimension ishad(nz+1)
      common/para/ishad,ity,ishmin,ishmax
      integer idp(8,8)
      data idp/1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1/
c
      ish=0
      do 101 i=1,nai
      if(iag(i).ne.3.and.iai(i).ge.1)print*,'iag,iai=',iag(i),iai(i)
      if(iai(i).ge.1)ish=iai(i)
  101 continue
c
      if(ish.gt.0)then
      do 102 i=1,nai
      if(iag(i).eq.5.and.iai(i).eq.-1)ish=0
 102  continue 
      end if
c
      if(ish.gt.0.and.ishad(ish).ne.0)then
      ishmin=min0(ish,ishmin)
      ishmax=max0(ish,ishmax)
        if(ity.lt.-2 .or. ity.gt.2) ity=2
c
c-- color full fill ----
      if(ity.eq.0) then
      call sfseti('TYPE OF FILL',ity)
      call sfsgfa(xcs,ycs,ncs,dst,100000,ind,100000
     *            ,abs(ishad(ish)))
c
c-- color line fill
      else if(ity.gt.0) then 
      call sfseti('TYPE OF FILL',ity)
      call sfseti('AN',45)
      call sfseti('CH',0)
      call sfseti('DO',0)
      if(ishad(ish).lt.0) call sfseti('DO',1)
      dsp=(0.004-0.0010)/15.
c     rsp=0.0045-float(ish)*dsp*3.5
      rsp=0.012-abs(float(ishad(ish))-1.)*dsp*3.5
      call sfsetp(idp)
      call sfsetr('SP',rsp)
      call sfsgfa(xcs,ycs,ncs,dst,100000,ind,100000
     *            ,abs(ishad(ish)))
c-- line fill
      else 
      call sfseti('TYPE OF FILL',ity)
      call sfseti('AN',45)
      call sfseti('CH',0)
      call sfseti('DO',0)
      if(ishad(ish).lt.0) call sfseti('DO',1)
      dsp=(0.004-0.0010)/15.
      rsp=0.012-abs(float(ishad(ish))-1.)*dsp*3.5
      call sfsetp(idp)
      call sfsetr('SP',rsp)
      call sfwrld(xcs,ycs,ncs-1,dst,100000,ind,100000)
      call sfnorm(xcs,ycs,ncs-1,dst,100000,ind,100000)
      end if
c
      if(ity.eq.-2.and.ishad(ish).gt.0)then
      call sfseti('AN',45+90)
      call sfseti('CH',0)
      call sfseti('DO',0)
c     if(ishad(ish).lt.0) call sfseti('DO',1)
      dsp=(0.004-0.0010)/15.
      rsp=0.012-abs(float(ishad(ish))-1.)*dsp*3.5
      call sfsetp(idp)
      call sfsetr('SP',rsp)
      call sfwrld(xcs,ycs,ncs-1,dst,100000,ind,100000)
      call sfnorm(xcs,ycs,ncs-1,dst,100000,ind,100000)
      end if
      endif
      return
      end
c
      subroutine param (ishad,ncl,ity)
c    
c  set up some parameters shared by "shader" and "plotbar"
c
      parameter (nz=40)
      dimension ishad(ncl)
      common /para/ishd(nz+1),iz,ishmin,ishmax
      data ishmin/99/,ishmax/0/ 
      iz=ity
      do i=1,nz+1
      ishd(i)=0
      end do
      do i=1,ncl+1
      ishd(i)=ishad(i)
      end do
      return 
      end
 
 
      subroutine plotbar (llb,rbar,barw,barl)
c
c     llb(ncl+2) : a character array providing a list of labels 
c                  for the bar
c     rbar         : specifying the position of the bar ( value must be
c                    great than or equal 1.0 )
c                    if 1.0 <= rbar < 2.0 -- below the plot     
c                    if 2.0 <= rbar < 3.0 -- above the plot
c                    if 3.0 <= rbar < 4.0 -- the right of the plot
c                    if 4.0 <= rbar < 5.0 -- the left of the plot
c  barw : the width of the bar
c  barl : the extending length of the bar
c      barw=0.3   ! you can change it (0.2 to 0.7)
c      barl=1.0   ! you can change it (0.4 to 1.0)
c     
c     provider : Chin-Tzu Fong
c                         
      parameter (nz=40)
      dimension ishad(nz+1)
      character*(*) llb(nz+2)
      common/para/ishad,ity,ishmin,ishmax
c
c
      ic=ishmax-ishmin+1
c
      if(ity .eq. 0)then
      do 10 i=1,nz+1
      ishad(i)=abs(ishad(i))
  10  continue
      end if
c
      call gslwsc(2.0)
c
c     call cpgetr('VPL',vl)
c     call cpgetr('VPR',vr)
c     call cpgetr('VPB',vb)
c     call cpgetr('VPT',vt)
      call getset(vl,vr,vb,vt,ul,ur,ub,ut,ll)
c
      if(rbar.ge.1.0.and.rbar.lt.2.0)then
      if(vb.lt.0.1)go to 99
      dx=(vb-0.02)-0.08
      vtz=(vb-0.02)-dx*(rbar-1.0)
      vbz=vtz-0.08
      vlz=vl+(1.-barl)*(vr-vl)
c
      call lblbar(0,vlz,vr,vbz,vtz,ic,1.0,barw,ishad(ishmin),
     *            ity,llb(ishmin),ic+1,1)
c   
      else if(rbar.ge.2.0.and.rbar.lt.3.0)then
      if(vt.gt.0.9)go to 99
      dx=0.92-(vt+0.02)
      vbz=(vt+0.02)+dx*(rbar-2.0) 
      vtz=vbz+0.08
      vlz=vl+(1.-barl)*(vr-vl)
c
      call lblbar(0,vlz,vr,vbz,vtz,ic,1.0,barw,ishad(ishmin),
     *            ity,llb(ishmin),ic+1,2)
c
      else if(rbar.ge.3.0.and.rbar.lt.4.0)then
      if(vr.gt.0.9)go to 99
      dx=0.92-(vr+0.02)
      vlz=(vr+0.02)+dx*(rbar-3.0)
      vrz=vlz+0.08
      vbz=vb+(barl)*(vt-vb)
c
      call lblbar(1,vlz,vrz,vb,vbz,ic,barw,1.0,ishad(ishmin),
     *            ity,llb(ishmin),ic+1,1)
c  
      else if(rbar.ge.4.0.and.rbar.lt.5.0)then
      if(vl.lt.0.1)go to 99
      dx=(vl-0.02)-0.08
      vrz=(vl-0.02)-dx*(rbar-4.0)
      vlz=vrz-0.08
      vbz=vb+(barl)*(vt-vb)
c
      call lblbar(1,vlz,vrz,vb,vbz,ic,barw,1.0,ishad(ishmin),
     *            ity,llb(ishmin),ic+1,2)
c  
      end if      
c
      ishmin=99   !reinitialize "ishmin" and "ishmax"
      ishmax=0
      go to 100
c    
  99  print*,'There is not enough space to include labelbar,'
      print*,'change the position of it.'
c
 100  continue
      return
      end      
c
      subroutine lbfill(ity,xcs,ycs,ncs,indx)
      dimension xcs(*),ycs(*)
      dimension dst(100000),ind(100000)
      integer idp(8,8)
      data idp/1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1,
     +         1,1,1,1,1,1,1,1/
c
c-- color line fill
      if(ity.gt.0) then 
      call sfseti('TYPE OF FILL',ity)
      call sfseti('AN',45)
      call sfseti('CH',0)
      call sfseti('DO',0)
      if(indx.lt.0) call sfseti('DO',1)
      dsp=(0.004-0.0010)/15.
c     rsp=0.0045-float(ish)*dsp*3.5
      rsp=0.012-abs(float(indx)-1.)*dsp*3.5
      call sfsetp(idp)
      call sfsetr('SP',rsp)
      call sfsgfa(xcs,ycs,ncs,dst,100000,ind,100000
     *            ,abs(indx))
c-- line fill
      else 
      call sfseti('TYPE OF FILL',ity)
      call sfseti('AN',45)
      call sfseti('CH',0)
      call sfseti('DO',0)
      if(indx.lt.0) call sfseti('DO',1)
      dsp=(0.004-0.0010)/15.
      rsp=0.012-abs(float(indx)-1.)*dsp*3.5
      call sfsetp(idp)
      call sfsetr('SP',rsp)
      call sfwrld(xcs,ycs,ncs-1,dst,100000,ind,100000)
      call sfnorm(xcs,ycs,ncs-1,dst,100000,ind,100000)
      end if
c
      if(ity.eq.-2.and.indx.gt.0)then
      call sfseti('AN',45+90)
      call sfseti('CH',0)
      call sfseti('DO',0)
c     if(indx.lt.0) call sfseti('DO',1)
      dsp=(0.004-0.0010)/15.
      rsp=0.012-abs(float(indx)-1.)*dsp*3.5
      call sfsetp(idp)
      call sfsetr('SP',rsp)
      call sfwrld(xcs,ycs,ncs-1,dst,100000,ind,100000)
      call sfnorm(xcs,ycs,ncs-1,dst,100000,ind,100000)
      end if
      return
      end


      subroutine bckclr(iline,itext,ibck)
      dimension xd(5),yd(5)
      data xd/0.,1.,1.,0.,0./
      data yd/0.,0.,1.,1.,0./
      call gsplci(iline)      ! line color
      call gstxci(itext)      ! text color 
      call gsfais(1)          ! set fill area style to be solid fill
      call gsfaci(ibck)        ! set fill area color index
      call gfa(5,xd,yd)       ! filled area
      call cpseti('LBC',ibck)    ! label box color index 
      return
      end 


      subroutine filland(xcs,ycs,ncs,iai,iag,nai)
      dimension xcs(*),ycs(*),iai(*),iag(*)
c
      iai1=-1
      iai5=-1
c
      do 101 i=1,nai
      if(iag(i).eq.1) iai1=iai(i)
      if(iag(i).eq.5) iai5=iai(i)
 101  continue
c
      if(iai1.gt.0.and.mapaci(iai1).ne.1)then   ! land area
c        if(iai5.ge.0)then      ! exclude the area overlapped with
         call gsfais(1)         ! label boxes     
         call gsfaci(0)
         call gfa(ncs-1,xcs,ycs)
c        endif
      endif
c
      return
      end 
      
      subroutine mapch1(vl,vr,vb,vt,ndot)
      call mappos(vl,vr,vb,vt)
      call mapset('MA',0.,0.,0.,0.)
      call mapsti('LA',0)
      call mapsti('DO',ndot)
      call mapsti('EL',0)
      call mapsti('PE',1)
      call mapstc('OU','CO')
      call mapsti('GR',0)
      call maproj('CE',0.,-180.,0.)
      call gslwsc(1.0)
      call mapdrw
      call cpseti('SET',0)  
      call cpseti('MAP',1)
      call cpsetr('XC1',0.)
      call cpsetr('XCM',360.)
      call cpsetr('YC1',-90.)
      call cpsetr('YCN',90.)
      return
      end


      subroutine maplbm1(key)
      dimension xlb(6),ylb(7)
      character*(*) key
      character xlb*4,ylb*3
      data xlb/' 30E',' 90E','150E','150W',' 90W',' 30W'/
      data ylb/'90S','60S','30S',' 0 ','30N','60N','90N'/
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      slw=6.25*(fr-fl)
      call gslwsc(slw)
      call perim(12,1,6,1)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      slw=3.75*(fr-fl)
      call gslwsc(slw)
      wlg1=0.01625*(fr-fl)
      wlg2=0.01875*(fr-fl)
      do i=1,6
      dx=(fr-fl)/12 
      x =(fl+dx)+(i-1)*2*dx  
      fbz=fb-0.025*(fr-fl) 
      call plchhq(x,fbz,xlb(i),wlg1,0.,0.)
      end do
      do i=1,7
      dy=(ft-fb)/6
      y =fb+(i-1)*dy
      flz=fl-0.0375*(fr-fl) 
      call plchhq(flz,y,ylb(i),wlg1,0.,0.)
      end do
      x2=fl+(fr-fl)/2.
      ftz=ft+0.05
      call plchhq(x2,ftz,key,wlg2,0.,0.)
      call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      return 
      end

      subroutine mapch2_n(vl,vr,vb,vt,ndot)
      call mappos(vl,vr,vb,vt)
      call mapsti('DO',ndot)
      call mapsti('LA',0)
      call mapsti('EL',1)
      call mapsti('PE',0)
      call mapstc('OU','CO')
      call mapsti('GR',30)
      call maproj('ST',90.,0.,120.)
c     call mapset('AN',90.,90.,90.,90.)
      call mapset('AN',92.5,92.5,92.5,92.5)
      call gslwsc(1.0)
      call mapdrw
      call cpseti('SET',0)
      call cpseti('MAP',1)
      call cpsetr('XC1',0.)
      call cpsetr('XCM',360.)
c     call cpsetr('YC1',0.)
      call cpsetr('YC1',-2.5)
      call cpsetr('YCN',90.)
      return
      end


      subroutine mapch2_s(vl,vr,vb,vt,ndot)
      call mappos(vl,vr,vb,vt)
      call mapsti('DO',ndot)
      call mapsti('LA',0)
      call mapsti('EL',1)
      call mapsti('PE',0)
      call mapstc('OU','CO')
      call mapsti('GR',30)
      call maproj('ST',-90.,0.,60.)
      call mapset('AN',92.5,92.5,92.5,92.5)
      call gslwsc(1.0)
      call mapdrw
      call cpseti('SET',0)
      call cpseti('MAP',1)
      call cpsetr('XC1',0.)
      call cpsetr('XCM',360.)
      call cpsetr('YC1',-90.)
      call cpsetr('YCN',2.5)
      return
      end


      subroutine maplbm2_n(key)
      dimension sx(12),xc(361),yc(361)
      character*(*) key
      character sx*4
      data sx/'120E','150E','180 ','150W','120W',' 90W',' 60W',' 30W',
     *        '  0 ',' 30E',' 60E',' 90E'/
      pi=4.0*atan(1.0)
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      r=amin1((fr-fl),(ft-fb))/2.
      opx=fl+(fr-fl)/2.
      opy=fb+(ft-fb)/2. 
      do i=1,361
      theta=(i-1)*1.
      xc(i)=opx+r*cos(theta*pi/180.)
      yc(i)=opy+r*sin(theta*pi/180.)
      end do
      slw=6.25*(fr-fl)
      call gslwsc(slw)
      call curved(xc,yc,361)
      slw=3.75*(fr-fl)
      call gslwsc(slw)
      wlg1=0.01625*(fr-fl)
      wlg2=0.01875*(fr-fl)
      do i=1,12
      theta=270.+(i-1)*30.
      if(theta.ge.360.)theta=theta-360.
      wcos=cos(theta*pi/180.)
      wsin=sin(theta*pi/180.)
      wf=(fr-fl)
      rr= (r+0.0425*wf)-abs(wsin)*0.0125*wf
      x= opx+rr*wcos
      y= opy+rr*wsin
      call plchhq(x,y,sx(i),wlg1,0.,0.)
      end do
      x2=fl+(fr-fl)/2.
      ftz=ft+0.05
      call plchhq(x2,ftz,key,wlg2,0.,0.)
      call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      return
      end


      subroutine maplbm2_s(key)
      dimension sx(12),xc(361),yc(361)
      character*(*) key
      character sx*4
      data sx/'120E',' 90E',' 60E','30E ',' 0  ',' 30W',' 60W',' 90W',
     *        '120W','150W','180 ','150E'/
      pi=4.0*atan(1.0)
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      r=amin1((fr-fl),(ft-fb))/2.
      opx=fl+(fr-fl)/2.
      opy=fb+(ft-fb)/2. 
      do i=1,361
      theta=(i-1)*1.
      xc(i)=opx+r*cos(theta*pi/180.)
      yc(i)=opy+r*sin(theta*pi/180.)
      end do
      slw=6.25*(fr-fl)
      call gslwsc(slw)
      call curved(xc,yc,361)
      slw=3.75*(fr-fl)
      call gslwsc(slw)
      wlg1=0.01625*(fr-fl)
      wlg2=0.01875*(fr-fl)
      do i=1,12
      theta=270.+(i-1)*30.
      if(theta.ge.360.)theta=theta-360.
      wcos=cos(theta*pi/180.)
      wsin=sin(theta*pi/180.)
      wf=(fr-fl)
      rr= (r+0.0425*wf)-abs(wsin)*0.0125*wf
      x= opx+rr*wcos
      y= opy+rr*wsin
      call plchhq(x,y,sx(i),wlg1,0.,0.)
      end do
      x2=fl+(fr-fl)/2.
      ftz=ft+0.05
      call plchhq(x2,ftz,key,wlg2,0.,0.)
      call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      return
      end

      subroutine dfclrs_clr
      dimension rgbv(3,15)
      dimension rgbv1(3,300)
c 15 basic color defined as "NCAR GRAPHICS GUIDE TO NEW UTILITYS" 
      data rgbv /1.00,1.00,1.00, !1 white 
     +           0.70,0.70,0.70, !2
     +           0.75,0.50,1.00, !3
     +           0.50,0.00,1.00, !4
     +           0.00,0.00,1.00, !5
     +           0.00,0.50,1.00, !6
     +           0.00,1.00,1.00, !7
     +           0.00,1.00,0.60, !8
     +           0.00,1.00,0.00, !9
     +           0.70,1.00,0.00, !10
     +           1.00,1.00,0.00, !11 
     +           1.00,0.75,0.00, !12
     +           1.00,0.38,0.38, !13
     +           1.00,0.00,0.38, !14
     +           1.00,0.00,0.00/ !15 
c
      call gscr(1,0,0.,0.,0.)
      do  i=1,15
      call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
      enddo
c
cseri-1  color(16-20)
c
      do k=16,25
      kc=k-15
      rgbv1(1,k)=0.1+0.04*kc
      rgbv1(2,k)= 0.
      rgbv1(3,k)=rgbv1(1,k)*2 
      enddo
c
cseri-2  color(21-25)
      do k=26,35
      kc=k-25
      rgbv1(1,k)=0.5 
      rgbv1(2,k)=0.1*kc 
      rgbv1(3,k)=1.0 
      enddo
c
cseri-3  color(26-30)
      do k=36,45
      kc=k-35
      rgbv1(1,k)=0.5+0.05*kc 
      rgbv1(2,k)=1. 
      rgbv1(3,k)=1.0
      enddo
c
cseri-4  color(31-35)
      do k=46,55
      kc=k-45
      rgbv1(1,k)=1.0
      rgbv1(2,k)=1.
      rgbv1(3,k)=1.-0.1*kc 
      enddo
c
cseri-5  color(36-40)
      do k=56,65
      kc=k-55
      rgbv1(1,k)=1.0
      rgbv1(2,k)=1.- 0.1*kc
      rgbv1(3,k)=0. 
      enddo
c
cseri-6  color(41-45) 
      do k=66,75
      kc=k-65
      rgbv1(1,k)=1.0 
      rgbv1(2,k)=1.0-0.1*kc 
      rgbv1(3,k)=1.0-0.1*kc
      enddo
c
cseri-13  color(76-80)
c
      do k=76,90
      kc=k-75
      rgbv1(1,k)=1.-0.05*kc
      rgbv1(2,k)=1. 
      rgbv1(3,k)=1.-0.05*kc
      enddo
c
cseri-14  color(81-85)
c
      do k=91,100
      kc=k-90
      rgbv1(1,k)=0.25-0.011*kc
      rgbv1(2,k)= rgbv1(1,k)*4.
      rgbv1(3,k)=rgbv1(1,k)
      enddo
c
      do 100 i=16,100
      call gscr(1,i,rgbv1(1,i),rgbv1(2,i),rgbv1(3,i))
  100 continue
c
      do  ii=101,115
	i=ii-100
      zzz=float(i-1)/14
      rgbv(1,ii)=1.0-zzz*.8  
      rgbv(2,ii)=1.0-zzz*.8  
      rgbv(3,ii)=1.0-zzz*.8  
      call gscr(1,ii,rgbv(1,ii),rgbv(2,ii),rgbv(3,ii))
      enddo

      do k=116,126
      kc=k-115
      rgbv1(1,k)=.08*kc               
      rgbv1(2,k)= .04*kc                
      rgbv1(3,k)=0.                
      enddo
c
      do i=116,126
      call gscr(1,i,rgbv1(1,i),rgbv1(2,i),rgbv1(3,i))
      enddo
c
      return
      end
	subroutine dfclrs
      dimension rgbv2(3,25)
      dimension rgbv(3,15)
      dimension rgbv1(3,300)

      dimension r(100),r_1(100)
c 15 basic color defined as "NCAR GRAPHICS GUIDE TO NEW UTILITYS"
      data rgbv /1.00,1.00,1.00, !1 white
     +           0.70,0.70,0.70, !2
     +           0.75,0.50,1.00, !3
     +           0.50,0.00,1.00, !4
     +           0.00,0.00,1.00, !5
     +           0.00,0.50,1.00, !6
     +           0.00,1.00,1.00, !7
     +           0.00,1.00,0.60, !8
     +           0.00,1.00,0.00, !9
     +           0.70,1.00,0.00, !10
     +           1.00,1.00,0.00, !11
     +           1.00,0.75,0.00, !12
     +           1.00,0.38,0.38, !13
     +           1.00,0.00,0.38, !14
     +           1.00,0.00,0.00/ !15
c
	do i=1,225
        call gscr(1,i,1.,1.,1.)
	enddo
      call gscr(1,0,0.,0.,0.)

      do  i=1,15
      call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
      enddo
	IWKID=1
      CALL GSCR(IWKID,16, 0.0, 0.9, 1.0)
      CALL GSCR(IWKID,17, 0.9, 0.25, 0.0)
      CALL GSCR(IWKID,18, 1.0, 0.0, 0.2)
      CALL GSCR(IWKID,19, 1.0, 0.65, 0.0)
      CALL GSCR(IWKID,20, 1.0, 1.0, 0.0)
      CALL GSCR(IWKID,21, 0.7, 1.0, 0.2)
      CALL GSCR(IWKID,22, 0.5, 1.0, 0.0)
      CALL GSCR(IWKID,23, 0.2, 1.0, 0.5)
      CALL GSCR(IWKID,24, 0.2, 0.8, 0.2)
      CALL GSCR(IWKID,25, 0.0, 0.75, 1.0)
      CALL GSCR(IWKID,26, 0.25, 0.45, 0.95)
      CALL GSCR(IWKID,27, 0.4, 0.35, 0.8)
      CALL GSCR(IWKID,28, 0.6, 0.0, 0.8)
      CALL GSCR(IWKID,29, 0.85, 0.45, 0.8)
      CALL GSCR(IWKID,30, 0.8, 0.8, 1.0)

      do  i=1,25
      zzz=float(i-1)/24
      rgbv2(1,i)=1. -zzz*.8
      rgbv2(2,i)=1. -zzz*.8
      rgbv2(3,i)=1. -zzz*.8
      call gscr(1,i+30,rgbv2(1,i),rgbv2(2,i),rgbv2(3,i))
      enddo
      call gscr(1,55,0.,0.,0.)                                 
		
	dr=31./255.
	r_0=252./255.	! r_0 is the begining of the color table
	n=int(r_0/dr)
	do i=1,n
	r(i)=(r_0-(i-1)*dr)
	r_1(i)=1.-r(i)
	enddo
	nnn=55

	r0=255./255.
	r1=0./255.
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r0,r_1(i),r1)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r(i),r0,r1)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r1,r0,r_1(i))    
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r1,r(i),r0)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r_1(i),r1,r0)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r0,r1,r(i))
	enddo
	nnn=nnn+1
      CALL GSCR(IWKID,nnn, 1.,1.,1.)       

c
c
c

cseri-1  color(16-20)
c
      do k=16,25,2
      kc=k-15
      rgbv1(1,k)=0.1+0.04*kc
      rgbv1(2,k)= 0.
      rgbv1(3,k)=rgbv1(1,k)*2
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo
c
cseri-2  color(21-25)
      do k=26,35,2
      kc=k-25
      rgbv1(1,k)=0.5
      rgbv1(2,k)=0.1*kc
      rgbv1(3,k)=1.0
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo
c
cseri-3  color(26-30)
      do k=36,45,2
      kc=k-35
      rgbv1(1,k)=0.5+0.05*kc
      rgbv1(2,k)=1.
      rgbv1(3,k)=1.0
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
	enddo
c
cseri-4  color(31-35)
      do k=46,55,2
      kc=k-45
      rgbv1(1,k)=1.0
      rgbv1(2,k)=1.
      rgbv1(3,k)=1.-0.1*kc
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo
c
cseri-5  color(36-40)
      do k=56,65,2
      kc=k-55
      rgbv1(1,k)=1.0
      rgbv1(2,k)=1.- 0.1*kc
      rgbv1(3,k)=0.
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo
c
cseri-6  color(41-45)
      do k=66,75,2
      kc=k-65
      rgbv1(1,k)=1.0
      rgbv1(2,k)=1.0-0.1*kc
      rgbv1(3,k)=1.0-0.1*kc
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
	enddo
c
cseri-13  color(76-80)
c
      do k=76,90,2
      kc=k-75
      rgbv1(1,k)=1.-0.05*kc
      rgbv1(2,k)=1.
      rgbv1(3,k)=1.-0.05*kc
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo
c
cseri-14  color(81-85)
c
      do k=91,100,2
      kc=k-90
      rgbv1(1,k)=0.25-0.011*kc
      rgbv1(2,k)= rgbv1(1,k)*4.
      rgbv1(3,k)=rgbv1(1,k)
	nnn=nnn+1
        CALL GSCR(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo


	nnn=nnn+1
      CALL GSCR(IWKID,nnn, 1.,1.,1.)       
	goto 91
	r0=215./255.
	r1=40./255.

	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r0,r_1(i),r1)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r(i),r0,r1)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r1,r0,r_1(i))    
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r1,r(i),r0)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r_1(i),r1,r0)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r0,r1,r(i))
	enddo
	nnn=nnn+1
      CALL GSCR(IWKID,nnn, 1.,1.,1.)       
91	continue
	r0=175./255.
	r1=80./255.

	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r0,r_1(i),r1)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r(i),r0,r1)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r1,r0,r_1(i))    
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r1,r(i),r0)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r_1(i),r1,r0)
	enddo
	do i=1,n
	nnn=nnn+1
        CALL GSCR(IWKID,nnn,r0,r1,r(i))
	enddo
	nnn=nnn+1
      CALL GSCR(IWKID,nnn, 1.,1.,1.)       
      CALL GSCR(IWKID,nnn, 1.,1.,1.)       
c
      do k=116,126
      kc=k-115
      rgbv1(1,k)=.08*kc
      rgbv1(2,k)= .04*kc
      rgbv1(3,k)=0.
	nnn=nnn+1
      call gscr(1,nnn,rgbv1(1,k),rgbv1(2,k),rgbv1(3,k))
      enddo
c

      return
      end
