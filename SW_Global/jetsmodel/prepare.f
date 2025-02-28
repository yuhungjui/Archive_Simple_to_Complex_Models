      subroutine prepare
c----------------------------------------------------------------
c  define constant and compute something need in the model and
c  read parameters which control the model run from the namelist
c
c  constants which define in this subroutine are following :
c
c  rad   : radius of earth
c  omega : angular velocity of earth
c  cp    :
c  rgas  : gas constant
c  capa  :
c  pi    :
c  radsq : rad**2
c  omega : angular velocity of earth
c  hmean : mean depth of model ( m )
c  taui  : initial time integration( hour )
c  taue  : end of model time integration( hour )
c  tauo  : time intervel of model output( hour )
c  dt    : time step of model intgration ( second )
c  linear : logical,whether model is linear or not.
c  baro   : logical,whether model is baro-tropic model or not
c  nx     : gaussian grid point number on e-w direction
c  my     : gaussian grid point number on s-n direction
c  jtrun  : triangular truncation wave number
c  mlmax  : number of spectrum coefficients
c  mlsort : index array of wave (  or p(m,n) where, m is horizontal
c           wave number and n is total wave number )
c  eps4   : array of n(n+1)/radsq
c  cim    : array of m
c  weight : array of gaussian weight
c  sinl   : sin of gaussian latitude
c  onocos : 1./cos(gsuaaian latitude)**2
c  cosl   : cos of gaussian latitude
c  cor    : corlioris parameter
c  poly   : legendre polynomial
c  dpoly  : derative of legendre polynomial with latitude
c  trigs  : use by fft package
c  ifax   : used by fft package
c  np     : maxmin allowd typhoon number
c  ntyp   : number of typhoon in domain
c  rlon   : lontitude position of typhoon
c  rlat   : latitude position of typhoon
c  vm     : maxmun wind of rinkin vortex
c  rvm    : radius of maxmun wind of vortex
c  rm2    : maxnun radius of typhoon can reach
c  bb     : shap of rinkin vortex
c  xshap  : the ratio of x-y of topo
c  topamp : the height of topo
c  topmax : the radius of x of topo
c  topang : the angle of rotation of topo
c  topshap: control shap gradient of topo
c-------------------------------------------------------------------------
c
      include '../include/param.h'
      include '../include/const.h'
c
      character yes*1,name1*16,name2*16
c
      namelist/modlst/hmean,taui,taue,tauo,dt,baro,linear
     *, rnkut4,domod,dohdi,freq,tfilter,hfilt,bb,rm2
     *, ktop,fricv,fricd,bogus,umean,sr,numtyp
     *, t1lat,t1lon,t1vm,t1rvm,t2lat,t2lon,t2vm,t2rvm
     *, baro,linear,         ! Order is not important.

      data hmean /1000./
      data taui/0.0/, taue/48.0/, tauo/12.0/, dt/600.0/
      data baro/.true./, linear/.true./
      data rnkut4/.true./, domod/.false./, dohdi/.false./
      data tfilter/0.15/, hfilt/1.0e16/, freq/0.5/
      data rm2/6000./, bb/0.63/
      data numtyp/2/
      data t1lat/15.0/, t1lon/140.0/,t1vm/40.0/, t1rvm/1200.0/
      data t2lat/15.0/, t2lon/140.0/,t2vm/40.0/, t2rvm/1200.0/
      data ktop/0/,fricv/1.0e-5/,fricd/1.0e-5/
      data name1/'typhoon1'/, name2/'typhoon2'/
      data bogus/.true./
      data umean/0./
C
      open(unit=1,file='./namlsts',form='formatted')
      read(1,modlst,end=40)
 40   print modlst
      print*,'                               '
c
      ntyp=numtyp
c
      typname(1)=name1
      rlon(1)=t1lon
      rlat(1)=t1lat
      vm(1)=t1vm
      rvm(1)=t1rvm
c
      if(ntyp.eq.2)then
      typname(2)=name2
      rlon(2)=t2lon
      rlat(2)=t2lat
      vm(2)=t2vm
      rvm(2)=t2rvm
      endif
c
      call cons
c
      hmean=hmean*grav
c
      nxmy8=nx*my*8

      open(11,file='../DATAOUT/U.OUT',access='direct'
     &, form='unformatted',recl=nxmy8,status='unknown')
      open(12,file='../DATAOUT/V.OUT',access='direct'
     &, form='unformatted',recl=nxmy8,status='unknown')
      open(13,file='../DATAOUT/VOR.OUT',access='direct'
     &, form='unformatted',recl=nxmy8,status='unknown')
      open(14,file='../DATAOUT/PHI.OUT',access='direct'
     &, form='unformatted',recl=nxmy8,status='unknown')
      open(15,file='../DATAOUT/PV.OUT',access='direct'
     &, form='unformatted',recl=nxmy8,status='unknown')

      return
      end
