c
      common/swconst/rad,radsq,omega,cp,rgas,pi,capa,grav
     *, hmean,taui,taue,tauo,dt,tau,mlsort(jtrun,jtrun)
     *, eps4(mlmax),cim(mlmax),weight(my),sinl(my),onocos(my)
     *, cosl(my),cor(my),poly(mlmax,my/2),dpoly(mlmax,my/2)
     *, msort(mlmax),lsort(mlmax),tfilter,hfilt,freq
     *, baro,linear,rnkut4,domod,dohdi
     *, idtg,vm(np),rvm(np),rm2,bb,rlon(np),rlat(np)
     *, center(5000,2,np),ix(np),jy(np),lbot,ltop,ntyp,typname(np)
     *, ktop,fricv,fricd
     *, bogus,umean,sr
     *, epsd16(mlmax),ran(4004001),its
c
      logical baro,linear,rnkut4,domod,dohdi,bogus,sr
      character filen*32,typname*16
      common /files/filen
