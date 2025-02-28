      program exampl7
c
      dimension zdat(23,14),rwrk(1000),iwrk(1000),iama(12000)
      dimension iasf(13)
      dimension scra(1000),ycra(1000)
      dimension iaia(10),igia(10)
c
      dimension lind(14)
      character*10 llbs(15)
c
      external colram
c
      data iasf /13*1/
c
      data lind / 2,3,4,5,6,7,8,9,10,11,12,13,14,15 /
c
      call opngks
c
      call gsclip (0)
c
      call gsasf (iasf)
c
      call gsfais (1)
c
      call dfclrs
c
      call gendat (zdat,23,23,14,20,20,-136.148,451.834)
c
      time=second(dumi)
c
      call cpsetr ('vpb - viewport bottom',.25)
c
      call cpseti ('nor - numeric omission flags',0)
c
      call cpseti ('cls - contour level selector',-13)
c
      call cprect (zdat,23,23,14,rwrk,1000,iwrk,1000)
c
      call arinam (iama,12000)
      call cpclam (zdat,rwrk,iwrk,iama)
c
      call arscam (iama,xcra,ycra,1000,iaia,igia,10,colram)
c
      call gsplci (0)
      call cpcldr (zdat,rwrk,iwrk)
      call gsplci (1)
c
      call cpgetr ('zmn',zmin)
      call cpgetr ('zmx',zmax)
c
      do 102 i=1,15
         call cpsetr ('zdv - z data value',
     +                zmin+real(i-1)*(zmax-zmin)/14.)
         call cpgetc ('zdv - z data value',llbs(i))
 102  continue
c
      call lbseti ('cbl - color of box lines',0)
      call lblbar (0,.05,.95,.15,.25,14,1.,.5,lind,0,llbs,15,1)
c
      call capsap ('example 7',time,iama,12000)
      call labtop ('example 7',.017)
      call bndary
c
      call frame
c
      call clsgks
c
      stop
c
      end
c
c======================================================
      subroutine colram (xcra,ycra,ncra,iaia,igia,naia)
c
      dimension xcra(*),ycra(*),iaia(*),igia(*)
c
      ifll=1
c
      do 101 i=1,naia
         if (iaia(i).lt.0) ifll=0
 101  continue
c
      if (ifll.ne.0) then
        ifll=0
        do 103 i=1,naia
          if (igia(i).eq.3) ifll=iaia(i)
 103    continue
        if (ifll.gt.0.and.ifll.lt.15) then
          call gsfaci (ifll+1)
          call gfa (ncra-1,xcra,ycra)
        endif
      endif
c
      return
c
      end
