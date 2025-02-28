      subroutine qmaxn3 (fld,t1,t2,i1,j1,k1,im,jm,lm)
c
c  subroutine to print max and min in 2-d layer (or a subarray)
c  within a 3-d field stored with n-s index slowest varying
c
c  ***input***
c
c  fld: input array
c  t1,t2: 2* 8 character caption
c  i1: starting index of first dimension (e-w)
c  j1: starting index of third dimension (n-s)
c  k1: starting index of second dimension (vertical)
c  im: first dimension
c  jm: third dimension
c  lm: second dimension
c
c *****************************************************************
c
cfj
      dimension fld(im,lm,jm)
      character*8 t1, t2
c
      xmin= 1.0e25
      xmax= -1.0e25
      imin=1
      jmin=1
      imax=1
      jmax=1
c
      do 10 j=j1,jm
      do 10 i=i1,im
      if (fld(i,k1,j).le.xmin) then
      xmin= fld(i,k1,j)
      jmin= j
      imin= i
      endif
      if (fld(i,k1,j).gt.xmax) then
      xmax= fld(i,k1,j)
      jmax= j
      imax= i
      endif
   10 continue
cfj
      print 9000, t1, t2
      print 8995, imax,jmax,xmax,imin,jmin,xmin
cfj
 9000 format (1h0,2a8)
 8995 format(' imax=',i4,' jmax=',i4,' xlarg=',g20.12
     *,' imin=',i4,' jmin=',i4,' xsmal=',g20.12)
c
      return
      end
