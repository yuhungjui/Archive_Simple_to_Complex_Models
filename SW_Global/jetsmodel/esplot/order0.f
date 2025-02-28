      program put_to_order

c--------------------------------------------------------------
c
c put the wave coefficients to their right position.
c
c written by LLMing
c 2005 Aug. 21
c
c--------------------------------------------------------------



      real enave(2080),enord(64,64),enplot(64,64)
      real temp,enmax

      open(45,file='enave0.dat',status='unknown',form='formatted') 
      open(46,file='enplot0.dat',status='unknown',form='formatted') 

      do i=1,2080,1
         read(45,*)enave(i)
      enddo

c-----------------put to order----------------------------------

      do i=1,63,1     ! 62
         enord(i,i)=enave(i)
      enddo
 
      do i=64,124,1    ! 60
         enord(i-63,i-61)=enave(i)    ! 2
      enddo

      do i=125,183,1     ! 58 
         enord(i-124,i-120)=enave(i)   ! 4
      enddo
 
      do i=184,240,1    ! 56
         enord(i-183,i-177)=enave(i)   ! 6
      enddo
 
      do i=241,295,1    ! 54
         enord(i-240,i-232)=enave(i)   ! 8
      enddo
 
      do i=296,348,1    ! 52
         enord(i-295,i-285)=enave(i)   ! 10
      enddo
 
      do i=349,399,1    ! 50
         enord(i-348,i-336)=enave(i)   ! 12
      enddo
 
      do i=400,448,1    ! 48
         enord(i-399,i-385)=enave(i)   ! 14
      enddo
    
c---------ok----------------------------------------- 
 
      do i=449,495,1    ! 46
         enord(i-448,i-432)=enave(i)   ! 16
      enddo
  
      do i=496,540,1    ! 44
         enord(i-495,i-477)=enave(i)   ! 18
      enddo
  
      do i=541,583,1    ! 42
         enord(i-540,i-520)=enave(i)   ! 20
      enddo
  
      do i=584,624,1    ! 40
         enord(i-583,i-561)=enave(i)   ! 22
      enddo
  
      do i=625,663,1    ! 38
         enord(i-624,i-600)=enave(i)   ! 24
      enddo
  
      do i=664,700,1    ! 36
         enord(i-663,i-637)=enave(i)   ! 26
      enddo
  
      do i=701,735,1    ! 34
         enord(i-700,i-672)=enave(i)   ! 28
      enddo

c-----------------------------------------------------------------      
  
      do i=736,768,1    ! 32
         enord(i-735,i-705)=enave(i)   ! 30
      enddo
  
      do i=769,799,1    ! 30
         enord(i-768,i-736)=enave(i)   ! 32
      enddo

      do i=800,828,1    ! 28
         enord(i-799,i-765)=enave(i)   ! 34
      enddo
  
      do i=829,855,1    ! 26
         enord(i-828,i-792)=enave(i)   ! 36
      enddo
  
      do i=856,880,1    ! 24
         enord(i-855,i-817)=enave(i)   ! 38
      enddo
c---------------------------------------------------------------
  
      do i=881,903,1    ! 22
         enord(i-880,i-840)=enave(i)   ! 40
      enddo
  
      do i=904,924,1    ! 20
         enord(i-903,i-861)=enave(i)   ! 42
      enddo
  
      do i=925,943,1    ! 18
         enord(i-924,i-880)=enave(i)   ! 44
      enddo
  
      do i=944,960,1    ! 16
         enord(i-943,i-897)=enave(i)   ! 46
      enddo
  
      do i=961,975,1    ! 14
         enord(i-960,i-912)=enave(i)   ! 48
      enddo
c----------------------------------------------------------------
  
      do i=976,988,1    ! 12
         enord(i-975,i-925)=enave(i)   ! 50
      enddo
  
      do i=989,999,1    ! 10
         enord(i-988,i-936)=enave(i)   ! 52
      enddo
  
      do i=1000,1008,1    ! 8
         enord(i-999,i-945)=enave(i)   ! 54
      enddo
  
      do i=1009,1015,1    ! 6
         enord(i-1008,i-952)=enave(i)   ! 56
      enddo
  
      do i=1016,1020,1    ! 4
         enord(i-1015,i-957)=enave(i)   ! 58
      enddo
  
      do i=1021,1023,1    ! 2
         enord(i-1020,i-960)=enave(i)   ! 60
      enddo
  
      do i=1024,1024,1    ! 0
         enord(i-1023,i-961)=enave(i)   ! 62
      enddo

c---------------------------------------------------------------------------------

c---------------------------------------------------------------------------------
  
      do i=1025,1087,1    ! 62
         enord(i-1024,i-1023)=enave(i)   ! 1
      enddo
  
      do i=1088,1148,1    ! 60
         enord(i-1087,i-1084)=enave(i)   ! 3
      enddo
  
      do i=1149,1207,1    ! 58
         enord(i-1148,i-1143)=enave(i)   ! 5
      enddo
  
      do i=1208,1264,1    ! 56
         enord(i-1207,i-1200)=enave(i)   ! 7
      enddo
  
      do i=1265,1319,1    ! 54
         enord(i-1264,i-1255)=enave(i)   ! 9
      enddo
c----------------------------------------------------------------------
  
      do i=1320,1372,1    ! 52
         enord(i-1319,i-1308)=enave(i)   ! 11
      enddo
  
      do i=1373,1423,1    ! 50
         enord(i-1372,i-1359)=enave(i)   ! 13
      enddo
  
      do i=1424,1472,1    ! 48
         enord(i-1423,i-1408)=enave(i)   ! 15
      enddo

      do i=1473,1519,1    ! 46
         enord(i-1472,i-1455)=enave(i)   ! 17
      enddo
  
      do i=1520,1564,1    ! 44
         enord(i-1519,i-1500)=enave(i)   ! 19
      enddo
c----------------------------------------------------------------------------
  
      do i=1565,1607,1    ! 42
         enord(i-1564,i-1543)=enave(i)   ! 21
      enddo
  
      do i=1608,1648,1    ! 40
         enord(i-1607,i-1584)=enave(i)   ! 23
      enddo
  
      do i=1649,1687,1    ! 38
         enord(i-1648,i-1623)=enave(i)   ! 25
      enddo
  
      do i=1688,1724,1    ! 36
         enord(i-1687,i-1660)=enave(i)   ! 27
      enddo
  
      do i=1725,1759,1    ! 34
         enord(i-1724,i-1695)=enave(i)   ! 29
      enddo
c---------------------------------------------------------------------------
  
      do i=1760,1792,1    ! 32
         enord(i-1759,i-1728)=enave(i)   ! 31
      enddo
  
      do i=1793,1823,1    ! 30
         enord(i-1792,i-1759)=enave(i)   ! 33
      enddo
  
      do i=1824,1852,1    ! 28
         enord(i-1823,i-1788)=enave(i)   ! 35
      enddo
  
      do i=1853,1879,1    ! 26
         enord(i-1852,i-1815)=enave(i)   ! 37
      enddo
  
      do i=1880,1904,1    ! 24
         enord(i-1879,i-1840)=enave(i)   ! 39
      enddo
  
c---------------------------------------------------------------------

      do i=1905,1927,1    ! 22
         enord(i-1904,i-1863)=enave(i)   ! 41
      enddo
  
      do i=1928,1948,1    ! 20
         enord(i-1927,i-1884)=enave(i)   ! 43
      enddo
  
      do i=1949,1967,1    ! 18
         enord(i-1948,i-1903)=enave(i)   ! 45
      enddo
  
      do i=1968,1984,1    ! 16
         enord(i-1967,i-1920)=enave(i)   ! 47
      enddo
  
      do i=1985,1999,1    ! 14
         enord(i-1984,i-1935)=enave(i)   ! 49
      enddo
c-----------------------------------------------------------------------
  
      do i=2000,2012,1    ! 12
         enord(i-1999,i-1948)=enave(i)   ! 51
      enddo
  
      do i=2013,2023,1    ! 10
         enord(i-2012,i-1959)=enave(i)   ! 53
      enddo
  
      do i=2024,2032,1    ! 8
         enord(i-2023,i-1968)=enave(i)   ! 55
      enddo
  
      do i=2033,2039,1    ! 6
         enord(i-2032,i-1975)=enave(i)   ! 57
      enddo
  
      do i=2040,2044,1    ! 4
         enord(i-2039,i-1980)=enave(i)   ! 59
      enddo
  
      do i=2045,2047,1    ! 2
         enord(i-2044,i-1983)=enave(i)   ! 61
      enddo
  
      do i=2048,2048,1    ! 0
         enord(i-2047,i-1984)=enave(i)   ! 63
      enddo
c-------------------------------------------------------------------------------
c----------------------------------------------------for n=64

      do i=2,64,1
         enord(i,64)=enave(i/2+2048)
      enddo

c--------------------------------------------------------------

      temp=enord(1,1)
      do i=1,64,1
         do j=1,64,1
            if (enord(i,j) .GT. temp) temp=enord(i,j)
         enddo
      enddo
      enmax=temp

      print*,'enmax=',enmax

      do i=1,64,1
         do j=1,64,1
            enplot(i,j)=enord(i,j)/enmax
         enddo
      enddo

      do j=1,64,1
         do i=1,64,1
            write(46,*)enplot(i,j)
         enddo
      enddo

      stop
      end
