      program onedenspc

      integer jtrun,mlmax
      parameter(jtrun=64,mlmax=2080)
      integer msort(mlmax),lsort(mlmax),mlsort(jtrun,jtrun)
      real enave(mlmax),oned(jtrun)

      open(35,file='enave0.dat')
      open(36,file='onedspc0.dat')
 
      do i=1,mlmax,1
         read(35,*)enave(i)
      enddo

      do i=1,jtrun,1
         oned(i)=0.
      enddo

      do i=1,mlmax,1  
         call sortml(jtrun,i,msort,lsort,mlsort)

         if (lsort(i) .eq. 1) then
              oned(1)=enave(i)      

         else if (lsort(i) .eq. 2) then
              oned(2)=oned(2)+enave(i)
         else if (lsort(i) .eq. 3) then
              oned(3)=oned(3)+enave(i)
         else if (lsort(i) .eq. 4) then
              oned(4)=oned(4)+enave(i)
         else if (lsort(i) .eq. 5) then
              oned(5)=oned(5)+enave(i)
         else if (lsort(i) .eq. 6) then
              oned(6)=oned(6)+enave(i)
         else if (lsort(i) .eq. 7) then
              oned(7)=oned(7)+enave(i)
         else if (lsort(i) .eq. 8) then
              oned(8)=oned(8)+enave(i)
         else if (lsort(i) .eq. 9) then
              oned(9)=oned(9)+enave(i)
         else if (lsort(i) .eq. 10) then
              oned(10)=oned(10)+enave(i)

         else if (lsort(i) .eq. 11) then
              oned(11)=oned(11)+enave(i)
         else if (lsort(i) .eq. 12) then
              oned(12)=oned(12)+enave(i)
         else if (lsort(i) .eq. 13) then
              oned(13)=oned(13)+enave(i)
         else if (lsort(i) .eq. 14) then
              oned(14)=oned(14)+enave(i)
         else if (lsort(i) .eq. 15) then
              oned(15)=oned(15)+enave(i)
         else if (lsort(i) .eq. 16) then
              oned(16)=oned(16)+enave(i)
         else if (lsort(i) .eq. 17) then
              oned(17)=oned(17)+enave(i)
         else if (lsort(i) .eq. 18) then
              oned(18)=oned(18)+enave(i)
         else if (lsort(i) .eq. 19) then
              oned(19)=oned(19)+enave(i)
         else if (lsort(i) .eq. 20) then
              oned(20)=oned(20)+enave(i)

         else if (lsort(i) .eq. 21) then
              oned(21)=oned(21)+enave(i)
         else if (lsort(i) .eq. 22) then
              oned(22)=oned(22)+enave(i)
         else if (lsort(i) .eq. 23) then
              oned(23)=oned(23)+enave(i)
         else if (lsort(i) .eq. 24) then
              oned(24)=oned(24)+enave(i)
         else if (lsort(i) .eq. 25) then
              oned(25)=oned(25)+enave(i)
         else if (lsort(i) .eq. 26) then
              oned(26)=oned(26)+enave(i)
         else if (lsort(i) .eq. 27) then
              oned(27)=oned(27)+enave(i)
         else if (lsort(i) .eq. 28) then
              oned(28)=oned(28)+enave(i)
         else if (lsort(i) .eq. 29) then
              oned(29)=oned(29)+enave(i)
         else if (lsort(i) .eq. 30) then
              oned(30)=oned(30)+enave(i)

         else if (lsort(i) .eq. 31) then
              oned(31)=oned(31)+enave(i)
         else if (lsort(i) .eq. 32) then
              oned(32)=oned(32)+enave(i)
         else if (lsort(i) .eq. 33) then
              oned(33)=oned(33)+enave(i)
         else if (lsort(i) .eq. 34) then
              oned(34)=oned(34)+enave(i)
         else if (lsort(i) .eq. 35) then
              oned(35)=oned(35)+enave(i)
         else if (lsort(i) .eq. 36) then
              oned(36)=oned(36)+enave(i)
         else if (lsort(i) .eq. 37) then
              oned(37)=oned(37)+enave(i)
         else if (lsort(i) .eq. 38) then
              oned(38)=oned(38)+enave(i)
         else if (lsort(i) .eq. 39) then
              oned(39)=oned(39)+enave(i)
         else if (lsort(i) .eq. 40) then
              oned(40)=oned(40)+enave(i)

         else if (lsort(i) .eq. 41) then
              oned(41)=oned(41)+enave(i)
         else if (lsort(i) .eq. 42) then
              oned(42)=oned(42)+enave(i)
         else if (lsort(i) .eq. 43) then
              oned(43)=oned(43)+enave(i)
         else if (lsort(i) .eq. 44) then
              oned(44)=oned(44)+enave(i)
         else if (lsort(i) .eq. 45) then
              oned(45)=oned(45)+enave(i)
         else if (lsort(i) .eq. 46) then
              oned(46)=oned(46)+enave(i)
         else if (lsort(i) .eq. 47) then
              oned(47)=oned(47)+enave(i)
         else if (lsort(i) .eq. 48) then
              oned(48)=oned(48)+enave(i)
         else if (lsort(i) .eq. 49) then
              oned(49)=oned(49)+enave(i)
         else if (lsort(i) .eq. 50) then
              oned(50)=oned(50)+enave(i)

         else if (lsort(i) .eq. 51) then
              oned(51)=oned(51)+enave(i)
         else if (lsort(i) .eq. 52) then
              oned(52)=oned(52)+enave(i)
         else if (lsort(i) .eq. 53) then
              oned(53)=oned(53)+enave(i)
         else if (lsort(i) .eq. 54) then
              oned(54)=oned(54)+enave(i)
         else if (lsort(i) .eq. 55) then
              oned(55)=oned(55)+enave(i)
         else if (lsort(i) .eq. 56) then
              oned(56)=oned(56)+enave(i)
         else if (lsort(i) .eq. 57) then
              oned(57)=oned(57)+enave(i)
         else if (lsort(i) .eq. 58) then
              oned(58)=oned(58)+enave(i)
         else if (lsort(i) .eq. 59) then
              oned(59)=oned(59)+enave(i)
         else if (lsort(i) .eq. 60) then
              oned(60)=oned(60)+enave(i)

         else if (lsort(i) .eq. 61) then
              oned(61)=oned(61)+enave(i)
         else if (lsort(i) .eq. 62) then
              oned(62)=oned(62)+enave(i)
         else if (lsort(i) .eq. 63) then
              oned(63)=oned(63)+enave(i)

         else
              oned(64)=oned(64)+enave(i)
         end if
      enddo

      do i=1,jtrun,1
         write(36,*)oned(i)
      enddo

      stop
      end
