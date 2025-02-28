      program mkrandom

      integer n
      parameter(n=4004001)
      real ran1(n),ran(n)

      open(21,file='1random.dat')
      open(22,file='random.dat')

      do i=1,n,1
         read(21,*)ran1(i)
         ran(i)=ran1(i)*2.-1.  
         write(22,*)ran(i)
      enddo

      stop
      end
