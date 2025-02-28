      real function second( )
      real etime, tarray(2)
c
c   Returns CPU time in seconds--UNIX version
c
      second = etime( tarray )
      second = tarray(1)
      return
      end
