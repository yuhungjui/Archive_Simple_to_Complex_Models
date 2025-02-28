      subroutine rhs(linear)
      logical linear
C
      if(linear)then
        call lrhs
      else
        call nrhs
      endif
c
      return
      end
