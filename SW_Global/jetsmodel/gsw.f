      program gsw
c
      include '../include/param.h'
      include '../include/const.h'
c
      call prepare
c
      call rdy2go

      if(rnkut4)then
        call rk4time
        print *,'call rk4time'
      else
        call lptime
        print *,'call lptime'
      endif

      stop
      end
