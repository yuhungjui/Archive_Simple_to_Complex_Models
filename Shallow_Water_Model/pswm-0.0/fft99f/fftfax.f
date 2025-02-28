      subroutine fftfax(n,ifax,trigs)
      integer           n, ifax(13)
      SINGLE_OR_DOUBLE  trigs(3*n/2+1)
c
c mode 3 is used for real/half-complex transforms.  It is possible
c to do complex/complex transforms with other values of mode, but
c documentation of the details were not available when this routine
c was written.
c
      integer nfax, mode
      save    mode
      data    mode /3/
c
      call fax (ifax, n, mode)
      nfax = ifax(1)
      if ( nfax.le.0 .or. ifax(nfax+1).gt.5 ) then
          print *, ' fftfax error:  invalid transform length n =', n
          print *, ' ============>  transforms will not be computed'
          ifax(1) = -99
          return
      end if
      call fftrig (trigs, n, mode)
      return
      end
