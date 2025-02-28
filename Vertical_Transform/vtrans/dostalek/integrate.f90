MODULE integrate
!===============================================================================
! integrate.f90		Jack Dostalek
!
! This module creates and explicit interface to the integration routines in
! intlib.f90.
!
! History
! 16 Sep 2008	Programming begun
!
!===============================================================================
!
INTERFACE
!
!-------------------------------------------------------------------------------
SUBROUTINE cubint(ftab,xtab,ntab,ia,ib,result,error)

USE kinds

IMPLICIT NONE

REAL(wrk), DIMENSION(ntab), INTENT(IN) :: ftab
REAL(wrk), DIMENSION(ntab), INTENT(IN) :: xtab
INTEGER, INTENT(IN) :: ntab
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ib
REAL(wrk), INTENT(OUT) :: result
REAL(wrk), INTENT(OUT) :: error

END SUBROUTINE cubint
!-------------------------------------------------------------------------------
SUBROUTINE simpsn(h,y,num,result)

USE kinds

IMPLICIT NONE

REAL(wrk), INTENT(IN) :: h
REAL(wrk), DIMENSION(num), INTENT(IN) :: y
INTEGER, INTENT(IN) :: num
REAL(wrk), INTENT(OUT) :: result

END SUBROUTINE simpsn
!-------------------------------------------------------------------------------
END INTERFACE

END MODULE integrate
