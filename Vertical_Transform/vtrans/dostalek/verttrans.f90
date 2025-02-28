MODULE verttrans
!===============================================================================
! verttrans.f90		Jack Dostalek
!
! This module creates an explicit interface to the vertical transform routines
! written by Scott Fulton.
!
! History:
! 20 Jun 2008 	Programming begun
!
!===============================================================================
!
INTERFACE
!
!-------------------------------------------------------------------------------
SUBROUTINE vtset(m,p,rt,n,w,work,con,c,tfor,tinv,ierr) 

USE kinds

IMPLICIT NONE

INTEGER, INTENT(IN) :: m
REAL(dbl), DIMENSION(m), INTENT(IN) :: p
REAL(dbl), DIMENSION(m,5), INTENT(INOUT) :: rt
INTEGER, INTENT(IN) :: n
REAL(dbl), DIMENSION(0:n,0:1), INTENT(OUT) :: w
REAL(dbl), DIMENSION(0:n,0:n), INTENT(OUT) :: work
REAL(dbl), DIMENSION(4), INTENT(OUT) :: con
REAL(dbl), DIMENSION(0:n), INTENT(OUT) :: c
REAL(dbl), DIMENSION(0:n,0:n), INTENT(OUT) :: tfor
REAL(dbl), DIMENSION(0:n,0:n), INTENT(OUT) :: tinv
INTEGER, INTENT(OUT) :: ierr

END SUBROUTINE vtset
!-------------------------------------------------------------------------------
SUBROUTINE veval(fvbf,md,m,p,con,n,w,fval)
USE kinds
IMPLICIT NONE

REAL(dbl), DIMENSION(0:n), INTENT(IN) :: fvbf
INTEGER, INTENT(IN) :: md
INTEGER, INTENT(IN) :: m
REAL(dbl), DIMENSION(m), INTENT(IN) :: p
REAL(dbl), DIMENSION(4), INTENT(IN) :: con
INTEGER, INTENT(IN) :: n
REAL(dbl), DIMENSION(0:n,2), INTENT(OUT) :: w
REAL(dbl), DIMENSION(m,2), INTENT(OUT) :: fval

END SUBROUTINE veval

!-------------------------------------------------------------------------------
FUNCTION rrgval(pval,con,m,rt,ierr) 
USE kinds
IMPLICIT NONE
REAL(dbl) :: rrgval
REAL(dbl), INTENT(IN) :: pval
REAL(dbl), DIMENSION(3), INTENT(IN) :: con
INTEGER, INTENT(IN) :: m
REAL(dbl), DIMENSION(m,5), INTENT(IN) :: rt
INTEGER, INTENT(OUT) :: ierr

END FUNCTION rrgval
!-------------------------------------------------------------------------------
SUBROUTINE csset(n,x,f,ibc1,ibcn,ierr)
USE kinds
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(dbl), DIMENSION(n), INTENT(IN) :: x
REAL(dbl), DIMENSION(n,4), INTENT(INOUT) :: f
INTEGER, INTENT(IN) :: ibc1
INTEGER, INTENT(IN) :: ibcn
INTEGER, INTENT(OUT) :: ierr

END SUBROUTINE csset
!-------------------------------------------------------------------------------
FUNCTION csval(xval,md,n,x,f,isw)
USE kinds
IMPLICIT NONE
REAL(dbl) :: csval
REAL(dbl), INTENT(IN) :: xval
INTEGER, INTENT(IN) :: md
INTEGER, INTENT(IN) :: n
REAL(dbl), DIMENSION(n), INTENT(IN) :: x
REAL(dbl), DIMENSION(n,4), INTENT(IN) :: f
INTEGER, INTENT(IN) :: isw

END FUNCTION csval

!-------------------------------------------------------------------------------
SUBROUTINE vproj(m,p,f,in,ip,ng,pg,wg,prom,n,con,w,fvbf,ierr)
USE kinds
IMPLICIT NONE
INTEGER, INTENT(IN) :: m
REAL(dbl), DIMENSION(m), INTENT(IN) :: p
REAL(dbl), DIMENSION(m,4), INTENT(INOUT) :: f
INTEGER, INTENT(IN) :: in
INTEGER, INTENT(IN) :: ip
INTEGER, INTENT(IN) :: ng
REAL(dbl), DIMENSION(ng), INTENT(INOUT) :: pg
REAL(dbl), DIMENSION(ng), INTENT(INOUT) :: wg
REAL(dbl), DIMENSION(0:n,0:n), INTENT(INOUT) :: prom
INTEGER, INTENT(IN) :: n
REAL(dbl), DIMENSION(3), INTENT(IN) :: con
REAL(dbl), DIMENSION(n+1), INTENT(IN) :: w
REAL(dbl), DIMENSION(0:n), INTENT(OUT) :: fvbf
INTEGER, INTENT(OUT) :: ierr

END SUBROUTINE vproj
!-------------------------------------------------------------------------------
SUBROUTINE vtran(f,n,tmat,g) 
USE kinds
IMPLICIT NONE
REAL(dbl), DIMENSION(0:n), INTENT(IN) :: f
INTEGER, INTENT(IN) :: n
REAL(dbl), DIMENSION(0:n,0:n), INTENT(IN) :: tmat
REAL(dbl), DIMENSION(0:n), INTENT(OUT) :: g

END SUBROUTINE vtran
!-------------------------------------------------------------------------------
SUBROUTINE gaussl(n,xa,xb,ab,wt,ierr)
USE kinds
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(dbl), INTENT(IN) :: xa
REAL(dbl), INTENT(IN) :: xb
REAL(dbl), DIMENSION(n), INTENT(OUT) :: ab
REAL(dbl), DIMENSION(n), INTENT(OUT) :: wt
INTEGER, INTENT(OUT) :: ierr
END SUBROUTINE gaussl
!-------------------------------------------------------------------------------
END INTERFACE

END MODULE verttrans
