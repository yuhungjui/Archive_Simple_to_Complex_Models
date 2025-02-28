MODULE vnm_transform
!===============================================================================
! vnm_transform.f90	
!
! This module "oversees" the vertical normal mode transform of the forcing for
! the global omega equation.
!
!===============================================================================
USE kinds
USE constants
USE integrate
USE glomega_grid
USE glomega_read
USE verttrans

IMPLICIT NONE
SAVE

REAL(sgl), DIMENSION(npmod) :: tbar1		!Base state temperature
REAL(dbl), DIMENSION(np) :: tbar		!Base state temperature
REAL(dbl), DIMENSION(np,5) :: rt
REAL(dbl), DIMENSION(4) :: con
REAL(dbl), DIMENSION(0:himode) :: cl
REAL(dbl), DIMENSION(0:himode,0:himode) :: tfor
REAL(dbl), DIMENSION(0:himode,0:himode) :: tinv
REAL(dbl), DIMENSION(0:himode) :: fvbf		!Vertical basis functions
REAL(dbl), DIMENSION(0:himode,0:himode) :: prom	!Projection matrix
REAL(dbl), DIMENSION(np,2) :: fval
INTEGER :: ierr
REAL(dbl), DIMENSION(np,himode+4) :: eigarr	!Holds vert. struct. info.
REAL(dbl), DIMENSION(np,0:himode) :: vlp	!Derivatives of eigenfunctions
REAL(dbl), DIMENSION(0:himode) :: epsl		!Lamb's parameter

! Variables for Gauss-Legendre quadrature
INTEGER, PARAMETER :: ngaus = 2*himode		!Number of gaussian levels
REAL(dbl), DIMENSION(ngaus) :: glp
REAL(dbl), DIMENSION(ngaus) :: glwt

!Working variables
INTEGER, PRIVATE :: i,j,k,l,m,n			!Indexing variables
REAL(dbl), PRIVATE, DIMENSION(:), ALLOCATABLE :: w1
REAL(dbl), PRIVATE, DIMENSION(:,:), ALLOCATABLE :: f2 
REAL(dbl), PRIVATE, DIMENSION(:,:), ALLOCATABLE :: g2 

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE vert_transform

! Jordan 1958 "hurricane season" sounding
tbar1(1)=-73.5_wrk+theta0  
tbar1(2)=-67.6_wrk+theta0 
tbar1(3)=-55.2_wrk+theta0 
tbar1(4)=-43.3_wrk+theta0 
tbar1(5)=-33.2_wrk+theta0 
tbar1(6)=-17.7_wrk+theta0 
tbar1(7)=-6.9_wrk+theta0  
tbar1(8)=8.6_wrk+theta0   
tbar1(9)=17.3_wrk+theta0  
tbar1(10)=21.4_wrk+theta0 
tbar1(11)=26.0_wrk+theta0 

! The basic state temperature, tbar1, is at the model levels. Use a cubic
! spline to interpolate to the computational vertical levels.
ALLOCATE(f2(npmod,4))
DO k=1,npmod
  f2(k,1)=DBLE(tbar1(k))
ENDDO 

CALL csset(npmod,DBLE(pmod),f2,3,3,ierr)
DO k=1,np
  tbar(k)=csval(DBLE(p(k)),0,npmod,DBLE(pmod),f2,k)
ENDDO
DEALLOCATE(f2)

! Compute Gauss-Legendre abscissas (pressures) and weights for integration 
! between ptop and pbot
CALL gaussl(ngaus,DBLE(ptop),DBLE(pbot),glp,glwt,ierr)

! Fill in first column of array rt, which contains the temperature times the gas
! constant
DO k=1,np
  rt(k,1)=DBLE(rd)*tbar(k)
  rt(k,2)=DBLE(zero)	
ENDDO  

! Set up vertical transform
ALLOCATE(f2(0:himode,0:1))
ALLOCATE(g2(0:himode,0:himode)) 
CALL vtset(np,DBLE(p),rt,himode,f2,g2,con,cl,tfor,tinv,ierr)

DEALLOCATE(g2)
DEALLOCATE(f2) 

! Array eigarr will contain some values relevant to the vertical transform.
! The values include pressure, temperature, stability, and the eigenfunctions.
! Pressure goes in first column and the base-state temperature in second column
DO k=1,np
  eigarr(k,1)=DBLE(p(k))
  eigarr(k,2)=tbar(k)
ENDDO 

! Next is static stability. I do not need to call csset, as vtset has set up the
! cubic spline interpolate for rt.
DO k=np,1,-1
  eigarr(k,3)=DBLE(kap)*csval(DBLE(p(k)),0,np,DBLE(p),rt,0) & 
             -DBLE(p(k))*csval(DBLE(p(k)),1,np,DBLE(p),rt,0)		          
  eigarr(k,3)=eigarr(k,3)/DBLE(rd)
ENDDO   

! Get eigenfunctions, derivatives of eigenfunctions with respect to pressure and 
! Lamb's parameter for each vertical mode  
DO l=0,himode
  ALLOCATE(w1(2*(himode+1)))
  ALLOCATE(f2(np,2))
  CALL veval(tinv(0,l),1,np,DBLE(p),con,himode,w1,f2)
  DEALLOCATE(w1)
  DO k=np,1,-1
    eigarr(k,4+l)=f2(k,1)
    vlp(k,l)=f2(k,2)
  ENDDO 
  DEALLOCATE(f2)
  epsl(l)=(DBLE(two*erot*erad)/cl(l))**2   
ENDDO 

WRITE(lulog,500) himode
500 FORMAT('Last vertical mode kept (himode) for subroutine vtset: ',I2)
WRITE(lulog,*)

END SUBROUTINE vert_transform

!-------------------------------------------------------------------------------        
END MODULE vnm_transform




