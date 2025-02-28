!===============================================================================
module pswm_vars
!-------------------------------------------------------------------------------
! Purpose:
!    Contains PSWM model variables (those which change during a run)
!
! See also:
!    pswm_pars  model parameters (user-specifiable constants)
!    pswm_cons  model constants  (not user-specifiable, constant)
!-------------------------------------------------------------------------------

use kinds

implicit none   ! all module variables untyped by default
public          ! all module variables and routines accessible by default
save            ! all module variables static

!-------------------------------------------------------------------------------
! model variables
!-------------------------------------------------------------------------------

! variables for time discretization
real    (kind=dbl_kind) :: time   ! model time
integer (kind=int_kind) :: nstep  ! current time step index
integer (kind=int_kind) :: now    ! index of current solution copy

real (kind=dbl_kind), allocatable, dimension(:,:,:) :: & 
   u, v, p, d, z, q       ! prognostic variables
real (kind=dbl_kind), allocatable, dimension(:,:,:) :: & 
   uf, vf, pf, df, zf, qf ! corresponding forcing (tendencies)

!-------------------------------------------------------------------------------
end module pswm_vars
!===============================================================================
