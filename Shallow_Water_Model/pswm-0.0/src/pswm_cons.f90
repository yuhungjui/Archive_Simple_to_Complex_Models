!===============================================================================
module pswm_cons
!-------------------------------------------------------------------------------
! Purpose:
!    Contains PSWM model constants (not user-specifiable but constant)
!
! See also:
!    pswm_pars  model parameters (user-specifiable constants)
!    pswm_vars  model variables  (change during a run)
!-------------------------------------------------------------------------------

use kinds

implicit none   ! all module variables untyped by default
public          ! all module variables and routines accessible by default

!-------------------------------------------------------------------------------
! model constants
!-------------------------------------------------------------------------------

character (len=4),    parameter :: model = "pswm"     ! model name
real (kind=dbl_kind), parameter :: vnum  = 0.9        ! model version
character (len=8),    parameter :: vdate = "08/01/07" ! version date

integer (kind=int_kind) :: & ! indices of terms in the time splitting:
   iterm_fast = 0 , & ! index of the fast terms (usually treated implicitly)
   iterm_slow = 0 , & ! index of the slow terms (usually treated explicitly)
   iterm_diss = 0     ! index of dissipation terms (explicitly, separately)

integer (kind=int_kind) :: &
   num_sol = 0, &! number of solution copies (set in pswm_setup)
   num_for = 0   ! number of forcing  copies (set in pswm_setup)

real (kind=dbl_kind), allocatable, dimension(:), save :: x, y ! coordinates

!-------------------------------------------------------------------------------
end module pswm_cons
!===============================================================================
