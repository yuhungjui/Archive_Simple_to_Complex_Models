!===============================================================================
module kinds
!-------------------------------------------------------------------------------
! Purpose: 
!    Defines standard kinds of variables.
!-------------------------------------------------------------------------------
integer, parameter :: int_kind  = kind(1),                 & ! default integer
                      log_kind  = kind(.true.),            & ! default logical
                      real_kind = selected_real_kind(6),   & ! default real
                      dbl_kind  = selected_real_kind(13)     ! default double
!-------------------------------------------------------------------------------
end module kinds
!===============================================================================
