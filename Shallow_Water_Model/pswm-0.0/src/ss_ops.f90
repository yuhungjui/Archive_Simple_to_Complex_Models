!===============================================================================
module ss_ops
!-------------------------------------------------------------------------------
! Purpose: 
!    Contains routines to compute operations in 2D Fourier spectral space.
!
! Description:
!    Fields are stored in 2D arrays of the form f(-1:nx,-1:ny).
!    When in physical space, values are explicitly periodic:
!
!       f(j,k):  value at (x,y) = (j*dx,k*dy)  j=-1,...,nx, k=-1,...,ny
!
!    The x-transform gives coefficients of Fourier modes l=0,...,mx:
!
!       f(2*l-1,k):  x-mode l (real part) at y = k*dy  k=-1,...,ny
!       f(2*l  ,k):  x-mode l (imag part) at y = k*dy  k=-1,...,ny
!
!    The y-transform then gives four real numbers which give the coefficients
!    of the Fourier modes  fhat_(l,m) for  l=0, ..., mx  as follows:
!       
!       for  m=0,...,my:  real part:  f(2*l-1,2*m-1) - f(2*l  ,2*m)
!                         imag part:  f(2*l  ,2*m-1) + f(2*l-1,2*m)
!       
!       for -m=0,...,my:  real part:  f(2*l-1,2*m-1) + f(2*l  ,2*m)
!                         imag part:  f(2*l  ,2*m-1) - f(2*l-1,2*m)
!
!    The values for  -l=0, ..., mx  could be obtained using the fact
!    that  fhat_(-l,-m)  is the complex conjugate of  fhat_(l,m).
!    Note that these complex numbers never need to be formed explicitly.
!
!    The transform lengths  nx  and  ny  are determined from the size 
!    of the array  f,  which must be dimensioned as  f(-1:nx,-1:ny).
!    To transform only part of the array, input that section only.
!
!    The spectral truncation  mx  and  my  can be specified in the 
!    argument of each routine, or set "globally" (for all routines)
!    by calling setup_ss_ops beforehand (values specified in the argument
!    list take precedence).  If not specified in either of these ways
!    (or specified nonpositive), the default spectral truncation is  
!    mx = nx/2  and  my = my/2  (larger values are flagged as errors).
!    To avoid aliasing in quadratic terms, use  nx > 3*mx, ny > 3*my.
!
!    Likewise, the domain length  xl  and  yl  needed for derivatives
!    can be specified in the argument of the derivative routines, or set 
!    globally by calling setup_ss_ops beforehand (the argument list
!    takes precedence).  If not specified in either of these ways
!    (or specified nonpositive), the default length is  xl = yl = 2*pi.
!
! Use: 
!    The routine  setup_ss_ops  can be used to set the spectral truncation
!    and domain lengths for all routines as described above.
!
!    The use of each public routine (see below) is described in its comments.
!    No initialization is required, but you may wish to call unset_fft when
!    done to release space allocated for FFTs (for a clean exit).
!
! Required routines:  fft99 and fftfax from Temperton FFT package fft99f
!
! Author: 
!    Scott R. Fulton
!
! History: 
!    Fortran77 version written in 1995 for explicit model swm/fg/rot.
!    Converted to Fortran95 in 2007 for the semi-implicit model pswm.
!    Added routines setup_fft and unset_fft to manage storage dynamically.
!-------------------------------------------------------------------------------

use kinds, only : int_kind, dbl_kind

implicit none
private         ! all entry points and variables are local by default

!-------------------------------------------------------------------------------
public           & ! routines to be called from elsewhere:
   setup_ss_ops, & ! transform to spectral space
   tospec,       & ! transform to spectral space
   tophys,       & ! transform to physical space
   xderiv,       & ! x-derivative (in spectral space)
   yderiv,       & ! y-derivative (in spectral space)
   nlap,         & ! negative Laplacian (in spectral space)
   mhsolve,      & ! modified Helmholtz solver (in spectral space)
   unset_fft       ! reclaims space allocated for FFTs (for a clean exit)
!-------------------------------------------------------------------------------

! local variables

integer (kind=int_kind) :: nx, ny              ! transform lengths
integer (kind=int_kind) :: lx, ly              ! local copy of 2*mx, 2*my
integer (kind=int_kind) :: l, m                ! Fourier mode indices
integer (kind=int_kind) :: nx_max=0, ny_max=0  ! max field size (for work only)
real (kind=dbl_kind), allocatable, save :: work(:,:) ! work space for fft99

! Note:  Scaling the last Fourier mode by two is needed in some contexts
! (e.g., Fourier psuedospectral multigrid methods); see details in:
!    Brandt, A., S. R. Fulton, and G. D. Taylor, 1985:
!    Improved spectral multigrid methods for periodic elliptic problems.
!    J. Comp. Phys., vol. 58, pp. 96-112.
! In most cases this is not needed (and may even degrade the performance).
! The following parameter may be set by calling  setup_ss_ops:
logical :: scale_last_mode = .false. ! scale modes  mx  and  my  by two?

! dynamic variables for storing factors for transforms (linked list)

!type, public :: fax_type
type :: fax_type ! stores all factors needed for one transform length:
   integer (kind=int_kind) :: n = 0    ! transform length for set of factors
   integer (kind=int_kind) :: ifax(13) ! integer factors (factorization of n)
   real    (kind=dbl_kind), pointer :: tfax(:) ! trig factors for this n
   type (fax_type), pointer ::  next => null() ! points to next set of factors
end type fax_type

type (fax_type), target  :: fax            ! head of linked list of factors
type (fax_type), pointer :: xfax => null() ! points to factors for x transform
type (fax_type), pointer :: yfax => null() ! points to factors for y transform
save fax, xfax, yfax

! handy constants

real (kind=dbl_kind), parameter :: zero = 0.0_dbl_kind
real (kind=dbl_kind), parameter :: one  = 1.0_dbl_kind
real (kind=dbl_kind), parameter :: two  = 2.0_dbl_kind

! default spectral truncation and domain lengths

integer (kind=int_kind) :: mx_def=0   , my_def=0    ! spectral truncation
real    (kind=dbl_kind) :: xl_def=zero, yl_def=zero ! domain lengths

contains

!===============================================================================
subroutine setup_ss_ops( mx, my, xl, yl, scale_mode_m )
!-------------------------------------------------------------------------------
! Purpose:  
!   Sets default spectral truncation and domain lengths for ss_ops routines
!
! Arguments:
   integer (kind=int_kind), intent(in), optional :: mx, my ! last modes to keep
   real    (kind=dbl_kind), intent(in), optional :: xl, yl ! domain lengths
   logical, intent(in), optional :: scale_mode_m ! scale modes mx, my by two?
!
! Note:
!   The values  mx, my, xl,  and  yl  may be overridden by values specified 
!   in the argument lists of the ss_ops routines (see the individual routines).
!   For explanation of scale_mode_m, see comments for scale_last_mode above.
!-------------------------------------------------------------------------------

   if ( present(mx) ) mx_def = mx
   if ( present(my) ) my_def = my
   if ( present(xl) ) xl_def = xl
   if ( present(yl) ) yl_def = yl
   if ( present(scale_mode_m) ) scale_last_mode = scale_mode_m

!-------------------------------------------------------------------------------
end subroutine setup_ss_ops
!===============================================================================


!===============================================================================
subroutine tospec( f, mx, my )
!-------------------------------------------------------------------------------
! Purpose:  
!   Transforms a field to spectral space
!
! Arguments:
   real (kind=dbl_kind), intent(inout) :: f(-1:,-1:)       ! field to transform
   integer (kind=int_kind), intent(in), optional :: mx, my ! last modes to keep
!
! Note:
!   If mx or my is supplied, modes greater than mx or my are set to zero.
!-------------------------------------------------------------------------------

   nx = ubound( f, 1 ); ny = ubound( f, 2 );
   lx = 2*mx_def; if ( present(mx) ) lx = 2*mx; if (lx<=0) lx = nx;
   ly = 2*my_def; if ( present(my) ) ly = 2*my; if (ly<=0) ly = ny;
   if (lx>nx) stop "tospec:  spectral truncation is too large in x"
   if (ly>ny) stop "tospec:  spectral truncation is too large in y"
   call setup_fft

! forward  x  transform

   call fft99( f, work, xfax%tfax, xfax%ifax, 1, nx+2, nx, ny+2, -1 )
   if ( lx<nx ) f(lx+1:nx,:) = zero
   if ( scale_last_mode ) f(lx-1:lx,:) = f(lx-1:lx,:)*two

! forward  y  transform

   call fft99( f, work, yfax%tfax, yfax%ifax, nx+2, 1, ny, lx+2, -1 )
   if ( ly<ny ) f(-1:lx,ly+1:ny) = zero
   if ( scale_last_mode ) f(-1:lx,ly-1:ly) = f(-1:lx,ly-1:ly)*two

!-------------------------------------------------------------------------------
end subroutine tospec
!===============================================================================


!===============================================================================
subroutine tophys( f, mx, my )
!-------------------------------------------------------------------------------
! Purpose:  
!   Transforms a field to physical space
!
! Arguments:
   real (kind=dbl_kind), intent(inout) :: f(-1:,-1:)       ! field to transform
   integer (kind=int_kind), intent(in), optional :: mx, my ! last modes to keep
!
! Note:
!   If mx or my is supplied, modes greater than mx or my are set to zero.
!-------------------------------------------------------------------------------

   nx = ubound( f, 1 ); ny = ubound( f, 2 );
   lx = 2*mx_def; if ( present(mx) ) lx = 2*mx; if (lx<=0) lx = nx;
   ly = 2*my_def; if ( present(my) ) ly = 2*my; if (ly<=0) ly = ny;
   if (lx>nx) stop "tophys:  spectral truncation is too large in x"
   if (ly>ny) stop "tophys:  spectral truncation is too large in y"
   call setup_fft

! inverse  y  transform

   if ( ly<ny ) f(-1:lx,ly+1:ny) = zero
   if ( scale_last_mode ) f(-1:lx,ly-1:ly) = f(-1:lx,ly-1:ly)/two
   call fft99( f, work, yfax%tfax, yfax%ifax, nx+2, 1, ny, lx+2, +1 )

! inverse  x  transform

   if ( lx<nx ) f(lx+1:nx,-1:ny) = zero
   if ( scale_last_mode ) f(lx-1:lx,-1:ny) = f(lx-1:lx,-1:ny)/two
   call fft99( f, work, xfax%tfax, xfax%ifax, 1, nx+2, nx, ny+2, +1 )

!-------------------------------------------------------------------------------
end subroutine tophys
!===============================================================================


!===============================================================================
subroutine xderiv( f, g, mx, my, xl )
!-------------------------------------------------------------------------------
! Purpose:  
!   Computes the x-derivative of a field in spectral space
!
! Arguments:
   real (kind=dbl_kind), intent(in)  :: f(-1:,-1:) ! input  array
   real (kind=dbl_kind), intent(out) :: g(-1:,-1:) ! output array
   integer (kind=int_kind), intent(in), optional :: mx, my ! last modes to keep
   real (kind=dbl_kind),    intent(in), optional :: xl     ! x-domain length
!
! Note:
!    If same array is passed to f and g then the output overwrites the input.
!    For default values of mx, my, and xl see comments above and setup_ss_ops.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: jr, ji      ! indices
   real (kind=dbl_kind)    :: fact, factx ! factors

! sizes and scale factor

   nx = ubound( f, 1 ); ny = ubound( f, 2 );
   if (ubound( g, 1 )/=nx .or. ubound( g, 2 )/=ny ) &
      stop "xderiv:  g  is wrong shape"
   lx = 2*mx_def; if ( present(mx) ) lx = 2*mx; if (lx<=0) lx = nx;
   ly = 2*my_def; if ( present(my) ) ly = 2*my; if (ly<=0) ly = ny;
   if (lx>nx) stop "xderiv:  spectral truncation is too large in x"
   if (ly>ny) stop "xderiv:  spectral truncation is too large in y"
   fact = xl_def; if ( present(xl) ) fact = xl; 
   factx = two*acos( -one ); if ( fact>0 ) factx = factx/fact
   if (.not.allocated(work)) call setup_fft

! compute the derivative (in spectral space)

   do l=0,lx/2
      jr = 2*l-1
      ji = 2*l
      fact = l*factx
      work( 1,-1:ly) =          f(jr,-1:ly) ! store real part (in case  g=f)
         g(jr,-1:ly) = -fact*   f(ji,-1:ly)
         g(ji,-1:ly) =  fact*work( 1,-1:ly)
   end do

!-------------------------------------------------------------------------------
end subroutine xderiv
!===============================================================================


!===============================================================================
subroutine yderiv( f, g, mx, my, yl )
!-------------------------------------------------------------------------------
! Purpose:  
!   Computes the y-derivative of a field in spectral space
!
! Arguments:
   real (kind=dbl_kind), intent(in)  :: f(-1:,-1:) ! input  array
   real (kind=dbl_kind), intent(out) :: g(-1:,-1:) ! output array
   integer (kind=int_kind), intent(in), optional :: mx, my ! last modes to keep
   real (kind=dbl_kind),    intent(in), optional :: yl     ! y-domain length
!
! Note:
!    If same array is passed to f and g then the output overwrites the input.
!    For default values of mx, my, and yl see comments above and setup_ss_ops.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: kr, ki      ! indices
   real (kind=dbl_kind)    :: fact, facty ! factors

! sizes and scale factor

   nx = ubound( f, 1 ); ny = ubound( f, 2 );
   if (ubound( g, 1 )/=nx .or. ubound( g, 2 )/=ny ) &
      stop "yderiv:  g  is wrong shape"
   lx = 2*mx_def; if ( present(mx) ) lx = 2*mx; if (lx<=0) lx = nx;
   ly = 2*my_def; if ( present(my) ) ly = 2*my; if (ly<=0) ly = ny;
   if (lx>nx) stop "yderiv:  spectral truncation is too large in x"
   if (ly>ny) stop "yderiv:  spectral truncation is too large in y"
   fact = yl_def; if ( present(yl) ) fact = yl; 
   facty = two*acos( -one ); if ( fact>0 ) facty = facty/fact
   if (.not.allocated(work)) call setup_fft

! compute the derivative (in spectral space)

   do m=0,ly/2
      kr = 2*m-1
      ki = 2*m
      fact = m*facty
      work(-1:lx, 1) =          f(-1:lx,kr) ! store real part (in case  g=f)
         g(-1:lx,kr) = -fact*   f(-1:lx,ki)
         g(-1:lx,ki) =  fact*work(-1:lx, 1)
   end do

!-------------------------------------------------------------------------------
end subroutine yderiv
!===============================================================================


!===============================================================================
subroutine nlap( f, g, mlap, mx, my, xl, yl )
!-------------------------------------------------------------------------------
! Purpose:  
!   Computes the *negative* Laplacian of a field in spectral space:  
!   g = (-del^2)^mlap (f)  [default:  mlap = 1  for -Laplacian]
!
! Arguments:
   real (kind=dbl_kind), intent(in)  :: f(-1:,-1:) ! input field
   real (kind=dbl_kind), intent(out) :: g(-1:,-1:) ! array for output
   integer (kind=int_kind), intent(in), optional :: mlap   ! operator exponent
   integer (kind=int_kind), intent(in), optional :: mx, my ! last modes to keep
   real (kind=dbl_kind),    intent(in), optional :: xl, yl ! domain size
!
! Note:
!    If same array is passed to f and g then the output overwrites the input.
!    If mlap is not supplied then operator is -del^2 (i.e., default  mlap=1).
!    For default values of mx, my, xl, yl see comments above and setup_ss_ops.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: m_lap ! operator exponent (local copy)
   real    (kind=dbl_kind) :: factx, facty, factx2, facty2, fact ! factors

! determine which operator to use

   m_lap = 1; if (present(mlap)) m_lap = mlap;
   if (m_lap<0 ) stop "mlap:  negative exponent on Laplacian"
   if (m_lap==0) then
      g = f; return
   end if

! sizes and scale factors

   nx = ubound( f, 1 ); ny = ubound( f, 2 );
   if (ubound( g, 1 )/=nx .or. ubound( g, 2 )/=ny ) &
      stop "nlap:  g  is wrong shape"
   lx = 2*mx_def; if ( present(mx) ) lx = 2*mx; if (lx<=0) lx = nx;
   ly = 2*my_def; if ( present(my) ) ly = 2*my; if (ly<=0) ly = ny;
   if (lx>nx) stop "nlap:  spectral truncation is too large in x"
   if (ly>ny) stop "nlap:  spectral truncation is too large in y"
   fact = xl_def; if ( present(xl) ) fact = xl; 
   factx = two*acos( -one ); if ( fact>0 ) factx = factx/fact
   fact = yl_def; if ( present(yl) ) fact = yl; 
   facty = two*acos( -one ); if ( fact>0 ) facty = facty/fact
   factx2 = factx*factx

! compute the operator in spectral space

   do m=0,ly/2
      facty2 = m*m*facty*facty
      do l=0,lx/2
         fact = (l*l*factx2 + facty2)**m_lap
         g(2*l-1,2*m-1) = fact*f(2*l-1,2*m-1)
         g(2*l  ,2*m-1) = fact*f(2*l  ,2*m-1)
         g(2*l-1,2*m  ) = fact*f(2*l-1,2*m  )
         g(2*l  ,2*m  ) = fact*f(2*l  ,2*m  )
      end do
   end do

!-------------------------------------------------------------------------------
end subroutine nlap
!===============================================================================


!===============================================================================
subroutine mhsolve( g, f, lambda, fmean, mx, my, xl, yl )
!-------------------------------------------------------------------------------
! Purpose:  
!   Solves the modified Helmholtz equation in spectral space:
!
!      lambda*f - del^2(f) = g
!
! Argument:
   real (kind=dbl_kind), intent(in)  :: g(-1:,-1:) ! right-hand side
   real (kind=dbl_kind), intent(out) :: f(-1:,-1:) ! output array
   real (kind=dbl_kind), intent(in),  optional :: lambda    ! Helmholtz const.
   real (kind=dbl_kind), intent(in),  optional :: fmean     ! mean value of f
   integer (kind=int_kind), intent(in), optional :: mx, my  ! last modes to keep
   real (kind=dbl_kind),    intent(in), optional :: xl, yl  ! domain size
!
! Note:
!    If same array is passed to f and g then the output overwrites the input.
!    If not supplied, the default values are  lambda = 0  and  fmean = 0.
!    For default values of mx, my, xl, yl see comments above and setup_ss_ops.
!-------------------------------------------------------------------------------

! local variables

   real    (kind=dbl_kind) :: & ! various real factors
      factx, facty, factx2, facty2, fact, fbar, lam

! sizes and scale factors

   nx = ubound( g, 1 ); ny = ubound( g, 2 );
   if (ubound( f, 1 )/=nx .or. ubound( f, 2 )/=ny ) &
      stop "mhsolve:  f  is wrong shape"
   lx = 2*mx_def; if ( present(mx) ) lx = 2*mx; if (lx<=0) lx = nx;
   ly = 2*my_def; if ( present(my) ) ly = 2*my; if (ly<=0) ly = ny;
   if (lx>nx) stop "mhsolve:  spectral truncation is too large in x"
   if (ly>ny) stop "mhsolve:  spectral truncation is too large in y"
   fact = xl_def; if ( present(xl) ) fact = xl; 
   factx = two*acos( -one ); if ( fact>0 ) factx = factx/fact
   fact = yl_def; if ( present(yl) ) fact = yl; 
   facty = two*acos( -one ); if ( fact>0 ) facty = facty/fact
   factx2 = factx*factx
   lam  = zero; if (present( lambda )) lam  = lambda
   fbar = zero; if (present(  fmean )) fbar = fmean

! compute the solution in spectral space

   do m=0,ly/2
      facty2 = m*m*facty*facty
      do l=0,lx/2
         if ( lam==zero .and. l==0 .and. m==0 ) then
            f(2*l-1,2*m-1) = fbar
            f(2*l  ,2*m-1) = zero
            f(2*l-1,2*m  ) = zero
            f(2*l  ,2*m  ) = zero
         else
            fact = one/(lam + l*l*factx2 + facty2)
            f(2*l-1,2*m-1) = fact*g(2*l-1,2*m-1)
            f(2*l  ,2*m-1) = fact*g(2*l  ,2*m-1)
            f(2*l-1,2*m  ) = fact*g(2*l-1,2*m  )
            f(2*l  ,2*m  ) = fact*g(2*l  ,2*m  )
         end if
      end do
   end do

!-------------------------------------------------------------------------------
end subroutine mhsolve
!===============================================================================


!===============================================================================
subroutine setup_fft
!-------------------------------------------------------------------------------
! Purpose:  
!   Sets up factors and work space needed for fft99
!
! Note:
!   You may release all space allocated by calling unset_fft (for a clean exit)
!-------------------------------------------------------------------------------

   if (nx<1 .or. ny<1) stop "setup_fft:  bad transform size nx and/or ny"

! get the factors needed for the FFTs

   call get_fax( nx, xfax )
   call get_fax( ny, yfax )

! allocate work space for the transforms and derivatives

   if (nx>nx_max .or. ny>ny_max) then
      if (allocated(work)) deallocate(work)
      nx_max = nx; ny_max = ny
      allocate(work(-1:nx_max,-1:ny_max))
   end if

!-------------------------------------------------------------------------------
end subroutine setup_fft
!===============================================================================


!===============================================================================
subroutine get_fax( n, fx )
!-------------------------------------------------------------------------------
! Purpose:  
!   Finds or computes factors for fft99 with transform length  n
!
! Arguments:
   integer (kind=int_kind), intent(in) :: n  ! transform length
   type (fax_type), pointer            :: fx ! points to factors for n
!
! On input:
!    n    specifies transform length
!    fx   is the pointer to be assigned (can point to previous factors)
!
! On return:
!    fx   points to the factors for n
!
! Note:
!   The factors computed are stored in memory allocated as needed.
!   To release this memory when no longer needed, call reset_fax.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind), parameter :: trace = 0 ! controls printed trace:
   ! trace=0 for no trace, trace=1 to flag new n, trace=2 for full information

   if (trace>=2) print *, "get_fax:  called with n", n 

   if (.not. associated(fx) .or. fx%n/=n) fx => fax ! start at head of list
   do while (fx%n/=n) ! find this n in the linked list and point to it
      if (trace>=2) print *, "fx%n=", fx%n 
      if (fx%n==0 ) then ! end of list without finding n--compute the factors
         if (trace>=1) print *, &
            "get_fax:  computing factors for transform length ", n
         fx%n = n
         allocate (fx%tfax(3*n/2+1))
         call fftfax( n, fx%ifax, fx%tfax )
         if ( fx%ifax(1)==-99 ) stop 'get_fax:  bad transform length  n'
      else 
         if (.not.associated(fx%next)) allocate(fx%next); 
         fx => fx%next
      end if
   end do

!-------------------------------------------------------------------------------
end subroutine get_fax
!===============================================================================


!===============================================================================
subroutine unset_fft
!-------------------------------------------------------------------------------
! Purpose:  
!   Reclaims memory allocated for FFTs (for clean exit from execution)
!-------------------------------------------------------------------------------

! local variable

   type (fax_type), pointer :: fx ! points to factors in the list

! reclaim the space for factors, starting from the head of the linked list

   fx => fax; call free_fax( fx )

! reclaim the work space

   if (allocated(work)) deallocate(work)

contains

   !----------------------------------------------------------------------------
   recursive subroutine free_fax( fx_point )
   ! reclaims node fx_point of linked list
   type (fax_type), pointer :: fx_point ! points to node to reclaim

   if (associated(fx_point%next)) then
      call free_fax(fx_point%next)
      deallocate (fx_point%next)
   end if
   if (associated(fx_point%tfax)) deallocate (fx_point%tfax)
   end subroutine free_fax
   !----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
end subroutine unset_fft
!===============================================================================


!-------------------------------------------------------------------------------
end module ss_ops
!===============================================================================
