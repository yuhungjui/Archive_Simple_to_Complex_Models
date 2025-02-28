!===============================================================================
module pswm_setup
!-------------------------------------------------------------------------------
! Purpose: 
!    Routines for setting up PSWM:  Periodic Shallow Water Model
!-------------------------------------------------------------------------------

use kinds, only : int_kind, dbl_kind
use pswm_pars
use pswm_cons
use pswm_vars
use pswm_terms
use sitpack_interface
use sitpack
use ss_ops

implicit none   ! all undeclared module variables are untyped
private         ! all entry points and variables are local by default

!-------------------------------------------------------------------------------
!public variables and routines to be called from the main program:

public          &! routines to be called from main program (in this order):
   get_runid,   &! to get the run identifier
   setup_time,  &! to set up the time scheme
   setup_space, &! to set up the space discretization
   pswm_start,  &! to initialize the model solution
   pswm_cleanup  ! to cleanup after completing the model run
!-------------------------------------------------------------------------------

contains

!===============================================================================
subroutine get_runid
!-------------------------------------------------------------------------------
! Purpose:  
!    Gets the run identifier  runid  from command line or file runid.txt
!-------------------------------------------------------------------------------

   integer (kind=int_kind), parameter :: iunit = 10  ! input unit
   logical :: exists ! for answer from file inquiry

   inquire ( file='runid.txt', exist=exists )
   if (command_argument_count()>0) then ! runid is first command line argument
      call get_command_argument( 1, runid )
   else if ( exists ) then
      open  ( unit=iunit, file='runid.txt', status='old' )
      read  ( unit=iunit, fmt='(a)' ) runid
      close ( unit=iunit )
   else
      runid = "pswm_run"
   end if

!-------------------------------------------------------------------------------
end subroutine get_runid
!===============================================================================


!===============================================================================
subroutine setup_time( slabel )
!-------------------------------------------------------------------------------
! Purpose:  
!    Sets up the time differencing scheme, based on model parameter itd:
!
!  explicit schemes (itd<0)       implicit schemes (itd>0)
!  itd=-1:  FOR (forward)         itd=+1:  FOR/BACK (forward/backward)
!  itd=-2:  LFG (leapfrog)        itd=+2:  TRAP/LFG (trapezoidal/leapfrog)
!  itd=-3:  AB3 (Adams-Bashforth) itd=+3:  GAM2/AB3 (semi-implicit AB3)
!
! Arguments:
   character (len=*), intent(out) :: slabel  ! label for time scheme
!-------------------------------------------------------------------------------

! local variables for time scheme

   integer (kind=int_kind), parameter :: max_terms = 3  ! max # of terms
   character (len=8)  :: scheme(max_terms) = " "        ! time scheme (methods)
   integer (kind=int_kind) :: num_terms                 ! number of terms
   integer (kind=int_kind) :: k                         ! term index
   
! specify the time splitting:  indices for implicit/explicit/dissipation terms

   iterm_fast = 1; iterm_slow = 2;
   
! choose the methods for the implicit and explicit terms

   select case (itd)
      case (-1) ! forward (explicit)
         scheme(iterm_fast) = "FOR"
         scheme(iterm_slow) = "FOR"
         num_sol = 1; num_for = 1;
         write (*,*) "setup_time warning:  forward scheme is unstable"
      case (-2) ! leapfrog (explicit)
         scheme(iterm_fast) = "LEAP"
         scheme(iterm_slow) = "LEAP"
         num_sol = 2; num_for = 2;
      case (-3) ! AB3 (explicit)
         scheme(iterm_fast) = "AB3"
         scheme(iterm_slow) = "AB3"
         num_sol = 1; num_for = 5;
      case (+1) ! forward/backward (semi-implicit)
         scheme(iterm_fast) = "BACK"
         scheme(iterm_slow) = "FOR"
         write (*,*) "setup_time warning:  forward scheme is unstable"
         num_sol = 1; num_for = 1;
      case (+2) ! trapezoidal/leapfrog (semi-implicit)
         scheme(iterm_fast) = "TRAP2"
         scheme(iterm_slow) = "LEAP"
         num_sol = 2; num_for = 2;
      case (+3) ! GAM2/AB3 (semi-implicit)
         scheme(iterm_fast) = "GAM2"
         scheme(iterm_slow) = "AB3"
         num_sol = 1; num_for = 4;
      case default
         write (*,*) "Unknown time discretization:  itd = ", itd
         stop "Unknown time discretization specified by parameter itd"
   end select

! choose the method for the dissipation (and sponge) terms

   if ( cdiss>zero .or. tdiss>zero .or. csponge>zero .or. tsponge>zero ) then
      iterm_diss = 3;
      scheme(iterm_diss) = "FOR"
   else
      iterm_diss = 0;
   end if

! set the scheme label for output

   num_terms = max( iterm_fast, iterm_slow, iterm_diss )
   slabel = trim(scheme(1))
   do k=2,num_terms
      slabel = trim(slabel)//"/"//trim(scheme(k))
   end do

! initialize the time differencing scheme

!  call sitpack_trace( "sitpack_step" )
   call sitpack_setp( "num_sol", num_sol )
   call sitpack_setp( "num_for", num_for )
   time  = zero
   now   = 0
   nstep = 0
   call sitpack_init( time, now )
   call sitpack_set_scheme( scheme )

!-------------------------------------------------------------------------------
end subroutine setup_time
!===============================================================================


!===============================================================================
subroutine setup_space
!-------------------------------------------------------------------------------
! Purpose:  
!    Sets up the space discretization scheme, based on model parameter ieq
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: j, k  ! grid point indices

! some sanity checks

   if (mx==0) mx = (nx-1)/3; if (my==0) my = (ny-1)/3
   if ( 1>mx .or. 2*mx>nx ) stop "setup_space:  mx or nx out of range"
   if ( 1>my .or. 2*my>ny ) stop "setup_space:  my or ny out of range"
   if ( 3*mx>nx .or. 3*my>ny ) write (*,*) "setup_space warning:  ", &
      "mx or my too large--allows aliasing in quadratic terms"
   if ( xl<=zero .or. yl<=zero ) stop "setup_space:  xl or yl not positive"

! set up the transform grid: active area is j=0,...,nx-1 and k=0,...,ny-1
! (note:  f(x(nx))=f(x(0)) by periodicity but leave x(nx) alone for output)

   allocate (x(-1:nx)); x=x0+(/(j,j=-1,nx)/)*(xl/nx); x(-1)=x(nx-1); !x(nx)=x(0)
   allocate (y(-1:ny)); y=y0+(/(k,k=-1,ny)/)*(yl/ny); y(-1)=y(ny-1); !y(ny)=y(0)

! set up the arrays for the solution

   allocate (u(-1:nx,-1:ny,num_sol)); allocate (v(-1:nx,-1:ny,num_sol))
   allocate (p(-1:nx,-1:ny,num_sol)); allocate (d(-1:nx,-1:ny,num_sol))
   allocate (z(-1:nx,-1:ny,num_sol)); allocate (q(-1:nx,-1:ny,num_sol))

! set up the arrays for the forcing

   select case (abs(ieq))
   case (1) ! momentum form
      allocate (uf(-1:nx,-1:ny,num_for)); allocate (vf(-1:nx,-1:ny,num_for));
   case (2) ! vorticity-divergence form
      allocate (zf(-1:nx,-1:ny,num_for)); allocate (df(-1:nx,-1:ny,num_for));
   case (3) ! potential vorticity-divergence form
      allocate (qf(-1:nx,-1:ny,num_for)); allocate (df(-1:nx,-1:ny,num_for));
   case default
      write (*,*) "unrecognized value of  ieq = ", ieq 
      stop        "unrecognized value of  ieq"
   end select
   allocate (pf(-1:nx,-1:ny,num_for));
   
! initialize stuff for computing the model terms

   call setup_terms

!-------------------------------------------------------------------------------
end subroutine setup_space
!===============================================================================


!===============================================================================
subroutine pswm_start
!-------------------------------------------------------------------------------
! Purpose:  
!    Initializes the model variables
!-------------------------------------------------------------------------------

   call put_model_solution( time, now )

!-------------------------------------------------------------------------------
end subroutine pswm_start
!===============================================================================


!===============================================================================
subroutine pswm_cleanup
!-------------------------------------------------------------------------------
! Purpose:  
!    Cleans up from model run:  release allocated space
!-------------------------------------------------------------------------------

   logical, parameter :: trace = .false. ! trace this routine

   if (trace) print *, "pswm_cleanup:  calling unset_fft"
   call unset_fft

   if (trace) print *, "pswm_cleanup:  calling unset_terms"
   call unset_terms

   if (trace) print *, "pswm_cleanup:  deallocating forcing"
   if (allocated(pf)) deallocate(pf)
   if (allocated(df)) deallocate(df)
   if (allocated(qf)) deallocate(qf)
   if (allocated(zf)) deallocate(zf)
   if (allocated(vf)) deallocate(vf)
   if (allocated(uf)) deallocate(uf)

   if (trace) print *, "pswm_cleanup:  deallocating solution"
   if (trace) print *, "pswm_cleanup:  deallocating q"
   if (allocated(q)) deallocate(q)
   if (trace) print *, "pswm_cleanup:  deallocating z"
   if (allocated(z)) deallocate(z)
   if (trace) print *, "pswm_cleanup:  deallocating d"
   if (allocated(d)) deallocate(d)
   if (trace) print *, "pswm_cleanup:  deallocating p"
   if (allocated(p)) deallocate(p)
   if (trace) print *, "pswm_cleanup:  deallocating v"
   if (allocated(v)) deallocate(v)
   if (trace) print *, "pswm_cleanup:  deallocating u"
   if (allocated(u)) deallocate(u)

   if (trace) print *, "pswm_cleanup:  deallocating coordinates"
   if (allocated(x)) deallocate(x); if (allocated( y)) deallocate( y)

!-------------------------------------------------------------------------------
end subroutine pswm_cleanup
!===============================================================================


!-------------------------------------------------------------------------------
end module pswm_setup
!===============================================================================
