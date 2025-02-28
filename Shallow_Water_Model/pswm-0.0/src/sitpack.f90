!===============================================================================
module sitpack
!-------------------------------------------------------------------------------
! Purpose: 
!    Model-independent routines for semi-implicit time integration.
!
! Overview: 
!    This module contains routines to control semi-implicit time integration.
!    You provide the code to compute the terms in the model equations, and
!    sitpack calls them as needed to step the model forward in time.  Your
!    code does not have to deal with time stepping at all; sitpack takes care
!    of all the bookkeeping.  Many standard methods are provided, and you
!    can also define your own custom semi-implicit scheme.
!    
!    Sitpack can control any model consisting of prognostic equations
!
!       dv/dt = F_1(v,w,t) + F_2(v,w,t) + ... + F_n(v,w,t)                 (1)
!
!    and diagnostic equations
!
!       C(v,w,t) = 0 [solvable for w in terms of v and t]                  (2)
!
!    The prognostic variables  v  and diagnostic variables  w  are vectors
!    (e.g., collections of values at gridpoints or spectral coefficients).
!    The terms  F_k  are to be treated by different methods for  k=1, ..., n.
!    Typically, one of these methods will be implicit and the others explicit.
!    You must decide which terms to include in which forcing term  F_k,  
!    and which time integration method to use for each  F_k.  Here we will
!    use the word "method" to denote a time integration method for a single
!    term, and the word "scheme" to denote a combination of such methods.
!
!    Any semi-implicit scheme generates an implicit problem
!
!       v - (tau_1*F_1(v,w,t) + ... + tau_n*F_n(v,w,t) = G,  C(v,w,t) = 0  (3)
!
!    to be solved at each time step.  Here  tau_k  is a fraction (or multiple)
!    of the time step (tau_k=0  iff method  k  is explicit),  v  and  w  are 
!    evaluated at time  t,  and  G  is a combination of values of  v  and  F_k
!    at various times.  You will supply a routine to solve this problem.
!
! Use:
!    You must provide storage for multiple copies of the solution  v
!    (indexed by integers  s=1,2,...,num_sol) and the forcing for  v  
!    (indexed by integers  f=1,2,...,num_for).  The diagnostic variables  w
!    may also be stored in the solution copies (see below for details).
!    How you actually store the variables is irrelevant, as long as you can
!    associate each copy of solution or forcing with an index  s  or  f.
!    You will probably want to access these variables through a module;
!    examples are provided in the test programs for sitpack.
!
!    The numbers of copies (num_sol and num_for) needed depend on the 
!    semi-implicit scheme; the most needed should be  num_sol = m+1  and 
!    num_for = m*n + 1,  where  m  is the number of steps in the scheme 
!    you are using.  If you set the internal parameters  num_sol  and  num_for
!    (via sitpack_setp), then  sitpack  checks each time a new copy is needed
!    to verify that it exists (otherwise, you can check this yourself).
!    If you want to find the minimum storage needed for a specific scheme, 
!    run the model once with lots of storage and note which copies are used.
!    Alternatively, you could allocate new copies dynamically during the run.
!
!    Sitpack keeps track of model time for you; if you want to keep track of
!    time separately, it would make sense to check that the two times agree.
!
!    You must also provide the following routines for these specified tasks:
!
!    *  subroutine compute_model_terms( k, s, t, f, c )
!       integer (kind=int_kind), intent(in) :: k ! index of term
!       integer (kind=int_kind), intent(in) :: s ! index of solution copy
!       real    (kind=dbl_kind), intent(in) :: t ! model time
!       integer (kind=int_kind), intent(in) :: f ! index of forcing copy
!       real    (kind=dbl_kind), intent(in) :: c ! forcing coefficient
!       ACTION: compute the terms  F_k  in (1) based on solution copy  s
!               at time  t  and (depending on the value of  c):
!               if c==0:  PUT the result     F_k  into forcing copy  f
!               if c/=0:  ADD the product  c*F_k    to forcing copy  f
!
!    *  subroutine compute_lincomb_sol( ns, s, cs, fout )
!       integer (kind=int_kind), intent(in) :: ns     ! # of solution copies
!       integer (kind=int_kind), intent(in) :: s(ns)  ! corresponding indices
!       real    (kind=dbl_kind), intent(in) :: cs(ns) ! corresponding constants
!       integer (kind=int_kind), intent(in) :: fout   ! index of output forcing
!       ACTION: compute the linear combination
!               forcing(fout) = sum(i=1:ns) cs(i)*solution(s(i))
!               (where "solution" here is only the prognostic variables  v)
!
!    *  subroutine compute_lincomb_for( nf, f, cf, fout )
!       integer (kind=int_kind), intent(in) :: nf     ! # of forcing copies
!       integer (kind=int_kind), intent(in) :: f(nf)  ! corresponding indices
!       real    (kind=dbl_kind), intent(in) :: cf(nf) ! corresponding constants
!       integer (kind=int_kind), intent(in) :: fout   ! index of output forcing
!       ACTION: compute the linear combination
!               forcing(fout) = forcing(fout) + sum(i=1:nf) cf(i)*forcing(f(i))
!               (note:  must ADD the linear combination to current contents)
!
!    *  subroutine solve_implicit_eqns( tau, t, f, s )
!       real    (kind=dbl_kind), intent(in) :: tau(n) ! scaled time step for (3)
!       real    (kind=dbl_kind), intent(in) :: t      ! model time
!       integer (kind=int_kind), intent(in) :: f      ! index of forcing copy
!       integer (kind=int_kind), intent(in) :: s      ! index of solution copy
!       ACTION: given the values  tau(1), ..., tau(n),  solve the implicit 
!               problem (3) at time  t  based on the forcing copy  f  and 
!               put the resulting solution into solution copy  s
!
!    *  subroutine put_model_solution( t, s )
!       real    (kind=dbl_kind), intent(in) :: t  ! initial model time
!       integer (kind=int_kind), intent(in) :: s  ! index of solution copy
!       ACTION: put the model solution at time  t  into solution copy  s
!       NOTE:   called only if  direct_start = .true. (see details below)
!               otherwise, this can simply be a dummy routine
!
!    There are two possibilities for dealing with diagnostic variables  w:
!    (a) compute them along with  v  in routine  solve_implicit_eqns 
!        and store in the same solution copy
!    (b) recompute them from  v  in routine  compute_model_terms 
!        each time they are needed
!    If  w  is involved in solving (3) then (a) would be more natural.
!    Also, (a) would probably be faster, but (b) might take less storage.
!
!    The routines listed above (and the storage they access) must be packaged
!    into the module "sitpack_interface" which is used by this module.
!
!    To run the model using  sitpack  you must do the following:
!
!         [call  sitpack_setp        to set internal parameters if needed]
!
!       1. call  sitpack_init        to initialize sitpack, AND (if desired):
!          call  sitpack_add_method  to define additional methods
!
!       2. call  sitpack_set_scheme  to select a method for each term, OR:
!          call  sitpack_set_custom  to define a custom time scheme directly
!
!       3. Initialize your model variables (you might call  put_model_solution)
!          [Note:  steps 2 and 3 may be reversed if desired.]
!
!       4. call  sitpack_step  (repeatedly) to step the model forward in time
!
!         [call  sitpack_restart     to restart the scheme if/when desired]
!
!    For specific details on usage see the comments in the individual routines.
!    Several other user-callable routines are also available for tracing,
!    debugging, etc., as listed in the "public" statement below.
!
!    Note that a multistep method with  m>1  steps requires some way to compute 
!    the first  m-1  values, referred to as the starting values.  Usually you 
!    will want to do this by initializing your m-step method by a sequence of 
!    methods using 1, 2, ..., m-1 steps.  If you use  sitpack_set_scheme,
!    default starting methods are used (you can change the defaults by calling
!    sitpack_add_method); if you use  sitpack_set_custom,  you specify these
!    directly, along with the "main" m-step scheme to be used for subsequent 
!    steps.  Alternatively, you may choose to set these starting values 
!    directly (e.g., from a known analytical solution or from a previous run).  
!    To do so, set the internal parameter  direct_start  to .true. (see below) 
!    and supply the subroutine put_model_solution (see above).
!
!    To allow restarting a model run without restarting the time scheme:
!
!         call  sitpack_write_dump  to write internal data to a restart file
!         call  sitpack_read_dump   to read this data from a restart file
!
!    (see further details in the comments of these routines).
!
! Internal parameters: 
!    The following parameters may be set by calling sitpack_setp:
!
!    name         type      default  meaning
!    ----         ----      -------  -------
!    num_sol      int_kind     0     number of solution copies available
!                                    (if zero then ignored and not checked)
!
!    num_for      int_kind     0     number of forcing copies available
!                                    (if zero then ignored and not checked)
!
!    iunit        int_kind     6     unit number for all output (6=stdout)
!                                    (messages, errors, tracing, etc.)
!                                    If  unit<0  then output is to abs(iunit)
!                                    and execution halts on errors or warnings.
!
!    direct_start logical   .false.  call put_model_solution for starting values
!
!    The generic interface sitpack_setp handles all types.  For example:  
!       call sitpack_setp( "iunit", 20_int_kind )
!       call sitpack_setp( "direct_start", .true. )
!
! Modules used: 
!    sitpack_interface  see description above

!    kinds              to define the following default kinds:
!                       int_kind  for default integers
!                       dbl_kind  for default reals
!
! Method: 
!    Sitpack uses a general combined linear multistep methods as detailed in:
!    Fulton, S. R., "Semi-implicit time differencing" (manuscript).
!
! Notes: 
!    This version is constructed to minimize first the number of times model
!    terms must be computed and second the number of terms which must be stored.
!    It may be possible (in a future version) to keep the first requirement 
!    but loosen the second in order to reduce the data movement.
!
! Author: 
!    Scott R. Fulton (fulton@clarkson.edu)
!
! Brief history: 
!    This package was conceived in the late 1970s (in Fortran 77), but the
!    original version was awkward and saw limited use.  A deeper understanding 
!    of combined linear multistep methods, the need to implement semi-implicit 
!    methods in a climate model, and the advent of new tools in Fortran 90
!    combined to prompt the creation of this new code from scratch in 2002.
!
! Future plans: 
!    Redo internal storage to allocate space as needed for the time scheme.
!    Implement "bootstrap start" with fractional time steps for accuracy.
!    Allow for Runge-Kutta methods and other combinations.
!    Permit groups of equations to be stepped successively.
!
! Finally: 
!    If you use a relaxation method to solve the implicit equations, you can
!    refer to the overall solution method as "sitpack and relax" (sorry!).
!-------------------------------------------------------------------------------

use kinds, only : int_kind, dbl_kind
!Note:  if not using module kinds, comment this out and uncomment this line:
!integer, parameter :: int_kind  = kind(1), dbl_kind  = selected_real_kind(13)

use sitpack_interface, only :  &! user-supplied routines (see above)
   compute_model_terms,    &
   compute_lincomb_sol,    &
   compute_lincomb_for,    &
   solve_implicit_eqns,    &
   put_model_solution

implicit none
private         ! make everything in this module invisible by default
save

public                  &! user-callable routines:
   sitpack_setp,        &! generic interface to set internal parameters
   sitpack_init,        &! initialize sitpack
   sitpack_add_method,  &! add a new method to those defined
   sitpack_set_scheme,  &! select the method to use for each term
   sitpack_set_custom,  &! define a custom time scheme directly
   sitpack_step,        &! do one time step
   sitpack_setp_int,    &! set an internal parameter (type int_kind)
   sitpack_setp_dbl,    &! set an internal parameter (type dbl_kind)
   sitpack_setp_log,    &! set an internal parameter (type logical)
   sitpack_trace,       &! turn on tracing options for these routines
   sitpack_details,     &! returns some internal details (for testing/debugging)
   sitpack_write_dump,  &! writes internal data to a restart file
   sitpack_read_dump     ! reads  internal data from a restart file

!===============================================================================
interface sitpack_setp
!-------------------------------------------------------------------------------
! Purpose: 
!    Generic interface to set internal parameters for  sitpack  (any type)
!
! Usage:
!       call sitpack_setp( param_name, value, iecho )
!
!    where:
!       param_name is the name of the parameter to set (a character string)
!       value      is the value to assign (corresponding type)
!       iecho      [optional] is the unit number to echo this assignment to
!                  (echo to stdout if  iecho  is supplied but not positive)
!
! Notes:
!    See comments above for list of available internal parameters and types.
!    See individual type routines listed below for details.
!-------------------------------------------------------------------------------
   module procedure sitpack_setp_int, sitpack_setp_dbl, sitpack_setp_log
end interface
!===============================================================================

! data for standard methods (with room for the user to add more)

integer, parameter :: max_method =  20 ! maximum number of standard methods
integer, parameter :: max_steps  =  10 ! max number of time steps  (m)
integer, parameter :: len_name   =   8 ! length of method name
type method_type
   character (len=len_name) :: name   ! method name
   integer (kind=int_kind)  :: m      ! number of steps
   real    (kind=dbl_kind)  :: coef_sol(0:max_steps) ! for solution
   real    (kind=dbl_kind)  :: coef_for(0:max_steps) ! for forcing
   integer (kind=int_kind)  :: istart ! index of starting method
end type method_type
type (method_type) :: method(max_method)      ! data for standard methods
integer (kind=int_kind) :: num_method =  0 ! number of standard methods defined

! parameters to specify maximum sizes needed (for storing coefficients)

integer, parameter :: max_terms = 4 ! max number of terms  (n)

! numbers specifying the current semi-implicit scheme

integer (kind=int_kind) :: num_steps = 0 ! number of time steps (m)
integer (kind=int_kind) :: num_terms = 0 ! number of terms (n)

! coefficients for the semi-implicit scheme:
! first  index  j=0,...,max_steps indexes time level  new-j 
! middle index  k=1,...,max_terms indexes terms
! last   index  l=1,...,max_steps indexes initialization steps
!               (l=1 for first step, l=2 for second, ..., l=m for remaining)

real (kind=dbl_kind) ::                        &!  coefficients for:
   c_sol(0:max_steps,          max_steps) = 0, &!  solution
   c_for(0:max_steps,max_terms,max_steps) = 0   !  forcing

character (len=len_name) :: methods(max_steps,max_terms) ! method labels
integer (kind=int_kind)  :: order(max_steps,max_terms)   ! accuracy

logical store_sol(max_steps,          max_steps) ! store the solution?
logical store_for(max_steps,max_terms,max_steps) ! store the forcing?

! "pointers" for the solution and forcing.  For each  j=0, ..., m:
!    p(j)     is the index of the solution/forcing copy at time level  new-j
!    p(j)=0   if that solution/forcing need not be stored
!    time(j)  is the corresponding model time

integer (kind=int_kind) ::              &!  "pointers" for:
   p_sol(0:max_steps)           = 0,    &!  solution
   p_for(0:max_steps,max_terms) = 0      !  forcing
real    (kind=dbl_kind) :: time(0:max_steps) ! corresponding model time

! variables for tracking the time steps, etc.
integer (kind=int_kind) :: nstep = 0   ! # of steps since start
integer (kind=int_kind) :: mstep = 0   ! # of steps since restart
real    (kind=dbl_kind) :: last_dt = 0 ! size of previous time step

! internal parameters (integer)
integer (kind=int_kind) :: &
   iunit   = 6, &! unit number for output (messages, errors, tracing, etc.)
   num_sol = 0, &! number of solution copies available (if zero, not checked)
   num_for = 0   ! number of forcing  copies available (if zero, not checked)

! internal parameters (logical)
logical :: direct_start = .false. ! starting values from put_model_solution?

! routine names to trace and corresponding trace options
integer, parameter      :: mtrace = 20 ! max number of routines to trace
integer (kind=int_kind) :: ntrace =  0 ! number of routines to trace
character (len=64)      :: trace_name(mtrace)   = " " ! routines to trace
integer (kind=int_kind) :: trace_option(mtrace) =  -1 ! trace options

! other variables, constants, etc.

character (len=128) :: message ! to hold trace and error messages
real    (kind=dbl_kind), parameter :: zero = 0.0_dbl_kind ! as it says...
real    (kind=dbl_kind), parameter :: one  = 1.0_dbl_kind ! as it says...

!-------------------------------------------------------------------------------

contains

!===============================================================================
subroutine sitpack_setp_int( param_name, value, iecho )
!-------------------------------------------------------------------------------
! Purpose: 
!    Sets internal parameters for  sitpack  (integer version)
!
! Arguments:
!
   character (len=*), intent(in) ::  param_name            ! name of parameter
   integer (kind=int_kind), intent(in) :: value            ! integer value
   integer (kind=int_kind), intent(in), optional :: iecho  ! unit number 
!
! Notes:
!    To echo the value set, supply the argument  iecho  (a positive number
!    specifies the unit number to echo to, otherwise echos to stdout).
!    See comments above for list of available internal parameters.
!    You can call the generic interface  sitpack_setp for any type.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: lecho ! unit number for echo

! echo value

   if ( present( iecho ) ) then
      lecho = iecho
      if ( lecho<=0 ) lecho = 6
      write (lecho,*) "sitpack_setp:  set parameter ", trim(param_name), &
                      " to ", value
   end if

! set the value

   select case (trim(param_name))

      case ("num_sol")
         num_sol = value
      case ("num_for")
         num_for = value
      case default
         call sitpack_message( "sitpack_setp", &
            "unknown parameter name = "//trim(param_name), 1 )

   end select 

!-------------------------------------------------------------------------------
end subroutine sitpack_setp_int
!===============================================================================


!===============================================================================
subroutine sitpack_setp_dbl( param_name, value, iecho )
!-------------------------------------------------------------------------------
! Purpose: 
!    Sets internal parameters for  sitpack  (real version)
!
! Arguments:
!
   character (len=*), intent(in) ::  param_name            ! name of parameter
   real    (kind=dbl_kind), intent(in) :: value            ! real value
   integer (kind=int_kind), intent(in), optional :: iecho  ! unit number 
!
! Notes:
!    To echo the value set, supply the argument  iecho  (a positive number
!    specifies the unit number to echo to, otherwise echos to stdout).
!    See comments above for list of available internal parameters.
!    You can call the generic interface  sitpack_setp for any type.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: lecho ! unit number for echo

! echo value

   if ( present( iecho ) ) then
      lecho = iecho
      if ( lecho<=0 ) lecho = 6
      write (lecho,*) "sitpack_setp:  set parameter ", trim(param_name), &
                      " to ", value
   end if

! set the value

   select case (trim(param_name))

      case default
         call sitpack_message( "sitpack_setp", &
            "unknown parameter name = "//trim(param_name), 1 )

   end select 

!-------------------------------------------------------------------------------
end subroutine sitpack_setp_dbl
!===============================================================================


!===============================================================================
subroutine sitpack_setp_log( param_name, value, iecho )
!-------------------------------------------------------------------------------
! Purpose: 
!    Sets internal parameters for  sitpack  (logical version)
!
! Arguments:
!
   character (len=*), intent(in) ::  param_name            ! name of parameter
   logical, intent(in) :: value                            ! logical value
   integer (kind=int_kind), intent(in), optional :: iecho  ! unit number 
!
! Notes:
!    To echo the value set, supply the argument  iecho  (a positive number
!    specifies the unit number to echo to, otherwise echos to stdout).
!    See comments above for list of available internal parameters.
!    You can call the generic interface  sitpack_setp for any type.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: lecho ! unit number for echo

! echo value

   if ( present( iecho ) ) then
      lecho = iecho
      if ( lecho<=0 ) lecho = 6
      write (lecho,*) "sitpack_setp:  set parameter ", trim(param_name), &
                      " to ", value
   end if

! set the value

   select case (trim(param_name))

      case ("direct_start")
         direct_start = value
      case default
         call sitpack_message( "sitpack_setp", &
            "unknown parameter name = "//trim(param_name), 1 )

   end select 

!-------------------------------------------------------------------------------
end subroutine sitpack_setp_log
!===============================================================================


!===============================================================================
subroutine sitpack_init( t, s )
!-------------------------------------------------------------------------------
! Purpose: 
!    Initializes  sitpack
!
! Calling order: 
!    First call  sitpack_init  to initialize sitpack.
!    Then call  sitpack_set_scheme  to select the method for each term
!    (or  call  sitpack_set_custom  to define a custom scheme directly).
!    Then initialize your model variables (in the solution copy indexed 
!    by  s)  and call  sitpack_step  (repeatedly) to run the model.
!
! Arguments: 
   real    (kind=dbl_kind), intent(in)    ::  t ! initial model time
   integer (kind=int_kind), intent(inout) ::  s ! index of solution copy
!
! Notes: 
!    The argument  s  specifies the index of for initial copy of the solution.
!    If  s=0  on input, that index is chosen internally and returned in  s.
!    If  s>0  on input, that index is used and not changed (which may be 
!    useful in restarting a model).  Either way, after calling  sitpack_init
!    you must put the initial values into the solution copy indexed by  s
!    before calling  sitpack_step  to step the solution forward in time, 
!
!    This routine defines all standard time differencing methods listed in
!    sitpack_add_method and makes them available for use in sitpack_set_scheme.
!    If you want to define additional methods (or redefine the standard ones),
!    call  sitpack_add_method  yourself (after calling  sitpack_init).
!-------------------------------------------------------------------------------

! define the standard methods (if not already defined)

   if ( get_method( "FOR"   )==0 )  call sitpack_add_method( "FOR" )
   if ( get_method( "AB2"   )==0 )  call sitpack_add_method( "AB2" )
   if ( get_method( "AB3"   )==0 )  call sitpack_add_method( "AB3" )
   if ( get_method( "AB4"   )==0 )  call sitpack_add_method( "AB4" )
   if ( get_method( "LEAP"  )==0 )  call sitpack_add_method( "LEAP" )
   if ( get_method( "BACK"  )==0 )  call sitpack_add_method( "BACK" )
   if ( get_method( "TRAP"  )==0 )  call sitpack_add_method( "TRAP" )
   if ( get_method( "TRAP2" )==0 )  call sitpack_add_method( "TRAP2" )
   if ( get_method( "AM3"   )==0 )  call sitpack_add_method( "AM3" )
   if ( get_method( "GAM1"  )==0 )  call sitpack_add_method( "GAM1" )
   if ( get_method( "GAM2"  )==0 )  call sitpack_add_method( "GAM2" )

! initialize time, step counts, and the current solution pointer

   time  = zero
   time(0) = t
   nstep = 0
   mstep = 0
   last_dt = zero
   p_sol = 0
   p_for = 0
   if ( s>0 ) then
      if ( num_sol>0 .and. s>num_sol ) then
         write (message,*) "input  s = ", s ," is too large"
         call sitpack_message( "sitpack_init", message, 2 )
      end if
      p_sol(0) = s
   else
      p_sol(0) = get_new_pointer( "solution", p_sol, num_sol )
      s = p_sol(0)
   end if

!-------------------------------------------------------------------------------
end subroutine sitpack_init
!===============================================================================


!===============================================================================
subroutine sitpack_add_method( name, coef_sol, coef_for, start, theta )
!-------------------------------------------------------------------------------
! Purpose: 
!    Defines a time differencing method
!
! Required Argument:
   character (len=*), intent(in) :: name ! name of method (up to 8 characters)
!
! Optional Arguments:
   real    (kind=dbl_kind), intent(in), optional :: coef_sol(0:) ! for solution
   real    (kind=dbl_kind), intent(in), optional :: coef_for(0:) ! for forcing
   character (len=*), intent(in), optional :: start  ! name of starting method
   real    (kind=dbl_kind), intent(in), optional :: theta ! coef of new forcing
!
! Details:
!    This routine defines either a STANDARD or a CUSTOM method as follows.
!
!    1. STANDARD METHODS:  
!
!      name   description                   type steps order start theta
!      ----   -----------                   ---- ----- ----- ----- -----
!      FOR    Forward (Euler)               exp    1     1
!      AB2    Second-order Adams-Bashforth  exp    2     2    FOR
!      AB3    Third-order  Adams-Bashforth  exp    3     3    AB2
!      AB4    Fourth-order Adams-Bashforth  exp    4     4    AB3
!      LEAP   Leapfrog (midpoint)           exp    2     2    FOR
!      BACK   Euler-backward                imp    1     1
!      TRAP   Trapezoidal                   imp    1     2
!      TRAP2  Two-step trapezoidal          imp    2     2*   BACK  0.5
!      AM3    Third-order Adams-Moulton     imp    2     3    TRAP
!      GAM1   Generalized Adams-Moulton     imp    1     2*         0.5
!      GAM2   Generalized Adams-Moulton     imp    2*    2*   BACK  1.25
!
!    Here:
!       name    is the method name (also used in sitpack_set_scheme)
!       type    is explict (exp) or implicit (imp)
!       steps   is one less than the number of time levels used
!       order   is the order of accuracy of the method
!       theta   is the default value used (if not supplied in argument list)
!       (*means with default theta)
!
!    Notes:
!       Each of these methods can be defined by calling sitpack_add_method
!       with its name only, e.g., call sitpack_add_method( "AB3" ).
!       If you supply  coef_sol  and/or  coef_for  it will be ignored.
!       All of these methods are predefined by sitpack_init, but you can 
!       change some details by calling  sitpack_add_method  as follows.
!
!       Methods with more than one step need a starting method with fewer steps.
!       To get something other than the default, specify the starting method,
!       e.g., call sitpack_add_method( "AB3", start="FOR" ).
!
!       Some implicit methods are actually one-parameter families indexed by
!       the parameter  theta  (weight of the forcing at the new time level).
!       To get something other than the default, specify the value you want,
!       e.g., call sitpack_add_method( "GAM1", theta=0.6_dbl_kind ).
!
!       GAM1:  theta=1.0 and theta=0.5 (default) give the BACK and TRAP 
!       methods, respectively.  The order is 1 unless theta=0.5 (default), 
!       in which case the order is 2.
!
!       TRAP2:  same as GAM1 except that it is spread over two steps so it
!       can be combined with the leapfrog method to yield the traditional 
!       semi-implicit scheme.
!
!       GAM2:  theta=0.5 and theta=5/12 give  TRAP and AM3, respectively.  
!       The order is 2 unless theta=5/12, in which case the order is 3.  
!       The scheme has 2 steps unless theta=0.5, in which case it has 1 step.  
!       The default value theta=1.25 gives a second-order two-step implicit 
!       method which is A-stable and damps all high frequencies.
!
!    2. CUSTOM METHODS:  
!
!    To define a new method not on the above list, you must supply the
!    following input arguments:
!
!       name, coeff_sol, coeff_for, and (if a multistep method) start
!
!    If you supply the input argument  theta,  it will be ignored.
!
!    The time integration method is defined for the differential equation
!    dv/dt = F(v)  by the difference equation:
!
!       (1/dt)*sum( coef_sol(j)*v(i-j), j=0,...,m ) = 
!              sum( coef_for(j)*F(i-j), j=0,...,m )
!
!    where  dt  is the time step and  v(i)  and  F(i)  are the solution and 
!    forcing, respectively, at time level  t(i) = t(0) + i*dt.  
!
!    In this routine, the number of steps  m  is taken as the largest integer
!    j  such that coef_sol(j) or coef_for(j) is nonzero.
!
!    If  m>1  (i.e., this is a multistep method) you must specify the method
!    to use to compute the starting values by supplying  start  (which must
!    identify a previously-defined method with less than  m  steps).
!
!    For examples, see the standard methods already defined in sitpack_init.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: i               ! index for loops
   integer (kind=int_kind) :: m               ! number of steps in method
   integer (kind=int_kind) :: imethod, istart ! indices of methods
   real    (kind=dbl_kind) :: cs(0:max_steps) ! local copy of coef_sol
   real    (kind=dbl_kind) :: cf(0:max_steps) ! local copy of coef_for
   real    (kind=dbl_kind) :: t, t0, t1, z    ! weights

! identify starting method (if specified)

   istart = 0
   if ( present( start ) ) then
      istart = get_method( start )
      if ( istart==0 ) then
         write (message,*) "specified starting method ", trim( start ), &
                           "is not defined"
         call sitpack_message( "sitpack_add_method", message, 1 )
         return
      end if
   end if

! set the coefficients

   cs = zero
   cf = zero
   z  = zero

   select case (trim(name))

      case ("FOR") ! Forward (Euler)
         cs(0:1) = (/ 1, -1 /)
         cf(0:1) = (/ 0,  1 /)

      case ("AB2") ! Second-order Adams-Bashforth
         cs(0:2) = (/ 1, -1, 0 /)
         cf(0:2) = (/ 0,  3,-1 /)*one/2
         if ( istart==0 ) istart  = get_method( "FOR" )

      case ("AB3") ! Third-order Adams-Bashforth
         cs(0:3) = (/ 1, -1,  0, 0 /)
         cf(0:3) = (/ 0, 23,-16, 5 /)*one/12
         istart  = get_method( "AB2" )

      case ("AB4") ! Fourth-order Adams-Bashforth
         cs(0:4) = (/ 1, -1,  0,  0,  0 /)
         cf(0:4) = (/ 0, 55,-59, 37, -9 /)*one/24
         if ( istart==0 ) istart  = get_method( "AB3" )

      case ("LEAP") ! Leapfrog (midpoint)
         cs(0:2) = (/ 1,  0, -1 /)*one/2
         cf(0:2) = (/ 0,  1,  0 /)
         if ( istart==0 ) istart  = get_method( "FOR" )

      case ("BACK") ! Backward (Euler-backward)
         cs(0:1) = (/ 1, -1 /)
         cf(0:1) = (/ 1,  0 /)

      case ("TRAP") ! Trapezoidal
         cs(0:1) = (/ 1, -1 /)
         cf(0:1) = (/ 1,  1 /)*one/2

      case ("TRAP2") ! Two-step trapezoidal
         t = 0.5
         if ( present( theta ) ) t = theta
         cs(0:2) = (/ 1,  0,  -1 /)*one/2
         cf(0:2) = (/ t,  z, 1-t /)
         if ( istart==0 ) istart  = get_method( "BACK" )

      case ("AM3") ! Third-order Adams-Moulton
         cs(0:2) = (/ 1, -1,  0 /)
         cf(0:2) = (/ 5,  8, -1 /)*one/12
         if ( istart==0 ) istart  = get_method( "TRAP" )

      case ("GAM1") ! Generalized Adams-Moulton (first order)
         t = 0.5
         if ( present( theta ) ) t = theta
         cs(0:1) = (/ 1,  -1 /)
         cf(0:1) = (/ t, 1-t /)

      case ("GAM2") ! Generalized Adams-Moulton (second order)
         t = 1.25
         if ( present( theta ) ) t = theta
         t0 = (3 - 4*t)/2
         t1 = t - 0.5
         if ( t==0.5_dbl_kind ) t1 = zero
         cs(0:2) = (/ 1, -1,  0 /)
         cf(0:2) = (/ t, t0, t1 /)
         if ( istart==0 ) istart  = get_method( "BACK" )

      case default ! user-defined method
         if ( present(coef_sol) .and. present(coef_for) ) then
            if ( (ubound(coef_sol,1)>max_steps .and. &
                  maxval(abs(coef_sol(max_steps+1:)))>0) .or. &
                 (ubound(coef_for,1)>max_steps .and. &
                  maxval(abs(coef_for(max_steps+1:)))>0) ) then
               write (message,*) "method ", trim(name), &
                           " has too many steps for current dimensions"
               call sitpack_message( "sitpack_add_method", message, 1 )
               return
            end if
            do i=0,min( max_steps, ubound(coef_sol,1) )
               cs(i) = coef_sol(i)
            end do
            do i=0,min( max_steps, ubound(coef_for,1) )
               cf(i) = coef_for(i)
            end do
         else
            write (message,*) "method ", trim(name), " is new--", &
                              "must specify  coef_sol  and  coef_for"
            call sitpack_message( "sitpack_add_method", message, 1 )
            return
         end if

   end select 

! determine the number of steps

   m = 0
   do i=0,max_steps
      if ( cs(i)/=zero .or. cf(i)/=zero) m = i
   end do
   if ( m == 0 ) then
      write (message,*) "method ", trim(name), " has no steps!"
      call sitpack_message( "sitpack_add_method", message, 1 )
      return
   end if

! identify the starting method (if any)

   if ( m>1 ) then
      if ( istart==0 .and. .not.present( start ) ) then
         write (message,*) "must specify  start  for a multistep method"
         call sitpack_message( "sitpack_add_method", message, 1 )
         return
      end if
      if ( method(istart)%m>=m ) then
         write (message,*) "starting method ", trim( start ), &
                           " has too many steps"
         call sitpack_message( "sitpack_add_method", message, 1 )
         return
      end if
   else
      if ( present( start ) ) then
         write (message,*) "method ", trim(name), " has one step--", &
                           "no starting method needed"
         call sitpack_message( "sitpack_add_method", message, 1 )
         return
      end if
   end if

! check to see whether the method is already defined

   imethod = get_method( name )
   if ( imethod==0 ) then
      if ( num_method>=max_method ) then
         write (message,*) "no room left to define method ", trim(name)
         call sitpack_message( "sitpack_add_method", message, 1 )
         return
      end if
      num_method = num_method + 1
      imethod = num_method
      !write (message,*) "defining method ", trim(name)
      !call sitpack_message( "sitpack_add_method", message )
   else
      write (message,*) "redefining method ", trim(name)
      call sitpack_message( "sitpack_add_method", message )
      if ( present( start ) ) then
         write (message,*) "specified starting method ", &
                           trim(method(istart)%name)
         call sitpack_message( "sitpack_add_method", message )
      end if
      if ( present( theta ) ) then
         write (message,*) "specified theta = ", theta
         call sitpack_message( "sitpack_add_method", message )
      end if
   end if

! store the values defining the method

   method(imethod)%name = name
   method(imethod)%m = m
   method(imethod)%coef_sol = cs
   method(imethod)%coef_for = cf
   method(imethod)%istart = istart

!-------------------------------------------------------------------------------
end subroutine sitpack_add_method
!===============================================================================


!===============================================================================
function get_method( name ) result (imethod)
!-------------------------------------------------------------------------------
! Purpose:
!    Looks up a method in the list of standard methods
!    (intended for use from within sitpack only--not a public routine)
!
! Argument:
   character (len=*), intent(in) :: name ! name of method to look for
!
! Result:
   integer (kind=int_kind) :: imethod ! index of method (0 if not found)
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: i ! loop index

   imethod = 0
   do i=1,num_method
      if ( name == method(i)%name ) then
         imethod = i
         return
      end if
   end do
   return

!-------------------------------------------------------------------------------
end function get_method
!===============================================================================


!===============================================================================
subroutine sitpack_set_scheme( scheme )
!-------------------------------------------------------------------------------
! Purpose: 
!    Select a semi-implicit time integration scheme (call sitpack_init first)
!
! Argument:
! 
   character (len=*), intent(in) :: scheme(:) 
!
!    scheme  is an array which specifies which method to use for each of the
!    forcing terms  k=1, ..., n  [see eq. (1) in the sitpack overview above].
!    Each value in  scheme  is the name of one of the standard methods 
!    defined in  sitpack_init (or defined by calling sitpack_add_method).
!    The number  n  of terms is taken as the index of the last non-blank value.
!
!    Examples:
!
!       scheme = (/ "TRAP2", "LEAP " /)  the usual trapezoidal/leapfrog scheme:
!           term 1:  TRAP2  trapezoidal (implicit)
!           term 2:  LEAP   leapfrog    (explicit)
!
!       scheme = (/ "GAM2", "AB3 ", "FOR " /) might be good for a climate model:
!           term 1:  GAM2 (implicit) [gravity wave terms]
!           term 2:  AB3  (explicit) [advective and Coriolis terms]
!           term 3:  FOR  (explicit) [dissipation, specified forcing, etc.]
!
!    Notes:
!       Any methods to be combined must have the same coefficients for the
!       solution (coef_sol in sitpack_init).  This will be checked internally.
!
!       Method names are case-sensitive and trailing blanks are ignored (but
!       may be needed in explicit array constructors as in the above examples).
!       Only the first 8 characters of the method name are significant.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: imethod ! index of method
   integer (kind=int_kind) :: k, l    ! loop indices
   real    (kind=dbl_kind) ::        &! copy of c_sol including each term
      c_tmp(0:max_steps,max_terms,max_steps) = 0 

! determine the number of terms and steps

   num_terms = size( scheme )
   do k=num_terms,1,-1
      if ( scheme(k)==" " ) num_terms = num_terms - 1
   end do
   if ( num_terms<=0 .or. num_terms>max_terms ) then
      write (message,*) "too many (or too few) terms:  num_terms = ", num_terms
      call sitpack_message( "sitpack_set_scheme", message, 1 )
      return
   end if
   num_steps = 0
   do k=1,num_terms
      imethod = get_method( scheme(k) )
      if ( imethod == 0 ) then
         write (message,*) "method for term k = ", k, " is undefined"
         call sitpack_message( "sitpack_set_scheme", message, 1 )
         return
      end if
      num_steps = max( num_steps, method(imethod)%m )
   end do

! initialize all coefficients to zero (i.e., not used)

   c_sol = 0
   c_for = 0
   methods = " "

! set up the method for term  k  (set starting methods recursively)

   c_tmp = 0
   do k=1,num_terms
      call set_std_method( num_steps, get_method( scheme(k)) )
   end do

! make sure the methods are compatible (solution coefficients must match)

   c_sol = 0
   do l=1,num_steps
      c_sol(:,l) = c_tmp(:,1,l)
      do k=2,num_terms
         if (maxval(abs(c_sol(:,l) - c_tmp(:,k,l))) > 5*epsilon(one)) then
            write (message,*) "scheme(1) = ", trim( scheme(1) )
            call sitpack_message( "sitpack_set_scheme", message, 1 )
            write (message,*) "scheme(k) = ", trim( scheme(k) )
            call sitpack_message( "sitpack_set_scheme", message, 1 )
            write (message,*) "methods not compatible for step  l = ", l
            call sitpack_message( "sitpack_set_scheme", message, 2 )
         end if
      end do ! loop on k
   end do ! loop on l

! determine what needs to be stored after each step

   call process_scheme( "sitpack_set_scheme" )

! trace output

   if (trace_routine( "sitpack_set_scheme" )>0) then
      call output_scheme
   end if

! quit if the scheme is inconsistent

   if (maxval(order(1:num_steps,1:num_terms))<=0) then
      call sitpack_message( "sitpack_set_scheme", "scheme is inconsistent", 2 )
   end if

! if the scheme was previously defined, reset everything to start fresh

   if (mstep>0) call sitpack_restart

! trace output

   if (trace_routine( "sitpack_set_scheme" )>0) then
      call sitpack_message( "sitpack_set_scheme", "set up scheme" )
   end if

contains ! internal (recursive) routine to set standard method parameters

   recursive subroutine set_std_method( l, imethod )
   integer (kind=int_kind), intent(in) :: l       ! index of time level
   integer (kind=int_kind), intent(in) :: imethod ! index of method

   integer (kind=int_kind) :: istart  ! method with which to start this one

   methods(l,k) = method(imethod)%name
   c_tmp(:,k,l) = method(imethod)%coef_sol
   c_for(:,k,l) = method(imethod)%coef_for
   if (l>1) then
      if ( l>method(imethod)%m ) then
         istart = imethod
      else
         istart = method(imethod)%istart
      end if
      call set_std_method( l-1, istart )
   end if
   end subroutine set_std_method

!-------------------------------------------------------------------------------
end subroutine sitpack_set_scheme
!===============================================================================


!===============================================================================
subroutine sitpack_set_custom( coef_sol, coef_for, names )
!-------------------------------------------------------------------------------
! Purpose: 
!    Define a custom semi-implicit scheme directly
!
! Required Arguments:
   real (kind=dbl_kind), intent(in) ::  &!
      coef_sol(0:,  :), &! coefficients for solution
      coef_for(0:,:,:)   ! coefficients for forcing
!
! Optional Argument:
   character (len=*), intent(in), optional :: names(:,:) ! method names
!                         (names(l,k) is the name for step l for term k)
!
! Notes:
!    The semi-implicit scheme is defined by the equation:
!
!       (1/dt)*sum( coef_sol(j,  l)*v(i-j),   j=0,...,m ) = 
!         sum( sum( coef_for(j,k,l)*F_k(i-j), j=0,...,m ), k=1,...,n )
!
!    where  v(i)  is the solution and  F_k(i)  is the forcing term  k  at
!    time level  t(i).  The schemes specified for l=1, ..., m-1 are used to 
!    compute the starting values for the m-step scheme l=m, which is used for
!    steps  m  and beyond (until next call to call sitpack_restart, if any).
!    Thus, scheme  l  may be at most an  l-step scheme (l=1, ..., m).
!
!    The values of  m  and  n  are taken as the largest values for which the
!    corresponding coefficients  coef_sol  and/or  coef_for  are nonzero.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: k, l ! indices (for terms and steps)
   integer (kind=int_kind) :: m, n ! short names for num_steps and num_terms

! initialize all coefficients to zero (i.e., not used)

   c_sol = 0
   c_for = 0
   methods = " "

! determine the number of terms

   n = 0
   do k=1,ubound( coef_for, dim=2 )
      if ( maxval(abs(coef_for(:,k,:)))/=zero ) n = k
   end do
   num_terms = n
   if ( num_terms<=0 .or. num_terms>max_terms ) then
      write (message,*) "too many (or too few) terms:  num_terms = ", num_terms
      call sitpack_message( "sitpack_set_custom", message, 1 )
      return
   end if

! determine the number of steps

   m = 0
   do l=1,min(ubound(coef_sol,1), ubound(coef_sol,2), &
              ubound(coef_for,1), ubound(coef_for,3) )
      if ( maxval(abs(coef_sol(l,  1:l)))/=zero ) m = l
      if ( maxval(abs(coef_for(l,:,1:l)))/=zero ) m = l
   end do
   num_steps = m
   if ( num_steps<=0 .or. num_steps>max_steps) then
      write (message,*) "too many (or too few) steps:  num_steps = ", num_steps
      call sitpack_message( "sitpack_set_custom", message, 1 )
      return
   end if

! check that the scheme for step  l  uses no more than  l  steps

   do l=1,num_steps
      if ((l<=ubound(coef_sol,dim=1).and.maxval(abs(coef_sol(l+1:,  l)))/=0)  &
      .or.(l<=ubound(coef_for,dim=1).and.maxval(abs(coef_for(l+1:,:,l)))/=0)) &
      then
         write (message,*) "step = ", l, " uses too many steps"
         call sitpack_message( "sitpack_set_custom", message, 1 )
         return
      end if
   end do

! set up the scheme as requested

   c_sol(0:m,    1:m) = coef_sol(0:m,    1:m)
   c_for(0:m,1:n,1:m) = coef_for(0:m,1:n,1:m)

! store the scheme labels if provided

   if (present(names)) then
      methods(1:num_steps,1:num_terms) = names(1:num_steps,1:num_terms) 
   end if

! determine what needs to be stored after each step

   call process_scheme( "sitpack_set_custom" )

! trace output

   if (trace_routine( "sitpack_set_custom" )>0) then
      call output_scheme
   end if

! quit if the scheme is inconsistent

   if (maxval(order(1:num_steps,1:num_terms))<=0) then
      call sitpack_message( "sitpack_set_custom", "scheme is inconsistent", 2 )
   end if

!-------------------------------------------------------------------------------
end subroutine sitpack_set_custom
!===============================================================================


!===============================================================================
subroutine process_scheme( called_from )
!-------------------------------------------------------------------------------
! Purpose: 
!    Processes parameters of semi-implicit scheme for use by sitpack_step
!
! Argument: 
   character (len=*), intent(in) :: called_from ! name of calling routine
!
! Note: 
!    Indexing here matches sitpack_step:  j=1  is current step
!-------------------------------------------------------------------------------

   integer (kind=int_kind) :: j, k, l, p          ! indices
   integer (kind=int_kind) :: j1, j2, jsol, j1max ! temps
   real    (kind=dbl_kind) :: ssum, fsum          ! temps
   real    (kind=dbl_kind) :: ztol                ! tolerance

! check that the scheme is workable

   if (count(c_sol(0,1:num_steps)==0)>0) then
      call sitpack_message( called_from, &
         "each scheme must involve the new time level", 2 )
   end if

! identify which solution and forcing copies must be stored after each step

   store_sol = .false.
   store_for = .false.
   do l=num_steps,1,-1

   ! first determine the last solution which is needed directly

       jsol = 0
       do j=1,num_steps
          if (c_sol(j,l)==zero) cycle
          jsol = j
       end do

   ! decide which forcing terms to store
   ! j1 is index of first forcing term to be used (and computed in this step)
   ! j2 is index of last  forcing term to be used (and usually discarded)

      j1max = 0
      do k=1,num_terms
         if (l==num_steps) then
             j1 = 0
             j2 = 0
             do j=1,num_steps
                if (c_for(j,k,l)==zero) cycle
                if (j1==0) j1 = j
                j2 = j
             end do
             if (j1<j2) store_for(j1:j2-1,k,l) = c_for(j1:j2-1,k,l)/=zero
             j1max = max( j1max, j1 )
         else
             do j=1,l
                if ((count(store_for(:,k,l+1))>0 .and. c_for(j+1,k,l+1)/=0) &
                    .or. store_for(j+1,k,l+1)) store_for(j,k,l) = .true.
             end do
         end if
      end do ! loop on  k

   ! store solution copies needed directly or for computing forcing 
      
      jsol = max( jsol, j1max )
      if (jsol>1) store_sol(1:jsol-1,l) = .true.
      if (l<num_steps) then
          do j=jsol,l
             if (c_sol(j+1,l+1)/=0 .or. store_sol(j+1,l+1)) &
                store_sol(j,l) = .true.
          end do
      end if
   end do ! loop on  l

! check the consistency and order of each scheme

   order = -1
   do l=1,num_steps

   ! set the zero tolerance relative to the size of the coefficients

      ztol = 10*epsilon(one)*maxval(abs(c_sol(0:l,l)))

   ! check the zero-order consistency condition

      ssum = sum(c_sol(0:l,l))
      if (abs(ssum)<ztol) order(l:,:) = 0

   ! check the first-order consistency condition

      ssum = dot_product( (/ (-(j-1)*one, j=0, l)/), c_sol(0:l,l) )
      do k=1,num_terms
         fsum = sum( c_for(0:l,k,l) )
         if (abs(ssum-fsum)<ztol .and. order(l,k)==0) order(l,k) = 1
      end do

   ! check the higher-order consistency conditions

      do p=2,l+1
         ssum = dot_product( (/ ((-j+1)**p*one, j=0, l)/), c_sol(0:l,l) )
         do k=1,num_terms
            fsum = p*dot_product( (/ ((-j+1)**(p-1)*one, j=0, l)/), &
                                  c_for(0:l,k,l) )
            if (abs(ssum-fsum)<ztol .and. order(l,k)==p-1) order(l,k) = p
         end do
      end do
       
   end do ! loop on  l

!-------------------------------------------------------------------------------
end subroutine process_scheme
!===============================================================================


!===============================================================================
subroutine output_scheme
!-------------------------------------------------------------------------------
! Purpose: 
!    Outputs parameters for a semi-implicit scheme
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: j, k, l ! indices
   integer (kind=int_kind) :: lunit   ! unit number for output
   character (len=128)     :: fmt     ! format for output
   character (len=1)       :: &! flag which copies are stored:
   flag_sol(0:max_steps,          max_steps), &! solution
   flag_for(0:max_steps,max_terms,max_steps)   ! forcing

! set the output unit from internal parameter  iunit

   lunit = abs(iunit)
   if ( lunit==0 ) lunit = 6

   flag_sol = " "
   where (store_sol) flag_sol(1:,:) = "*"
   flag_sol(0,:) = "*"
   flag_for = " "
   where (store_for) flag_for(1:,:,:) = "*"

   write (lunit,*)
   write (lunit,*) "Parameters defining time differencing scheme:"
   write (lunit,*) "num_steps =", num_steps
   write (lunit,*) "num_terms =", num_terms
   do l=1,num_steps
      fmt = "(//' step  l =',i2,':  coefficients for time level  new-j')"
      write (lunit,fmt) l
      fmt = "(/' step  solution   ',8(a8:3x))"
      write (lunit,fmt) (methods(l,k),k=1,num_terms)
      do j=0,num_steps
         fmt = "(i3,2x,9(f10.6,a1))"
         write (lunit,fmt) j, c_sol(j,l), flag_sol(j,l), &
                           (c_for(j,k,l), flag_for(j,k,l), k=1,num_terms)
      end do
      fmt = "('          order:',9(i6:5x))"
      write (lunit,fmt) (order(l,k), k=1, num_terms)
   end do
   write (lunit,*)
   write (lunit,*) "(* means store this copy at end of step)"

! output the orders for each step and term

   write (lunit,*)
   write (lunit,*) "Order of accuracy by step and term:"
   fmt = "(/'  step',8('  term',i2:))"
   write (lunit,fmt) (k, k=1,num_terms)
   do l=1,num_steps
      fmt =  "(i5,8(i6,2x))"
      write (lunit,fmt) l, (order(l,k), k=1,num_terms)
   end do

!-------------------------------------------------------------------------------
end subroutine output_scheme
!===============================================================================


!===============================================================================
function get_new_pointer( which, pointers, available ) result (p_new)
!-------------------------------------------------------------------------------
! Purpose:
!    gets a pointer to an unused copy of solution or forcing
!
! Arguments:
!
   character (len=*),       intent(in) :: which       ! label for errors
   integer (kind=int_kind), intent(in) :: pointers(:) ! those used so far
   integer (kind=int_kind), intent(in) :: available   ! copies available
!
! Result:
   integer (kind=int_kind) :: p_new ! new pointer
!-------------------------------------------------------------------------------

! local variable

   integer (kind=int_kind) :: p_max ! largest pointer value to allow

   if ( available>0 ) then
      p_max = available
   else
      p_max = maxval( pointers ) + 1
   end if
   do p_new=1,p_max
      if (count(pointers==p_new)==0 ) return
   end do

! error exit:  couldn't find a copy not already in use

   write (message,*) "storage for ", which ," is insufficient"
   call sitpack_message( "get_new_pointer", message, 2 )

!-------------------------------------------------------------------------------
end function get_new_pointer
!===============================================================================


!===============================================================================
subroutine sitpack_step( dt, s, nstep_si, mstep_si, time_si )
!-------------------------------------------------------------------------------
! Purpose: 
!    Executes one semi-implicit time step
!
! Required Arguments:
   real    (kind=dbl_kind), intent(in)  :: dt ! time step (must be nonzero)

! Optional Arguments:
   integer (kind=int_kind), intent(out), optional :: &!
       s  ! index of current solution copy
   integer (kind=int_kind), intent(out), optional :: &!
      nstep_si ! number of steps since start
   integer (kind=int_kind), intent(out), optional :: &!
      mstep_si ! number of steps since restart
   real    (kind=dbl_kind), intent(out), optional :: &!
      time_si  ! model time
!-------------------------------------------------------------------------------

! scaled coefficients for the current step

   real (kind=dbl_kind) ::      &!  
      cs(max_steps),            &!  solution coefficients
      cf(max_steps,max_terms)    !  forcing  coefficients

! keep track of which solution and forcing copies to keep for next step

   logical ::                           &!  
      keep_s(max_steps),                &!  keep solution copy?
      keep_f(max_steps,max_terms)        !  keep forcing  copy?

! keep track of which forcing terms we need and which we have

   logical ::                           &!
      needed(max_steps,max_terms),      &!  forcing term needed?
      stored(max_steps,max_terms)        !  forcing term stored?

! local variables for stacking "pointers" and coefficients

   integer (kind=int_kind) :: n_stack
   integer (kind=int_kind), dimension(max_steps*max_terms) :: p_stack
   real    (kind=dbl_kind), dimension(max_steps*max_terms) :: c_stack

! other local variables and parameters

   integer (kind=int_kind) :: j, k, l       ! indices
   real    (kind=dbl_kind) :: tau(max_terms)! scaled time step
   logical ::           &! set trace_sitpack_step to sum of values:
      trace_beg,        &!  1: message at beginning of step
      trace_end,        &!  2: message at the end   of step
      point_beg,        &!  4: pointers at beginning of step
      point_mid,        &!  8: pointers at middle    of step (after cycle)
      point_end          ! 16: pointers at the end   of step

! trace output

   j = trace_routine( "sitpack_step" )
   trace_beg = mod(j   ,2)==1
   trace_end = mod(j/ 2,2)==1
   point_beg = mod(j/ 4,2)==1
   point_mid = mod(j/ 8,2)==1
   point_end = mod(j/16,2)==1
   if (trace_beg) then
      write (message,*) "starting  time step  nstep =", nstep+1
      call sitpack_message( "sitpack_step", message, header=2 )
   end if

! sanity checks

   if (dt==zero) then
      call sitpack_message( "sitpack_step", "time step must be nonzero", 2 )
   end if
   if (last_dt/=zero .and. dt/=last_dt) then
      call sitpack_message( "sitpack_step", "time step changed", 1 )
      call sitpack_restart
   end if
   last_dt = dt
   if (maxval(abs(c_sol))==zero .or. maxval(abs(c_for))==zero) then
      call sitpack_message( "sitpack_step", "time scheme not defined yet", 2 )
   end if
   if (p_sol(0)==0) then
      call sitpack_message( "sitpack_step", "sitpack not initialized yet", 2 )
   end if

! choose the scheme for the appropriate step (initialization or otherwise)

   mstep = mstep+1
   l = min( mstep, num_steps )          ! scheme for current step
   cs = -c_sol(1:,  l)/c_sol(0,l)       ! scaled coefficients for solution
   cf =  c_for(1:,:,l)*dt/c_sol(0,l)    ! scaled coefficients for solution
   tau = c_for(0,:,l)*dt/c_sol(0,l)     ! scaled time step for implicit problem
   keep_s = store_sol(:,  l)            ! which solution copies do we keep?
   keep_f = store_for(:,:,l)            ! which forcing  copies do we keep?
   if (point_beg) then
      call output_pointers(l,"sitpack_step:  pointers at beginning of step")
   end if

! cycle the "pointers" for this step (get all previously used information)

   p_sol(0:l)   = cshift( p_sol(0:l)  , shift=-1 )
   time(0:l)    = cshift( time(0:l)   , shift=-1 )
   p_for(0:l,:) = cshift( p_for(0:l,:), shift=-1, dim = 1 )
   time(0) = time(1) + dt
   nstep   = nstep + 1
   if (point_mid) then
      call output_pointers(l,"sitpack_step:  pointers after cycling")
   end if

! set "pointer" to space for RHS of implicit problem

   if (p_for(0,1)==0) p_for(0,1) = &
      get_new_pointer( "forcing" , pack(p_for,p_for/=0), num_for )
   if (point_mid) then
      call output_pointers(l,"sitpack_step:  set pointer for RHS")
   end if

! set up the linear combination of previous solutions 

   n_stack = count( cs(1:)/=0 )
   p_stack = 0
   c_stack = 0
   p_stack(1:n_stack) = pack( p_sol(1:), mask = cs(1:)/=0 )
   c_stack(1:n_stack) = pack( cs(1:)   , mask = cs(1:)/=0 )
   call compute_lincomb_sol( n_stack, p_stack, c_stack, p_for(0,1) )

! decide which forcing terms we need (to use) and which we have

   needed = cf/=zero
   stored = p_for(1:,:)/=0

! initialize the stack of forcing terms

   n_stack = 0
   p_stack = 0
   c_stack = 0

! compute or look up the forcing terms (add in those which will not be kept)

   do k=1,num_terms
   do j=num_steps,1,-1  ! reverse order so copies freed can be reused
      if (keep_f(j,k) .and. .not.stored(j,k)) then ! compute and store it
         if (p_for(j,k)==0) p_for(j,k) = &
            get_new_pointer( "forcing" , pack(p_for,p_for/=0), num_for )
         call compute_model_terms( k, p_sol(j), time(j), p_for(j,k), zero )
         stored(j,k) = .true.
      end if
      if (needed(j,k)) then
         if (keep_f(j,k)) then ! stack term to be used in linear combination
            n_stack          = n_stack + 1
            p_stack(n_stack) = p_for(j,k)
            c_stack(n_stack) = cf(j,k)
         else if (stored(j,k)) then ! use this term and free the space
            call compute_lincomb_for( 1, p_for(j,k), cf(j,k), p_for(0,1) )
            p_for(j,k) = 0
         else  ! use this term on the fly
            call compute_model_terms(k, p_sol(j), time(j), p_for(0,1), cf(j,k))
         end if
      end if
   end do ! loop on time level (j)
   end do ! loop on term (k)

! add in the remaining forcing terms (those which will be kept)

   if (n_stack>0) &
      call compute_lincomb_for( n_stack, p_stack, c_stack, p_for(0,1) )

! free space and (re)set "pointer" to space for new solution 

   where (.not.keep_s) p_sol(1:) = 0
   if (p_sol(0)==0) then
      p_sol(0) = get_new_pointer( "solution", p_sol, num_sol )
   end if
   if (point_mid) then
      call output_pointers(l,"sitpack_step:  set pointer for new solution")
   end if
   where (.not.keep_f) p_for(1:,:) = 0

! solve the implicit problem (or get starting values directly)

   if ( l<num_steps .and. direct_start ) then
      call put_model_solution( time(0), p_sol(0) )
   else
      call solve_implicit_eqns( tau(1:num_terms),time(0),p_for(0,1),p_sol(0) )
   end if

! free up the space for the RHS 

   p_for(0,1) = 0

! return the information about the new solution

   if (present(s))        s        = p_sol(0)
   if (present(nstep_si)) nstep_si = nstep
   if (present(mstep_si)) mstep_si = mstep
   if (present(time_si))  time_si  = time(0)

! trace output

   if (point_end) then
      call output_pointers(l,"sitpack_step:  pointers at end of step")
   end if
   if (trace_end) then
      write (message,*) "completed time step  nstep =", nstep, &
                        "  time =", time(0)
      call sitpack_message( "sitpack_step", message, footer=2 )
   end if

! all done

!-------------------------------------------------------------------------------
end subroutine sitpack_step
!===============================================================================


!===============================================================================
subroutine sitpack_restart
!-------------------------------------------------------------------------------
! Purpose: 
!    Restarts the semi-implicit scheme at the current time
!-------------------------------------------------------------------------------

   integer (kind=int_kind) :: s ! pointer to current solution

   if ( p_sol(0) == 0 ) then
      write (message,*) "must initialize scheme first"
      call sitpack_message( "sitpack_restart", message, 1 )
      return
   end if
   s = p_sol(0)
   p_sol = 0
   p_for = 0
   p_sol(0) = s
   mstep = 0
   last_dt = zero

! trace output

   if (trace_routine( "sitpack_restart")>0) then
      write (message,*) "restarted scheme at time =", time(0), &
                        "  nstep =", nstep
      call sitpack_message( "sitpack_restart", message )
   end if

!-------------------------------------------------------------------------------
end subroutine sitpack_restart
!===============================================================================


!===============================================================================
subroutine output_pointers( l, mess )
!-------------------------------------------------------------------------------
! Purpose: 
!    Outputs pointers (to stdin) to help trace/debug the semi-implicit code
!
! Arguments:
   integer (kind=int_kind), intent(in) :: l     ! index of current step
   character (len=*),       intent(in) :: mess  ! label for header
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: j, k ! indices
   integer (kind=int_kind) :: lunit ! unit number for output
   character (len=128) :: fmt ! format for output

! set the output unit from internal parameter  iunit

   lunit = abs(iunit)
   if ( lunit==0 ) lunit = 6

   write (lunit,*) mess
   write (lunit,*)
   write (lunit,*)  "                 solution  forcing pointers"
   fmt = "(' step     time    pointers  ',8(a8:2x))"
   write (lunit,fmt) (methods(l,k),k=1,num_terms)
   do j=0,num_steps
      fmt = "(i4,2x,es10.2,i7,2x,8(i5:5x))"
      write (*,fmt) j, time(j), p_sol(j), (p_for(j,k), k=1,num_terms)
   end do

!-------------------------------------------------------------------------------
end subroutine output_pointers
!===============================================================================


!===============================================================================
subroutine sitpack_trace( routine_name, itrace )
!-------------------------------------------------------------------------------
! Purpose: 
!    Turns on tracing options for routines in this module
!
! Arguments:
   character (len=*), intent(in) :: routine_name           ! routine to trace
   integer (kind=int_kind), intent(in), optional :: itrace ! select trace option
!
! Notes:
!    The name  routine_name  must be lowercase (unrecognized names are ignored).
!
!    Positive  itrace  values specify trace options; see code of individual 
!    routines for details (not all routines may have tracing options).
!    If  itrace  is omitted, the default value  itrace=1  is used.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: iroutine ! index of routine name in list
   integer (kind=int_kind) :: ltrace   ! local copy of itrace

! look up  routine_name  in the list of routines to trace

   iroutine = trace_routine( routine_name )

! set the value of  itrace  if not supplied

   ltrace = 1
   if (present(itrace)) ltrace = itrace
   if ( iroutine>=0 ) then
      trace_option(iroutine) = ltrace
   else
      if (ntrace>=mtrace) then
         call sitpack_message( "sitpack_trace", "too many routines to trace" )
         return
      end if
      ntrace = ntrace+1
      trace_name(ntrace)   = routine_name
      trace_option(ntrace) = max( ltrace, 0 )
   end if

!-------------------------------------------------------------------------------
end subroutine sitpack_trace
!===============================================================================


!===============================================================================
function trace_routine( routine_name ) result (itrace)
!-------------------------------------------------------------------------------
! Purpose: 
!    Internal routine to look up trace options set by user via  sitpack_trace
!
! Argument:
   character (len=*), intent(in) :: routine_name ! routine to trace
!
! Result:
   integer (kind=int_kind) :: itrace ! specifies the trace option set:
!    itrace>0  trace the routine (meaning depends on which routine)
!    itrace=0  no trace
!    itrace<0  no trace option set for this routine
!-------------------------------------------------------------------------------

! local variable

   integer (kind=int_kind) :: iroutine ! index of routine name in list

! see if routine_name is already on the list

   do iroutine=1,ntrace
      if ( trace_name(iroutine) == routine_name ) then
         itrace = trace_option(iroutine)
         return
      end if
   end do
   itrace = -1
   return

!-------------------------------------------------------------------------------
end function trace_routine
!===============================================================================


!===============================================================================
subroutine sitpack_write_dump( idump, iformat )
!-------------------------------------------------------------------------------
! Purpose: 
!    Dumps internal data to a restart file
!
! Arguments:
   integer (kind=int_kind), intent(in) :: &
      idump,   & ! unit number to dump to
      iformat    ! nonzero to write formatted dump, zero for binary
!
! Usage Notes:
!    This routine writes all relevant internal data from sitpack to a file.
!    You must also write all model variables (all solution and forcing copies).
!
!    To use this data to continue a model run from the dump, simply call
!    sitpack_read_dump (and restore all of your model variables and forcing 
!    terms as they were at the dump time), and then continue by calling
!    sitpack_step exactly as you were doing.  You do not need to call
!    sitpack_init unless you want to change the time scheme.
!-------------------------------------------------------------------------------

   if ( iformat==0 ) then

!    write binary dump

      write (unit=idump) num_sol, num_for, iunit, direct_start
      write (unit=idump) num_steps, num_terms, c_sol, c_for
      write (unit=idump) nstep, mstep, time, last_dt, p_sol, p_for

   else

!    write formatted dump

      write (unit=idump,fmt="('-----Dump of data from sitpack-----')")
      write (unit=idump,fmt=*) num_sol, num_for, iunit, direct_start
      write (unit=idump,fmt=*) num_steps, num_terms, c_sol, c_for
      write (unit=idump,fmt=*) nstep, mstep, time, last_dt, p_sol, p_for
      write (unit=idump,fmt="('------end of data from sitpack-----')")

   end if

!-------------------------------------------------------------------------------
end subroutine sitpack_write_dump
!===============================================================================


!===============================================================================
subroutine sitpack_read_dump( idump, iformat )
!-------------------------------------------------------------------------------
! Purpose: 
!    Reads internal data from a restart file written by sitpack_write_dump
!
! Arguments:
   integer (kind=int_kind), intent(in) :: &
      idump,   & ! unit number to read from
      iformat    ! nonzero to read formatted dump, zero for binary
!-------------------------------------------------------------------------------

   if ( iformat==0 ) then

!    read binary dump

      read (unit=idump) num_sol, num_for, iunit, direct_start
      read (unit=idump) num_steps, num_terms, c_sol, c_for
      read (unit=idump) nstep, mstep, time, last_dt, p_sol, p_for

   else

!    read formatted dump

      read (unit=idump,fmt=*)
      read (unit=idump,fmt=*) num_sol, num_for, iunit, direct_start
      read (unit=idump,fmt=*) num_steps, num_terms, c_sol, c_for
      read (unit=idump,fmt=*) nstep, mstep, time, last_dt, p_sol, p_for
      read (unit=idump,fmt=*)

   end if

! reset the remaining internal parameters needed for this scheme

   call process_scheme( "sitpack_read_dump" )

!-------------------------------------------------------------------------------
end subroutine sitpack_read_dump
!===============================================================================



!===============================================================================
subroutine sitpack_details( ps, pf, t )
!-------------------------------------------------------------------------------
! Purpose: 
!    Returns internal details (pointers and times) for debugging
!
! Arguments:
   integer (kind=int_kind), intent(out) :: &
      ps(0:),   &! pointers to the solution (index j   for time level)
      pf(0:,:)   ! pointers to the forcing  (index j,k for time level, term)
   real    (kind=dbl_kind), intent(out) :: t(0:) ! corresponding model times
!-------------------------------------------------------------------------------

   ps(0:num_steps)           = p_sol(0:num_steps)
   pf(0:num_steps,num_terms) = p_for(0:num_steps,num_terms)
    t(0:num_steps)           =  time(0:num_steps)

!-------------------------------------------------------------------------------
end subroutine sitpack_details
!===============================================================================


!===============================================================================
subroutine sitpack_message( called_from, message, level, header, footer )
!-------------------------------------------------------------------------------
! Purpose: 
!    Message and error handler for sitpack
!
! Arguments:
   character (len=*), intent(in) :: called_from   ! calling routine
   character (len=*), intent(in) :: message       ! single line to output
   integer (kind=int_kind), intent(in), optional :: level
!    level = 0:  message [default]
!    level = 1:  warning
!    level = 2:  error [halts execution]
   integer (kind=int_kind), intent(in), optional :: header, footer
!    =0:  insert formfeed (Ctrl-L) before/after message
!    >0:  number of blank lines to insert before/after message
!
! Notes:
!    Errors are written to unit number  iunit  (can be reset by sitpack_setp).
!    It is up to you to ensure that this unit is open for output.
!    If  iunit<0, abs(iunit) is used and execution quits if level>0.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: lunit       ! unit number for output
   integer (kind=int_kind) :: loop        ! loop index
   integer (kind=int_kind) :: error_level ! loop index
   character (len=48) :: error_type       ! type of error

! set the output unit from internal parameter  iunit

   lunit = abs(iunit)
   if ( lunit==0 ) lunit = 6

! set the error type

   error_type = " "
!  if ( present(level) .and. level==1 ) error_type = "Warning in"
!  if ( present(level) .and. level==2 ) error_type = "Error in"
   if ( present(level) ) then
      error_level = level
      if (level==1 ) error_type = "Warning in"
      if (level==2 ) error_type = "Error in"
   else
      error_level = 0
      error_type = " "
   end if

! output the header (if present)

   if ( present(header) ) then
      if ( header>0 ) then
         do loop=1,header
            write(lunit,*) " "
         end do
      else if ( header==0 ) then
         write(lunit,*) ""
      end if
   end if

! output the message

   if ( error_level>0 ) then
      write(lunit,*) " "
      write(lunit,*) trim(error_type)//" "//trim(called_from)//":"
      write(lunit,*) trim(message)
      write(lunit,*) " "
   else if ( called_from/=" " ) then
      write(lunit,*) trim(called_from)//":  "//trim(message)
   else
      write(lunit,*) trim(message)
   end if

! output the footer (if present)

   if ( present(footer) ) then
      if ( footer>0 ) then
         do loop=1,footer
            write(lunit,*) " "
         end do
      else if ( footer==0 ) then
         write(lunit,*) ""
      end if
   end if

! halt execution (if appropriate)

   if ( error_level>=2 .or. (error_level>0 .and. iunit<0) ) &
      stop "Fortran stop in sitpack"

!-------------------------------------------------------------------------------
end subroutine sitpack_message
!===============================================================================


!-------------------------------------------------------------------------------
end module sitpack
!===============================================================================
