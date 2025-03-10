!===============================================================================
module pswm_pars
!-------------------------------------------------------------------------------
! Purpose:
!    Contains PSWM model parameters (user-specifiable constants)
!    Also contains routines for reading and writing them.
!
! Notes:
!    The model parameters are listed and explained below (with defaults).
!    They may be read by the routine "read_params"; if any one is not
!    present in the input file it retains the default value assigned below.
!    Values may be written to a file or stdout by the routine "write_params".
!
! See also:
!    pswm_cons  model constants  (not user-specifiable, constant)
!    pswm_vars  model variables  (change during a run)
!    pswm_data  initial condition and forcing parameters (problem-dependent)
!-------------------------------------------------------------------------------

use kinds

implicit none   ! all module variables untyped by default
public          ! all module variables and routines accessible by default
save            ! all module variables static

! handy constants
real (kind=dbl_kind), parameter :: zero = 0.0_dbl_kind
real (kind=dbl_kind), parameter :: one  = 1.0_dbl_kind
real (kind=dbl_kind), parameter :: hour = 3600_dbl_kind
real (kind=dbl_kind), parameter :: day  = 24*hour
real (kind=dbl_kind), parameter :: km   = 1000_dbl_kind

!-------------------------------------------------------------------------------
! User-definable model parameters (with explanation and default values)
!-------------------------------------------------------------------------------

! Problem domain and physical constants:
! -------------------------------------
real (kind=dbl_kind) ::  & ! Problem domain (in meters):
   x0 = zero, y0 = zero, & ! origin (lower-left corner)
   xl = zero, yl = zero    ! domain length in  x  and  y
! if xl<=zero then it is reset to 2*abs(x0) (gives centered domain if x0<0)
! (same for yl in terms of y0)

real (kind=dbl_kind) :: & ! Reference phase speed:
   cval = one             ! c = sqrt( phibar )  (meters/second)

real (kind=dbl_kind) :: & ! Coriolis parameter (specify one or the other):
   fcor = zero,         & ! Coriolis parameter (1/seconds)
   dlat = zero            ! latitude (degrees)
! Note:  Specify either  fcor  or  dlat  (the other will be set accordingly)

! Model equations:
! ---------------
integer (kind=int_kind) :: ieq = 2 ! Specifies the form of the model equations:
!  ieq   predict   name
!  ---   -------   ----
!   1    u, v, p   momentum
!   2    z, d, p   vorticity/divergence
!   3    q, d, p   potential vorticity/divergence
! Note:  If  ieq<0  then the corresponding linear form will be used instead.

integer (kind=int_kind) :: icond = 4 ! Specifies the type of initialization:
!  icond=1:  initial  u, v, p  specified
!  icond=2:  initial  z, d, p  specified
!  icond=3:  initial  q, d, p  specified
!            Note:  the above match the corresponding values of abs(ieq)
!  icond=4:  initial  z  specified
!            (all other fields will be obtained from nonlinear balance)
!  Otherwise, all initial fields are assumed to be zero.
!  Note:  The initial fields are specified in the module pswm_data

integer (kind=int_kind) :: iforce = 0 ! Specifies the type of forcing:
!  iforce=0:  zero forcing
!  iforce=1:  nonseparable (functions of x, y, and t)
!  iforce=2:  separable (function of (x,y) times function of t)
!  Note:  The forcing is specified in the module pswm_data

! Dissipation (see also Sponge):
! -----------------------------
! The dissipation operator in each predictive equation is -cdiss*(-del^2)^idiss.
! You can specify  idiss  and either  cdiss  or  tdiss  as follows:

integer (kind=int_kind) :: idiss = 1 ! Type of disspation:
!  idiss>0:  (-del^2)^idiss [1 for del^2, 2 for del^4, ...]
!  idiss=0:  Rayleigh friction (operator is simply  -cdiss* )

real (kind=dbl_kind) :: & ! Amount of dissipation:
   cdiss = zero, & ! coefficient of dissipation (specified directly)
   tdiss = zero    ! corresponding e-folding time (in HOURS) for shortest wave
! Note:  Specify either  cdiss  or  tdiss  (the other will be set accordingly)
!        (if both are zero then there is no dissipation)

! Sponge (see also Dissipation):
! -----------------------------
! The model can include a "sponge" layer around the "boundaries" of the domain
! to mitigate the effects of using periodicity for non-periodic problems.
! The sponge layer is a space-dependent dissipation term added to each
! predictive equation; the form of each term is the same as the dissipation
! above, i.e., the negative Laplacian raised to the exponent  isponge.
! The coefficient is specified by either  csponge  or  tsponge, and then
! multiplied by a space-dependent factor which has the value one at the 
! domain boundaries and tends towards zero with length scale  lsponge
! (for details of the space dependence see routine sponge in pswm_terms).

integer (kind=int_kind) :: isponge = 0 ! Type of sponge:
!  isponge>0:  (-del^2)^isponge [1 for del^2, 2 for del^4, ...]
!  isponge=0:  Rayleigh friction (operator is simply  -csponge* )

integer (kind=int_kind) :: fsponge = 2 ! Sponge function:
!  fsponge=1:  linear "ramp" function (to zero at  lsponge)
!  fsponge=2:  cubic Hermite spline (to zero at  lsponge, smoother)
!  fsponge=3:  "bump" function (to zero at  lsponge, fully smooth)
!  fsponge=4:  Gaussian (e-folding scale  lsponge)
! Each function is one at the domain boundary and decays with scale  lsponge

real (kind=dbl_kind) :: lsponge = -0.1 ! Length scale of sponge:
!  lsponge>0  use length scale  lsponge (see setup_sponge in pswm_terms)
!  lsponge<0  use length scale  abs(lsponge)*(xl,yl) in (x,y), respectively
!  lsponge=0  there is no sponge

real (kind=dbl_kind) :: & ! Amplitude of sponge:
   csponge = zero, & ! coefficient of sponge (specified directly)
   tsponge = zero    ! corresponding e-folding time (in HOURS) for shortest wave
! Note:  Specify either csponge or tsponge (the other will be set accordingly)
!        (if both are zero then there is no sponge)

! Discretization:
! --------------
integer (kind=int_kind) :: & ! transform grid size
   nx = 64, & ! number of grid intervals in  x  on transform grid
   ny = 64    ! number of grid intervals in  y  on transform grid
integer (kind=int_kind) :: & ! spectral truncation (limited by nmax)
   mx =  0, & ! spectral truncation in  x  (index of last  x  mode carried)
   my =  0    ! spectral truncation in  y  (index of last  y  mode carried)
! If  mx  or  my  is zero, it will be reset to the largest value which avoids
! aliasing in the quadratic terms, i.e., such that  nx>3*mx  and  ny>3*my.
! Some typical values:
!   nx=ny:   32   48   64   96  128  192  256  384  512
!   mx=my:   10   15   21   31   42   63   85  127  170

integer (kind=int_kind) :: itd = 3 ! specifies the time discretization scheme
! Note:  see routine setup_time for choices of time discretization scheme

real (kind=dbl_kind) :: dt = one !  time step for execution (in seconds)

! EH TEMPORARY FIX FOR TIME DEPENDENT HEATING
!
!real (kind=dbl_kind) :: nnn = one 

! Output specification:
! --------------------
character (len=80) :: runid ="pswm_run" ! run identifier (used for file names)
! Note:  runid is NOT read by read_params (since it's used for file names).
! Use get_runid to set it from command line or read it from file runid.txt.

real (kind=dbl_kind) :: tstop = -one ! model time to end run (in seconds)
! (if negative, stop after abs(tstop) time steps)

real (kind=dbl_kind) :: & ! output time steps (interpretation below):
   dtout0 = -one,       & ! scalar diagnostics 
   dtout2 = zero,       & ! 2D field output
   dtoutr = zero          ! restart file
! These are interpreted as follows:
!    dtout>0  output every  freq  hours
!    dtout<0  output every -freq  time steps
!    dtout=0  no output

integer (kind=int_kind) :: & ! indices of output points (on transform grid):
   j1=0, j2=0, jinc=1,     & ! output at x(j) = x0 + j*xl/nx, (j=j1,j2,jinc)
   k1=0, k2=0, kinc=1        ! output at y(k) = y0 + k*yl/ny, (k=k1,k2,kinc)
!  If j1=j2=k1=k2=0 then the whole domain is output:  j=0,nx  and k=0,ny

integer (kind=int_kind) :: ifield = 31 ! specifies which field(s) to output:
!   1: u           x-component of velocity
!   2: v           y-component of velocity
!   4: p           scaled geopotential deviation:  p = grav*(h - href)/cval
!   8: d = delta   divergence
!  16: z = zeta    relative  vorticity
!  32: q           potential vorticity:  q = (fcor + zeta)/(1 + p/cval)
!  Set ifield to the sum of the desired values

integer (kind=int_kind) :: nsigd = 4 ! number of significant digits for output:
!  Fields are output one value per line in the format 1pew.d, with width w=d+8 
!  and d=nsigd-1 digits after the decimal point (see also mtout for details).

! Other parameters:
! ----------------
integer (kind=int_kind) :: iblow = 1 ! when to check:
!  iblow = 0:  never check
!  iblow = 1:  check every time step (not too expensive)
!  iblow = 2:  check only when scalar diagnostics are computed

real (kind=dbl_kind) :: blowup = 10*one ! condition to check
!  If blowup>0, then at each step the condition  c+p<blowup*c  is checked and 
!  the run is halted if it fails.  In the nonlinear cases the run is halted
!  if the solution bottoms out (c+p<=0), regardless of the value of blowup.

logical :: & ! logical parameters:
   output_ic   = .true. ! output initial conditions (if freq/=0)?

contains

!===============================================================================
subroutine read_params( in_unit, echo_unit )
!-------------------------------------------------------------------------------
! Purpose:
!    Reads the model parameters (from a file)
!
! Arguments:
   integer (kind=int_kind), intent(in), optional :: &
      in_unit, & ! unit to read params from [default:  read from stdin]
      echo_unit  ! unit to echo file to (6 for stdout) [default:  no echo]
!
! Notes:
!    The units in_unit and echo_unit (if specified) must be open.
!
!    The input is read one line at a time from the specified unit (or stdin).
!    Each input line is scanned for all of the model parameters listed in this 
!    module.  They may be in any order but must appear one per line in the 
!    following format:
!
!       name value [optional comment]
!       # lines starting with "#" are comments and may appear anywhere
!       ! lines starting with "!" are comments and may appear anywhere
!
!      The input is read until the end of file is reached or a line starting
!      with the word "end" is encountered.  This allows you to put parameters
!      to be read by other routines (e.g., read_icpars) in a single file.
!      For example:
!
!         # input file for model run
!         # model parameters
!         dt     3600      time step for model run [in seconds]
!         end of model parameters
!         # initial condition parameters
!         rmax    50       radius of maximum wind (to be read by read_icpars)
!         end of initial condition parameters
!
!    Separate names and values by whitespace; put multiword strings in quotes.
!    Any parameter not read in retains the default value assigned above.
!
!    Length and time units for input values can be specified by lines starting 
!    with "lunits" (for length units:  m  or km) or "tunits" (for time values:  
!    sec, min, hour, or day).  These units are used for subsequent values of
!    model parameters representing lengths and times only (not combinations).
!    For example, the file:
!
!       lunits    km
!       xl        2000
!       tunits    min
!       dt        2
!       tunits    day
!       tstop     10
!
!    specifies a ten-day model run with a time step of two minutes using 
!    the domain length 2000 km in x.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: iunit  ! local unit number for input
   character (len=40)      :: vname  ! for variable name
   character (len=132)     :: line   ! to hold one input line
   character (len=1)       :: cunits ! for one-character unit identifier
   integer (kind=int_kind) :: ios    ! iostatus
   real (kind=dbl_kind) :: &         ! length and time scales:
      lunits = one,        &         ! length units (in meters)
      tunits = one                   ! time   units (in seconds)
   logical             :: echo_input ! echo the input file?
   real (kind=dbl_kind), parameter :: omega = 7.292e-05_dbl_kind

! process the arguments input to this routine

   if (present( in_unit )) then
      iunit = in_unit
   else
      iunit = 5
   end if
   echo_input = present( echo_unit )

! read the lines, scanning for the model parameters

   readloop: do
      read ( unit=iunit, fmt="(a)", iostat=ios ) line
      if ( ios<0 ) then
         exit readloop
      else if ( ios>0 ) then
         write (*,*) " read_params:  cannot read input line"
         close ( unit=iunit )
         stop "read_params:  cannot read input line"
      else if ( echo_input ) then
         write (echo_unit,*) trim( line )
      end if
      if ( len_trim(line)==0 .or. line(1:1)=="#" .or. line(1:1)=="!" ) cycle
      read ( unit=line, fmt=* ) vname
      select case (vname)
         case ("end")
            return
! Problem domain and physical constants:
         case ("x0")
            read ( unit=line, fmt=* ) vname, x0
            x0 = lunits*x0
            if (xl<=zero) xl = 2*abs(x0)
         case ("xl")
            read ( unit=line, fmt=* ) vname, xl
            xl = lunits*xl
            if (xl<=zero) xl = 2*abs(x0)
         case ("y0")
            read ( unit=line, fmt=* ) vname, y0
            y0 = lunits*y0
            if (yl<=zero) yl = 2*abs(y0)
         case ("yl")
            read ( unit=line, fmt=* ) vname, yl
            yl = lunits*yl
            if (yl<=zero) yl = 2*abs(y0)
         case ("cval")
            read ( unit=line, fmt=* ) vname, cval
         case ("dlat")
            read ( unit=line, fmt=* ) vname, dlat
            fcor = 2*omega*sin( dlat*acos(-one)/180 )
         case ("fcor")
            read ( unit=line, fmt=* ) vname, fcor
            dlat = asin( fcor/2*omega )*180/acos(-one)
! Model equations:
         case ("ieq")
            read ( unit=line, fmt=* ) vname, ieq
         case ("icond")
            read ( unit=line, fmt=* ) vname, icond
         case ("iforce")
            read ( unit=line, fmt=* ) vname, iforce
! Dissipation
         case ("idiss")
            read ( unit=line, fmt=* ) vname, idiss
         case ("cdiss")
            read ( unit=line, fmt=* ) vname, cdiss
         case ("tdiss")
            read ( unit=line, fmt=* ) vname, tdiss
! Sponge
         case ("isponge")
            read ( unit=line, fmt=* ) vname, isponge
         case ("fsponge")
            read ( unit=line, fmt=* ) vname, fsponge
         case ("lsponge")
            read ( unit=line, fmt=* ) vname, lsponge
         case ("csponge")
            read ( unit=line, fmt=* ) vname, csponge
         case ("tsponge")
            read ( unit=line, fmt=* ) vname, tsponge
! Discretization
         case ("nx")
            read ( unit=line, fmt=* ) vname, nx
         case ("ny")
            read ( unit=line, fmt=* ) vname, ny
         case ("mx")
            read ( unit=line, fmt=* ) vname, mx
         case ("my")
            read ( unit=line, fmt=* ) vname, my
         case ("itd")
            read ( unit=line, fmt=* ) vname, itd
         case ("dt")
            read ( unit=line, fmt=* ) vname, dt
            dt = tunits*dt
! Output specification
         case ("tstop")
            read ( unit=line, fmt=* ) vname, tstop
            if ( tstop>zero ) tstop = tunits*tstop
         case ("dtout0")
            read ( unit=line, fmt=* ) vname, dtout0
            if ( dtout0>zero ) dtout0 = tunits*dtout0
         case ("dtout2")
            read ( unit=line, fmt=* ) vname, dtout2
            if ( dtout2>zero ) dtout2 = tunits*dtout2
         case ("dtoutr")
            read ( unit=line, fmt=* ) vname, dtoutr
            if ( dtoutr>zero ) dtoutr = tunits*dtoutr
         case ("j1")
            read ( unit=line, fmt=* ) vname, j1
         case ("j2")
            read ( unit=line, fmt=* ) vname, j2
         case ("jinc")
            read ( unit=line, fmt=* ) vname, jinc
         case ("k1")
            read ( unit=line, fmt=* ) vname, k1
         case ("k2")
            read ( unit=line, fmt=* ) vname, k2
         case ("kinc")
            read ( unit=line, fmt=* ) vname, kinc
         case ("ifield")
            read ( unit=line, fmt=* ) vname, ifield
         case ("nsigd")
            read ( unit=line, fmt=* ) vname, nsigd
! Units for time and length (for scaling input values only):
         case ("tunits")
          read ( unit=line, fmt=* ) vname, cunits
          if ( cunits.eq.'m' ) then
              tunits = 60*one
          else if ( cunits.eq.'h' ) then
              tunits = 60*60*one
          else if ( cunits.eq.'d' ) then
              tunits = 24*60*60*one
          else
              tunits = one
          end if
         case ("lunits")
          read ( unit=line, fmt=* ) vname, cunits
          if ( cunits.eq.'k' ) then
              lunits = 1000*one
          else
              lunits = one
          end if
! Other parameters
         case ("iblow")
            read ( unit=line, fmt=* ) vname, iblow
         case ("blowup")
            read ( unit=line, fmt=* ) vname, blowup
         case ("output_ic")
            read ( unit=line, fmt=* ) vname, output_ic
         case default
            write (*,*) " read_params:  unrecognized variable on line:"
            write (*,*) line
      end select
   end do readloop

!-------------------------------------------------------------------------------
end subroutine read_params
!===============================================================================


!===============================================================================
subroutine write_params( out_unit )
!-------------------------------------------------------------------------------
! Purpose:
!    Writes the values of all parameters to a text file
!
! Arguments:
   integer (kind=int_kind), intent(in), optional :: &
      out_unit ! unit number for output [default:  6 (stdout)]
!
! Notes:
!    For output to a file, you must open it first (and close it later).
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: iunit ! unit number for output

! write the output

   if ( present( out_unit ) ) then
      iunit = out_unit
   else
      iunit = 6
   end if

   write (iunit,*) "--------------------------------------------------------"
   write (iunit,*) "PSWM parameters:"
! Problem domain and physical constants:
   write (iunit,*) 'x0      = ', x0/km, " km"
   write (iunit,*) 'xl      = ', xl/km, " km"
   write (iunit,*) 'y0      = ', y0/km, " km"
   write (iunit,*) 'yl      = ', yl/km, " km"
   write (iunit,*) 'cval    = ', cval, " m/s"
   write (iunit,*) 'dlat    = ', dlat, " degrees"
   write (iunit,*) 'fcor    = ', fcor, " 1/s"
! Model equations:
   write (iunit,*) 'ieq     = ', ieq
   write (iunit,*) 'icond   = ', icond
   write (iunit,*) 'iforce  = ', iforce
! Dissipation:
   write (iunit,*) 'idiss   = ', idiss
   write (iunit,*) 'cdiss   = ', cdiss
   write (iunit,*) 'tdiss   = ', tdiss, " hours"
! Sponge:
   write (iunit,*) 'isponge = ', isponge
   write (iunit,*) 'fsponge = ', fsponge
   write (iunit,*) 'lsponge = ', lsponge
   write (iunit,*) 'csponge = ', csponge
   write (iunit,*) 'tsponge = ', tsponge, " hours"
! Discretization:
   write (iunit,*) 'nx      = ', nx
   write (iunit,*) 'ny      = ', ny
   write (iunit,*) 'mx      = ', mx
   write (iunit,*) 'my      = ', my
   write (iunit,*) 'itd     = ', itd
   write (iunit,*) 'dt      = ', dt, " seconds"
! Output specification:
   write (iunit,*) 'runid   = ', trim(runid)
   write (iunit,*) 'tstop   = ', tstop
   write (iunit,*) 'dtout0  = ', dtout0
   write (iunit,*) 'dtout2  = ', dtout2
   write (iunit,*) 'j1      = ', j1
   write (iunit,*) 'j2      = ', j2
   write (iunit,*) 'jinc    = ', jinc
   write (iunit,*) 'k1      = ', k1
   write (iunit,*) 'k2      = ', k2
   write (iunit,*) 'kinc    = ', kinc
   write (iunit,*) 'ifield  = ', ifield
   write (iunit,*) 'nsigd   = ', nsigd
! Other parameters:
   write (iunit,*) 'iblow   = ', iblow
   write (iunit,*) 'blowup  = ', blowup
   write (iunit,*) 'ouput_ic= ', output_ic
   write (iunit,*) "--------------------------------------------------------"

!-------------------------------------------------------------------------------
end subroutine write_params
!===============================================================================

!-------------------------------------------------------------------------------
end module pswm_pars
!===============================================================================
