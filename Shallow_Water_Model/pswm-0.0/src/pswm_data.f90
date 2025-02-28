!===============================================================================
module pswm_data
!-------------------------------------------------------------------------------
! Purpose:
!    Contains data for a model run:  initial conditions and forcing
!
! Description:
!    This module is the place where the data for a model run is defined.
!    The model calls the following routines:
!
!       put_initial  to put the initial values into the predicted variables
!       add_forcing  to compute the forcing of predicted variables
!       read_icpars  to read  initial condition parameters from a file
!       write_icpars to write initial condition parameters to a file
!
!    See the comments in each routine below for its calling sequence.
!    If you keep the interface to these routines the same, you should not
!    have to modify any of the rest of the model code to run your problem.
!
! See also:
!    pswm_terms  (for routines which call these routines)
!-------------------------------------------------------------------------------

use kinds
use pswm_pars
use pswm_cons
use ss_ops

implicit none ! all module variables untyped by default
private       ! all module variables and routines hidden by default

public          &! routines called from other places
   put_initial, &! specifies initial fields
   add_forcing, &! adds specified forcing
   add_forcing2, &! adds specified forcing two variables
   read_icpars, &! reads  initial condition parameters from a file
   write_icpars  ! writes initial condition parameters  to  a file

!-------------------------------------------------------------------------------
! Initial condition parameters (read by read_icpars, written by write_icpars):
!-------------------------------------------------------------------------------
! These parameters are stored as module variables in this module so they can 
! be accessed by these all routines in this module.  Any additional parameters 
! you need could be added similarly.  If you don't want to read these from a 
! file, you could declare them as "public" and then set them directly in your 
! main program (you'd need to include the statement "use pswm_data" there).

! parameters for initial Gaussian vortex (when  icond==1):
real (kind=dbl_kind) :: rmax = 50*km ! radius of maximum wind
real (kind=dbl_kind) :: vmax = 10    ! value  of maxium wind
real (kind=dbl_kind) :: imbalance = 0.01_dbl_kind  ! for initial imbalance
! parameters for initial PV ring (when  icond==4):
real (kind=dbl_kind) :: & ! radii which define the PV ring:
   r1 =  30*km, r2 =  45*km, r3 =  50*km, r4 =  65*km
real (kind=dbl_kind) :: & ! radii which define the mass sink:
   rm1 =  30*km, rm2 =  45*km, rm3 =  50*km, rm4 =  65*km
real (kind=dbl_kind) :: & ! relative vorticities (inner, outer):
   z1 = 4.1825e-04, z2 = 3.3450e-03 
real (kind=dbl_kind) :: & ! mass sink magnitudes (inner, outer):
   m1 = 4.1825e-04, m2 = 3.3450e-03 
real (kind=dbl_kind) :: pamp = 1.0e-02 ! amplitude relative to z2
integer (kind=int_kind) :: & ! min, max wavenumbers to perturb:
   minwav = 3, maxwav = 3
! parameters for elliptical vortex (when  icond==4):
real (kind=dbl_kind) :: & ! radii which define the PV ring:
   ri =  30*km, ro =  60*km, ze =  3.0e-03, ec =  0.7, lam = 2.0
! parameters for linearly increasing vortex (when  icond==4):
real (kind=dbl_kind) :: & ! radii which define the PV ring:
   r3l =  30*km, r4l =  60*km, zp1 =  3.0e-03
! parameters for both initial conditions:
real (kind=dbl_kind) :: xc = zero, yc = zero ! center of vortex/ring
!-------------------------------------------------------------------------------

! other module variables

real (kind=dbl_kind) :: xx, yy   ! coordinates reduced by periodicity
real (kind=dbl_kind) :: r, phi   ! polar coordinates relative to the center
logical, parameter   :: trace = .false. ! trace routines in this module?

contains

!===============================================================================
subroutine put_initial( which, f )
!-------------------------------------------------------------------------------
! Purpose:
!    Puts an initial field into a model variable (in physical space)
!
! Arguments:
   character (len=*), intent(in) :: which ! specifies which field to initialize
   real (kind=dbl_kind), dimension(-1:,-1:), intent(out) :: f ! array for field
!-------------------------------------------------------------------------------
! This version:  distinguishes between two different initial conditions,
! controlled by the model parameter  icond  (must be either 1 or 4 here).
!-------------------------------------------------------------------------------

   integer (kind=int_kind) :: j, k            ! indices for gridpoints
   real (kind=dbl_kind) :: tmp                ! EH: a temporary variable

   if (trace) write (*,*) "put_initial:  initializing ", which

   select case (which)

      case ("u","v","p") ! initial u,v,p specified (when icond=1)
         do k=0,ny-1
         do j=0,nx-1
            f(j,k) = uvp0fun( which, x(j), y(k) )
         end do
         end do

      case ("z") ! initial relative vorticity specified (when icond=4)
         do k=0,ny-1
         do j=0,nx-1
            f(j,k) = z0fun( x(j), y(k) )
!             f(j,k) = ellz0fun( x(j), y(k) )
!             f(j,k) = z0funlin( x(j), y(k) )
         end do
         end do

!         ! EH temporary code to read in vorticity from a file   
!         ! overwrite f(j,k) in previous step      
!
!         if (cval .ge. 200) then
!         open ( unit=15, file='e7o20hires_ed.init', status='old' )        
!         do k=0,ny-1
!         do j=0,nx-1
!            read(15,*) tmp
!            f(j,k) = tmp
!         end do
!         end do
!         close(15)
!         end if

      case default ! all other initial fields vanish
         f = zero

   end select

! make the initial field explicitly periodic

   f(-1,0:ny-1) = f(nx-1,0:ny-1)
   f(nx,0:ny-1) = f(   0,0:ny-1)
   f(:,-1) = f(:,ny-1)
   f(:,ny) = f(:,   0)

   if (trace) write (*,*) "... done initializing ", which

!-------------------------------------------------------------------------------
end subroutine put_initial
!===============================================================================


!===============================================================================
function uvp0fun( which, x, y )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes initial velocity and momentum
!
! Arguments:
   character (len=*),    intent(in) :: which ! which field to initialize
   real (kind=dbl_kind), intent(in) :: x, y  ! coordinates of point
   real (kind=dbl_kind) ::           uvp0fun ! returns value
!-------------------------------------------------------------------------------
! This version:  quasi-balanced initial vortex (called only if  icond==1)
!
! The scaled geopotential deviation is the Gaussian function
!
!    p(r) = -(1+imbalance)*pamp*exp( (1-(r/rmax)^2)/2 )
!
! where  r = sqrt( (x-xc)^2 + (y-yc)^2 )  is the distance from the center at
! (xc,yc).  The velocity is a circular vortex with tangential component
!
!    v(r) = vmax*(r/rmax)*exp( (1-(r/rmax)^2)/2 )
!
! which has its maximum value  v=vmax  at radius  r=rmax.  The amplitude of the
! Gaussian is  pamp = fcor*rmax*vmax/cval,  which gives geostrophic balance in
! the linear case; this is modified by the parameter  imbalance  which allows
! for unbalanced initial conditions (to generate some gravity-inertia waves).
! Default values of the parameters  xc, yc, rmax, vmax,  and  imbalance  are
! set above; these may be changed by reading them from a file by read_icpars.
!-------------------------------------------------------------------------------

! local variables:

   real (kind=dbl_kind) :: efact ! exponential factor in  v  and  p

! compute the location in polar coordinates (relative to the center)

   xx = x0 + modulo( x-xc-x0, xl ); yy = y0 + modulo( y-yc-y0, yl )
   r   = sqrt( xx**2 + yy**2 ); phi = atan2( yy, xx )

! radial structure of initial vorticity and perturbation

   efact = exp( (one - (r/rmax)**2)/2 )
   select case (which)

      case ("u") ! initial velocity component  u
         uvp0fun = -vmax*(r/rmax)*efact*sin(phi)

      case ("v") ! initial velocity component  v
         uvp0fun =  vmax*(r/rmax)*efact*cos(phi)

      case ("p") ! initial scaled geopotential deviation
         uvp0fun = -(one+imbalance)*fcor*rmax*vmax*efact/cval

      case default ! all other initial fields vanish
         stop "uvp0fun:  called with bad argument 'which'"

   end select

!-------------------------------------------------------------------------------
end function uvp0fun
!===============================================================================


!===============================================================================
function z0fun( x, y )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes initial vorticity at a point
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: x, y ! coordinates of point
   real (kind=dbl_kind) ::            z0fun ! initial vorticity
!-------------------------------------------------------------------------------
! This version:  PV ring  (called only if  icond==4)
!
! The relative vorticity  z  at radius r from the center (xc,yc) is:
!           ( z1                0 <= r <= r1
!           ( (transition)     r1 <  r <  r2
!    z(r) = ( z2               r2 <= r <= r3
!           ( (transition)     r3 <  r <  r4
!           ( 0                r4 <  r
! where the transitions are continuous and smooth (cubic Hermite splines).
! The outer band is perturbed with amplitude  amp*z2  in each of azimuthal
! wavenumbers  m = minwav, ..., maxwav  with random phase for each  m.  
! Default values of the parameters  xc, yc, r1, r2, r3, r4, z1, z2, minwav,
! maxwav,  and  pamp  are set above; these may be changed by reading them 
! from a file by read_icpars.
!-------------------------------------------------------------------------------

! local variables:

   real (kind=dbl_kind) :: zr       ! radial part of z0fun
   real (kind=dbl_kind) :: pr, pphi ! radial, azimuthal parts of perturbation
   real (kind=dbl_kind) :: phase(minwav:maxwav) ! phases of perturbation
   integer (kind=int_kind) :: i, m  ! indices
   logical              :: first = .true. ! first call? (if so, initialize)

! compute the location in polar coordinates (relative to the center)

   xx = x0 + modulo( x-xc-x0, xl ); yy = y0 + modulo( y-yc-y0, yl )
   r   = sqrt( xx**2 + yy**2 ); phi = atan2( yy, xx )

! radial structure of initial vorticity and perturbation

   if ( r<=r1 ) then
      zr = z1
      pr = zero
   else if ( r<r2 ) then
      zr = z1*shape_fun((r-r1)/(r2-r1)) + z2*shape_fun((r2-r)/(r2-r1))
      pr =                                z2*shape_fun((r2-r)/(r2-r1))
   else if ( r<=r3 ) then
      zr = z2
      pr = z2
   else if ( r<r4 ) then
      zr = z2*shape_fun((r-r3)/(r4-r3))
      pr = z2*shape_fun((r-r3)/(r4-r3))
   else
      zr = zero
      pr = zero
   end if

! phase of perturbation

   if ( first ) then
      call random_seed(size=m) ! returns number of digits used in seed
      call random_seed(put=(/(i,i=1,m)/)) ! set seed for repeatability
      call random_number( phase )
      phase = 2*acos( -one )*phase
      first = .false.
   end if

! azimuthal structure of perturbation

   pphi = zero
   do m=minwav,maxwav
      pphi = pphi + cos( m*(phi-phase(m)) )
   end do

! put the pieces together to form the initial vorticity

   z0fun = zr + pamp*pr*pphi

contains

!===============================================================================
function shape_fun( s )
!-------------------------------------------------------------------------------
! Purpose:
!    Internal function for cubic Hermite shape function (for transitions)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: s ! value in the interval [0,1]
   real (kind=dbl_kind) :: shape_fun ! array for field
!-------------------------------------------------------------------------------
   shape_fun = (1-s)**2*(1+2*s)
!-------------------------------------------------------------------------------
end function shape_fun
!===============================================================================
!-------------------------------------------------------------------------------
end function z0fun
!===============================================================================

!===============================================================================
function ellz0fun( x, y )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes initial vorticity at a point
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: x, y ! coordinates of point
   real (kind=dbl_kind) ::            ellz0fun ! initial vorticity
!-------------------------------------------------------------------------------
! This version:  elliptical vortex  (called only if  icond==4)
!
! The relative vorticity  z  at radius r from the center (xc,yc) is:
!           
!           (       1          0  <  r <  ri * a(phi)
!    z(r,phi)   1 - f(r')      ri*a(phi) <= r <= ro*a(phi)
!           (       0          ro*a(phi) <  r 
!           
!  This is the elliptical vortex of Guinn (1992)
!-------------------------------------------------------------------------------

! local variables:

   real (kind=dbl_kind) :: zr       ! radial part of z0fun
   real (kind=dbl_kind) :: rp, fl, al   
   integer (kind=int_kind) :: i, m  ! indices
   logical              :: first = .true. ! first call? (if so, initialize)

! compute the location in polar coordinates (relative to the center)

   xx = x0 + modulo( x-xc-x0, xl ); yy = y0 + modulo( y-yc-y0, yl )
   r   = sqrt( xx**2 + yy**2 ); phi = atan2( yy, xx )

! radial structure of initial vorticity and perturbation

   al = ((one - ec**2) / (one - ec**2*(cos( phi ))**2) )**0.5
   rp = (r - ri * al) / (ro * al - ri * al)
   fl = exp( -(lam/rp)*exp(one/(rp - one)) ) 

   if ( r<=ri * al ) then
      zr = ze
   else if ( r< ro * al ) then
      zr = ze * ( one - fl )
   else
      zr = zero
   end if

! put the pieces together to form the initial vorticity

   ellz0fun = zr 

!-------------------------------------------------------------------------------
end function ellz0fun
!===============================================================================

!===============================================================================
function z0funlin( x, y )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes initial vorticity at a point (linearly increasing)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: x, y ! coordinates of point
   real (kind=dbl_kind) ::            z0funlin ! initial vorticity
!-------------------------------------------------------------------------------
! This version:  PV ring  (called only if  icond==4)
!
! The relative vorticity  z  at radius r from the center (xc,yc) is:
!           ( z1                0 <= r <= r1
!           ( (transition)     r1 <  r <  r2
!    z(r) = ( z2               r2 <= r <= r3
!           ( (transition)     r3 <  r <  r4
!           ( 0                r4 <  r
! where the transitions are continuous and smooth (cubic Hermite splines).
! The outer band is perturbed with amplitude  amp*z2  in each of azimuthal
! wavenumbers  m = minwav, ..., maxwav  with random phase for each  m.  
! Default values of the parameters  xc, yc, r1, r2, r3, r4, z1, z2, minwav,
! maxwav,  and  pamp  are set above; these may be changed by reading them 
! from a file by read_icpars.
!-------------------------------------------------------------------------------

! local variables:

   real (kind=dbl_kind) :: zr       ! radial part of z0fun
   real (kind=dbl_kind) :: slope    ! vorticity slope 
   real (kind=dbl_kind) :: pr, pphi ! radial, azimuthal parts of perturbation
   real (kind=dbl_kind) :: phase(minwav:maxwav) ! phases of perturbation
   integer (kind=int_kind) :: i, m  ! indices
   logical              :: first = .true. ! first call? (if so, initialize)

! compute the location in polar coordinates (relative to the center)

   xx = x0 + modulo( x-xc-x0, xl ); yy = y0 + modulo( y-yc-y0, yl )
   r   = sqrt( xx**2 + yy**2 ); phi = atan2( yy, xx )
   slope = 2.0*zp1/(r3l+r4l)

! radial structure of initial vorticity and perturbation

   if ( r<=r3 ) then
      zr = slope * r
   else if ( r<r4 ) then
      zr = zp1*shape_fun((r-r3l)/(r4-r3l))
   else
      zr = zero
   end if

   z0funlin = zr 

contains

!===============================================================================
function shape_fun( s )
!-------------------------------------------------------------------------------
! Purpose:
!    Internal function for cubic Hermite shape function (for transitions)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: s ! value in the interval [0,1]
   real (kind=dbl_kind) :: shape_fun ! array for field
!-------------------------------------------------------------------------------
   shape_fun = (1-s)**2*(1+2*s)
!-------------------------------------------------------------------------------
end function shape_fun
!===============================================================================
!-------------------------------------------------------------------------------
end function z0funlin
!===============================================================================


!===============================================================================
subroutine add_forcing( which, cforce, f )
!-------------------------------------------------------------------------------
! Purpose:
!    Adds the specified forcing term to a model variable (in spectral space)
!    eh (mod 9-19-07): change to physical space for ease of making ring
!
! Arguments:
   character (len=*),    intent(in) :: which ! specifies which field to force
   real (kind=dbl_kind), intent(in) :: cforce ! multiplier for forcing
   real (kind=dbl_kind), dimension(-1:,-1:), intent(inout) :: f ! forcing term
   real (kind=dbl_kind)    :: ffac  ! forcing factor h --> p  
   real (kind=dbl_kind)    :: zetamax
   integer (kind=int_kind) :: j, k      
!
! Note (carefully!):
!    The field must be multiplied by cforce (or by 1 if cforce=0) before adding.
!    This interpretation of cforce is different than in compute_model_terms,
!    since this routine must add the forcing to terms already computed.
!-------------------------------------------------------------------------------
! This version:  dummy version intended as a hook for specifying forcing later.
! You could use the parameter iforce (pswm_pars) to the flag separable case.
!-------------------------------------------------------------------------------

! dummy forcing:  replace "zero" with your forcing term (based on "which")

   call tophys( f )

   select case (which)
      case ("p")
      ! mass sink ring
      do k=0,ny-1
      do j=0,nx-1
         ffac = cval + f(j,k)
!         ffac = cval
         if (cforce==0) then
            ! we add 0.5 as a test to stop Gibbs phenomena  
            f(j,k) = f(j,k) - ffac*msfun( x(j), y(k) )
         else
            f(j,k) = f(j,k) - cforce*ffac*msfun( x(j), y(k) )
         end if
      end do
      end do

      case default ! all forcing terms are zero in this version
      if (cforce==0) then
         f = f + zero
      else
         f = f + cforce*zero
      end if 
   end select

   ! make the forcing field explicitly periodic

   f(-1,0:ny-1) = f(nx-1,0:ny-1)
   f(nx,0:ny-1) = f(   0,0:ny-1)
   f(:,-1) = f(:,ny-1)
   f(:,ny) = f(:,   0)

   call tospec( f )

!-------------------------------------------------------------------------------
end subroutine add_forcing
!===============================================================================

!===============================================================================
subroutine add_forcing2( which, cforce, f, g )
!-------------------------------------------------------------------------------
! Purpose:
!    Adds the specified forcing term to a model variable (in spectral space)
!    eh (mod 9-19-07): change to physical space for ease of making ring
!
! Arguments:
   character (len=*),    intent(in) :: which ! specifies which field to force
   real (kind=dbl_kind), intent(in) :: cforce ! multiplier for forcing
   real (kind=dbl_kind), dimension(-1:,-1:), intent(inout) :: f ! forcing term
   real (kind=dbl_kind), dimension(-1:,-1:), intent(inout) :: g ! forcing term
   real (kind=dbl_kind)    :: ffac  ! forcing factor h --> p  
   real (kind=dbl_kind)    :: zetamax
   integer (kind=int_kind) :: j, k      
!
! Note (carefully!):
!    The field must be multiplied by cforce (or by 1 if cforce=0) before adding.
!    This interpretation of cforce is different than in compute_model_terms,
!    since this routine must add the forcing to terms already computed.
!-------------------------------------------------------------------------------
! This version:  dummy version intended as a hook for specifying forcing later.
! You could use the parameter iforce (pswm_pars) to the flag separable case.
!-------------------------------------------------------------------------------

! dummy forcing:  replace "zero" with your forcing term (based on "which")

   call tophys( f )
   call tophys( g )

   select case (which)

     case ("p")
!     THIN 60 m/s
!      zetamax = 0.003175
!     THIN 100 m/s
!     zetamax = 0.007619
!     THIN 150 m/s
!     zetamax = 0.011429
!     MIDDLE 60 m/s
!      zetamax = 0.001852
!     MIDDLE 100 m/s
      zetamax = 0.004444
!     MIDDLE 150 m/s
!     zetamax = 0.006667
!     THICK 60 m/s
!      zetamax = 0.001481
!     THICK 100 m/s
!      zetamax = 0.003556 
!     THICK 150 m/s
!      zetamax = 0.005333     
      ! mass sink ring logistic
      do k=0,ny-1
      do j=0,nx-1
!         ffac = cval + f(j,k)
         ffac = cval
         if (g(j,k) .ge. zetamax) g(j,k) = zetamax
         if (cforce==0) then
            ! we add 0.5 as a test to stop Gibbs phenomena  
!   vorticity forcing
            f(j,k) = f(j,k) - (1.0 - (fcor+g(j,k))/(fcor+zetamax))*ffac* &
                              msfun( x(j), y(k) )
!   PV forcing, pass PV in here
!            f(j,k) = f(j,k) - (1.0 - (g(j,k))/(fcor+zetamax))*ffac* &
!                              msfun( x(j), y(k) )
!   new formulation with heating proportional to vorticity
!            f(j,k) = f(j,k) - (g(j,k))/(fcor+zetamax)* &
!                              (1.0 - (g(j,k))/(fcor+zetamax))*ffac* &
!                              msfun( x(j), y(k) )
         else
            f(j,k) = f(j,k) - cforce*(1.0 - (fcor+g(j,k))/(fcor+zetamax))*ffac* &
                              msfun( x(j), y(k) )
!            f(j,k) = f(j,k) - cforce*(1.0 - (g(j,k))/(fcor+zetamax))*ffac* &
!                              msfun( x(j), y(k) )
!   new formulation with heating proportional to vorticity
!            f(j,k) = f(j,k) - cforce*(g(j,k))/(fcor+zetamax)* & 
!                              (1.0 - (g(j,k))/(fcor+zetamax))*ffac* &
!                              msfun( x(j), y(k) ) 
         end if
      end do
      end do

      case default ! all forcing terms are zero in this version
      if (cforce==0) then
         f = f + zero
      else
         f = f + cforce*zero
      end if 
   end select

   ! make the forcing field explicitly periodic

   f(-1,0:ny-1) = f(nx-1,0:ny-1)
   f(nx,0:ny-1) = f(   0,0:ny-1)
   f(:,-1) = f(:,ny-1)
   f(:,ny) = f(:,   0)

   call tospec( f )
   call tospec( g )   

!-------------------------------------------------------------------------------
end subroutine add_forcing2
!===============================================================================


!===============================================================================
function msfun( x, y )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes mass sink at a point
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: x, y ! coordinates of point
!   real (kind=dbl_kind), intent(in)  :: nnn
   real (kind=dbl_kind) ::            msfun ! mass sink
!-------------------------------------------------------------------------------
! This version:  Mass sink ring  (called only if  icond==4)
!
! The mass sink ms at radius r from the center (xc,yc) is:
!            ( 0                0 <= r <= r1
!            ( (transition)     r1 <  r <  r2
!    ms(r) = ( m2               r2 <= r <= r3
!            ( (transition)     r3 <  r <  r4
!            ( 0                r4 <  r
! where the transitions are continuous and smooth (cubic Hermite splines).
! The outer band is perturbed with amplitude  amp*z2  in each of azimuthal
! wavenumbers  m = minwav, ..., maxwav  with random phase for each  m.  
! Default values of the parameters  xc, yc, r1, r2, r3, r4, z1, z2, minwav,
! maxwav,  and  pamp  are set above; these may be changed by reading them 
! from a file by read_icpars.
!-------------------------------------------------------------------------------

! local variables:

   real (kind=dbl_kind)    :: ms       ! radial part of z0fun
   real (kind=dbl_kind)    :: tdotmax  ! max theta-dot heating rate
   real (kind=dbl_kind)    :: mm2       ! mass sink peak magnitude 
   real (kind=dbl_kind)    :: mm1       ! mass sink peak magnitude 

! compute the location in polar coordinates (relative to the center)

   xx = x0 + modulo( x-xc-x0, xl ); yy = y0 + modulo( y-yc-y0, yl )
   r   = sqrt( xx**2 + yy**2 ); phi = atan2( yy, xx )

! definte mass sink magnitude m2

   tdotmax = (m2 * 28935.2) / (r3**2 - r2**2)
!   tdotmax = 3617500.0 / (r3**2 - r2**2)   ! K/s max heating possible
   mm2 = tdotmax / 30.0  
   mm1 = zero   ! automatically zero for now

! radial structure of initial vorticity and perturbation

   if ( r<=rm1 ) then
      ms = mm1  
   else if ( r<rm2 ) then
      ms = mm1*shape_fun((r-rm1)/(rm2-rm1))+mm2*shape_fun((rm2-r)/(rm2-rm1)) 
   else if ( r<=rm3 ) then
      ms = mm2
   else if ( r<rm4 ) then
      ms = mm2*shape_fun((r-rm3)/(rm4-rm3))
   else
      ms = zero
   end if 

! put the pieces together to form the initial vorticity

   msfun = ms

contains

!===============================================================================
function shape_fun( s )
!-------------------------------------------------------------------------------
! Purpose:
!    Internal function for cubic Hermite shape function (for transitions)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: s ! value in the interval [0,1]
   real (kind=dbl_kind) :: shape_fun ! array for field
!-------------------------------------------------------------------------------
   shape_fun = (1-s)**2*(1+2*s)
!-------------------------------------------------------------------------------
end function shape_fun
!===============================================================================

!-------------------------------------------------------------------------------
end function msfun
!===============================================================================



!===============================================================================
subroutine read_icpars( in_unit, echo_unit )
!-------------------------------------------------------------------------------
! Purpose:
!    Reads the initial condition parameters (from a file)
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
!    Each input line is scanned for all of the initial condition  parameters 
!    listed in this module.  They may be in any order but must appear one per 
!    line in the following format:
!
!       name value [optional comment]
!       # lines starting with "#" are comments and may appear anywhere
!       ! lines starting with "!" are comments and may appear anywhere
!
!      The input is read until the end of file is reached or a line starting
!      with the word "end" is encountered.  For example:
!
!         # initial condition parameters
!         dt     3600      time step for model run [in seconds]
!         end of initial condition parameters
!
!    Separate names and values by whitespace; put multiword strings in quotes.
!    Any parameter not read in retains the default value assigned above.
!
!    Length and time units for input values can be specified by lines starting 
!    with "lunits" (for length units:  km or m [default]) or "zunits" (for 
!    vorticity values:  f or 1/s [default]).  These units are used for 
!    subsequent values of initial condition parameters representing lengths 
!    and vorticities only (not combinations).  For example, the file:
!
!       lunits    km
!       xc        300
!       yc        300
!       zunits    f
!       z1        0.1
!       z2        5
!
!    specifies a vortex centered at (300km,300km) with vorticity (PV) levels
!    z1 = 0.1*f  and  z2 = 5*f  (where  f  is the Coriolis parameter).
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: iunit  ! local unit number for input
   character (len=40)      :: vname  ! for variable name
   character (len=132)     :: line   ! to hold one input line
   character (len=1)       :: cunits ! for one-character unit identifier
   integer (kind=int_kind) :: ios    ! iostatus
   real (kind=dbl_kind) :: &         ! length and time scales:
      lunits = one,        &         ! length units (in meters)
      zunits = one                   ! vorticity units (in inverse seconds)
   logical             :: echo_input ! echo the input file?

! process the arguments input to this routine

   if (present( in_unit )) then
      iunit = in_unit
   else
      iunit = 5
   end if
   echo_input = present( echo_unit )

! read the lines, scanning for the initial condition parameters

   readloop: do
      read ( unit=iunit, fmt="(a)", iostat=ios ) line
      if ( ios<0 ) then
         exit readloop
      else if ( ios>0 ) then
         write (*,*) " read_icpars:  cannot read input line"
         close ( unit=iunit )
         stop "read_icpars:  cannot read input line"
      else if ( echo_input ) then
         write (echo_unit,*) trim( line )
      end if
      if ( len_trim(line)==0 .or. line(1:1)=="#" .or. line(1:1)=="!" ) cycle
      read ( unit=line, fmt=* ) vname
      select case (vname)
         case ("end")
            return
! Parameters for initial Gaussian vortex (when  icond==1):
         case ("rmax")
            read ( unit=line, fmt=* ) vname, rmax
            rmax = lunits*rmax
         case ("vmax")
            read ( unit=line, fmt=* ) vname, vmax
         case ("imbalance")
            read ( unit=line, fmt=* ) vname, imbalance
! Parameters for initial PV ring (when  icond==4):
         case ("r1")
            read ( unit=line, fmt=* ) vname, r1
            r1 = lunits*r1
         case ("r2")
            read ( unit=line, fmt=* ) vname, r2
            r2 = lunits*r2
         case ("r3")
            read ( unit=line, fmt=* ) vname, r3
            r3 = lunits*r3
         case ("r4")
            read ( unit=line, fmt=* ) vname, r4
            r4 = lunits*r4
         case ("z1")
            read ( unit=line, fmt=* ) vname, z1
            z1 = zunits*z1
         case ("z2")
            read ( unit=line, fmt=* ) vname, z2
            z2 = zunits*z2
         case ("pamp")
            read ( unit=line, fmt=* ) vname, pamp
         case ("minwav")
            read ( unit=line, fmt=* ) vname, minwav
         case ("maxwav")
            read ( unit=line, fmt=* ) vname, maxwav
! Parameters for elliptical vortex (when  icond==4):
         case ("ze")
            read ( unit=line, fmt=* ) vname, ze
            ze = zunits*ze
         case ("ri")
            read ( unit=line, fmt=* ) vname, ri
            ri = lunits*ri
         case ("ro")
            read ( unit=line, fmt=* ) vname, ro
            ro = lunits*ro
         case ("ec")
            read ( unit=line, fmt=* ) vname, ec
         case ("lam")
            read ( unit=line, fmt=* ) vname, lam
! Parameters for linear vortex (when  icond==4):
         case ("r3l")
            read ( unit=line, fmt=* ) vname, r3l
            r3l = lunits*r3l
         case ("r4l")
            read ( unit=line, fmt=* ) vname, r4l
            r4l = lunits*r4l
         case ("zp1")
            read ( unit=line, fmt=* ) vname, zp1
            zp1 = zunits*zp1
! Parameters for initial mass sink ring :
         case ("rm1")
            read ( unit=line, fmt=* ) vname, rm1
            rm1 = lunits*rm1
         case ("rm2")
            read ( unit=line, fmt=* ) vname, rm2
            rm2 = lunits*rm2
         case ("rm3")
            read ( unit=line, fmt=* ) vname, rm3
            rm3 = lunits*rm3
         case ("rm4")
            read ( unit=line, fmt=* ) vname, rm4
            rm4 = lunits*rm4
         case ("m1")
            read ( unit=line, fmt=* ) vname, m1
         case ("m2")
            read ( unit=line, fmt=* ) vname, m2
! Parameters for both initial conditions:
         case ("xc")
            read ( unit=line, fmt=* ) vname, xc
            xc = lunits*xc
         case ("yc")
            read ( unit=line, fmt=* ) vname, yc
            yc = lunits*yc
! Units for length and vorticity (for scaling input values only):
         case ("zunits")
            read ( unit=line, fmt=* ) vname, cunits
            if ( cunits.eq.'f' ) then
               zunits = fcor
            else
               zunits = one
            end if
         case ("lunits")
            read ( unit=line, fmt=* ) vname, cunits
            if ( cunits.eq.'k' ) then
               lunits = 1000*one
            else
               lunits = one
            end if
         case default
            write (*,*) " read_icpars:  unrecognized variable on line:"
            write (*,*) line
      end select
   end do readloop

!-------------------------------------------------------------------------------
end subroutine read_icpars
!===============================================================================


!===============================================================================
subroutine write_icpars( out_unit )
!-------------------------------------------------------------------------------
! Purpose:
!    Writes the values of all initia condition parameters to a text file
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
   if (icond==1) then
      write (iunit,*) "Initial condition parameters (Gaussian vortex):"
      write (iunit,*) 'xc        = ', xc/km, " km"
      write (iunit,*) 'yc        = ', yc/km, " km"
      write (iunit,*) 'rmax      = ', rmax/km, " km"
      write (iunit,*) 'vmax      = ', vmax, " m/s"
      write (iunit,*) 'imbalance = ', imbalance
   else if (icond==4) then
      write (iunit,*) "Initial condition parameters (PV ring):"
      write (iunit,*) 'xc        = ', xc/km, " km"
      write (iunit,*) 'yc        = ', yc/km, " km"
      write (iunit,*) 'r1        = ', r1/km, " km"
      write (iunit,*) 'r2        = ', r2/km, " km"
      write (iunit,*) 'r3        = ', r3/km, " km"
      write (iunit,*) 'r4        = ', r4/km, " km"
      write (iunit,*) 'z1        = ', z1, " 1/s = ", z1/fcor, "*f"
      write (iunit,*) 'z2        = ', z2, " 1/s = ", z2/fcor, "*f"
      write (iunit,*) 'pamp      = ', pamp
      write (iunit,*) 'minwav    = ', minwav
      write (iunit,*) 'maxwav    = ', maxwav
      write (iunit,*) "Initial condition parameters (elliptical vortex):"
      write (iunit,*) 'ze        = ', ze, " 1/s"
      write (iunit,*) 'ri        = ', ri/km, " km"
      write (iunit,*) 'ro        = ', ro/km, " km"
      write (iunit,*) 'ec        = ', ec 
      write (iunit,*) 'lam        = ', lam
      write (iunit,*) "Initial condition parameters (linear vortex):"
      write (iunit,*) 'r3l        = ', r3l/km, " km"
      write (iunit,*) 'r4l        = ', r4l/km, " km"
      write (iunit,*) 'zp1        = ', zp1, " 1/s = ", z1/fcor, "*f"
      write (iunit,*) "Initial condition parameters (mass sink):"
      write (iunit,*) 'rm1        = ', rm1/km, " km"
      write (iunit,*) 'rm2        = ', rm2/km, " km"
      write (iunit,*) 'rm3        = ', rm3/km, " km"
      write (iunit,*) 'rm4        = ', rm4/km, " km"
      write (iunit,*) 'm1        = ', m1, " K/day = "
      write (iunit,*) 'm2        = ', m2, " K/day = "
   end if
   write (iunit,*) "--------------------------------------------------------"

!-------------------------------------------------------------------------------
end subroutine write_icpars
!===============================================================================


!-------------------------------------------------------------------------------
end module pswm_data
!===============================================================================
