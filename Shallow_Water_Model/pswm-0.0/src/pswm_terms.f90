!===============================================================================
module pswm_terms
!-------------------------------------------------------------------------------
! Purpose:
!    Primary code to compute the terms for PSWM:  Periodic Shallow Water Model
!
! Note:
!    This module does NOT use the module pswm_vars where the model variables
!    and forcing are stored at multiple time levels as arranged by sitpack.
!    Thus, routines in this module know nothing about multiple time levels:
!    all solution and forcing arrays are passed to these routines at one
!    time level only (through the argument list).
!
! See also:
!    sitpack_interface  (for principal routines which call these routines)
!-------------------------------------------------------------------------------

use kinds
use pswm_pars
use pswm_cons
use ss_ops
use pswm_data

implicit none   ! all module variables untyped by default
private         ! all module variables and routines hidden by default

public          &! routines called from other places (mainly sitpack_interface)
   setup_terms,         &! initializations needed for this module
   initialize_soln,     &! put initial solution into the model
   compute_terms_uv,    &! model terms (momentum form)
   compute_terms_zd,    &! model terms (vorticity/divergence form)
   compute_terms_qd,    &! model terms (potential vorticity/divergence form)
   solve_implicit_uv,   &! solve implicit problem (momentum form)
   solve_implicit_zd,   &! solve implicit problem (v/d form)
   solve_implicit_qd,   &! solve implicit problem (PV/d form)
   compute_diagnostics, &! compute scalar diagnostics
   unset_terms           ! to release allocated memory (for clean exit)

! module variables

real (kind=dbl_kind), dimension(:,:), allocatable, save :: & ! arrays for terms
      w1, w2, w3, w4, w5, w6, sponge
real (kind=dbl_kind) :: lambda ! constant in modified Helmholtz problem
real (kind=dbl_kind) :: tau_c  ! constant tau*cval for implicit problem
logical :: linear              ! use linearized equations? (set from ieq)

! internal parameters (not user-accessible--set at compile time only)
logical :: trace       = .false. ! trace routines in this module?
logical :: recompute_p = .false. ! recompute p from continuity for mass consv.

contains

!===============================================================================
subroutine setup_terms
!-------------------------------------------------------------------------------
! Purpose:
!    Performs initializations needed for this module
!-------------------------------------------------------------------------------

   real (kind=dbl_kind) :: pi, kfact ! 3.14... and the wavenumber factor

   if (trace) write (*,*) "Running setup_terms..."

! allocate the workspace needed for various model terms

   if (allocated(w1)) deallocate(w1); allocate(w1(-1:nx,-1:ny)); w1 = zero
   if (allocated(w2)) deallocate(w2); allocate(w2(-1:nx,-1:ny)); w2 = zero
   if (allocated(w3)) deallocate(w3); allocate(w3(-1:nx,-1:ny)); w3 = zero
   if (allocated(w4)) deallocate(w4); allocate(w4(-1:nx,-1:ny)); w4 = zero
   if (allocated(w5)) deallocate(w5); allocate(w5(-1:nx,-1:ny)); w5 = zero
   if (allocated(w6)) deallocate(w6); allocate(w6(-1:nx,-1:ny)); w6 = zero

! store the global parameters needed for spectral space operations

   call setup_ss_ops( mx, my, xl, yl )

! compute the dissipation coefficient  cdiss  and e-folding time  tdiss

   if ( idiss>0 ) then
!     pi = acos( -one ); kfact = ((2*pi*mx/xl)**2 + (2*pi*my/yl)**2)**idiss
      pi = acos( -one ); kfact = 2*pi*max( mx/xl, my/yl )
      kfact = (kfact**2)**idiss
   else if ( idiss==0 ) then
      kfact = one
   else
      write (*,*) "setup_terms:  idiss = ", idiss, "<0 not allowed"
      stop        "setup_terms:  idiss<0  not allowed"
   end if
   if ( tdiss>zero ) then
      cdiss = one/(kfact*tdiss*3600)
   else if ( cdiss>0 ) then
      tdiss = (one/(cdiss*kfact))/3600
   else
      cdiss = zero; tdiss = zero
   end if

! compute the sponge coefficient  csponge  and e-folding time  tsponge

   if ( isponge>0 ) then
!     pi = acos( -one ); kfact = ((2*pi*mx/xl)**2 + (2*pi*my/yl)**2)**isponge
      pi = acos( -one ); kfact = 2*pi*max( mx/xl, my/yl )
      kfact = (kfact**2)**isponge
   else if ( isponge==0 ) then
      kfact = one
   else
      write (*,*) "setup_terms:  isponge = ", isponge, "<0 not allowed"
      stop        "setup_terms:  isponge<0  not allowed"
   end if
   if ( tsponge>zero ) then
      csponge = one/(kfact*tsponge*3600)
   else if ( csponge>0 ) then
      tsponge = (one/(csponge*kfact))/3600
   else
      csponge = zero; tsponge = zero
   end if

! set up the sponge amplitude function (if needed)

   if ( lsponge==zero ) then
      csponge = zero; tsponge = zero
   end if
   if (csponge>zero ) call setup_sponge
   if (trace) write (*,*) "... finished setup_terms."

!-------------------------------------------------------------------------------
end subroutine setup_terms
!===============================================================================


!===============================================================================
subroutine setup_sponge
!-------------------------------------------------------------------------------
! Purpose:
!    Defines the amplitude function (in physical space) for the sponge layer
!
! Description:
!    The sponge function is defined on the transform grid and has values in
!    the range [0,1]; it multiplies the sponge coefficient term.  At each point
!    the value is based on the normalized distance  r  which is defined as the
!    distance to the nearest boundary, normalized by the corresponding length
!    scale as specified by the model paramter  lsponge  (see module pswm_pars).
!    Various sponge functions are included below and selected by  fsponge.
!    For each the value is 1 at the boundary (r=0) and tails off with length
!    scale 1 (some vanish beyond r=1, others decay with that scaling in r).
!-------------------------------------------------------------------------------

   real (kind=dbl_kind) :: xs, ys  ! sponge length scales in  x  and  y
   real (kind=dbl_kind) :: r       ! normalized length to nearest boundary
   integer (kind=int_kind) :: j, k ! indices for physical space points

   if (trace) write (*,*) "Running setup_sponge..."

! set the sponge length scales

   if ( lsponge>zero ) then
      xs = lsponge; ys = lsponge
   else if ( lsponge<zero ) then
      xs = abs(lsponge)*xl; ys = abs(lsponge)*yl
   else
      return
   end if

! allocate the space needed

   if ( allocated(sponge) ) deallocate(sponge); allocate(sponge(-1:nx,-1:ny)) 

! evaluate the sponge function based on normalized distance to nearest boundary

   do k=0,ny-1
   do j=0,nx-1
      r = min( (x(j)-x0)/xs, (x0+xl-x(j))/xs, (y(k)-y0)/ys, (y0+yl-y(k))/ys)
      select case (fsponge)

         case (1) ! linear "ramp" function:  goes to zero at r=1
            if ( zero<=r .and. r<one ) then
               sponge(j,k) = 1-r
            else
               sponge(j,k) = zero
            end if

         case (2) ! cubic Hermite:  goes to zero at r=1
            if ( zero<=r .and. r<one ) then
               sponge(j,k) = (1+2*r)*(1-r)**2 
            else
               sponge(j,k) = zero
            end if

         case (3) ! Bump function:  goes to zero at r=1 and smooth
            if ( zero==r ) then
               sponge(j,k) = one
            else if ( r<one ) then
               sponge(j,k) = exp( r*r/(r*r-one) )
            else
               sponge(j,k) = zero
            end if

         case (4) ! Gaussian:  e-folding length 1 
            sponge(j,k) = exp( -r*r )

         case default
            write (*,*) "setup_sponge:  fsponge = ", fsponge, " not allowed"
            stop        "setup_sponge:  fsponge value not allowed"

      end select
   end do
   end do

! enforce explicit periodicity

   sponge(-1,0:ny-1) = sponge(nx-1,0:ny-1)
   sponge(nx,0:ny-1) = sponge(   0,0:ny-1)
   sponge(:,-1) = sponge(:,ny-1)
   sponge(:,ny) = sponge(:,   0)

   if (trace) write (*,*) "... finished setup_sponge."

!-------------------------------------------------------------------------------
end subroutine setup_sponge
!===============================================================================


!===============================================================================
subroutine initialize_soln( u, v, p, d, z, q )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes the model terms (momentum form)
!
! Arguments:
   real (kind=dbl_kind), dimension(:,:), intent(out) :: &
      u, v, p, d, z, q ! model fields (at initial time level)
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running initialize_soln..."
   linear = ieq<0
   select case (icond)

      case (1) ! direct initialization -- momentum form
         call put_initial( 'u', u ); call tospec( u )
         call put_initial( 'v', v ); call tospec( v )
         call put_initial( 'p', p ); 
!         write(*,*) "bef spec", p(10,10); 
         call tospec( p ); 
!         write(*,*) "after spec", p(10,10);
         call xderiv( u, w1 ) ! u_x
         call yderiv( u, w2 ) ! u_y
         call xderiv( v, w3 ) ! v_x
         call yderiv( v, w4 ) ! v_y
         d = w1 + w4; z = w3 - w2; q = (fcor + z)/(one + p/cval)
!         stop

      case (2) ! direct initialization -- vorticity/divergence form
         call put_initial( 'z', z ); call tospec( z )
         call put_initial( 'd', d ); call tospec( d )
         call put_initial( 'p', p ); call tospec( p )
         call mhsolve( d, w1 ) ! -chi
         call mhsolve( z, w3 ) ! -psi
         call yderiv( w1, w2 ) ! -chi_y
         call xderiv( w1, w1 ) ! -chi_x
         call yderiv( w3, w4 ) ! -psi_y
         call xderiv( w3, w3 ) ! -psi_x
         u = -w1 + w4; v = -w2 - w3; q = (fcor + z)/(one + p/cval)

      case (3) ! direct initialization -- potential vorticity/divergence form
         call put_initial( 'q', q ); call tospec( q )
         call put_initial( 'd', d ); call tospec( d )
         call put_initial( 'p', p ); call tospec( p )
         z = (one + p/cval)*q - fcor
         call mhsolve( d, w1 ) ! -chi
         call mhsolve( z, w3 ) ! -psi
         call yderiv( w1, w2 ) ! -chi_y
         call xderiv( w1, w1 ) ! -chi_x
         call yderiv( w3, w4 ) ! -psi_y
         call xderiv( w3, w3 ) ! -psi_x
         u = -w1 + w4; v = -w2 - w3

      case (4) ! initialization using nonlinear balance (nondivergent)
         d = zero
         call put_initial( 'z', z ); call tospec( z )
         call mhsolve( z, w1 ) ! -psi
         call yderiv( w1, u ); ! -psi_y =  u
         call xderiv( w1, v ); ! -psi_x = -v
         call xderiv(  v, w1 ) ! -psi_xx
         call yderiv(  v, w2 ) ! -psi_xy
         call yderiv(  u, w3 ) ! -psi_yy
         v = -v
         call tophys( w1 ); call tophys( w2 ); call tophys( w3 )
         w1 = w1*w3 - w2*w2 ! psi_xx * psi_yy - (psi_xy)^2
         call tospec( w1 ); p = -(fcor*z + 2*w1)/cval; call mhsolve( p, p )
         call tophys( z); call tophys(p)
         q = (fcor + z)/(one + p/cval)
         call tospec( z); call tospec(p); call tospec(q)

   case default ! all initial fields vanish
      u = zero; v = zero; p = zero; d = zero; z = zero; q = zero

   end select
   if (trace) write (*,*) "... finished initialize_soln."

!-------------------------------------------------------------------------------
end subroutine initialize_soln
!===============================================================================


!===============================================================================
subroutine compute_terms_uv( kterm, cforce, u, v, p, d, z, uf, vf, pf )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes the model terms (momentum form)
!
! Arguments:
   integer (kind=int_kind), intent(in) :: kterm  ! specifies terms to compute
   real    (kind=dbl_kind), intent(in) :: cforce ! forcing coefficient
   real    (kind=dbl_kind), dimension(:,:), intent(in) :: &
      u, v, p, d, z ! input fields (at one time level)
   real    (kind=dbl_kind), dimension(:,:), intent(inout) :: &
      uf, vf, pf ! corresponding forcing from explicit terms
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running compute_terms_uv..."

! decide which terms we're computing

   linear = ieq<0

!------------------------------------------------------------------
   if (kterm==iterm_fast) then ! fast terms:  gravity-wave (linear)
!------------------------------------------------------------------

      call xderiv( p, w1 ) ! p_x
      call yderiv( p, w2 ) ! p_y
      if (cforce==zero) then
         uf = -cval*w1
         vf = -cval*w2
         pf = -cval*d
      else
         uf = uf - (cforce*cval)*w1
         vf = vf - (cforce*cval)*w2
         pf = pf - (cforce*cval)*d
      end if

!------------------------------------------------------------------------------
   else if (kterm==iterm_slow) then ! slow terms:  Coriolis, nonlinear, forcing
!------------------------------------------------------------------------------

      if ( linear ) then
         if (cforce==zero) then
            uf = +fcor*v
            vf = -fcor*u
            pf =  zero
         else
            uf = +(cforce*fcor)*v
            vf = -(cforce*fcor)*u
         end if
      else
!        transform to physical space:  eta, p, u, v
         w1 = z; w3 = p; w5 = u; w6 = v
         call tophys( w1 ); w1 = fcor + w1 ! eta
         call tophys( w3 ); if (iblow==1) call check_for_blowup( w3 )
         call tophys( w5 )
         call tophys( w6 )

!        compute nonlinear terms in physical space:
         w2 = w6*w1; w1 = w5*w1 ! v*eta, u*eta
         w4 = w6*w3; w3 = w5*w3 ! v*p  , u*p
         w5 = (w5*w5 + w6*w6)/2 ! K (kinetic energy)

!        transform nonlinear terms back to physical space
         call tospec( w1 ) ! u*eta
         call tospec( w2 ) ! v*eta
         call tospec( w3 ) ! u*p
         call tospec( w4 ) ! u*p
         call tospec( w5 ) ! K

!        forcing for u equation:  F_u + v*eta - K_x
!        forcing for v equation:  F_u - u*eta - K_y
!        forcing for p equation:  F_p - (u*p)_x - (v*p)_y
         call yderiv( w5, w6 ) ! K_y
         call xderiv( w5, w5 ) ! K_x
         call xderiv( w3, w3 ) ! (u*p)_x
         call yderiv( w4, w4 ) ! (v*p)_y
         if (cforce==zero) then
            uf =   w2 - w5
            vf = -(w1 + w6)
            pf = -(w3 + w4) 
         else
            uf = uf + cforce*(w2 - w5)
            vf = vf - cforce*(w1 + w6)
            pf = pf - cforce*(w3 + w4)
         end if
      end if

! add the specified forcing (if any)

      if (iforce/=0) then
         call add_forcing( 'u', cforce, uf )
         call add_forcing( 'v', cforce, vf )
         call add_forcing( 'p', cforce, pf )
      end if

!------------------------------------------------------
   else if (kterm==iterm_diss) then ! dissipation terms
!------------------------------------------------------

      if (cdiss/=zero) then
         call nlap( u, w1, idiss )
         call nlap( v, w2, idiss )
         call nlap( p, w3, idiss )
      end if
      if (csponge/=zero) then
         call nlap( u, w4, isponge )
         call nlap( v, w5, isponge )
         call nlap( p, w6, isponge )
         call tophys( w4 ); w4 = sponge*w4; call tospec( w4 )
         call tophys( w5 ); w5 = sponge*w5; call tospec( w5 )
         call tophys( w6 ); w6 = sponge*w6; call tospec( w6 )
      end if
      if (cforce==zero) then
         uf = -(cdiss*w1 + csponge*w4)
         vf = -(cdiss*w2 + csponge*w5)
         pf = -(cdiss*w3 + csponge*w6)
      else
         uf = uf - cforce*(cdiss*w1 + csponge*w4)
         vf = vf - cforce*(cdiss*w2 + csponge*w5)
         pf = pf - cforce*(cdiss*w3 + csponge*w6)
      end if

!---------------------------------------
   else
      write (*,*) "compute_terms_uv:  kterm = ", kterm, " out of range"
      stop        "compute_terms_uv:  kterm out of range"
   end if
   if (trace) write (*,*) "... finished compute_terms_uv."

!-------------------------------------------------------------------------------
end subroutine compute_terms_uv
!===============================================================================


!===============================================================================
subroutine compute_terms_zd( kterm, cforce, z, d, p, u, v, zf, df, pf )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes the model terms (vorticity/divergence form)
!
! Arguments:
   integer (kind=int_kind), intent(in) :: kterm  ! specifies terms to compute
   real    (kind=dbl_kind), intent(in) :: cforce ! forcing coefficient
   real    (kind=dbl_kind), dimension(:,:), intent(in) :: &
      z, d, p, u, v ! input fields (at one time level)
   real    (kind=dbl_kind), dimension(:,:), intent(inout) :: &
      zf, df, pf ! corresponding forcing from explicit terms
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running compute_terms_zd..."

! decide which terms we're computing

   linear = ieq<0

!------------------------------------------------------------------
   if (kterm==iterm_fast) then ! fast terms:  gravity-wave (linear)
!------------------------------------------------------------------

      call nlap( p, w1 ) ! -del^2(p)
      if (cforce==zero) then
         df =  cval*w1
         pf = -cval*d
      else
         df = df + (cforce*cval)*w1
         pf = pf - (cforce*cval)*d
      end if

!------------------------------------------------------------------------------
   else if (kterm==iterm_slow) then ! slow terms:  Coriolis, nonlinear, forcing
!------------------------------------------------------------------------------

      if ( linear ) then
         if (cforce==zero) then
            zf = -fcor*d ! F_z - div(f*vector_v) = F_z - f*d
            df =  fcor*z ! F_d + rot(f*vector_v) = F_d + f*z
            pf =  zero   ! F_p
         else
            zf = zf - (cforce*fcor)*d
            df = df + (cforce*fcor)*z
         end if
      else
!        transform to physical space:  eta, p, u, v
         w1 = z; w3 = p; w5 = u; w6 = v
         call tophys( w1 ); w1 = fcor + w1 ! eta
         call tophys( w3 ); if (iblow==1) call check_for_blowup( w3 )
         call tophys( w5 )
         call tophys( w6 )

!        compute nonlinear terms in physical space:
         w2 = w6*w1; w1 = w5*w1 ! v*eta, u*eta
         w4 = w6*w3; w3 = w5*w3 ! v*p  , u*p
         w5 = (w5*w5 + w6*w6)/2 ! K (kinetic energy)

!        transform nonlinear terms back to physical space
         call tospec( w1 ) ! u*eta
         call tospec( w2 ) ! v*eta
         call tospec( w3 ) ! u*p
         call tospec( w4 ) ! u*p
         call tospec( w5 ) ! K

!        forcing for p equation:  F_p - div(p*vector_v)
         call xderiv( w3, w3 ) ! (u*p)_x
         call yderiv( w4, w4 ) ! (v*p)_y
         if (cforce==zero) then
            pf = -(w3 + w4) 
         else
            pf = pf - cforce*(w3 + w4)
         end if

!        forcing for z equation:  F_z - (u*eta)_x - (v*eta)_y
!        forcing for d equation:  F_d + (v*eta)_x - (u*eta)_y - del^2(K)
         call yderiv( w1, w3 ) ! (u*eta)_y
         call xderiv( w1, w1 ) ! (u*eta)_x
         w4 = w2;                          !  v*eta
         call yderiv( w2, w4 ) ! (v*eta)_y
         call xderiv( w2, w2 ) ! (v*eta)_x
         call   nlap( w5, w5 ) ! -del^2(K)
         if (cforce==zero) then
            zf = -(w1 + w4)
            df =  (w2 - w3) + w5
         else
            zf = zf - cforce*(w1 + w4)
            df = df + cforce*((w2 - w3) + w5)
         end if
      end if

! add the specified forcing (if any)

! temp code PV
      call tophys( z ); call tophys( p)
      w6 = (fcor + z)/(one + p/cval)
      call tospec( z ); call tospec( p); call tospec( w6)

      if (iforce/=0) then
! Comment in for normal forcing
!         call add_forcing( 'z', cforce, zf )
!         call add_forcing( 'd', cforce, df )
!         call add_forcing( 'p', cforce, pf )
! Trick to do logistic forcing
         call add_forcing2( 'p', cforce, pf, w6 )
      end if

!------------------------------------------------------
   else if (kterm==iterm_diss) then ! dissipation terms
!------------------------------------------------------

      if (cdiss/=zero) then
         call nlap( z, w1, idiss )
         call nlap( d, w2, idiss )
         call nlap( p, w3, idiss )
      end if
      if (csponge/=zero) then
         call nlap( z, w4, isponge )
         call nlap( d, w5, isponge )
         call nlap( p, w6, isponge )
         call tophys( w4 ); w4 = sponge*w4; call tospec( w4 )
         call tophys( w5 ); w5 = sponge*w5; call tospec( w5 )
         call tophys( w6 ); w6 = sponge*w6; call tospec( w6 )
      end if
      if (cforce==zero) then
         zf = -(cdiss*w1 + csponge*w4 + 0.000001*z)
         df = -(cdiss*w2 + csponge*w5 + 0.000001*d)
!         zf = -(cdiss*w1 + csponge*w4)
!         df = -(cdiss*w2 + csponge*w5)
         pf = -(cdiss*w3 + csponge*w6)
      else
         zf = zf - cforce*(cdiss*w1 + csponge*w4 + 0.000001*z)
         df = df - cforce*(cdiss*w2 + csponge*w5 + 0.000001*d)
!         zf = zf - cforce*(cdiss*w1 + csponge*w4)
!         df = df - cforce*(cdiss*w2 + csponge*w5)
         pf = pf - cforce*(cdiss*w3 + csponge*w6)
      end if

!---------------------------------------
   else
      write (*,*) "compute_terms_zd:  kterm = ", kterm, " out of range"
      stop        "compute_terms_zd:  kterm out of range"
   end if
   if (trace) write (*,*) "... finished compute_terms_zd."

!-------------------------------------------------------------------------------
end subroutine compute_terms_zd
!===============================================================================


!===============================================================================
subroutine compute_terms_qd( kterm, cforce, q, d, p, u, v, z, qf, df, pf )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes the model terms (potential vorticity/divergence form)
!
! Arguments:
   integer (kind=int_kind), intent(in) :: kterm  ! specifies terms to compute
   real    (kind=dbl_kind), intent(in) :: cforce ! forcing coefficient
   real    (kind=dbl_kind), dimension(:,:), intent(in) :: &
      q, d, p, u, v, z ! input fields (at one time level)
   real    (kind=dbl_kind), dimension(:,:), intent(inout) :: &
      qf, df, pf ! corresponding forcing from explicit terms
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running compute_terms_qd..."

! decide which terms we're computing

   linear = ieq<0

!------------------------------------------------------------------
   if (kterm==iterm_fast) then ! fast terms:  gravity-wave (linear)
!------------------------------------------------------------------

      call nlap( p, w1 ) ! -del^2(p)
      if (cforce==zero) then
         qf = zero
         df =  cval*w1
         pf = -cval*d
      else
         df = df + (cforce*cval)*w1
         pf = pf - (cforce*cval)*d
      end if

!------------------------------------------------------------------------------
   else if (kterm==iterm_slow) then ! slow terms:  Coriolis, nonlinear, forcing
!------------------------------------------------------------------------------

      if ( linear ) then
         if (cforce==zero) then
            qf = zero
            df = fcor*z
            pf = zero
         else
            df = df + (cforce*fcor)*z
         end if
      else
!        transform u and v to physical space:
         w5 = u; call tophys( w5 )
         w6 = v; call tophys( w6 )

!        forcing for q equation:  F_q - (u*q_x + v*q_y)
         call xderiv( q, w1 ) ! q_x
         call yderiv( q, w2 ) ! q_y
         call tophys( w1 );
         call tophys( w2 )
         w1 = w5*w1; call tospec( w1 ) ! u*q_x
         w2 = w6*w2; call tospec( w2 ) ! v*q_x
         w1 = w1 + w2 ! u*q_x + v*q_y

!        forcing for d equation:  F_d + (v*eta)_x - (u*eta)_y - del^2(K)
         w2 = z; call tophys( w2 ); w2 = fcor + w2 ! eta
         w3 = w6*w2; w2 = w5*w2 ! v*eta, u*eta
         w4 = (w5*w5 + w6*w6)/2 ! K (kinetic energy)
         call tospec( w2 )
         call tospec( w3 )
         call tospec( w4 )
         call xderiv( w3, w3 ) ! (v*eta)_x
         call yderiv( w2, w2 ) ! (u*eta)_y
         call   nlap( w4, w4 ) ! -del^2(K)
         w2 = w3 - w2 + w4 ! F_z + (v*eta)_x - (u*eta)_y - del^2(K)

!        forcing for p equation:  F_p - (u*p)_x - (v*p)_y
         w3 = p; call tophys( w3 ); if (iblow==1) call check_for_blowup( w3 )
         w4 = w6*w3; w3 = w5*w3 ! v*p  , u*p
         call tospec( w3 )
         call tospec( w4 )
         call xderiv( w3, w3 ) ! (u*p)_x
         call yderiv( w4, w4 ) ! (v*p)_y

!        store or add the forcing:
         if (cforce==zero) then
            qf = -w1
            df =  w2
            pf = -(w3 + w4)
         else
            qf = qf - cforce*w1
            df = df + cforce*w2
            pf = pf - cforce*(w3 + w4)
         end if
      end if

! add the specified forcing (if any)

      if (iforce/=0) then
         call add_forcing( 'q', cforce, qf )
         call add_forcing( 'd', cforce, df )
         call add_forcing( 'p', cforce, pf )
      end if

!------------------------------------------------------
   else if (kterm==iterm_diss) then ! dissipation terms
!------------------------------------------------------

      if (cdiss/=zero) then
         call nlap( q, w1, idiss )
         call nlap( d, w2, idiss )
         call nlap( p, w3, idiss )
      end if
      if (csponge/=zero) then
         call nlap( q, w4, isponge )
         call nlap( d, w5, isponge )
         call nlap( p, w6, isponge )
         call tophys( w4 ); w4 = sponge*w4; call tospec( w4 )
         call tophys( w5 ); w5 = sponge*w5; call tospec( w5 )
         call tophys( w6 ); w6 = sponge*w6; call tospec( w6 )
      end if
      if (cforce==zero) then
         qf = -(cdiss*w1 + csponge*w4)
         df = -(cdiss*w2 + csponge*w5)
         pf = -(cdiss*w3 + csponge*w6)
      else
         qf = qf - cforce*(cdiss*w1 + csponge*w4)
         df = df - cforce*(cdiss*w2 + csponge*w5)
         pf = pf - cforce*(cdiss*w3 + csponge*w6)
      end if

!---------------------------------------
   else
      write (*,*) "compute_terms_qd:  kterm = ", kterm, " out of range"
      stop        "compute_terms_qd:  kterm out of range"
   end if
   if (trace) write (*,*) "... finished compute_terms_qd."

!-------------------------------------------------------------------------------
end subroutine compute_terms_qd
!===============================================================================


!===============================================================================
subroutine solve_implicit_uv( tau, uf, vf, pf, u, v, p, d, z, q )
!-------------------------------------------------------------------------------
! Purpose:
!    Solves the implicit problem (momentum form)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: tau ! scaled time step
   real (kind=dbl_kind), dimension(:,:), intent(in) :: &
      uf, vf, pf ! forcing for u, v, and p equations
   real (kind=dbl_kind), dimension(:,:), intent(out) :: &
      u, v, p, d, z, q ! output fields
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running solve_implicit_uv..."
   if (tau==zero) then ! fully explicit case
      u = uf; v = vf; p = pf;
   else

!   form the right-hand side MH equation:  G = P - tau*c*(U_x + V_y)

      tau_c = tau*cval; lambda = one/tau_c**2
      call xderiv( uf, w1 ) ! U_x
      call yderiv( vf, w2 ) ! V_y
      p = lambda*(pf - tau_c*(w1 + w2)) ! G

!   solve the modified Helmholtz (MH) equation for p

      call mhsolve( p, p, lambda )

!   compute the corresponding velocity field

      call xderiv( p, w1 ) ! p_x
      call yderiv( p, w2 ) ! p_y
      u = uf - tau_c*w1; v = vf - tau_c*w2

!   recompute  p  from the continuity equation to ensure conservation of mass

      if ( recompute_p ) then
         call xderiv( u, w1 ) ! u_x
         call yderiv( v, w2 ) ! v_y
         p = pf - tau_c*(w1 + w2)
      end if
   end if

! compute the corresponding divergence, vorticity, and potential vorticity

   call xderiv( u, w1 ) ! u_x
   call yderiv( u, w2 ) ! u_y
   call xderiv( v, w3 ) ! v_x
   call yderiv( v, w4 ) ! v_y
   d = w1 + w4; z = w3 - w2; q = (fcor + z)/(one + p/cval)
   if (trace) write (*,*) "... finished solve_implicit_uv."

!-------------------------------------------------------------------------------
end subroutine solve_implicit_uv
!===============================================================================


!===============================================================================
subroutine solve_implicit_zd( tau, zf, df, pf, u, v, p, d, z, q )
!-------------------------------------------------------------------------------
! Purpose:
!    Solves the implicit problem (vorticity/divergence form)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: tau ! scaled time step
   real (kind=dbl_kind), dimension(:,:), intent(in) :: &
      zf, df, pf ! forcing for z, d, and p equations
   real (kind=dbl_kind), dimension(:,:), intent(out) :: &
      u, v, p, d, z, q ! output fields
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running solve_implicit_zd..."

! vorticity is treated explicitly

   z = zf

   if (tau==zero) then ! fully explicit case
      d = df; p = pf;
   else

!   form the right-hand side MH equation:  G = P - tau*c*D

      tau_c = tau*cval; lambda = one/tau_c**2
      p = lambda*(pf - tau_c*df) ! G

!   solve the modified Helmholtz (MH) equation for p

      call mhsolve( p, p, lambda )

!   compute the corresponding divergence field

      call nlap( p, w1 ) ! -del^2(p)
      d = df + tau_c*w1

!   recompute  p  from the continuity equation to ensure conservation of mass

      if ( recompute_p ) p = pf - tau_c*d
   end if

! compute the corresponding velocity and potential vorticity

   call mhsolve( d, w1 ) ! -chi
   call mhsolve( z, w3 ) ! -psi
   call yderiv( w1, w2 ) ! -chi_y
   call xderiv( w1, w1 ) ! -chi_x
   call yderiv( w3, w4 ) ! -psi_y
   call xderiv( w3, w3 ) ! -psi_x
   u = -w1 + w4; v = -w2 - w3; q = (fcor + z)/(one + p/cval)
   call tophys( z); call tophys( p)
   q = (fcor + z)/(one + p/cval)
   call tospec( z); call tospec( p); call tospec( q)
   if (trace) write (*,*) "... finished solve_implicit_zd."

!-------------------------------------------------------------------------------
end subroutine solve_implicit_zd
!===============================================================================


!===============================================================================
subroutine solve_implicit_qd( tau, qf, df, pf, u, v, p, d, z, q )
!-------------------------------------------------------------------------------
! Purpose:
!    Solves the implicit problem (potential vorticity/divergence form)
!
! Arguments:
   real (kind=dbl_kind), intent(in) :: tau ! scaled time step
   real (kind=dbl_kind), dimension(:,:), intent(in) :: &
      qf, df, pf ! forcing for q, d, and p equations
   real (kind=dbl_kind), dimension(:,:), intent(out) :: &
      u, v, p, d, z, q ! output fields
!-------------------------------------------------------------------------------

   if (trace) write (*,*) "Running solve_implicit_qd..."

! potential vorticity is treated explicitly

   q = qf

   if (tau==zero) then ! fully explicit case
      d = df; p = pf;
   else

!   form the right-hand side MH equation:  G = P - tau*c*D

      tau_c = tau*cval; lambda = one/tau_c**2
      p = lambda*(pf - tau_c*df) ! G

!   solve the modified Helmholtz (MH) equation for p

      call mhsolve( p, p, lambda )

!   compute the corresponding divergence field

      call nlap( p, w1 ) ! -del^2(p)
      d = df + tau_c*w1

!   recompute  p  from the continuity equation to ensure conservation of mass

      if ( recompute_p ) p = pf - tau_c*d
   end if

! compute the corresponding vorticity and velocity

   z = (one + p/cval)*q - fcor
   call mhsolve( d, w1 ) ! -chi
   call mhsolve( z, w3 ) ! -psi
   call yderiv( w1, w2 ) ! -chi_y
   call xderiv( w1, w1 ) ! -chi_x
   call yderiv( w3, w4 ) ! -psi_y
   call xderiv( w3, w3 ) ! -psi_x
   u = -w1 + w4; v = -w2 - w3
   if (trace) write (*,*) "... finished solve_implicit_qd."

!-------------------------------------------------------------------------------
end subroutine solve_implicit_qd
!===============================================================================


!===============================================================================
subroutine compute_diagnostics( u, v, p, d, z, q, diagnostics )
!-------------------------------------------------------------------------------
! Purpose:
!    Computes the diagnostics (scalars associated with model variables)
!
! Arguments:
   real (kind=dbl_kind), dimension(:,:), intent(in) :: &
      u, v, p, z, d, q ! input fields (at one time level)
   real (kind=dbl_kind), intent(out) :: diagnostics(3,10) ! diagnostics:
!    for each of the fields indexed by  ifld  as follows:
!                  ifld= 1:  u
!                  ifld= 2:  v
!                  ifld= 3:  p
!                  ifld= 4:  d = delta (divergence)
!                  ifld= 5:  z = zeta      (relative  vorticity)
!                  ifld= 6:  q = e/(1+p/cval) (potential vorticity)
!                  ifld= 7:  e = eta = fcor+z (absolute  vorticity)
!                  ifld= 8:  Z = z*z/2     (enstrophy)
!                  ifld= 9:  K = (u*u+v*v)*(cval+p)/2 (kinetic energy)
!                  ifld=10:  A = cval*p*p/2 (available potential energy)
!    diagnostics(1,ifld):  minimum pointwise value
!    diagnostics(2,ifld):  maximum pointwise value
!    diagnostics(3,ifld):  domain average value
!
!   Note:
!      The domain average of a field is its integral over the model domain
!      divided by the area of the model domain, and the energies are missing
!      the factor  cval/g  (since  g  isn't needed elsewhere in the model).
!-------------------------------------------------------------------------------

   integer (kind=int_kind) :: ifld  ! field index

   if (trace) write (*,*) "Running compute_diagnostics..."

! transform the model variables to physical space

   w1 = u; call tophys( w1 )
   w2 = v; call tophys( w2 )
   w3 = p; call tophys( w3 )
   w4 = d; call tophys( w4 )
   w5 = z; call tophys( w5 )
   w6 = q; call tophys( w6 )

! diagnostics for u:

   ifld = 1;
   diagnostics(1,ifld) = minval( w1 ); diagnostics(2,ifld) = maxval( w1 );
   diagnostics(3,ifld) = sum( w1(1:,1:) )/(xl*yl);

! diagnostics for v:

   ifld = 2;
   diagnostics(1,ifld) = minval( w2 ); diagnostics(2,ifld) = maxval( w2 );
   diagnostics(3,ifld) = sum( w2(1:,1:) )/(xl*yl);

! diagnostics for p:

   ifld = 3; if ( iblow==2) call check_for_blowup( w3 )
   diagnostics(1,ifld) = minval( w3 ); diagnostics(2,ifld) = maxval( w3 );
   diagnostics(3,ifld) = sum( w3(1:,1:) )/(xl*yl);

! diagnostics for d:

   ifld = 4;
   diagnostics(1,ifld) = minval( w4 ); diagnostics(2,ifld) = maxval( w4 );
   diagnostics(3,ifld) = sum( w4(1:,1:) )/(xl*yl);

! diagnostics for z:

   ifld = 5;
   diagnostics(1,ifld) = minval( w5 ); diagnostics(2,ifld) = maxval( w5 );
   diagnostics(3,ifld) = sum( w5(1:,1:) )/(xl*yl);

! diagnostics for q:

   ifld = 6;
   diagnostics(1,ifld) = minval( w6 ); diagnostics(2,ifld) = maxval( w6 );
   diagnostics(3,ifld) = sum( w6(1:,1:) )/(xl*yl);

! diagnostics for e = eta (absolute vorticity)

   ifld = 7;
   diagnostics(:,ifld) = diagnostics(:,5) + fcor

! diagnostics for Z (enstrophy)

   ifld = 8; w5 = w5*w5/2
   diagnostics(1,ifld) = minval( w5 ); diagnostics(2,ifld) = maxval( w5 );
   diagnostics(3,ifld) = sum( w5(1:,1:) )/(xl*yl);

! diagnostics for K (kinetic energy)

   ifld = 9; w1 = (w1*w1 + w2*w2)*(cval + w3)/2
   diagnostics(1,ifld) = minval( w1 ); diagnostics(2,ifld) = maxval( w1 );
   diagnostics(3,ifld) = sum( w1(1:,1:) )/(xl*yl);

! diagnostics for A (available potential energy)

   ifld = 10; w3 = cval*w3*w3/2
   diagnostics(1,ifld) = minval( w3 ); diagnostics(2,ifld) = maxval( w3 );
   diagnostics(3,ifld) = sum( w3(1:,1:) )/(xl*yl);

   if (trace) write (*,*) "... finished compute_diagnostics."

!-------------------------------------------------------------------------------
end subroutine compute_diagnostics
!===============================================================================


!===============================================================================
subroutine check_for_blowup( p )
!-------------------------------------------------------------------------------
! Purpose:
!    Internal routine to check for bottom out or blowup
!
! Argument:
   real (kind=dbl_kind), intent(in) :: p(:,:) ! scaled deviation phi (physical)
!-------------------------------------------------------------------------------

! check for  p  too small or too large

   if (cval+minval(p)<=zero) then 
      write (*,*) "Solution has bottomed out--quitting"
      stop
   else if (blowup>zero .and. cval+maxval(p)>blowup*cval) then
      write (*,*) "Solution has blown up--quitting"
      stop
   end if

!-------------------------------------------------------------------------------
end subroutine check_for_blowup
!===============================================================================


!===============================================================================
subroutine unset_terms
!-------------------------------------------------------------------------------
! Purpose:
!    Releases the memory allocated by setup_terms
!-------------------------------------------------------------------------------

! release the workspace

   if (allocated(sponge)) deallocate(sponge)
   if (allocated(w6)) deallocate(w6)
   if (allocated(w5)) deallocate(w5)
   if (allocated(w4)) deallocate(w4)
   if (allocated(w3)) deallocate(w3)
   if (allocated(w2)) deallocate(w2)
   if (allocated(w1)) deallocate(w1)

!-------------------------------------------------------------------------------
end subroutine unset_terms
!===============================================================================


!-------------------------------------------------------------------------------
end module pswm_terms
!===============================================================================
