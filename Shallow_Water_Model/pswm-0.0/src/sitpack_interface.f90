!===============================================================================
module sitpack_interface
!-------------------------------------------------------------------------------
! Purpose: 
!    Interface between the semi-implicit time differencing module sitpack
!    and the model code.  Includes the routines called by sitpack to
!    compute and combine the model terms and solve the implicit equations.  
!
! Note: 
!    These routines are mainly an interface to the routines in the module
!    pswm_terms, which compute the terms in the model equations but know 
!    nothing about the time integration or multiple time levels:  they access
!    the solution and forcing as two-dimensional arrays at one time level
!    passed through the argument list rather than from the module pswm_vars.
!
! See also:  
!    sitpack    (for details of the semi-implicit method)
!    pswm_terms (contains the code for computing the model terms)
!-------------------------------------------------------------------------------

use kinds
use pswm_pars
use pswm_cons
use pswm_vars
use pswm_terms

implicit none
private         ! all entry points and variables are local by default

!-------------------------------------------------------------------------------

public          &! routines defining the interface to module  sitpack:
   compute_model_terms, &
   compute_lincomb_sol, &
   compute_lincomb_for, &
   solve_implicit_eqns, &
   put_model_solution

logical, parameter :: trace = .false. ! trace routines in this module?

contains

!===============================================================================
subroutine compute_model_terms( kterm, spoint, sitime, fpoint, cforce )
!-------------------------------------------------------------------------------
! Purpose:  
!    Compute model terms  F_kterm  based on solution  spoint  at time  sitime
!       If cforce==0:  PUT the result          F_kterm  into forcing  fpoint
!       If cforce/=0:  ADD the product  cforce*F_kterm    to forcing  fpoint
!
! Arguments:
   integer (kind=int_kind), intent(in) :: kterm  ! index of term
   integer (kind=int_kind), intent(in) :: spoint ! index of solution copy
   real    (kind=dbl_kind), intent(in) :: sitime ! model time from sitpack
   integer (kind=int_kind), intent(in) :: fpoint ! index of forcing copy
   real    (kind=dbl_kind), intent(in) :: cforce ! forcing coefficient
!
! Note:
!    This interface is required by the module sitpack
!-------------------------------------------------------------------------------

   if (trace) print *, "Starting compute_model_terms ..."
   if (trace) print *, "spoint = ", spoint, "  fpoint = ", fpoint

! compute the desired terms (according to the selected equation form)

   select case (abs(ieq))
      case (1) ! momentum form
         call compute_terms_uv( kterm,  cforce, &
                      u(:,:,spoint),  v(:,:,spoint),  p(:,:,spoint),  &
                      d(:,:,spoint),  z(:,:,spoint),                  &
                     uf(:,:,fpoint), vf(:,:,fpoint), pf(:,:,fpoint) )
      case (2) ! vorticity/divergence form
         call compute_terms_zd( kterm,  cforce, &
                      z(:,:,spoint),  d(:,:,spoint),  p(:,:,spoint),  &
                      u(:,:,spoint),  v(:,:,spoint),                  &
                     zf(:,:,fpoint), df(:,:,fpoint), pf(:,:,fpoint) )
      case (3) ! potential vorticity/divergence form
         call compute_terms_qd( kterm,  cforce, &
                      q(:,:,spoint),  d(:,:,spoint),  p(:,:,spoint),  &
                      u(:,:,spoint),  v(:,:,spoint),  z(:,:,spoint),  &
                     qf(:,:,fpoint), df(:,:,fpoint), pf(:,:,fpoint) )
      case default
         write (*,*) "compute_model_terms:  ieq=", ieq, " out of range"
         stop        "compute_model_terms:  ieq out of range"
   end select

   if (trace) print *, "... finished compute_model_terms."

!-------------------------------------------------------------------------------
end subroutine compute_model_terms
!===============================================================================


!===============================================================================
subroutine compute_lincomb_sol( nsoln, spoint, csoln, fout )
!-------------------------------------------------------------------------------
! Purpose:  
!    Compute the linear combination
!       forcing(fout) = sum(i=1:nsoln) csoln(i)*solution(spoint(i))
!    (where "solution" here is only the prognostic variables)
! Arguments:
   integer (kind=int_kind), intent(in) :: nsoln         ! # of solution copies
   integer (kind=int_kind), intent(in) :: spoint(nsoln) ! matching indices
   real    (kind=dbl_kind), intent(in) :: csoln(nsoln)  ! matching constants
   integer (kind=int_kind), intent(in) :: fout          ! index output forcing
! Note:
!    This interface is required by the module sitpack
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: i ! loop index

   if (trace) print *, "Starting compute_lincomb_sol ..."
   if (trace) print *, "spoint = ", spoint, "  fout = ", fout

! compute the desired linear combination

   select case (abs(ieq))

      case (1) ! momentum form
         uf(:,:,fout) = csoln(1)*u(:,:,spoint(1))
         vf(:,:,fout) = csoln(1)*v(:,:,spoint(1))
         pf(:,:,fout) = csoln(1)*p(:,:,spoint(1))
         do i=2,nsoln
            uf(:,:,fout) = uf(:,:,fout) + csoln(i)*u(:,:,spoint(i))
            vf(:,:,fout) = vf(:,:,fout) + csoln(i)*v(:,:,spoint(i))
            pf(:,:,fout) = pf(:,:,fout) + csoln(i)*p(:,:,spoint(i))
         end do

      case (2) ! vorticity/divergence form
         zf(:,:,fout) = csoln(1)*z(:,:,spoint(1))
         df(:,:,fout) = csoln(1)*d(:,:,spoint(1))
         pf(:,:,fout) = csoln(1)*p(:,:,spoint(1))
         do i=2,nsoln
            zf(:,:,fout) = zf(:,:,fout) + csoln(i)*z(:,:,spoint(i))
            df(:,:,fout) = df(:,:,fout) + csoln(i)*d(:,:,spoint(i))
            pf(:,:,fout) = pf(:,:,fout) + csoln(i)*p(:,:,spoint(i))
         end do

      case (3) ! potential vorticity/divergence form
         qf(:,:,fout) = csoln(1)*q(:,:,spoint(1))
         df(:,:,fout) = csoln(1)*d(:,:,spoint(1))
         pf(:,:,fout) = csoln(1)*p(:,:,spoint(1))
         do i=2,nsoln
            qf(:,:,fout) = qf(:,:,fout) + csoln(i)*q(:,:,spoint(i))
            df(:,:,fout) = df(:,:,fout) + csoln(i)*d(:,:,spoint(i))
            pf(:,:,fout) = pf(:,:,fout) + csoln(i)*p(:,:,spoint(i))
         end do

      case default
         write (*,*) "compute_lincomb_sol:  ieq=", ieq, " out of range"
         stop        "compute_lincomb_sol:  ieq out of range"
   end select

   if (trace) print *, "... finished compute_lincomb_sol."

!-------------------------------------------------------------------------------
end subroutine compute_lincomb_sol
!===============================================================================


!===============================================================================
subroutine compute_lincomb_for( nforce, fpoint, cforce, fout )
!-------------------------------------------------------------------------------
! Purpose:  
!    Compute the linear combination
!       forcing(fout) = forcing(fout) 
!                       + sum(i=1:nforce) cforce(i)*forcing(fpoint(i))
!    (note:  must ADD the linear combination to current contents)
! Arguments:
   integer (kind=int_kind), intent(in) :: nforce         ! # of forcing copies
   integer (kind=int_kind), intent(in) :: fpoint(nforce) ! matching indices
   real    (kind=dbl_kind), intent(in) :: cforce(nforce) ! matching constants
   integer (kind=int_kind), intent(in) :: fout   ! index of output forcing
! Note:
!    This interface is required by the module sitpack
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: i ! loop index

   if (trace) print *, "Starting compute_lincomb_for ..."
   if (trace) print *, "fpoint = ", fpoint, "  fout = ", fout

! compute the desired linear combination

   select case (abs(ieq))

      case (1) ! momentum form
         do i=1,nforce
            uf(:,:,fout) = uf(:,:,fout) + cforce(i)*uf(:,:,fpoint(i))
            vf(:,:,fout) = vf(:,:,fout) + cforce(i)*vf(:,:,fpoint(i))
            pf(:,:,fout) = pf(:,:,fout) + cforce(i)*pf(:,:,fpoint(i))
         end do

      case (2) ! vorticity/divergence form
         do i=1,nforce
            zf(:,:,fout) = zf(:,:,fout) + cforce(i)*zf(:,:,fpoint(i))
            df(:,:,fout) = df(:,:,fout) + cforce(i)*df(:,:,fpoint(i))
            pf(:,:,fout) = pf(:,:,fout) + cforce(i)*pf(:,:,fpoint(i))
         end do

      case (3) ! potential vorticity/divergence form
         do i=1,nforce
            qf(:,:,fout) = qf(:,:,fout) + cforce(i)*qf(:,:,fpoint(i))
            df(:,:,fout) = df(:,:,fout) + cforce(i)*df(:,:,fpoint(i))
            pf(:,:,fout) = pf(:,:,fout) + cforce(i)*pf(:,:,fpoint(i))
         end do

      case default
         write (*,*) "compute_lincomb_for:  ieq=", ieq, " out of range"
         stop        "compute_lincomb_for:  ieq out of range"
   end select

   if (trace) print *, "... finished compute_lincomb_for."

!-------------------------------------------------------------------------------
end subroutine compute_lincomb_for
!===============================================================================


!===============================================================================
subroutine solve_implicit_eqns( tau, sitime, fpoint, spoint )
!-------------------------------------------------------------------------------
! Purpose:  
!    Given  tau,  solve the implicit problem at time  t  based on the forcing
!    copy fpoint and put the resulting solution in solution copy spoint
! Arguments:
   real    (kind=dbl_kind), intent(in) :: tau(:) ! scaled time step
   real    (kind=dbl_kind), intent(in) :: sitime ! model time from sitpack
   integer (kind=int_kind), intent(in) :: fpoint ! index of forcing copy
   integer (kind=int_kind), intent(in) :: spoint ! index of solution copy
! Notes:
!    Solves  sol(spoint) - tau(1)*implicit(sol(spoint)) = for(fpoint)
!
!    This interface routine is included here simply so that all routines 
!    called from module  sitpack  are in this module, and all details 
!    of the actual calculations are in the module  pswm_terms.
!
!    This interface is required by the module sitpack
!-------------------------------------------------------------------------------

   if (trace) print *, "Starting solve_implicit_eqns ..."
   if (trace) print *, "fpoint = ", fpoint, "  spoint = ", spoint

! sanity checks

   if ( count(tau/=0)>1 ) then
      print *, "tau = ", tau
      stop "solve_implicit_eqns:  multiple tau values"
   else if ( count(tau/=0)==1 .and. maxloc(tau,1,tau/=0)>1 ) then
      print *, "tau = ", tau
      stop "solve_implicit_eqns:  wrong term is implicit--shouldn't happen!"
   end if

! solve the implicit equations (according to the selected equation form)

   select case (abs(ieq))
      case (1) ! momentum form
         call solve_implicit_uv( tau(1), &
                     uf(:,:,fpoint), vf(:,:,fpoint), pf(:,:,fpoint), &
                      u(:,:,spoint),  v(:,:,spoint),  p(:,:,spoint), &
                      d(:,:,spoint),  z(:,:,spoint),  q(:,:,spoint) )
      case (2) ! vorticity/divergence form
         call solve_implicit_zd( tau(1), &
                     zf(:,:,fpoint), df(:,:,fpoint), pf(:,:,fpoint), &
                      u(:,:,spoint),  v(:,:,spoint),  p(:,:,spoint), &
                      d(:,:,spoint),  z(:,:,spoint),  q(:,:,spoint) )
      case (3) ! potential vorticity/divergence form
         call solve_implicit_qd( tau(1), &
                     qf(:,:,fpoint), df(:,:,fpoint), pf(:,:,fpoint), &
                      u(:,:,spoint),  v(:,:,spoint),  p(:,:,spoint), &
                      d(:,:,spoint),  z(:,:,spoint),  q(:,:,spoint) )
      case default
         write (*,*) "solve_implicit_eqns:  ieq=", ieq, " out of range"
         stop        "solve_implicit_eqns:  ieq out of range"
   end select

   if (trace) print *, "... finished solve_implicit_eqns."

!-------------------------------------------------------------------------------
end subroutine solve_implicit_eqns
!===============================================================================


!===============================================================================
subroutine put_model_solution( sitime, spoint )
!-------------------------------------------------------------------------------
! Purpose:  
!    Given the initial model time  t,  put the intitial solution into copy  s
!
! Arguments:
   real    (kind=dbl_kind), intent(in) :: sitime ! model time from sitpack
   integer (kind=int_kind), intent(in) :: spoint ! index of solution copy
!
! Note:
!    This interface is required by the module sitpack (for direct start only).
!    Here we use it for initializing the model variables; this "extra layer"
!    helps to keep the scope of the variables well-defined (see above).
!-------------------------------------------------------------------------------

   if (trace) print *, "Starting put_model_solution ..."
   if (trace) print *, "spoint = ", spoint

! sanity check

   if ( sitime/=zero ) stop "put_model_solution:  called with time>0"

! initialize the model variables (according to the selected method)

   call initialize_soln( u(:,:,spoint), v(:,:,spoint), p(:,:,spoint), &
                         d(:,:,spoint), z(:,:,spoint), q(:,:,spoint) )

   if (trace) print *, "... finished put_model_solution."

!-------------------------------------------------------------------------------
end subroutine put_model_solution
!===============================================================================


!-------------------------------------------------------------------------------
end module sitpack_interface
!===============================================================================
