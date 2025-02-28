!===============================================================================
module pswm_output
!-------------------------------------------------------------------------------
! Purpose:
!    Output routines for pswm:  Periodic Shallow Water Model
!
! Use:
!    First call   setup_output  to initialize stuff needed in this module.
!    Then call    write_output  at each time step--it does the rest.
!    You can call unset_output  at the end to release memory (for clean exit).
!
! Note:
!    What is output (and when) is controlled by various model parameters.
!    Scalar diagnostics are output for more variables than are available 
!    for fields--see data statement below for available variables.
!
! See also:
!    module pswm_pars       (for details of model parameters controlling output)
!    output routines--below (for details of file names and their contents)
!-------------------------------------------------------------------------------

use kinds
use pswm_pars
use pswm_cons
use pswm_vars
use pswm_terms
use ss_ops

implicit none   ! all module variables untyped by default
private         ! all module variables and routines hidden by default

public          &! routines to be called from main program or elsewhere:
   setup_output,   &! initialization routine
   write_output,   &! front-end to determine what and when to write
   unset_output     ! releases allocated memory for clean exit when done

! module variables

character (len=80)      :: fname ! for file names
real (kind=dbl_kind)    :: & ! time of last output of each type
   time_out0 = zero, time_out2 = zero, time_outr = zero 
integer (kind=int_kind) :: ounit = 30 ! output unit number
real (kind=dbl_kind), allocatable, save :: fout(:,:) ! array for output field
integer (kind=int_kind) :: ifld ! index for output fields
integer (kind=int_kind), save :: & ! indices for output points
   ja, jb, jc, ka, kb, kc, nxout, nyout
integer (kind=int_kind), parameter :: &
   nfld =  6, & ! number of fields available for fields output
   mfld = 10    ! number of fields available for scalar output
real (kind=dbl_kind) :: diagnostics(3,mfld) ! scalar diagnostic values
character (len= 1) :: lfld(mfld)       ! one-character field labels
character (len=26) :: labelf(mfld)*26  ! longer field labels
data  ( lfld(ifld), labelf(ifld), ifld=1,mfld ) / &
      'u', 'x-component of velocity   ',&
      'v', 'y-component of velocity   ',&
      'p', 'geopotential deviation    ',& ! p = g*(h - href)/cval
      'd', 'divergence                ',& ! d = delta
      'z', 'relative  vorticity       ',& ! z = zeta
      'q', 'potential vorticity       ',& ! q = eta/(1 + p/cval)
      'e', 'absolute  vorticity       ',& ! e = eta = fcor + zeta
      'Z', 'enstrophy                 ',& ! Z = z*z/2
      'K', 'kinetic energy            ',& ! K = (u*u + v*v)*(cval + p)/2
      'A', 'available potential energy' / ! A = cval*p*p/2

! Note:
!    The domain average of a field is its integral over the model domain
!    divided by the area of the model domain, and the energies are missing
!    the factor  cval/g  (since g isn't needed elsewhere in the model).

contains

!===============================================================================
subroutine setup_output
!-------------------------------------------------------------------------------
! Purpose:
!    Initialization routine for output
!-------------------------------------------------------------------------------

! determine the physical space points for output

   ja = max( j1, 0 ); jb = min( j2, nx ); jc = max( jinc, 1 );
   ka = max( k1, 0 ); kb = min( k2, nx ); kc = max( kinc, 1 );
   if ( ja==0 .and. jb==0 .and. ka==0 .and. kb==0 ) then
      ja = 0; jb = nx; ka = 0; kb = ny
   end if
   nxout = (jb-ja)/jc; nyout = (kb-ka)/kc;

! allocate space for the output grid

   if (allocated(fout)) deallocate( fout ); allocate( fout(-1:nx,-1:ny) )

!-------------------------------------------------------------------------------
end subroutine setup_output
!===============================================================================


!===============================================================================
subroutine write_output
!-------------------------------------------------------------------------------
! Purpose:
!    Output front end:  determines what to write and when
!-------------------------------------------------------------------------------

! write output of each kind--if the time is right

   if ( time>zero ) then
      if ( output_time( dtout0, time_out0 ) ) call output_scalars
      if ( output_time( dtout2, time_out2 ) ) call output_fields
      if ( output_time( dtoutr, time_outr ) ) call output_restart
   else if ( output_ic ) then ! initial output
      if ( dtout0/=zero ) call output_scalars
      if ( dtout2/=zero ) call output_fields
   end if

!-------------------------------------------------------------------------------
end subroutine write_output
!===============================================================================


!===============================================================================
subroutine output_scalars
!-------------------------------------------------------------------------------
! Purpose:
!    Writes output of scalars to the file runid.out (appends after first call)
!
! Description:
!    The first time this routine is called, it creates the text file runid.out
!    and writes a header containing the model version and values of  c  and  f.
!    Then the scalar diagnostics are appended to this file at this time and 
!    at subsequent output times.
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: i              ! loop index
   logical                 :: first = .true. ! flags the first call

! output formats

   character (len=132) :: vdfmt, fcfmt, sdfmt 
   vdfmt = "(a,' version',f5.2,' (',a,')')"
   fcfmt = "('fcor=',1p,e20.12,' beta=',e20.12,' cval=',e20.12)"
   sdfmt = "(/'step',i6,':  time=',1p,e11.5,&
         &'       minimum        maximum        average',/(a1,':  ',a,3e15.5))"

! compute the scalar diagnostics

   call compute_diagnostics( u(:,:,now), v(:,:,now), p(:,:,now), &
                             d(:,:,now), z(:,:,now), q(:,:,now), diagnostics )

! write the scalar diagnostics to the output file

   fname = trim(runid)//'.out'
   if ( first ) then
       open ( unit=ounit, file=fname, status='new' )
       write (ounit,vdfmt) model, vnum, vdate
       write (ounit,fcfmt) fcor, zero, cval
       first = .false.
   else
       open ( unit=ounit, file=fname, status='old', position='append' )
   end if
   write (ounit,sdfmt) nstep, time, &
      (lfld(ifld), labelf(ifld), (diagnostics(i,ifld), i=1,3), ifld=1,mfld)
   close ( unit=ounit )

!-------------------------------------------------------------------------------
end subroutine output_scalars
!===============================================================================


!===============================================================================
subroutine output_fields
!-------------------------------------------------------------------------------
! Purpose:
!    Writes output of 2D fields to a file
!
! Description:
!    Each time this routine is called, a separate file is written for each
!    field requested by the parameter  ifield.  The filename is runid_f.tnn, 
!    where:
!       runid  is the run identifier
!       f      identifies the field as specified by  ifield:  u, v, etc.
!       n      is the sequence number of the file:  00, 01, 02, ...
!    For example:  with  runid='pstest' and ifield=20 the following
!    files are generated:
!       first  call:  pstest_p.t00, pstest_z.t00
!       second call:  pstest_p.t01, pstest_z.t01
!       (etc.)
!
!    Each file has a five-line header:
!       PSWM version ?.?? (mm/dd/yy)
!       field ? for run runid
!       nxout, nyout
!       xaout, xbout, yaout, ybout
!       fcor, beta, cval, time [with beta=0.0 for compatibility with cswm]
!    followed by the field values (one per line) in the order:
!       ((field(j,k) = field(xout(j), yout(k)), j=0,nxout), k=0,nyout)
!    where
!       xout(j) = xaout + j*(xbout-xaout)/nxout, j=0,...,nxout
!       yout(j) = yaout + k*(ybout-yaout)/nyout, k=0,...,nyout
!       The format of the field values is controlled by the parameter nsigd:
!          nsigd>0  formatted output, nsigd digits in E format
!          nsigd=0  binary output
!    For more details see the comments describing  nsigd  in pswm_pars.
!
! Note:
!    This output is intended to be the same as that used in cswm.
!    Known differences:  
!       - potential vorticity is available for output
!       - output points must be a subset of the transform grid
!-------------------------------------------------------------------------------

! local variables

   integer (kind=int_kind) :: indout = 0 ! index of output call/file
   character (len=3)       :: ext        ! file extension
   logical                 :: binary     ! binary output format?
   integer (kind=int_kind) :: j, k, ifld ! loop indices

! output formats

   character (len=132) :: fmt1, fmt2, fmt3, fmt4, fmtn
   fmt1 = "(a,' version',f5.2,' (',a,')')"
   fmt2 = "('field ',a1,' for run ',a)"
   fmt3 = "(2i20)"
   fmt4 = "(1p,4e20.12)"

! determine the type of output requested and set the format accordingly

   binary = nsigd<=0
   if ( .not.binary ) then
      if ( nsigd<=2 ) then
         write (fmtn,'(''(1pe'',i1,''.'',i1,'')'')') nsigd+7, nsigd-1
      else if ( nsigd<=10 ) then
         write (fmtn,'(''(1pe'',i2,''.'',i1,'')'')') nsigd+7, nsigd-1
      else
         write (fmtn,'(''(1pe'',i2,''.'',i2,'')'')') nsigd+7, nsigd-1
      end if
   end if

! output each field requested

   field_loop:  do ifld=1,nfld
      if ( mod( ifield/2**(ifld-1), 2 )/=1 ) cycle

!   transform the output field to the physical space grid

      select case (lfld(ifld))
         case ("u")
            fout = u(:,:,now)
         case ("v")
            fout = v(:,:,now)
         case ("p")
            fout = p(:,:,now)
         case ("d")
            fout = d(:,:,now)
         case ("z")
            fout = z(:,:,now)
         case ("q")
            fout = q(:,:,now)
         case default
            write (*,*) "field ", lfld(ifld), " not available for output"
      end select
      call tophys( fout )

!   write the output file

      if ( indout<100 ) then
         write (ext,"('t',i2.2)") indout
      else
         write (ext,"(i3.3)") indout
      end if
      fname = trim(runid)//'_'//lfld(ifld)//'.'//ext
      if ( binary ) then
         open ( unit=ounit, file=fname, status='new', form='unformatted' )
         write (ounit) model, vnum, vdate
         write (ounit) lfld(ifld), trim(runid)
         write (ounit) nxout, nyout
         write (ounit) x(ja), x(jb), y(ka), y(kb)
         write (ounit) fcor, zero, cval, time
         write (ounit) ((fout(j,k), j=ja,jb,jc), k=ka,kb,kc)
      else
         open ( unit=ounit, file=fname, status='new' )
         write (ounit,fmt1) model, vnum, vdate
         write (ounit,fmt2) lfld(ifld), trim(runid)
         write (ounit,fmt3) nxout, nyout
         write (ounit,fmt4) x(ja), x(jb), y(ka), y(kb)
         write (ounit,fmt4) fcor, zero, cval, time
         write (ounit,fmtn) ((fout(j,k), j=ja,jb,jc), k=ka,kb,kc)
      end if
      close (ounit)

   end do field_loop

!   increment the output index 

   indout = indout + 1

!-------------------------------------------------------------------------------
end subroutine output_fields
!===============================================================================


!===============================================================================
subroutine output_restart
!-------------------------------------------------------------------------------
! Purpose:
!    Writes output of everything needed to restart the model
!
! Description:
!    This is (obviously) a dummy routine--you can write your own if you wish.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
end subroutine output_restart
!===============================================================================



!===============================================================================
subroutine unset_output
!-------------------------------------------------------------------------------
! Purpose:
!    Releases allocated space (for clean exit when done)
!-------------------------------------------------------------------------------

   if (allocated(fout)) deallocate( fout )

!-------------------------------------------------------------------------------
end subroutine unset_output
!===============================================================================


!===============================================================================
function output_time( dtout, time_out ) result (output_now)
!-------------------------------------------------------------------------------
! Purpose:
!    Internal function to decide whether it's time to produce output
!
! Arguments:
   real (kind=dbl_kind), intent(in)    :: dtout    ! time between output
   real (kind=dbl_kind), intent(inout) :: time_out ! time of last output
!
! Result:
   logical :: output_now ! time for output?
!
! Side Effects:
!    time_out is incremented (and is reset to zero if there is output)
!-------------------------------------------------------------------------------

! local variable

   integer (kind=int_kind) :: nstep_out ! number of steps between output

! check for output

   if ( dtout>zero ) then
      output_now = abs((time-time_out)-dtout)<0.75*dt
   else if ( nint(dtout)<zero ) then
      nstep_out = abs( nint( dtout ) )
      if ( nstep_out>0 .and. mod(nstep,nstep_out).eq.0 ) &
         output_now = .true.
   else
      output_now = .false.
   end if

! prepare time_out for the next step

   if ( output_now ) time_out = time

!-------------------------------------------------------------------------------
end function output_time
!===============================================================================


!-------------------------------------------------------------------------------
end module pswm_output
!===============================================================================
