!===============================================================================
program pswm_main
!-------------------------------------------------------------------------------
! Purpose: 
!
!    Main program to run the Periodic Shallow Water Model
!
! Description: 
!
!    First sets the runid (from the command line or the file runid.txt).
!    Then reads the file <runid>.dat:  first reads model parameters using
!    read_params,  and then the data parameters using  rdatap.  
!    Echos all input parameters to stdout and runs the model.
!
!    This program unit controls the time stepping via calls to sitpack.
!
! See also:
!    pswm_pars     for descriptions of model parameters
!    pswm_output   for description of model output
! 
! Author:
!    Scott R. Fulton
!
! Revison history:
!    07/23/07:  Original version (based in part on cswm and sitpack/test2)
!-------------------------------------------------------------------------------

   use kinds
   use sitpack_interface
   use sitpack
   use pswm_pars
   use pswm_vars
   use pswm_setup
   use pswm_data
   use pswm_output

   implicit none

! local variables

   integer (kind=int_kind), parameter :: iunit  = 10  ! input  unit
   integer (kind=int_kind), parameter :: ounit  =  6  ! output unit
   integer (kind=int_kind)            :: nsteps =  0  ! number of steps
   logical, parameter :: trace = .false. ! turns on/off tracing of time steps
   character (len=80) :: slabel = ""     ! time scheme label

! determine the runid

   call get_runid

! read the model parameters and initial condition parameters

   open ( unit=iunit, file=trim(runid)//'.dat', status='old' )
   call read_params( iunit )
   call read_icpars( iunit )
   close ( unit=iunit )

! set up the time discretization

   call setup_time( slabel )
   write (ounit,*) "--------------------------------------------------------"
   write (ounit,*) "Time discretization parameters:"
   write (ounit,*) "time scheme:  ", trim(slabel)
   write (ounit,*) "time step   dt = ", dt, " seconds"
   write (ounit,*) "--------------------------------------------------------"

! set up the space discretization

   call setup_space
   write (ounit,*) "--------------------------------------------------------"
   write (ounit,*) "Space discretization parameters:"
   write (ounit,*) "Transform length (physical grid):  nx = ", nx, " ny = ", ny
   write (ounit,*) "Spectral truncation  (last mode):  mx = ", mx, " my = ", my
   write (ounit,*) "--------------------------------------------------------"

! write the model parameters to the output file

   call write_params( ounit )
   call write_icpars( ounit )

! initialize the model variables

   call pswm_start

! produce the initial output

   call setup_output
   call write_output

! integrate the model forward in time
!   nnn = zero

   if ( tstop>0 ) then
      nsteps = nint( tstop/dt )
   else
      nsteps = nint( abs( tstop ) )
   end if
   time_loop: do while (nstep<nsteps)

!    do the time step

      call sitpack_step( dt, now, nstep_si=nstep, time_si=time )
      if ( trace ) print "(' finished step',i6,' time = ',f8.2, ' hours')", &
         nstep, time/hour

!    write the output (what and when is controlled by various model parameters)

      call write_output

!      nnn = (nnn + 1.0)*dt
   end do time_loop

! clean up

   call pswm_cleanup
   call unset_output
   stop

!-------------------------------------------------------------------------------
end program pswm_main
!===============================================================================
