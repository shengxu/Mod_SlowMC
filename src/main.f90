program main
!==============================================================================!
!> @mainpage SlowMC: Slowing Down Monte Carlo 
!>
!> @section Overview
!>
!> This program solves the slowing down neutron transport equation in either
!> infinite medium or effective two-region collision probability theory. It 
!> models parts of the same physics performed by the NJOY data processing code.
!> This code is for strictly academic purposes and allows the user to see the
!> relative impact of physics in the generation of multigroup cross sections
!> and on flux spectra. This code currently uses the following external 
!> libraries:
!>  - HDF5 v1.8.# 
!>
!> The package HDF5 can be downloaded from http://www.hdfgroup.org/HDF5/ 
!>
!> @section Compiling
!>
!> Compiling is as easy as running the Makefile with:
!>
!> @verbatim
!>   make xml-fortran
!>   make slowmc
!> @endverbatim
!>
!> @section Running
!>
!> To run SlowMC, execute the following:
!>
!> @verbatim
!>   slowmc
!> @endverbatim
!>
!==============================================================================!


  implicit none

  ! Added by S. Xu for debug  
!  open(10, file='Debug_inf.out', status='unknown')
!  open(11, file='Tally_inf.out', status='unknown')
!  open(997, file="./sampled_xs/capt.out", status='unknown')
!  open(998, file="./sampled_xs/scat.out", status='unknown')
!  open(999, file="./sampled_xs/fiss.out", status='unknown')

  ! initialize problem
  call initialize()

  ! run problem
  call run_problem()

  ! finalize problem 
  call finalize()
  
  ! Added by S. Xu for debug 
!  close(10)
!  close(11)
!  close(997)
!  close(998)
!  close(999)

  ! terminate program
  stop

contains

!===============================================================================
! INTIALIZE
!> @brief high level routine for intializing problem 
!===============================================================================

  subroutine initialize()

    use hdf5
    use global
    use execute,   only: allocate_problem
    use materials, only: compute_macro_xs, doppler_broaden
    use output,    only: print_heading
    use timing,    only: timer_start,timer_stop
    use random_lcg, only: initialize_prng
    use input,     only: read_input

#ifdef MPI
    use mpi
#endif

    ! local variables
    integer :: error ! hdf5 error
    real(8) :: rn    ! initial random number

!    ! begin timer
!    call timer_start(time_init)

    ! Setup MPI
    call setup_mpi()

      ! initialize the fortran hdf5 interface
      call h5open_f(error)

    if (master) then


      ! print heading information
      call print_heading()

!output debug information
#ifdef DEBUG
    open(404, file='MPI_DEBUG.txt', status='unknown')
#endif

    end if

    ! read input
    call read_input()

!    ! initalize random number generator
!    rn = rand(seed)

    ! Initialize random number generator
    call initialize_prng()

! commented out by S. Xu (Apr. 2012)
!    ! precompute macroscopic cross section of materials
!    call compute_macro_cross_sections()

!    ! end timer
!    call timer_stop(time_init)

  end subroutine initialize


!===============================================================================
! SETUP_MPI initilizes the Message Passing Interface (MPI) and determines the
! number of processors the problem is being run with as well as the rank of each
! processor.
!===============================================================================

  subroutine setup_mpi()

    use global,    only: rank,n_procs,master,mpi_err

#ifdef MPI
    use mpi
#endif

#ifdef MPI
!    integer        :: bank_blocks(5) ! Count for each datatype
!    integer        :: bank_types(5)  ! Datatypes
!    integer(MPI_ADDRESS_KIND) :: bank_disp(5)   ! Displacements
!    type(Bank)     :: b

!    mpi_enabled = .true.

    ! Initialize MPI
    call MPI_INIT(mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       print *, "Failed to initialize MPI."
       stop
    end if

    ! Determine number of processors
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       print*, "Could not determine number of processors."
       stop
    end if

    ! Determine rank of each processor
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       print *, "Could not determine MPI rank."
       stop
    end if

    ! Determine master
    if (rank == 0) then
       master = .true.
    else
       master = .false.
    end if

!    ! Determine displacements for MPI_BANK type
!    call MPI_GET_ADDRESS(b % id,  bank_disp(1), mpi_err)
!    call MPI_GET_ADDRESS(b % wgt, bank_disp(2), mpi_err)
!    call MPI_GET_ADDRESS(b % xyz, bank_disp(3), mpi_err)
!    call MPI_GET_ADDRESS(b % uvw, bank_disp(4), mpi_err)
!    call MPI_GET_ADDRESS(b % E,   bank_disp(5), mpi_err)

!    ! Adjust displacements 
!    bank_disp = bank_disp - bank_disp(1)
!    
!    ! Define MPI_BANK for fission sites
!    bank_blocks = (/ 1, 1, 3, 3, 1 /)
!    bank_types = (/ MPI_INTEGER8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8 /)
!    call MPI_TYPE_CREATE_STRUCT(5, bank_blocks, bank_disp, & 
!         bank_types, MPI_BANK, mpi_err)
!    call MPI_TYPE_COMMIT(MPI_BANK, mpi_err)

#else
    ! if no MPI, set processor to master
!    mpi_enabled = .false.
    rank = 0
    n_procs = 1
    master = .true.
#endif

  end subroutine setup_mpi

!===============================================================================
! RUN_PROBLEM
!> @brief main routine for executing the transport calculation
!===============================================================================

  subroutine run_problem()

    use global
    use materials, only: doppler_broaden, compute_macro_xs
    use tally,     only: bank_tallies
    use particle,  only: init_particle
    use physics,   only: sample_source,perform_physics !,get_eidx !commented out by S. Xu
    use timing,    only: timer_start,timer_stop
    use random_lcg, only: set_particle_seed

#ifdef MPI
    use mpi
#endif
!    use on_the_fly_xs_gen, only: doppler_broaden

    ! local variables
    integer(8) :: i  ! iteration counter

!    ! begin timer
!    call timer_start(time_run)

! every processor should know the number of history by reading it from either command line or input file
!#ifdef MPI
!    if (master) then
!      call MPI_BCAST(nhistories, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_err)
!    end if
!#endif

    ! begin loop over histories
    do i = rank + 1, nhistories, n_procs

      ! set the seed for random number generator for each particle history
      call set_particle_seed(i)

      ! intialize history
      call init_particle(neut)

      ! sample source energy
      call sample_source()

      ! begin transport of neutron
      do while (neut%alive)

        ! check for energy cutoff
        if (neut%E < emin) neut%E = 1.1e-11_8

! Commented out by S. Xu (Apr. 2012)
!        ! call index routine for first tally
!        eidx = get_eidx(neut%E)


!-----------------------------------------------------------------------
! Added by S. Xu (APr. 2012)
        ! doppler broadening
        call doppler_broaden()

        ! compute macroscopic xs
        call compute_macro_xs()
!-----------------------------------------------------------------------

        ! perform physics and also records collision tally
        call perform_physics()

      end do

      ! neutron is dead if out of transport loop (ecut or absorb) --> bank tally
      call bank_tallies()

!      if (master) then
!        ! print update to user
!        if (mod(i,nhistories/10) <= n_procs-1) then
!          write(*,'(/A,1X,I0,1X,A)') 'Simulated',i,'neutrons...'
!        end if
!      end if

    end do

!    ! end timer
!    call timer_stop(time_run)


  end subroutine run_problem

!===============================================================================
! FINALIZE
!> @brief routine that finalizes the problem
!===============================================================================

  subroutine finalize()

    use global
    use execute, only: deallocate_problem
    use tally,  only: reduce_tallies, finalize_tallies
    use hdf5
    use output, only: write_output

#ifdef MPI
    use mpi
#endif

    ! local variables
    integer :: error ! hdf5 error

    call reduce_tallies()


    if (master) then

      ! calculate statistics on tallies
      call finalize_tallies()

!      ! compute resonance integral
!      if (res_intg) then
!        call compute_res_intg()
!      end if

      ! write output
      call write_output()


! close debug output
#ifdef DEBUG
    close(404)
#endif
    end if 

      ! close the fortran interface
      call h5close_f(error)

    ! deallocate problem
    call deallocate_problem()

#ifdef MPI
    ! shut down MPI
    call MPI_Finalize(mpi_err)
#endif

  end subroutine finalize

end program main
