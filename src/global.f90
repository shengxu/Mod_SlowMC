!==============================================================================!
! MODULE: global 
!
!> @author Bryan Herman
!>
!> @brief Contains all of the global variables
!==============================================================================!

module global

  use tally_header, only: tally_type
  use material_header, only: n_source_type, thermal_type, iso_type, material_type
  use particle,  only: particle_type
  use timing,    only: Timer
  use constants

  implicit none
  save

  ! version information
  integer :: VERSION_MAJOR   = 0
  integer :: VERSION_MINOR   = 1
  integer :: VERSION_RELEASE = 1

  ! list all types
  type(particle_type)              :: neut
  type(material_type), allocatable :: mat(:)
  type(tally_type), allocatable    :: tal(:)
  type(tally_type), allocatable    :: reduced_tal(:)  ! store reduced tally information across different processors, used for final results

  ! list history input information
  integer(8) :: nhistories=10000
  integer :: seed=1
  integer :: source_type=1
  ! number of samples per xs for broadening
  integer :: sample_per_xs=1
  real(8) :: randerror=0.05

  ! list global vars that are set during run
! Commented out by S. Xu
!  integer :: eidx        ! energy index for cross sections
  integer :: n_tallies   ! number of tallies
  integer :: n_materials ! n_materials
  integer :: res_iso     ! resonant isotope id in material 1
  real(8) :: Dancoff     ! lattice Dancoff factor (C)
  real(8) :: radius      ! radius of fuel pin

  real(8)      :: T   ! temperature
  real(8) :: kT != 8.6173324e-5_8*300*1.0e-6_8

  ! set nu value
  real(8) :: nubar = 2.455_8

  ! set max and min energy
  real(8) :: emin = 1e-11_8
  real(8) :: emax = 20.0_8
  real(8) :: vmin = sqrt(2._8*1e-11_8/M_NEUT)
  real(8) :: vmax = sqrt(2._8*20._8/M_NEUT)


  ! timers
  type(Timer) :: time_init
  type(Timer) :: time_run

  ! analog counters for k-inf
  integer :: n_abs=0
  integer :: n_fiss=0
  real(8) :: ana_kinf_mean = 0.0_8
  real(8) :: ana_kinf_std  = 0.0_8
  real(8) :: col_kinf_mean = 0.0_8
  real(8) :: col_kinf_std  = 0.0_8
  integer :: reduced_n_abs=0
  integer :: reduced_n_fiss=0

!  !added by S. Xu (May. 2012)
!  real(8) :: col_kinf_mean_2 = 0.0_8
!  real(8) :: col_kinf_std_2  = 0.0_8
!  real(8) :: reaction_tally(3)

  character(len=255)  :: filename    
  character(len=255)  :: output_filename
  character(len=255)  :: output_path

!  For resonance integral tally
  logical :: res_intg = .false.   ! indicator for resonance integral
  integer :: res_intg_inf(2) = 0 ! store the index for micro_capture and flux tally
  real(8), allocatable :: res_intg_rec(:,:) ! to store resonance integral (with variance)
  integer :: res_intg_fd = 511   ! file descriptor for resonance integral

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

  integer :: n_procs     ! number of processes
  integer :: rank        ! rank of process
  logical :: master      ! master process?
!  logical :: mpi_enabled ! is MPI in use and initialized?
  integer :: mpi_err     ! MPI error code
!  integer :: MPI_BANK    ! MPI datatype for fission bank

!  ! No reduction at end of batch
!  logical :: no_reduce = .false.

end module global
