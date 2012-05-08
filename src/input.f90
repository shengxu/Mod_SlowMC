!==============================================================================!
! MODULE: input
!
!> @author Bryan Herman
!>
!> @brief Handles reading in the input xml file and intializing global vars
!==============================================================================!

module input

  implicit none
  private
  public read_input

contains

!===============================================================================
! READ_INPUT
!> @brief Reads the input xml file and sets global variables
!===============================================================================

  subroutine read_input
    
    ! added by S. Xu
    use parameters,        only: emin,emax

    use global,           only: nhistories,seed,source_type,mat, & !emin,emax,     & ! commented out by S. Xu
   &                            allocate_problem,tal,n_tallies,n_materials,    &
   &                            res_iso,Dancoff,radius, in_out_filename
    use materials,        only: setup_material,load_source,load_isotope
    use tally,            only: set_user_tally,set_spectrum_tally,             &
   &                            set_kinf_tally
    use xml_data_input_t
! Added by S. Xu (Apr. 2012)
    use constants, only: PI, K_BOLTZMANN, M_NUCLEON
    use parameters, only: T


    ! local variables
    logical                        :: file_exists  ! see if file exists
    character(len=255)             :: filename     ! filename to open
    real(8)                        :: N            ! temp number dens
    real(8)                        :: A            ! temp atomic weight
    real(8)                        :: vol          ! volume of region
    character(len=255)             :: path         ! path to isotope file
    character(len=255)             :: name         ! name of isotope
    logical                        :: thermal      ! contains thermal lib
    integer                        :: i            ! iteration counter
    integer                        :: j            ! iteration counter
    integer                        :: nisotopes    ! number of isotopes in mat
    integer                        :: react_type   ! reaction type
    integer                        :: isotope=0    ! isotope for micro mult
    integer                        :: region=0     ! region to tally micros
    real(8), allocatable           :: Ebins(:)     ! tally energy bins
    integer                        :: nbins        ! number of energy bins
    logical                        :: dv           ! divide tally by volumes 
    ! added by S. Xu
    logical               :: doppler ! indicte whether doppler braodening is performed
    integer               :: n_arg   ! number of input argument
    character(len=255)    :: arg_str ! input arguments
    character             :: nhist_1digt, nhist_power ! convert nhistories to string          
    logical               :: nhist_from_cmdline ! check whether nhistories is input from command line

! Modified by S. Xu (May 2012)
!    ! check for input file
!    filename = "input.xml"
    filename = "input.xml"
    nhist_from_cmdline = .false.
    i = 1
    do 
      call get_command_argument(i, arg_str)
      if (len_trim(arg_str) == 0) exit
      select case(i)
        case(1)
          filename = trim(arg_str)
        case(2)
          read(arg_str, *) nhistories
          nhist_from_cmdline = .true.
        case default
          print *, 'Input argument should be less than or equal to 2'
      end select
      i = i + 1
    end do
      
!    n_arg = command_argument_count()
!    if (n_arg == 0) then
!      filename = "input.xml"
!    else if(n_arg == 1) then
!      call get_command_argument(1, filename)
!    else if(n_arg == 2) then
!      call get_command_argument(2, arg_str)
!      read(arg_str, *) nhistories
!      nhist_from_cmdline = .true.
!    else
!      write(*,*) "Too many input argument!"
!    end if
    
    write(*,*) filename
    

    inquire(FILE=trim(filename), EXIST=file_exists)
    if (.not. file_exists) then
      write(*,*) 'Cannot read input file!'
      stop
    else

      ! tell user
      write(*,'(A/)') "Reading INPUT XML file..."

    end if

    ! read in input file
    call read_xml_file_input_t(trim(filename))

    ! read in settings
    if (.not. nhist_from_cmdline) then
      nhistories = settings_%histories
    end if

    write(*,*) "number of history:", nhistories

    write(nhist_power, '(i1)') int(log10(dble(nhistories)))
    write(nhist_1digt, '(i1)') int(nhistories/10**log10(dble(nhistories)))
!print *, nhist_1digt, nhist_power
    in_out_filename = filename(1:(len(trim(filename))-4))//'_'//nhist_1digt//'e'//nhist_power

!    write(*,*) in_out_filename

    seed = settings_%seed
    source_type = settings_%source_type
    
    ! Added by S. Xu (Apr. 2012) for doppler broadening
    T = settings_%temperature

    ! get size of materials
    n_materials = size(materials_%material)

    ! get size of tallies
    if (.not.associated(tallies_%tally)) then
      n_tallies = 2
    else
      n_tallies = size(tallies_%tally) + 2
    end if

    ! allocate problem
    call allocate_problem()

    ! begin loop around materials
    do i = 1,n_materials

      ! get number of isotopes and volume
      nisotopes = size(materials_%material(i)%nuclides)
      vol = materials_%material(i)%V

      ! set homogeneous volume
      if (trim(materials_%material(i)%type)=='homogeneous') vol = 1.0_8

      ! set up the material object
      call setup_material(mat(i),emin,emax,nisotopes,vol)

      ! begin loop over isotope materials
      do j = 1,mat(i)%nisotopes

        ! check volumes and number densities
        if (trim(materials_%material(i)%type)=='homogeneous') then

          ! set volume to 1 and don't adjust n dens
          N = materials_%material(i)%nuclides(j)%N

        else if (trim(materials_%material(i)%type)=='fuel') then

          ! don't adjust n dens
          N = materials_%material(i)%nuclides(j)%N

          ! check volume
          if (abs(vol - 0.0_8) < 1e-10_8) then
            write(*,*) 'Please enter a physical fuel volume!'
            stop
          end if

        else

          ! check volume
          if (abs(vol - 0.0_8) < 1e-10_8) then
            write(*,*) 'Please enter a physical moderator volume!'
            stop
          end if

          ! adjust number density by volume weighting
          N = materials_%material(i)%nuclides(j)%N*                            &
         &   (materials_%material(i)%nuclides(j)%V/vol)

        end if

        ! extract other info
        A = materials_%material(i)%nuclides(j)%A
        path = materials_%material(i)%nuclides(j)%path
        thermal = materials_%material(i)%nuclides(j)%thermal
        name = materials_%material(i)%nuclides(j)%name 
        ! Added by S. Xu
        doppler = materials_%material(i)%nuclides(j)%doppler 

        ! load the isotope into memory
        call load_isotope(mat(i),N,A,path,thermal,name, doppler)

        ! check for resonant isotope in material 1
        if (trim(materials_%material(i)%type)=='fuel' .and.                    &
       &    trim(settings_%res_iso) == trim(name)) then

          ! get Dancoff factor and resonant isotope
          res_iso = j
          Dancoff = settings_%Dancoff
          radius = settings_%radius

        end if

      end do

    end do

    ! begin loop over tallies
    do i = 1,n_tallies-2

      ! check for divide by volume
      dv = .false.
      if (tallies_%tally(i)%dv) dv = .true.

      ! set reaction type
      select case(trim(tallies_%tally(i)%type))
        case('flux')
          react_type = 0
        case('total')
          react_type = 1
        case('absorption')
          react_type = 2
        case('scattering')
          react_type = 3
        case('nufission')
          react_type = 4
        case('fission')
          react_type = 5
        case('diffusion')
          react_type = 6
        case('transport')
          react_type = 7
        case('micro_capture')
          react_type = 8
          isotope = tallies_%tally(i)%isotope
          region = tallies_%tally(i)%region
        case DEFAULT
          react_type = 0
      end select

      ! preallocate Ebins
      nbins = size(tallies_%tally(i)%Ebins) - 1
      if(.not. allocated(Ebins)) allocate(Ebins(nbins+1))

      ! set Ebins
      Ebins = tallies_%tally(i)%Ebins

      ! set up user-defined tallies
      call set_user_tally(tal(i),Ebins,nbins,react_type,isotope,region,  &
     &                    n_materials,dv)

      ! deallocate Ebins
      if(allocated(Ebins)) deallocate(Ebins)

    end do

    ! set up spectrum tally
    call set_spectrum_tally(tal(n_tallies-1),emax,emin,n_materials)

    ! set up k_inf tally
    call set_kinf_tally(tal(n_tallies),emax,emin,n_materials)

    ! load the source
    call load_source(mat(1),source_type,settings_%source_path)

  end subroutine read_input

end module input
