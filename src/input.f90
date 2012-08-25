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
    
    use global
    use execute,          only: allocate_problem
    use materials,        only: setup_material,load_source,load_isotope
    use tally,            only: set_user_tally,set_spectrum_tally, set_kinf_tally
    use xml_data_input_t

!    use constants, only: PI, K_BOLTZMANN, M_NUCLEON


    ! local variables
    logical                        :: file_exists  ! see if file exists
    character(len=255)             :: path    
    real(8)                        :: N            ! temp number dens
    real(8)                        :: A            ! temp atomic weight
    real(8)                        :: vol          ! volume of region
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
    logical               :: samplexs_from_cmdline ! check whether sample_per_xs is input from command line
    character(len=3)      :: samplexs_str
! Modified by S. Xu (May 2012)
!    ! check for input file
!    filename = "input.xml"
    filename = "input.xml"
    nhist_from_cmdline = .false.
    samplexs_from_cmdline = .false.
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
        case(3)
          read(arg_str, *) sample_per_xs
          samplexs_from_cmdline = .true.
        case default
          print *, 'Input argument should be less than or equal to 3'
      end select
      i = i + 1
    end do

    if (master) then
      write(*,'("Input file: ", A/)') trim(filename)
    end if
    
    output_filename = filename(1:(len(trim(filename))-4))

    inquire(FILE=trim(filename), EXIST=file_exists)
    if (.not. file_exists) then
      if (master) then
        write(*,*) 'Cannot read input file!'
      end if
      stop
    else
  
      ! tell user
      if (master) then
        write(*,'(A/)') "Reading INPUT XML file..."
      end if

    end if

    ! set to default values
    settings_%histories = nhistories
    settings_%seed = seed
    settings_%source_type = source_type
    settings_%output_path = ''

    ! read in input file
    call read_xml_file_input_t(trim(filename))

    ! read in settings
    if (.not. nhist_from_cmdline) then
      nhistories = settings_%histories
    end if

    if (master) then
      write(*,'(A,i15/)') "number of history: ", nhistories

    ! get the neutron histories to output file name
      write(nhist_power, '(i1)') int(log10(dble(nhistories)))
      write(nhist_1digt, '(i1)') int(nhistories/10**log10(dble(nhistories)))
    end if
!print *, nhist_1digt, nhist_power

    if (nhist_from_cmdline) then
      output_filename = trim(output_filename)//'_'//nhist_1digt//'e'//nhist_power
    end if

    if ( (.not. samplexs_from_cmdline) .and. (settings_%sample_per_xs/= 0) ) then
      sample_per_xs = settings_%sample_per_xs
    end if

    if (master) then
      write(*,'(A,i15/)') "number of sample per xs data: ", sample_per_xs

      write(samplexs_str,'(i3.3)') sample_per_xs
    end if

!write(*, *) samplexs_str

    if (samplexs_from_cmdline) then
      output_filename = trim(output_filename)//'_'//samplexs_str
    end if

    seed = settings_%seed
    source_type = settings_%source_type
    
    ! check output directory (not sure if there is better approach)
!print *, 'try another:',trim(settings_%output_path)//'/'
    output_path = './'//trim(settings_%output_path)//'/'
    if (master) then
      write(*, '(A, A/)') 'output_path = ', trim(output_path)

      inquire(FILE=trim(output_path), EXIST=file_exists)
      if (.not. file_exists) then
        call system('mkdir '//trim(output_path))
      end if

    end if
    
    
    ! Added by S. Xu (Apr. 2012) for doppler broadening
    T = settings_%temperature

!    ! for resonance integral
!    res_intg = settings_%res_intg

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
      

      ! set homogeneous volume
      if (trim(materials_%material(i)%type)=='homogeneous') then
        vol = 1.0_8
      else
        vol = materials_%material(i)%V
      end if

      ! set up the material object
      call setup_material(mat(i), emin, emax, nisotopes, vol)

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
            if (master) then
              write(*,*) 'Please enter a physical fuel volume!'
            end if
            stop
          end if

        else

          ! check volume
          if (abs(vol - 0.0_8) < 1e-10_8) then
            if (master) then
              write(*,*) 'Please enter a physical moderator volume!'
            end if
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

!          ! for resonance integral
!          if (res_intg) then
!            res_intg_inf(2) = i
!          end if

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

!          ! for resonance integral
!          if (res_intg) then
!            res_intg_inf(1) = i
!          end if

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

      ! set up user-defined tallies
      call set_user_tally(reduced_tal(i),Ebins,nbins,react_type,isotope,region,  &
     &                    n_materials,dv)

      ! deallocate Ebins
      if(allocated(Ebins)) deallocate(Ebins)

    end do

    ! load the source
    call load_source(mat(1),source_type,settings_%source_path)

    ! set up spectrum tally
    call set_spectrum_tally(tal(n_tallies-1),emax,emin,n_materials)
    ! set up reduced_tally
    call set_spectrum_tally(reduced_tal(n_tallies-1),emax,emin,n_materials)

    ! set up k_inf tally
    call set_kinf_tally(tal(n_tallies),emax,emin,n_materials)
    ! set up reduced_tally
    call set_kinf_tally(reduced_tal(n_tallies),emax,emin,n_materials)

  end subroutine read_input

end module input
