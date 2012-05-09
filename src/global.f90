!==============================================================================!
! MODULE: global 
!
!> @author Bryan Herman
!>
!> @brief Contains all of the global variables
!==============================================================================!

module global

  use materials, only: material_type
  use particle,  only: particle_type
  use tally,     only: tally_type
  use timing,    only: Timer

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

  ! list history input information
  integer :: nhistories
  integer :: seed
  integer :: source_type
  ! number of samples per xs for broadening
  integer :: sample_per_xs=1

  ! list global vars that are set during run
! Commented out by S. Xu
!  integer :: eidx        ! energy index for cross sections
  integer :: n_tallies   ! number of tallies
  integer :: n_materials ! n_materials
  integer :: res_iso     ! resonant isotope id in material 1
  real(8) :: Dancoff     ! lattice Dancoff factor (C)
  real(8) :: radius      ! radius of fuel pin
!  real(8) :: T           ! temperature, added by S. Xu
  logical :: res_intg    ! indicator for resonance integral
  integer :: res_intg_inf(2)=0 ! store the index for micro_capture and flux tally
  real(8), allocatable :: res_intg_rec(:,:) ! to store resonance integral (with variance)

! commented out by S. Xu (for compile dependency problem)
!  ! set max and min energy
!  real(8) :: emin = 1e-11_8
!  real(8) :: emax = 20.0_8

  ! kT value base on 300K
  real(8) :: kT != 8.6173324e-5_8*300*1.0e-6_8

  ! set nu value
  real(8) :: nubar = 2.455_8

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

!  !added by S. Xu (May. 2012)
!  real(8) :: col_kinf_mean_2 = 0.0_8
!  real(8) :: col_kinf_std_2  = 0.0_8
!  real(8) :: reaction_tally(3)

  character(len=255)  :: output_filename
  character(len=255)  :: output_path

contains

!==============================================================================
! ALLOCATE_PROBLEM
!> @brief allocates global variables for calculation
!==============================================================================

  subroutine allocate_problem()

    ! formal variables

    ! allocate tallies
    if (.not.allocated(tal)) allocate(tal(n_tallies))
    if (.not.allocated(mat)) allocate(mat(n_materials))

  end subroutine allocate_problem

!===============================================================================
! DEALLOCATE_PROBLEM
!> @brief deallocates global variables
!===============================================================================

  subroutine deallocate_problem()

    use materials, only: deallocate_material
    use tally,     only: deallocate_tally

    ! local variables
    integer :: i  ! loop counter

    ! deallocate within materials
    do i = 1,n_materials

      ! deallocate material
      call deallocate_material(mat(i))

    end do

    ! deallocate material variable
    if (allocated(mat)) deallocate(mat)

    ! deallocate within tallies
    do i = 1,n_tallies

      ! deallocate tally
      call deallocate_tally(tal(i))

    end do

    ! deallocate tally variable
    if (allocated(tal)) deallocate(tal)

  end subroutine deallocate_problem


! Added by S. Xu
!===============================================================================
! doppler_broaden
!> @brief routine routine to do the on the fly doppler broadening
!===============================================================================

  subroutine doppler_broaden()

    use materials, only: doppler_broaden_xs

    ! local variables
    integer :: i  ! loop counter

    ! begin loop over materals
    do i = 1,n_materials

      ! call routine to compute xs
      call doppler_broaden_xs(mat(i), neut%E, sample_per_xs)

    end do

  end subroutine doppler_broaden


!===============================================================================
! COMPUTE_MACRO_CROSS_SECTIONS
!> @brief routine that handles the call to compute macro cross sections
!===============================================================================

  subroutine compute_macro_cross_sections()

    use materials, only: compute_macroxs

    ! local variables
    integer :: i  ! loop counter

    ! begin loop over materals
    do i = 1,n_materials

      ! call routine to compute xs
      call compute_macroxs(mat(i))

    end do

  end subroutine compute_macro_cross_sections


!===============================================================================
! ADD_TO_TALLIES
!> @brief routine that adds temporary value to tallies
!===============================================================================

  subroutine add_to_tallies()

    use tally, only: add_to_tally

    ! local variables
    integer :: i            ! loop counter
    real(8) :: fact = 1.0_8 ! multiplier factor
    real(8) :: totxs        ! total macroscopic xs of material
    real(8) :: mubar        ! average cosine scattering angle

! Commented out by S. Xu (Apr. 2012)
!    ! compute macroscopic cross section
!    totxs = sum(mat(neut%region)%totalxs(eidx,:))

! Added by S. Xu
    totxs = sum(mat(neut%region)%xs_Total_brdn)

    ! begin loop over tallies
    do i = 1,n_tallies

      ! set multiplier
      select case(tal(i)%react_type)

        ! flux only
        case(0)
          fact = 1.0_8

        ! total 
        case(1)
          fact = totxs

        ! absorption
        case(2)
! Modified by S. Xu
!          fact = sum(mat(neut%region)%absorxs(eidx,:))
          fact = sum(mat(neut%region)%xs_Absb_brdn)

        ! scattering
        case(3)
! Modified by S. Xu
!          fact = sum(mat(neut%region)%scattxs(eidx,:))
          fact = sum(mat(neut%region)%xs_Scat_brdn)

        ! nufission
        case(4)
! Modified by S. Xu
!          fact = nubar*sum(mat(neut%region)%fissixs(eidx,:))
          fact = nubar*sum(mat(neut%region)%xs_Fiss_brdn)

        ! fission
        case(5)
! Modified by S. Xu
!          fact = sum(mat(neut%region)%fissixs(eidx,:))
          fact = sum(mat(neut%region)%xs_Fiss_brdn)

        ! diffusion coefficient
        case(6)
! Modified by S. Xu
!          fact = 1._8/(3._8*sum(mat(neut%region)%transxs(eidx,:)))
          fact = 1._8/(3._8*sum(mat(neut%region)%xs_Tran_brdn))

        ! transport 
        case(7)
! Modified by S. Xu
!          fact = sum(mat(neut%region)%transxs(eidx,:))
          fact = sum(mat(neut%region)%xs_Tran_brdn)

        ! micro capture
        case(8)
          if (neut%region == tal(i)%region) then
! Modified by S. Xu
!            fact = mat(neut%region)%isotopes(tal(i)%isotope)%xs_capt(eidx)
            fact = mat(neut%region)%isotopes(tal(i)%isotope)%xs_capt_brdn

          else
            cycle
          end if 

        ! default is flux tally
        case DEFAULT
          fact = 1.0_8

      end select

      ! call routine to add tally
      call add_to_tally(tal(i),fact,totxs,neut%E,neut%region)

    end do

  end subroutine add_to_tallies

!===============================================================================
! BANK_TALLIES
!> @brief routine that record temporary history information in tallies 
!===============================================================================

  subroutine bank_tallies()

    use tally, only: bank_tally

    ! local variables
    integer :: i  ! loop counter

    ! begin loop over tallies
    do i = 1,n_tallies

      ! call routine to bank tally
      call bank_tally(tal(i))

    end do

  end subroutine bank_tallies

!===============================================================================
! FINALIZE_TALLIES
!> @brief routine that calls another routine to compute tally statistics
!===============================================================================

  subroutine finalize_tallies()

    use tally, only: calculate_statistics

    ! local variables
    integer :: i ! loop counter
    integer :: j ! loop counter

    ! begin loop over tallies
    do i = 1,n_tallies

      ! call routine to compute statistics
      call calculate_statistics(tal(i),nhistories)

      ! normalize by volumes and histories if flux tally
      if (tal(i)%flux_tally .or. tal(i)%dv) then
        do j = 1,n_materials
          tal(i)%mean(:,j) = tal(i)%mean(:,j) / mat(j)%vol
        end do
      end if

    end do

    ! compute k_inf
    ana_kinf_mean = dble(n_fiss)*nubar/dble(nhistories)
    ana_kinf_std  = nubar*sqrt((dble(n_fiss)/dble(nhistories) -                &
   &                (dble(n_fiss)/dble(nhistories))**2)/dble(nhistories-1))
    col_kinf_mean = sum(tal(n_tallies)%mean)
    col_kinf_std  = sum(tal(n_tallies)%std)

! !Added by S. Xu (May. 2012)
!    col_kinf_mean_2 = sum(tal(6)%mean)
!    col_kinf_std_2  = sum(tal(6)%std)
!
!    reaction_tally(1) = sum(tal(2)%mean)
!    reaction_tally(2) = sum(tal(3)%mean)
!    reaction_tally(3) = sum(tal(4)%mean)

  end subroutine finalize_tallies


!===============================================================================
! COMPUTE_RES_INTG
!> @brief routine that computes resonance integrals
!===============================================================================

  subroutine compute_res_intg()

    use tally, only: calculate_statistics

    ! local variables
    integer :: i ! loop counter
    integer :: this_region
    integer :: n_bin

    if (tal(res_intg_inf(1))%nbins == tal(res_intg_inf(2))%nbins) then

      n_bin = tal(res_intg_inf(1))%nbins
      allocate(res_intg_rec(n_bin, 2))
      this_region = tal(res_intg_inf(1))%region

      res_intg_rec(:,1) = log(tal(res_intg_inf(1))%E(2:(n_bin+1))/tal(res_intg_inf(1))%E(1:n_bin)) * &
   &                    tal(res_intg_inf(1))%mean(:,this_region)/tal(res_intg_inf(2))%mean(:,this_region)
      res_intg_rec(:,2) = res_intg_rec(:,1)* &
   &          (tal(res_intg_inf(1))%std(:,this_region)/tal(res_intg_inf(1))%mean(:,this_region)+  &
   &           tal(res_intg_inf(2))%std(:,this_region)/tal(res_intg_inf(2))%mean(:,this_region))
      
    else

      write(*,*) "For computing resonance integral: Energy bins of flux and micro capture are not the same!"   
      stop

    end if

  end subroutine compute_res_intg


!===============================================================================
! WRITE_RES_INTG
!> @brief routine to output resonance integral
!===============================================================================

  subroutine write_res_intg(file_unit)

    integer          :: file_unit
    integer          :: i

      do i = 1,size(res_intg_rec,1)
        write(file_unit,'(es9.2e2, " --", es9.2e2, 5x, f10.5, 1x, "+/-", 1x, f10.5)') &
            & tal(res_intg_inf(1))%E(i),tal(res_intg_inf(1))%E(i+1),res_intg_rec(i,1),res_intg_rec(i,2)
      end do

  end subroutine write_res_intg

end module global
