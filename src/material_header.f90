!==============================================================================!
! MODULE: material_header
!
!> @author Bryan Herman
!>
!> @brief Contains information about the isotopics of problem
!==============================================================================!

module material_header

  implicit none

  type :: n_source_type

    real(8), allocatable :: E(:)      ! energy range for fission source
    real(8)              :: cdf_width ! width of cdf bins from 0 to 1

  end type n_source_type

  type :: thermal_type

    integer               :: kTsize    ! size of kT vector
    integer               :: cdfsize   ! size of cdf
    real(8), allocatable  :: kTvec(:)  ! vector of kT values
    real(8), allocatable  :: Erat(:,:) ! energy
    real(8)               :: cdf_width ! width of cdf interval from 0 to 1 
    
  end type thermal_type

  type :: iso_type

    real(8)               :: N           ! number density
    real(8)               :: A           ! atomic weight
    real(8)               :: alpha       ! (A-1)^2/(A+1)^2
    real(8)               :: mubar       ! average cosine scattering angle
    real(8), allocatable  :: xs_capt(:)  ! capture micro xs
    real(8), allocatable  :: xs_scat(:)  ! scattering micro xs
    real(8), allocatable  :: xs_fiss(:)  ! fission micro xs
    character(len=255)    :: name        ! name of isotope

    logical               :: thermal     ! thermal scatterer
    type(thermal_type)    :: thermal_lib ! thermal library

    
    ! added by S. Xu (Apr. 2012)
    real(8)                     :: alpha_MB      ! for doppler braodening
    real(8)                     :: xs_capt_brdn  ! capture xs after doppler broadening
    real(8)                     :: xs_scat_brdn  ! scattering xs after doppler broadening
    real(8)                     :: xs_fiss_brdn  ! fission xs after doppler broadening
    real(8), allocatable  :: engy_capt(:)  ! energy grid for capture xs
    real(8), allocatable  :: engy_scat(:)  ! energy grid for scattering xs
    real(8), allocatable  :: engy_fiss(:)  ! energy grid for fission xs
    integer               :: capt_size     ! size of xs_capt
    integer               :: scat_size     ! size of xs_scat
    integer               :: fiss_size     ! size of xs_fiss
    logical               :: doppler ! indicte whether doppler braodening is performed

  end type iso_type

  type, public ::  material_type

    type(n_source_type)           :: source        ! the source of neutrons
    type(iso_type), allocatable :: isotopes(:)   ! 1-D array of isotopes in mat
    integer                     :: nisotopes     ! number of isotopes in mat
    integer                     :: curr_iso      ! the current isotope
    integer                     :: npts          ! number of points in energy
    real(8)                     :: E_width       ! width of energy interval
    real(8)                     :: E_min         ! min energy
    real(8)                     :: E_max         ! max energy
    real(8)                     :: vol           ! volume of region
    real(8), allocatable        :: totalxs(:,:)  ! array of macroscopic tot xs
    real(8), allocatable        :: scattxs(:,:)  ! array of macroscopic scat xs
    real(8), allocatable        :: absorxs(:,:)  ! array of macroscopic abs xs
    real(8), allocatable        :: captuxs(:,:)  ! array of macroscopic capt xs
    real(8), allocatable        :: fissixs(:,:)  ! array of macroscopic fiss xs
    real(8), allocatable        :: transxs(:,:)  ! array of macro transport xs

    ! added by S. Xu (Apr. 2012)
    real(8), allocatable        :: xs_Total_brdn(:)  ! array of total macro xs for each isotope after doppler broadening
    real(8), allocatable        :: xs_Capt_brdn(:)  ! macro capture xs after doppler broadening
    real(8), allocatable        :: xs_Scat_brdn(:)  ! macro scattering xs after doppler broadening
    real(8), allocatable        :: xs_Fiss_brdn(:)  ! macro fission xs after doppler broadening
    real(8), allocatable        :: xs_Absb_brdn(:)  ! macro absorb xs after doppler broadening
    real(8), allocatable        :: xs_Tran_brdn(:)  ! macro transport xs after doppler broadening

  end type material_type

end module material_header
