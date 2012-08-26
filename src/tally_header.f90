module tally_header


  implicit none

  type, public :: tally_type

    real(8), allocatable :: E(:)      ! user defined energy structure
    real(8), allocatable :: val(:,:)    ! the temporary value
    real(8), allocatable :: sum0(:,:)    ! the sum for the mean and var
    real(8), allocatable :: sum0_sq(:,:) ! the sum for the variable
    real(8), allocatable :: mean(:,:)   ! mean of tallies
    real(8), allocatable :: std(:,:)    ! standard deviation of tallies
    logical :: flux_tally = .false.   ! is this the flux tally
    integer :: nbins                  ! number of tally regions
    real(8) :: width                  ! the uniform width
    real(8) :: emax                   ! max e
    real(8) :: emin                   ! min e
    integer :: react_type             ! reaction type id
    integer :: isotope                ! isotope number
    integer :: region                 ! region for micros
    logical :: dv = .false.           ! divide by volume of regions

  end type tally_type

end module tally_header
