module parameters
use constants, only: M_NEUT
    implicit none

    real(8)      :: T   ! temperature

    ! set max and min energy
    real(8) :: emin = 1e-11_8
    real(8) :: emax = 20.0_8
    real(8) :: vmin = sqrt(2._8*1e-11_8/M_NEUT)
    real(8) :: vmax = sqrt(2._8*20._8/M_NEUT)

end module parameters
