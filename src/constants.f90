module constants

    implicit none

    real(8), parameter:: ZERO = 0.0_8
    real(8), parameter:: ONE = 1.0_8
    real(8), parameter:: PI = 3.1415926535898_8
    real(8), parameter:: K_BOLTZMANN = 8.617342e-11_8      ! unit: MeV/K
    real(8), parameter:: M_NEUT = 939.565378_8/(3.0e+8)**2  !unit: MeV/((m/s)^2)
    real(8), parameter:: M_NUCLEON = 931.494061_8/(3.0e+8)**2  !unit: MeV/((m/s)^2)

end module constants
