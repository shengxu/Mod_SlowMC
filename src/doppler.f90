module doppler

  use constants, only: ZERO, ONE, PI, K_BOLTZMANN
  use BiSearch  ! binary search for the interval
!  use proj_record    ! count the time that each while loop is executed

  implicit none

  real(8), parameter :: sqrt_pi_inv = ONE / sqrt(PI)

contains

!===============================================================================
! BROADEN takes a microscopic cross section at a temperature T_1 and Doppler
! broadens it to a higher temperature T_2 based on a method originally developed
! by Cullen and Weisbin (see "Exact Doppler Broadening of Tabulated Cross
! Sections," Nucl. Sci. Eng. 60, 199-229 (1976)). The only difference here is
! the F functions are evaluated based on complementary error functions rather
! than error functions as is done in the BROADR module of NJOY.
!===============================================================================

  subroutine broaden(energy, xs, A_target, T, energypt, sigmapt)

    real(8), dimension(:), intent(in)  :: energy   ! energy grid
    real(8), dimension(:), intent(in)  :: xs       ! unbroadened cross section
    real(8), intent(in)  :: A_target    ! mass number of target
    real(8), intent(in)  :: T           ! temperature (difference)
    real(8), intent(in)  :: energypt     ! input energy where broadened xs is to be evaluated
    ! real(8), intent(out) :: sigmaNew(:) ! broadened cross section
    real(8), intent(out) :: sigmapt    ! broadened cross section at one point

    real(8)              :: intg_width = 4.0_8  ! parameter for integration width

    ! Added by S. Xu (Mar. 2012)
    integer              :: energyintv ! the interval that the input energypt falls in
    logical              :: stat      ! whether the enrgypt is at a grid point or not
    real(8)              :: xspt

    integer              :: k     ! loop indices
    integer              :: n        ! number of energy points
    real(8)              :: F_a(0:4) ! F(a) functions as per C&W
    real(8)              :: F_b(0:4) ! F(b) functions as per C&W
    real(8)              :: H(0:4)   ! H functions as per C&W
    ! real(8), allocatable :: x(:)     ! proportional to relative velocity
    real(8)              :: y        ! proportional to neutron velocity
    real(8)              :: y_sq     ! y**2
    real(8)              :: y_inv    ! 1/y
    real(8)              :: y_inv_sq ! 1/y**2
    real(8)              :: alpha    ! constant equal to A/kT
    real(8)              :: slope    ! slope of xs between adjacent points
    real(8)              :: Ak, Bk   ! coefficients at each point
    real(8)              :: a, b     ! values of x(k)-y and x(k+1)-y
    real(8)              :: xk, xkp1 ! values of x(k) and x(k+1), added by S. Xu(Mar. 2012)



    ! Determine alpha parameter -- have to convert k to MeV/K
    alpha = A_target/(K_BOLTZMANN * T)

    ! Allocate memory for x and assign values
    n = size(energy)

!    allocate(x(n))
!    x = sqrt(alpha * energy)

       sigmapt    = ZERO
       ! x = sqrt(alpha *energypt)
       y        = sqrt(alpha *energypt)
       y_sq     = y*y
       y_inv    = ONE / y
       y_inv_sq = y_inv / y

!       do k=1,3
!          count_0(k) = 0
!       end do

!        print *, energy
       ! find the interval that energypt falls in
       call binarySearch(energy, energypt, energyintv, stat)
!       write(11, *), n, energypt, energyintv, stat

       ! =======================================================================
       ! EVALUATE FIRST TERM FROM x(k) - y = 0 to -4
       k = energyintv
       a = ZERO
       call calculate_F(F_a, ZERO)

       if (k<1) then  ! energypt below the lowest energy grid
          ! Since x = 0, this implies that a = -y
          F_b = F_a
          a = -y

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          H = F_a - F_b

          ! Add contribution to broadened cross section, assume 1/v shape
          xspt = sqrt(energy(1)/energypt)*xs(1)
          sigmapt = sigmapt + xspt*y*(y_inv_sq*H(1) + y_inv*H(0))

       else if (k>=n) then   ! energypt above the highest energy grid
          ! Assume xs constant here, therefore, slope is zero.
          F_b = F_a
          xk = sqrt(alpha*energy(k))
          a = xk - y
          call calculate_F(F_a, a)
          H = F_a - F_b

          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          sigmapt = sigmapt + Ak*xs(n)


       else if (.not. stat) then  ! energypt falls in the energy range, but not on a grid point

          F_b = F_a
          xkp1 = sqrt(alpha*energy(k+1))
          xk = sqrt(alpha*energy(k))
          a = xk - y
          ! Calculate F and H functions
          call calculate_F(F_a, a)
          H = F_a - F_b

          ! Calculate A(k), B(k), and slope terms
          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)

          slope = (xs(k+1)-xs(k))/(xkp1**2-xk**2)
          ! interpolate to find the xs corresponding to energypt
          xspt = xs(k) + slope*(y_sq-xk**2)

          ! Add contribution to broadened cross section
          sigmapt = sigmapt + Ak*(xs(k) - slope*xk**2) + slope*Bk

       else     !energy falls on a grid point (except for the highest energy grid)

          xk = sqrt(alpha*energy(k))

       end if

!       write(11, *) 'after first if ', sigmapt

       do while (a >= -intg_width .and. k > 1)

          if (energy(k-1) == energy(k)) then
              k = k - 1
              cycle
          end if

          ! Move to next point
          F_b = F_a
          xkp1 = xk
          k = k - 1
          xk = sqrt(alpha*energy(k))
          a = xk - y

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          H = F_a - F_b

          ! Calculate A(k), B(k), and slope terms
          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)
          slope = (xs(k+1) - xs(k)) / (xkp1**2 - xk**2)

          ! Add contribution to broadened cross section
          sigmapt = sigmapt + Ak*(xs(k) - slope*xk**2) + slope*Bk

!          count_0(1) = count_0(1) + 1
       end do

!       write(11, *) 'after first do loop ', sigmapt
       ! =======================================================================
       ! EXTEND CROSS SECTION TO 0 ASSUMING 1/V SHAPE

!       if (k == 1 .and. a >= -4.0) then
       if (k == 1) then
          ! Since x = 0, this implies that a = -y
          F_b = F_a
          a = -y
          xkp1 = xk

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          H = F_a - F_b

          ! Add contribution to broadened cross section
          sigmapt = sigmapt + xs(k)*xkp1*(y_inv_sq*H(1) + y_inv*H(0))
       end if

!       write(11, *) 'after 1/V ', sigmapt
       ! =======================================================================
       ! EVALUATE FIRST TERM FROM x(k) - y = 0 to +4

       k = energyintv
       b = ZERO
       call calculate_F(F_b, ZERO)

       if (k<1) then  ! energypt below the lowest energy grid

          F_a = F_b
          xkp1 = sqrt(alpha*energy(1))
          b = xkp1 - y

          ! Calculate F and H functions
          call calculate_F(F_b, b)
          H = F_a - F_b

          ! Add contribution to broadened cross section
          ! xspt already evaluated above
          ! xspt = sqrt(energy(1)/energypt)*xs(1)
          sigmapt = sigmapt + xspt*y*(y_inv_sq*H(1) + y_inv*H(0))

          k = k + 1

       else if (k>=n) then   ! energypt above the highest energy grid
          ! Assume xs constant here, therefore, slope is zero.
          a = ZERO
          call calculate_F(F_a, a)

          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          sigmapt = sigmapt + Ak*xs(n)

          k = k + 1

       else if (.not. stat) then  ! energypt falls in the energy range, but not on a grid point

          F_a = F_b
          xk = sqrt(alpha*energy(k))
          xkp1 = sqrt(alpha*energy(k+1))
          b = xkp1 - y
          ! Calculate F and H functions
          call calculate_F(F_b, b)
          H = F_a - F_b

          ! Calculate A(k), B(k), and slope terms
          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)

          slope = (xs(k+1)-xs(k))/(xkp1**2-xk**2)
          ! interpolate to find the xs corresponding to energypt
          ! xspt already evaluated above
          ! xspt = xs(k) + slope*(y_sq-xk**2)

          ! Add contribution to broadened cross section
          sigmapt = sigmapt + Ak*(xspt - slope*y_sq) + slope*Bk

          k = k + 1

       else     !energy falls on a grid point (except for the highest energy grid)

          xkp1 = sqrt(alpha*energy(k))

       end if

!       write(11, *) 'after second if ', sigmapt

       do while (b <= intg_width .and. k < n)

          if (energy(k+1) == energy(k)) then
              k = k + 1
              cycle
          end if

          ! Move to next point
          F_a = F_b
          xk = xkp1
          xkp1 = sqrt(alpha*energy(k+1))
          b = xkp1 - y

          ! Calculate F and H functions
          call calculate_F(F_b, b)
          H = F_a - F_b

          ! Calculate A(k), B(k), and slope terms
          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)
          slope = (xs(k+1) - xs(k)) / (xkp1**2 - xk**2)

          ! Add contribution to broadened cross section
          sigmapt = sigmapt + Ak*(xs(k) - slope*xkp1**2) + slope*Bk

          k = k + 1

!          count_0(2) = count_0(2) + 1
       end do

!       write(11, *) 'after second do ', sigmapt

       ! =======================================================================
       ! EXTEND CROSS SECTION TO INFINITY ASSUMING CONSTANT SHAPE

       if (k == n .and. b <= 4.0) then
          ! Calculate F function at last energy point
          xk = sqrt(alpha*energy(k))
          a = xk - y
          call calculate_F(F_a, a)

          ! Add contribution to broadened cross section
          sigmapt = sigmapt + xs(k) * (y_inv_sq*F_a(2) + 2.0*y_inv*F_a(1) + F_a(0))
       end if

!       write(11, *) 'after constant xs ', sigmapt

       ! =======================================================================
       ! EVALUATE SECOND TERM FROM x(k) + y = 0 to +4

       if (y <= intg_width) then
          ! Swap signs on y
          y = -y
          y_inv = -y_inv
          k = 0
          xkp1 = sqrt(alpha*energy(k+1))

          ! Calculate a and b based on 0 and x(1)
          a = -y
          b = xkp1 - y

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          call calculate_F(F_b, b)
          H = F_a - F_b

          ! Add contribution to broadened cross section
          sigmapt = sigmapt - xs(k) * xkp1 * (y_inv_sq*H(1) + y_inv*H(0))

!          write(11, *) 'after i/V for second term ', sigmapt

          k = k + 1

          ! Now progress forward doing the remainder of the second term
          do while (b <= intg_width)

             if (energy(k-1) == energy(k)) then
                k = k + 1
                cycle
             end if

             ! Move to next point
             F_a = F_b
             xk = xkp1
             xkp1 = sqrt(alpha*energy(k+1))
             b = xkp1 - y

             ! Calculate F and H functions
             call calculate_F(F_b, b)
             H = F_a - F_b

             ! Calculate A(k), B(k), and slope terms
             Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
             Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)
             slope = (xs(k+1) - xs(k)) / (xkp1**2 - xk**2)

             ! Add contribution to broadened cross section
             sigmapt = sigmapt - Ak*(xs(k) - slope*xk**2) - slope*Bk

             k = k + 1

!             count_0(3) = count_0(3) + 1
          end do

!          write(11, *) 'after do loop for second term ', sigmapt
       end if

  end subroutine broaden

!===============================================================================
! CALCULATE_F evaluates the function:
!
!    F(n,a) = 1/sqrt(pi)*int(z^n*exp(-z^2), z = a to infinity)
!
! The five values returned in a vector correspond to the integral for n = 0
! through 4. These functions are called over and over during the Doppler
! broadening routine.
!===============================================================================

  subroutine calculate_F(F, a)

    real(8), intent(inout) :: F(0:4)
    real(8), intent(in)    :: a

    F(0) = 0.5*erfc(a)
    F(1) = 0.5*sqrt_pi_inv*exp(-a*a)
    F(2) = 0.5*F(0) + a*F(1)
    F(3) = F(1)*(1.0 + a*a)
    F(4) = 0.75*F(0) + F(1)*a*(1.5 + a*a)

  end subroutine calculate_F

!#ifdef NO_F2008
!===============================================================================
! ERFC computes the complementary error function of x
!===============================================================================

!  function erfc(x) result(y)
!
!    real(8), intent(in) :: x
!    real(8)             :: y
!
!    real(8) :: a1 =  0.254829592_8
!    real(8) :: a2 = -0.284496736_8
!    real(8) :: a3 =  1.421413741_8
!    real(8) :: a4 = -1.453152027_8
!    real(8) :: a5 =  1.061405429_8
!    real(8) :: p  =  0.3275911_8
!    real(8) :: t
!
!    ! Abramowitz and Stegun formula 7.1.26
!    t = 1.0_8/(1.0_8 + p*abs(x))
!    y = (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
!
!    ! Account for negative values of x
!    y = sign(y,x)
!
!  end function erfc
!#endif

end module doppler
