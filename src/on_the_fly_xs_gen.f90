!==============================================================================!
! MODULE: on_the_fly_xs_gen
!
!> @author S. Xu
!>
!> @brief Provide on the fly xs generetion
!==============================================================================!

module on_the_fly_xs_gen

  implicit none
  private
  public :: energy_doppler_broadened


contains

!===============================================================================
! doppler_broaden
!> @brief routine routine to do the on the fly doppler broadening
!===============================================================================
!  subroutine doppler_broaden()
!
!    use global, only: mat, n_materials, neut
!
!    ! local variables
!    integer :: i, j  ! loop counter
! !    type(iso_type), pointer :: iso ! pointer to current isotope
!
!    ! begin loop over materals
!    do j = 1,n_materials
!
!      do i = 1, mat(j)%nisotopes
!
!print *, 'in doppler_broaden'
!print *, j, i
!print *, mat(j)%isotopes(i)%engy_capt
!
!        call xs_gen(neut%E, mat(j)%isotopes(i)%engy_capt, mat(j)%isotopes(i)%xs_capt, &
!             & mat(j)%isotopes(i)%alpha_MB, mat(j)%isotopes(i)%xs_capt_brdn)
!        call xs_gen(neut%E, mat(j)%isotopes(i)%engy_scat, mat(j)%isotopes(i)%xs_scat, &
!             & mat(j)%isotopes(i)%alpha_MB, mat(j)%isotopes(i)%xs_scat_brdn)
!        call xs_gen(neut%E, mat(j)%isotopes(i)%engy_fiss, mat(j)%isotopes(i)%xs_fiss, &
!             & mat(j)%isotopes(i)%alpha_MB, mat(j)%isotopes(i)%xs_fiss_brdn)
!
!      end do
!
!    end do
!
!  end subroutine doppler_broaden


!===============================================================================
! xs_gen
!> @brief routine to generate on the fly xs
!===============================================================================

  subroutine energy_doppler_broadened(v, alpha_MB, v_brdn)
    
    use constants, only: M_NEUT, PI
    use global,    only: vmin

    ! local variables
    real(8), intent(in)        :: v         ! incident neutron velocity
    real(8), intent(in)        :: alpha_MB  ! for doppler broadening
    real(8),  intent(out)       :: v_brdn     ! relative velocity
    real(8)  :: v_target  ! target velocity
!    real(8)  :: theta     ! angle between two velocities
    real(8)  :: mu         ! cosine of the angle between two velocities
    real(8)  :: rn        ! random number
    real(8)  :: v_most_prob ! most probable velocity
    real(8)  :: uplimit   ! uplimit for rejection sampling
    
    ! sample relative
    rn = rand()
!    theta = PI*rn
    mu = 2._8*rn-1._8
    
    ! sample target velocity
    v_most_prob = 1._8/sqrt(alpha_MB)
    uplimit = 1.0001_8*pdf_MB(v_most_prob, alpha_MB)

    do while (.true.)
      rn = rand()
      v_target = 20.0_8*v_most_prob*rn  ! sample a velocity interval of [0, 20]*v_most_probable
      rn = rand()
      if (uplimit*rn < pdf_MB(v_target, alpha_MB)) exit
    end do

    v_brdn = sqrt(v**2+v_target**2-2._8*v*v_target*mu)
    if (v_brdn < vmin) then
      v_brdn = vmin
    end if


  end subroutine energy_doppler_broadened


!  subroutine energy_doppler_broadened(E, alpha_MB, E_brdn)
!
!    use constants, only: M_NEUT, PI
!    use parameters, only: emin
!
!    ! local variables
!    real(8), intent(in)        :: E         ! incident neutron energy
!    real(8), intent(in)        :: alpha_MB  ! for doppler broadening
!    real(8),  intent(out)       :: E_brdn     ! neutron energy after broadening
!    real(8)  :: v_inc     ! incident neutron velocity
!    real(8)  :: v_target  ! target velocity
! !   real(8)  :: v_rel     ! relative velocity in the direction of neutron velocity
! !   real(8)  :: theta     ! angle between two velocities
!    real(8)  :: mu         ! cosine of the angle between two velocities
!    real(8)  :: rn        ! random number
!    real(8)  :: v_most_prob ! most probable velocity
!    real(8)  :: uplimit   ! uplimit for rejection sampling
!
!
!    v_inc = sqrt(2*E/M_NEUT)
!
!    ! sample relative
!    rn = rand()
! !   theta = PI*rn
!    mu = 2._8*rn-1._8
!
!    ! sample target velocity
!    v_most_prob = 1._8/sqrt(alpha_MB)
!    uplimit = 1.0001_8*pdf_MB(v_most_prob, alpha_MB)
!
!    do while (.true.)
!      rn = rand()
!      v_target = 20.0_8*v_most_prob*rn  ! sample a velocity interval of [0, 20]*v_most_probable
!      rn = rand()
!      if (uplimit*rn < pdf_MB(v_target, alpha_MB)) exit
!    end do
!
!    E_brdn = 0.5_8*M_NEUT*((v_inc-v_target*mu)**2.+v_target**2.*(1._8-mu**2.))
!    if (E_brdn < emin) then
!      E_brdn = emin
!    end if
!
!
!  end subroutine energy_doppler_broadened

!===============================================================================
! pdf_MB
!> @brief Calculate the M-B probabilty density function
!===============================================================================
  function pdf_MB(v, alpha_MB) result(prob_v)

    use constants, only: PI

    real(8)  :: alpha_MB
    real(8)  :: v
    real(8)  :: prob_v

    prob_v = 4._8/sqrt(PI)*alpha_MB**(1.5_8)*v**2._8*exp(-alpha_MB*v**2)

  end function pdf_MB

end module on_the_fly_xs_gen
