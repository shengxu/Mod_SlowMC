module execute
  
  use global

  implicit none

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

    if (.not.allocated(reduced_tal)) allocate(reduced_tal(n_tallies))

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

      call deallocate_tally(reduced_tal(i))

    end do

    ! deallocate tally variable
    if (allocated(tal)) deallocate(tal)

  end subroutine deallocate_problem


!===============================================================================
! COMPUTE_RES_INTG
!> @brief routine that computes resonance integrals
!===============================================================================

!  subroutine compute_res_intg()

!    use tally, only: calculate_statistics

!    ! local variables
!    integer :: i ! loop counter
!    integer :: this_region
!    integer :: n_bin

!    if (tal(res_intg_inf(1))%nbins == tal(res_intg_inf(2))%nbins) then

!      n_bin = tal(res_intg_inf(1))%nbins
!      allocate(res_intg_rec(n_bin, 2))
!      this_region = tal(res_intg_inf(1))%region

!      res_intg_rec(:,1) = log(tal(res_intg_inf(1))%E(2:(n_bin+1))/tal(res_intg_inf(1))%E(1:n_bin)) * &
!   &                    tal(res_intg_inf(1))%mean(:,this_region)/tal(res_intg_inf(2))%mean(:,this_region)
!      res_intg_rec(:,2) = res_intg_rec(:,1)* &
!   &          (tal(res_intg_inf(1))%std(:,this_region)/tal(res_intg_inf(1))%mean(:,this_region)+  &
!   &           tal(res_intg_inf(2))%std(:,this_region)/tal(res_intg_inf(2))%mean(:,this_region))
!      
!    else

!      write(*,*) "For computing resonance integral: Energy bins of flux and micro capture are not the same!"   
!      stop

!    end if

!  end subroutine compute_res_intg


!===============================================================================
! WRITE_RES_INTG
!> @brief routine to output resonance integral
!===============================================================================

!  subroutine write_res_intg(file_unit)

!    integer          :: file_unit
!    integer          :: i

!      do i = 1,size(res_intg_rec,1)
!        write(file_unit,'(es9.2e2, " --", es9.2e2, 5x, f10.5, 1x, "+/-", 1x, f10.5)') &
!            & tal(res_intg_inf(1))%E(i),tal(res_intg_inf(1))%E(i+1),res_intg_rec(i,1),res_intg_rec(i,2)
!      end do

!  end subroutine write_res_intg

end module execute
