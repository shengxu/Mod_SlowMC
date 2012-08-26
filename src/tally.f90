!==============================================================================!
! MODULE: tally 
!
!> @author Bryan Herman
!>
!> @brief Contains information about tallying quantities 
!==============================================================================!

module tally

!  use tally_header, only: tally_type
  use global

  implicit none

contains

!===============================================================================
! SET_USER_TALLY
!> @brief routine to intialize user-defined tallies
!===============================================================================

  subroutine set_user_tally(this,Ebins,nbins,react_type,isotope,region,n_materials,&
                            dv)

    ! formal variables
    type(tally_type) :: this        ! a tally
    integer          :: nbins       ! size of energy bins
    integer          :: react_type  ! reaction type
    integer          :: isotope     ! isotope for multiplier
    integer          :: region      ! region filter for isotope
    integer          :: n_materials ! number of material tally regions
    real(8)          :: Ebins(nbins+1)  ! vector of energy bins
    logical          :: dv          ! divide by volume   
 
    ! preallocate user-defined energy structure
    if (.not.allocated(this%E)) allocate(this%E(nbins+1))

    ! set divide by volume
    this%dv = dv

    ! set energy structure
    this%E = Ebins

    ! set number of energy bins
    this%nbins = nbins

    ! set reaction type
    this%react_type = react_type

    ! set isotope
    this%isotope = isotope

    ! set region
    this%region = region

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(nbins,n_materials))
    if(.not.allocated(this%sum0)) allocate(this%sum0(nbins,n_materials))
    if(.not.allocated(this%sum0_sq)) allocate(this%sum0_sq(nbins,n_materials))

    ! preallocate mean and stdev
    if (.not.allocated(this%mean)) allocate(this%mean(nbins,n_materials))
    if (.not.allocated(this%std))  allocate(this%std(nbins,n_materials))

    ! zero out tallies
    this%val = 0.0_8
    this%sum0 = 0.0_8
    this%sum0_sq = 0.0_8

  end subroutine set_user_tally

!===============================================================================
! SET_SPECTRUM_TALLY
!> @brief routine to initialize all tallies 
!===============================================================================

  subroutine set_spectrum_tally(this,emax,emin,n_materials)

    ! formal variables
    type(tally_type) :: this         ! a tally
    integer          :: n_materials  ! number of materials
    real(8)          :: emax         ! max e
    real(8)          :: emin         ! min e

    ! set up automatic flux tally
    this%flux_tally = .true.
    this%nbins = 5000
    this%emax = emax
    this%emin = emin
    this%width = (log10(emax) - log10(emin))/dble(this%nbins)

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(5000,n_materials))
    if(.not.allocated(this%sum0)) allocate(this%sum0(5000,n_materials))
    if(.not.allocated(this%sum0_sq)) allocate(this%sum0_sq(5000,n_materials))

    ! preallocate mean and stdev
    if (.not.allocated(this%mean)) allocate(this%mean(5000,n_materials))
    if (.not.allocated(this%std))  allocate(this%std(5000,n_materials))

    ! zero out tallies
    this%val = 0.0_8
    this%sum0 = 0.0_8
    this%sum0_sq = 0.0_8

  end subroutine set_spectrum_tally

!===============================================================================
! SET_KINF_TALLY
!> @brief routine to initiaize kinf nu-fission tally
!===============================================================================

  subroutine set_kinf_tally(this,emax,emin,n_materials)

    ! formal variables
    type(tally_type) :: this         ! a tally
    integer          :: n_materials  ! number of materials
    real(8)          :: emax         ! max e
    real(8)          :: emin         ! min e

    ! preallocate user-defined energy structure
    if (.not.allocated(this%E)) allocate(this%E(2))

    ! set energy structure
    this%E(1) = emin
    this%E(2) = emax
    this%nbins = 1

    ! set reaction type
    this%react_type = 4

    ! preallocate vectors
    if(.not.allocated(this%val)) allocate(this%val(1,n_materials))
    if(.not.allocated(this%sum0)) allocate(this%sum0(1,n_materials))
    if(.not.allocated(this%sum0_sq)) allocate(this%sum0_sq(1,n_materials))

    ! preallocate mean and stdev
    if (.not.allocated(this%mean)) allocate(this%mean(1,n_materials))
    if (.not.allocated(this%std))  allocate(this%std(1,n_materials))

    ! zero out tallies
    this%val = 0.0_8
    this%sum0 = 0.0_8
    this%sum0_sq = 0.0_8

  end subroutine set_kinf_tally


!===============================================================================
! ADD_TO_TALLIES
!> @brief routine that adds temporary value to tallies
!===============================================================================

  subroutine add_to_tallies()


    ! local variables
    integer :: i            ! loop counter
    real(8) :: fact = 1.0_8 ! multiplier factor
    real(8) :: totxs        ! total macroscopic xs of material
    real(8) :: mubar        ! average cosine scattering angle

!    ! compute macroscopic cross section
!    totxs = sum(mat(neut%region)%totalxs(eidx,:))

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
!          fact = sum(mat(neut%region)%absorxs(eidx,:))
          fact = sum(mat(neut%region)%xs_Absb_brdn)

        ! scattering
        case(3)
!          fact = sum(mat(neut%region)%scattxs(eidx,:))
          fact = sum(mat(neut%region)%xs_Scat_brdn)

        ! nufission
        case(4)
!          fact = nubar*sum(mat(neut%region)%fissixs(eidx,:))
          fact = nubar*sum(mat(neut%region)%xs_Fiss_brdn)

        ! fission
        case(5)
!          fact = sum(mat(neut%region)%fissixs(eidx,:))
          fact = sum(mat(neut%region)%xs_Fiss_brdn)

        ! diffusion coefficient
        case(6)
!          fact = 1._8/(3._8*sum(mat(neut%region)%transxs(eidx,:)))
          fact = 1._8/(3._8*sum(mat(neut%region)%xs_Tran_brdn))

        ! transport 
        case(7)
!          fact = sum(mat(neut%region)%transxs(eidx,:))
          fact = sum(mat(neut%region)%xs_Tran_brdn)

        ! micro capture
        case(8)
          if (neut%region == tal(i)%region) then
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
! ADD_TO_TALLY
!> @brief routine to add quantities during transport of a particle 
!===============================================================================

  subroutine add_to_tally(this,fact,totxs,E,region)


    ! formal variables
    type(tally_type) :: this     ! a tally 
    integer          :: region   ! region id
    real(8)          :: fact     ! multiplier for tally
    real(8)          :: totxs    ! totalxs
    real(8)          :: E        ! neutron energy

    ! local variables
    integer :: i      ! iteration counter
    integer :: idx=0  ! index in tally grid
    logical :: stat_bisch ! indicator for binary search

    ! use uniform grid sampling if flux tally
    if (this%flux_tally) then

      ! calculate index
      idx = ceiling((log10(E) - log10(this%emin))/this%width)

    else

      ! check for output bounds
      if (E < minval(this%E) .or. E > maxval(this%E)) then
        print *,'FATAL: Energy out of tally bounds, check input bounds!'
        stop
      end if


      ! begin loop around energy vector to get index
      do i = 1,size(this%E)
        if (E < this%E(i)) then
          idx = i - 1
          exit
        end if
      end do

!    !use binary search instead (Added by S. Xu, Apr. 2012)
!    call binarySearch(this%E, E, idx, stat_bisch)

! Added by S. Xu
! write(11, '(e19.8e3, 3x, i3)') E, this%react_type
    end if

    ! add to tally
    if (idx /= 0) this%val(idx,region) = this%val(idx,region) + fact/totxs

  end subroutine add_to_tally

!===============================================================================
! BANK_TALLIES
!> @brief routine that record temporary history information in tallies 
!===============================================================================

  subroutine bank_tallies()

!    ! formal variables
!    type(tally_type) :: this ! a tally

    integer  :: i

    do i = 1, n_tallies
        ! record to sums
        tal(i)%sum0    = tal(i)%sum0    + tal(i)%val
        tal(i)%sum0_sq = tal(i)%sum0_sq + tal(i)%val**2

!#ifdef DEBUG
!    if (master) then
!        write(404, *) tal(i)%val
!        write(404, *) tal(i)%sum0
!    end if
!#endif

        ! zero out temp value
        tal(i)%val = 0.0_8

    end do

  end subroutine bank_tallies


!===============================================================================
! REDUCE_TALLIES
!> @brief routine to reduce tally across different processors
!===============================================================================

  subroutine reduce_tallies()


#ifdef MPI
    use mpi
#endif

    ! local variable
    integer  :: i
    integer  :: data_count

#ifdef MPI    

#ifdef DEBUG
  write(404, *) 'n_tallies = ', n_tallies, ' from processor ', rank
#endif
    do i = 1, n_tallies
        data_count = n_materials*tal(i)%nbins

!#ifdef DEBUG
!  write(*, *) 'i = ', i, 'n_materials = ', n_materials, 'tal(i)%nbins = ', tal(i)%nbins, "from processor ", rank
!  write(*, *) 'i = ', i, 'data_count = ', data_count, "from processor ", rank
!#endif

#ifdef DEBUG
  write(404, *) 'tal(', i, ')%sum0 : ', tal(i)%sum0,  "from processor ", rank
  write(404, *) 'reduced_tal(', i, ')%sum0: ', reduced_tal(i)%sum0, ' from processor ', rank
#endif

        call MPI_REDUCE(tal(i)%sum0, reduced_tal(i)%sum0, data_count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

    if (mpi_err /= MPI_SUCCESS) then
       print *, "Failed to reduce tal(", i, ")%sum0."
       stop
    end if

        call MPI_REDUCE(tal(i)%sum0_sq, reduced_tal(i)%sum0_sq, data_count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

    if (mpi_err /= MPI_SUCCESS) then
       print *, "Failed to reduce tal(", i, ")%sum0_sq."
       stop
    end if



    end do

    call MPI_REDUCE(n_abs, reduced_n_abs, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_REDUCE(n_fiss, reduced_n_fiss, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

#else
    do i = 1, n_tallies
        reduced_tal(i)%sum0 = tal(i)%sum0
        reduced_tal(i)%sum0_sq = tal(i)%sum0_sq
    end do

    reduced_n_abs = n_abs
    reduced_n_fiss = n_fiss
#endif

#ifdef DEBUG

!#ifdef MPI
!  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!#endif

  write(404, *) 'tal(1)%sum0 : ', tal(1)%sum0,  "from processor ", rank
  if (master) then

      write(404, *) 'reduced_tal(1)%sum0: ', reduced_tal(1)%sum0, ' from processor ', rank

      write(404, *) 'n_abs = ', n_abs, ' reduced_n_abs = ', reduced_n_abs, "from processor ", rank
  end if
#endif

  end subroutine reduce_tallies


!===============================================================================
! FINALIZE_TALLIES
!> @brief routine that calls another routine to compute tally statistics
!===============================================================================

  subroutine finalize_tallies()

    ! local variables
    integer :: i ! loop counter
    integer :: j ! loop counter


!#ifdef DEBUG
!  write(404, *) ' nhistories = ', nhistories, ' n_tallies = ', n_tallies, "from processor ", rank
!#endif

    ! begin loop over tallies
    do i = 1, 1

!#ifdef DEBUG
!  write(404, *) 'reduced_tal(', i, ')%sum0 = ', reduced_tal(i)%sum0, ' nhistories = ', nhistories, "from processor ", rank
!#endif
      ! call routine to compute statistics
      call calculate_statistics(reduced_tal(i),nhistories)

      ! normalize by volumes and histories if flux tally
      if (reduced_tal(i)%flux_tally .or. reduced_tal(i)%dv) then
        do j = 1,n_materials
          reduced_tal(i)%mean(:,j) = reduced_tal(i)%mean(:,j) / mat(j)%vol
        end do
      end if

    end do

    ! compute k_inf
    ana_kinf_mean = dble(reduced_n_fiss)*nubar/dble(nhistories)
    ana_kinf_std  = nubar*sqrt((dble(reduced_n_fiss)/dble(nhistories) -                &
   &                (dble(reduced_n_fiss)/dble(nhistories))**2)/dble(nhistories-1))
    col_kinf_mean = sum(reduced_tal(n_tallies)%mean)
    col_kinf_std  = sum(reduced_tal(n_tallies)%std)

! !Added by S. Xu (May. 2012)
!    col_kinf_mean_2 = sum(tal(6)%mean)
!    col_kinf_std_2  = sum(tal(6)%std)
!
!    reaction_tally(1) = sum(tal(2)%mean)
!    reaction_tally(2) = sum(tal(3)%mean)
!    reaction_tally(3) = sum(tal(4)%mean)

  end subroutine finalize_tallies

!===============================================================================
! CALCULATE_STATISTICS
!> @brief routine to compute mean and standard deviation of tallies
!===============================================================================

  subroutine calculate_statistics(this,n)

    ! formal variables
    type(tally_type) :: this ! a tally
    integer(8)          :: n    ! number of histories run

!#ifdef DEBUG
!  write(404, *) 'sum0 = ', this%sum0, ' nhistories = ', n, "from processor ", rank
!#endif

    ! compute mean
    this%mean =  this%sum0 / dble(n)

    ! compute standard deviation
    this%std = sqrt((this%sum0_sq/dble(n) - this%mean**2)/dble(n-1))

!#ifdef DEBUG
!  write(404, *) this%mean, this%std, "from processor ", rank
!#endif

  end subroutine calculate_statistics


!===============================================================================
! DEALLOCATE_TALLY
!> @brief routine to deallocate tally types
!===============================================================================

  subroutine deallocate_tally(this)

    ! formal variables
    type(tally_type) :: this ! a tally

    ! deallocate all
    if (allocated(this%E)) deallocate(this%E)
    if (allocated(this%val)) deallocate(this%val)
    if (allocated(this%sum0)) deallocate(this%sum0)
    if (allocated(this%sum0_sq)) deallocate(this%sum0_sq)

  end subroutine deallocate_tally

end module tally
