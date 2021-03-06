!==============================================================================!
! MODULE: materials
!
!> @author Bryan Herman
!>
!> @brief Contains information about the isotopics of problem
!==============================================================================!

module materials

!  use material_header, only: n_source_type, thermal_type, iso_type, material_type
  use global

  implicit none

contains

!===============================================================================
! SET_UP_MATERIALS
!> @brief routine that initializes the materials
!===============================================================================

  subroutine setup_material(this,emin,emax,nisotopes,vol)

    ! formal variables
    type(material_type) :: this      ! a material
    real(8)             :: emin      ! minimum energy to consider
    real(8)             :: emax      ! maximum energy to consider
    real(8)             :: vol       ! volume of material
    integer             :: nisotopes ! number of isotopes

    ! set number of isotopes
    this%nisotopes = nisotopes

    ! set volume
    this%vol = vol

    ! allocate isotopes array
    if (.not. allocated(this%isotopes)) allocate(this%isotopes(this%nisotopes))


    ! set up current isotope index
    this%curr_iso = 1

    ! set energy bounds
    this%E_min = emin
    this%E_max = emax

  end subroutine setup_material

!===============================================================================
! LOAD_ISOTOPE
!> @brief routine that loads isotope properties, xs, etc. into memory
!===============================================================================

  subroutine load_isotope(this,N,A,path,thermal,name, doppler, xsranderror)

    use hdf5

    ! formal variables
    type(material_type),target :: this    ! a material
    real(8)                    :: N       ! number density
    real(8)                    :: A       ! atomic weight
    character(len=255)         :: path    ! path to isotope
    character(len=255)         :: name    ! name of isotope
    logical                    :: thermal ! contains a thermal lib
    ! Added by S. Xu
    logical                    :: doppler !indicte whether doppler braodening is performed
    logical                    :: xsranderror ! indicte whether to use xs with random error

    ! local variables
    integer                        :: error        ! hdf5 error 
    integer(HID_T)                 :: hdf5_file    ! hdf5 file id
    integer(HID_T)                 :: dataset_id   ! hdf5 dataset id
    integer(HSIZE_T), dimension(1) :: dim1         ! dimension of hdf5 var
    integer(HSIZE_T), dimension(2) :: dim2         ! dimension of hdf5 var
    ! arrray size of capture, fission and scattering xs
    integer                        :: capt_size, fiss_size, scat_size
    type(thermal_type), pointer    :: therm

    integer                        :: i

    if (master) then
        ! display to user
        write(*,*) 'Loading isotope: ',trim(name)
    end if

    ! set parameters
    this%isotopes(this%curr_iso)%N = N
    this%isotopes(this%curr_iso)%A = A
    this%isotopes(this%curr_iso)%mubar = 2._8/(3._8*A)
    this%isotopes(this%curr_iso)%alpha = ((A-1._8)/(A+1._8))**2
    this%isotopes(this%curr_iso)%thermal = thermal
    this%isotopes(this%curr_iso)%name = name

    ! Added by S. Xu (Apr. 2012) to calculate the alpha factor for doppler broadening
    this%isotopes(this%curr_iso)%alpha_MB = A*M_NUCLEON/2._8/K_BOLTZMANN/T
    this%isotopes(this%curr_iso)%doppler = doppler
    this%isotopes(this%curr_iso)%xsranderror = xsranderror
    this%isotopes(this%curr_iso)%randerror = randerror

    ! open up hdf5 file
    call h5fopen_f(trim(path),H5F_ACC_RDWR_F,hdf5_file,error)

! Commented out by S. Xu (Apr. 2012)
!    ! read size of vector
!    call h5dopen_f(hdf5_file,"/vecsize",dataset_id,error)
!    dim1 = (/1/)
!    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,vecsize,dim1,error)
!    call h5dclose_f(dataset_id,error)


    dim1 = (/1/)
    call h5dopen_f(hdf5_file,"/capt_size",dataset_id,error)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,capt_size,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/fiss_size",dataset_id,error)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,fiss_size,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/scat_size",dataset_id,error)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,scat_size,dim1,error)
    call h5dclose_f(dataset_id,error)

! Added by S. Xu
    this%isotopes(this%curr_iso)%capt_size = capt_size
    this%isotopes(this%curr_iso)%scat_size = scat_size
    this%isotopes(this%curr_iso)%fiss_size = fiss_size

    ! added by S. Xu (Apr. 2012)
    ! allocate all energy grid vectors
    if (.not.allocated(this%isotopes(this%curr_iso)%engy_capt))  &
       & allocate(this%isotopes(this%curr_iso)%engy_capt(capt_size))
    if (.not.allocated(this%isotopes(this%curr_iso)%engy_scat))  &
       & allocate(this%isotopes(this%curr_iso)%engy_scat(scat_size))
    if (.not.allocated(this%isotopes(this%curr_iso)%engy_fiss))  &
       & allocate(this%isotopes(this%curr_iso)%engy_fiss(fiss_size))

    ! allocate all xs vectors
    if (.not.allocated(this%isotopes(this%curr_iso)%xs_scat))  &
       & allocate(this%isotopes(this%curr_iso)%xs_scat(scat_size))
    if (.not.allocated(this%isotopes(this%curr_iso)%xs_capt))  &
       & allocate(this%isotopes(this%curr_iso)%xs_capt(capt_size))
    if (.not.allocated(this%isotopes(this%curr_iso)%xs_fiss))  &
       & allocate(this%isotopes(this%curr_iso)%xs_fiss(fiss_size))

! Commented out by S. Xu (Apr. 2012)
!    ! keep the size
!    this%isotopes(this%curr_iso)%npts = vecsize

!print *, this%curr_iso
    ! added by S. Xu (Apr. 2012)
    ! zero out energy grid vectors
    this%isotopes(this%curr_iso)%engy_scat = 0.0_8
    this%isotopes(this%curr_iso)%engy_capt = 0.0_8
    this%isotopes(this%curr_iso)%engy_fiss = 0.0_8

    ! zero out xs vectors
    this%isotopes(this%curr_iso)%xs_scat = 0.0_8
    this%isotopes(this%curr_iso)%xs_capt = 0.0_8
    this%isotopes(this%curr_iso)%xs_fiss = 0.0_8

    ! added by S. Xu (Apr. 2012)
    ! read in energy grid
    call h5dopen_f(hdf5_file,"/engy_scat",dataset_id,error)
    dim1 = (/scat_size/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%engy_scat,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/engy_capt",dataset_id,error)
    dim1 = (/capt_size/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%engy_capt,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/engy_fiss",dataset_id,error)
    dim1 = (/fiss_size/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%engy_fiss,dim1,error)
    call h5dclose_f(dataset_id,error)


    ! read in xs
    call h5dopen_f(hdf5_file,"/xs_scat",dataset_id,error)
    dim1 = (/scat_size/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%xs_scat,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/xs_capt",dataset_id,error)
    dim1 = (/capt_size/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%xs_capt,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/xs_fiss",dataset_id,error)
    dim1 = (/fiss_size/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%xs_fiss,dim1,error)
    call h5dclose_f(dataset_id,error)
    
!    if (trim(this%isotopes(this%curr_iso)%name)=='U-238') then
!      open(12, file='test_capt_U238', status='unknown')
!      do i = 1,fiss_size
!        write(12, '(e19.8e3, 3x, e19.8e3)') this%isotopes(this%curr_iso)%engy_fiss(i), &
!        & this%isotopes(this%curr_iso)%xs_fiss(i)
!      end do
!      close(12)
!    end if

! Commented out by S. Xu (Apr. 2012)
!    ! get energy interval width
!    call h5dopen_f(hdf5_file,"/E_width",dataset_id,error)
!    dim1 = (/1/)
!    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%isotopes(this%curr_iso)%E_width,dim1,error)
!    call h5dclose_f(dataset_id,error)

    ! check for thermal scattering kernel and load that
    if (this%isotopes(this%curr_iso)%thermal) then

      ! set pointer
      therm => this%isotopes(this%curr_iso)%thermal_lib

      ! load sizes
      call h5dopen_f(hdf5_file,"/kTsize",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,therm%kTsize,dim1,error)
      call h5dclose_f(dataset_id,error)
      call h5dopen_f(hdf5_file,"/cdfsize",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,therm%cdfsize,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! read in cdf width
      call h5dopen_f(hdf5_file,"/cdf_width",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,therm%cdf_width,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! preallocate vectors
      if(.not.allocated(therm%kTvec)) allocate(therm%kTvec(therm%kTsize))
      if(.not.allocated(therm%Erat)) allocate(therm%Erat(therm%cdfsize,therm%kTsize))

      ! read in vectors
      call h5dopen_f(hdf5_file,"/kT",dataset_id,error)
      dim1 = (/therm%kTsize/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,therm%kTvec,dim1,error)
      call h5dclose_f(dataset_id,error)
      call h5dopen_f(hdf5_file,"/Erat",dataset_id,error)
      dim2 = (/therm%cdfsize,therm%kTsize/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,therm%Erat,dim2,error)
      call h5dclose_f(dataset_id,error)

    end if

    ! close hdf5 file
    call h5fclose_f(hdf5_file,error)

! Commented out by S. Xu (Apr. 2012)
    ! increment isotope counter
    this%curr_iso = this%curr_iso + 1

!#ifdef DEBUG
!  print *, 'from processor ', rank
!  print *, "isotope ", this%curr_iso, "from processor ", rank
!  print *, "capture xs data : "
!  print *, this%isotopes(this%curr_iso-1)%xs_capt(1:20)
!  print *, this%isotopes(this%curr_iso-1)%engy_capt(1:20)
!#endif

  end subroutine load_isotope

!===============================================================================
! LOAD_SOURCE
!> @brief routine to load fission source into memory
!===============================================================================

  subroutine load_source(this,source_type,source_path)

    use hdf5 

    ! formal variables
    type(material_type) :: this        ! a material
    integer             :: source_type ! 0 - fixed, 1 - fission
    character(len=255)  :: source_path ! path to source file

    ! local variables
    integer                        :: error        ! hdf5 error 
    integer(HID_T)                 :: hdf5_file    ! hdf5 file id
    integer(HID_T)                 :: dataset_id   ! hdf5 dataset id
    integer(HSIZE_T), dimension(1) :: dim1         ! dimension of hdf5 var
    integer                        :: vecsize      ! vector size for fission


    ! check for fission source
    if (source_type == 1) then

      ! open the fission source file
      call h5fopen_f(trim(source_path),H5F_ACC_RDWR_F,hdf5_file,error)

      ! open dataset and read in vector size
      call h5dopen_f(hdf5_file,"/vecsize",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,vecsize,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! open dataset and read in width of cdf interval
      call h5dopen_f(hdf5_file,"/cdf_width",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%source%cdf_width,dim1,  &
     &               error)
      call h5dclose_f(dataset_id,error)

      ! preallocate vectors in source object
      if(.not.allocated(this%source%E)) allocate(this%source%E(vecsize))

      ! open dataset and read in energy vector
      call h5dopen_f(hdf5_file,"/E",dataset_id,error)
      dim1 = (/vecsize/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%source%E,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! close the file
      call h5fclose_f(hdf5_file,error)

    end if

  end subroutine load_source


!===============================================================================
! doppler_broaden
!> @brief routine to do the on the fly doppler broadening
!===============================================================================
  subroutine doppler_broaden()

    use on_the_fly_xs_gen, only: energy_doppler_broadened
    use LinearInterpolation, only: LinInterp
    use random_lcg, only: prn
!    use constants, only: M_NEUT

    ! local variables
    real(8) :: E_brdn
    integer :: i, j, k  ! loop counter
    real(8) :: v
    real(8) :: v_brdn
    real(8) :: xs_capt_tmp,xs_scat_tmp,xs_fiss_tmp
    real(8) :: rn  !random number
!    real    :: dE ! for debug

    do k = 1, n_materials

        do i = 1, mat(k)%nisotopes

          if (mat(k)%isotopes(i)%doppler) then

            v = sqrt(2._8*neut%E/M_NEUT)
            mat(k)%isotopes(i)%xs_capt_brdn = 0.0_8
            mat(k)%isotopes(i)%xs_scat_brdn = 0.0_8
            mat(k)%isotopes(i)%xs_fiss_brdn = 0.0_8

            do j=1,sample_per_xs
            ! sample the relative kinetic energy
              call energy_doppler_broadened(v, mat(k)%isotopes(i)%alpha_MB, v_brdn)
              E_brdn = 0.5_8*M_NEUT*v_brdn**2

              ! broaden xs
              call LinInterp(mat(k)%isotopes(i)%engy_capt, mat(k)%isotopes(i)%xs_capt, E_brdn, xs_capt_tmp)
              mat(k)%isotopes(i)%xs_capt_brdn = mat(k)%isotopes(i)%xs_capt_brdn+v_brdn/v*xs_capt_tmp

              call LinInterp(mat(k)%isotopes(i)%engy_scat, mat(k)%isotopes(i)%xs_scat, E_brdn, xs_scat_tmp)
              mat(k)%isotopes(i)%xs_scat_brdn = mat(k)%isotopes(i)%xs_scat_brdn+v_brdn/v*xs_scat_tmp

              call LinInterp(mat(k)%isotopes(i)%engy_fiss, mat(k)%isotopes(i)%xs_fiss, E_brdn, xs_fiss_tmp)
              mat(k)%isotopes(i)%xs_fiss_brdn = mat(k)%isotopes(i)%xs_fiss_brdn+v_brdn/v*xs_fiss_tmp
            end do

            mat(k)%isotopes(i)%xs_capt_brdn = mat(k)%isotopes(i)%xs_capt_brdn/dble(sample_per_xs)
            mat(k)%isotopes(i)%xs_scat_brdn = mat(k)%isotopes(i)%xs_scat_brdn/dble(sample_per_xs)
            mat(k)%isotopes(i)%xs_fiss_brdn = mat(k)%isotopes(i)%xs_fiss_brdn/dble(sample_per_xs)

    !        write(997, '(es19.8e3, 3x, es19.8e3)')  neut%E, mat(k)%isotopes(i)%xs_capt_brdn
    !        write(998, '(es19.8e3, 3x, es19.8e3)')  neut%E, mat(k)%isotopes(i)%xs_scat_brdn
    !        write(999, '(es19.8e3, 3x, es19.8e3)')  neut%E, mat(k)%isotopes(i)%xs_fiss_brdn

          else if (mat(k)%isotopes(i)%xsranderror) then

            ! find xs
            call LinInterp(mat(k)%isotopes(i)%engy_capt, mat(k)%isotopes(i)%xs_capt, neut%E, mat(k)%isotopes(i)%xs_capt_brdn)
            rn = prn()
            mat(k)%isotopes(i)%xs_capt_brdn = mat(k)%isotopes(i)%xs_capt_brdn*(1 + randerror*(2*rn - 1))
            
            call LinInterp(mat(k)%isotopes(i)%engy_scat, mat(k)%isotopes(i)%xs_scat, neut%E, mat(k)%isotopes(i)%xs_scat_brdn)
            rn = prn()
            mat(k)%isotopes(i)%xs_scat_brdn = mat(k)%isotopes(i)%xs_scat_brdn*(1 + randerror*(2*rn - 1))
      
            call LinInterp(mat(k)%isotopes(i)%engy_fiss, mat(k)%isotopes(i)%xs_fiss, neut%E, mat(k)%isotopes(i)%xs_fiss_brdn)
            rn = prn()
            mat(k)%isotopes(i)%xs_fiss_brdn = mat(k)%isotopes(i)%xs_fiss_brdn*(1 + randerror*(2*rn - 1))

          else

            ! find xs
            call LinInterp(mat(k)%isotopes(i)%engy_capt, mat(k)%isotopes(i)%xs_capt, neut%E, mat(k)%isotopes(i)%xs_capt_brdn)

            call LinInterp(mat(k)%isotopes(i)%engy_scat, mat(k)%isotopes(i)%xs_scat, neut%E, mat(k)%isotopes(i)%xs_scat_brdn)
      
            call LinInterp(mat(k)%isotopes(i)%engy_fiss, mat(k)%isotopes(i)%xs_fiss, neut%E, mat(k)%isotopes(i)%xs_fiss_brdn)

          end if

    ! ! Added by S. Xu for Debug
    ! if (trim(mat(k)%isotopes(i)%name) == 'U-238' .and. E<=10.e-6 .and. E>=1.e-6) then
    !       dE = sqrt(2._8*M_NEUT*E/mat(k)%isotopes(i)%alpha_MB)
    !       write(10, '(2e19.8e3, 2x, f7.3, 3e19.8e3)') E, E_brdn, (E_brdn-E)/dE, &
    !    & mat(k)%isotopes(i)%xs_capt_brdn, mat(k)%isotopes(i)%xs_scat_brdn, mat(k)%isotopes(i)%xs_fiss_brdn
    ! end if

        end do

    end do

  end subroutine doppler_broaden

!===============================================================================
! COMPUTE_MACRO_XS
!> @brief routine to pre-compute macroscopic cross sections 
!===============================================================================

  subroutine compute_macro_xs()

    ! local variables
    integer                 :: i, j  ! loop counter

    do j = 1, n_materials

        ! allocate xs arrays
        if (.not.allocated(mat(j)%xs_Total_brdn)) allocate(mat(j)%xs_Total_brdn(mat(j)%nisotopes))
        if (.not.allocated(mat(j)%xs_Capt_brdn)) allocate(mat(j)%xs_Capt_brdn(mat(j)%nisotopes))
        if (.not.allocated(mat(j)%xs_Scat_brdn)) allocate(mat(j)%xs_Scat_brdn(mat(j)%nisotopes))
        if (.not.allocated(mat(j)%xs_Fiss_brdn)) allocate(mat(j)%xs_Fiss_brdn(mat(j)%nisotopes))
        if (.not.allocated(mat(j)%xs_Absb_brdn)) allocate(mat(j)%xs_Absb_brdn(mat(j)%nisotopes))
        if (.not.allocated(mat(j)%xs_Tran_brdn)) allocate(mat(j)%xs_Tran_brdn(mat(j)%nisotopes))

        ! zero out total xs
        mat(j)%xs_Total_brdn = 0.0_8

        ! begin loop over isotopes
        do i = 1,mat(j)%nisotopes

          ! multiply microscopic cross section by number density and append
          mat(j)%xs_Capt_brdn(i) = mat(j)%isotopes(i)%N*(mat(j)%isotopes(i)%xs_capt_brdn)
          mat(j)%xs_Fiss_brdn(i) = mat(j)%isotopes(i)%N*(mat(j)%isotopes(i)%xs_fiss_brdn)
          mat(j)%xs_Scat_brdn(i) = mat(j)%isotopes(i)%N*(mat(j)%isotopes(i)%xs_scat_brdn)
          mat(j)%xs_Absb_brdn(i) = mat(j)%xs_Capt_brdn(i)+mat(j)%xs_Fiss_brdn(i)
          mat(j)%xs_Total_brdn(i) = mat(j)%xs_Absb_brdn(i)+mat(j)%xs_Scat_brdn(i)
          mat(j)%xs_Tran_brdn(i) = mat(j)%xs_Total_brdn(i) - mat(j)%isotopes(i)%mubar*mat(j)%xs_Scat_brdn(i)

        end do

    end do

  end subroutine compute_macro_xs


!===============================================================================
! DEALLOCATE_MATERIAL
!> @brief routine to deallocate a material
!===============================================================================

  subroutine deallocate_material(this)

    ! formal variables
    type(material_type) :: this ! a material

    ! local variables
    integer :: i ! loop counter

    ! deallocate source information
    if (allocated(this%source%E)) deallocate(this%source%E)

    ! begin loop over isotopes for deallocation
    do i = 1,this%nisotopes

      ! deallocate thermal library
      if (allocated(this%isotopes(i)%thermal_lib%kTvec)) deallocate            &
     &             (this%isotopes(i)%thermal_lib%kTvec)
      if (allocated(this%isotopes(i)%thermal_lib%Erat)) deallocate             &
     &             (this%isotopes(i)%thermal_lib%Erat)

      ! deallocate xs
      if (allocated(this%isotopes(i)%xs_scat)) deallocate                      &
     &             (this%isotopes(i)%xs_scat)
      if (allocated(this%isotopes(i)%xs_capt)) deallocate                      &
     &             (this%isotopes(i)%xs_capt)
      if (allocated(this%isotopes(i)%xs_fiss)) deallocate                      &
     &             (this%isotopes(i)%xs_fiss)

     ! Added by S. Xu (Apr. 2012)
     ! deallocate xs
      if (allocated(this%isotopes(i)%engy_scat)) deallocate                      &
     &             (this%isotopes(i)%engy_scat)
      if (allocated(this%isotopes(i)%engy_capt)) deallocate                      &
     &             (this%isotopes(i)%engy_capt)
      if (allocated(this%isotopes(i)%engy_fiss)) deallocate                      &
     &             (this%isotopes(i)%engy_fiss)

    end do

    ! deallocate isotopes
    if (allocated(this%isotopes)) deallocate(this%isotopes) 

    ! deallocate macro xs
    if (allocated(this%totalxs)) deallocate(this%totalxs)
    if (allocated(this%scattxs)) deallocate(this%scattxs)
    if (allocated(this%absorxs)) deallocate(this%absorxs)
    if (allocated(this%captuxs)) deallocate(this%captuxs)
    if (allocated(this%fissixs)) deallocate(this%fissixs)
    if (allocated(this%transxs)) deallocate(this%transxs)

! Added by S. Xu (Apr. 2012)
    if (allocated(this%xs_Total_brdn)) deallocate(this%xs_Total_brdn)
    if (allocated(this%xs_Capt_brdn)) deallocate(this%xs_Capt_brdn)
    if (allocated(this%xs_Scat_brdn)) deallocate(this%xs_Scat_brdn)
    if (allocated(this%xs_Fiss_brdn)) deallocate(this%xs_Fiss_brdn)
    if (allocated(this%xs_Absb_brdn)) deallocate(this%xs_Absb_brdn)
    if (allocated(this%xs_Tran_brdn)) deallocate(this%xs_Tran_brdn)

  end subroutine deallocate_material

end module materials
