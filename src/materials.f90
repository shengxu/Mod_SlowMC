!==============================================================================!
! MODULE: materials
!
!> @author Bryan Herman
!>
!> @brief Contains information about the isotopics of problem
!==============================================================================!

module materials

  implicit none
  private
  public :: setup_material,load_source,load_isotope,compute_macroxs,           &
 &          deallocate_material, doppler_broaden_xs  ! last added by S. Xu

  type :: source_type

    real(8), allocatable :: E(:)      ! energy range for fission source
    real(8)              :: cdf_width ! width of cdf bins from 0 to 1

  end type source_type

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

    type(source_type)           :: source        ! the source of neutrons
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

  subroutine load_isotope(this,N,A,path,thermal,name, doppler)

    use hdf5
! Added by S. Xu (Apr. 2012)
    use constants, only: PI, K_BOLTZMANN, M_NUCLEON
    use parameters, only: T


    ! formal variables
    type(material_type),target :: this    ! a material
    real(8)                    :: N       ! number density
    real(8)                    :: A       ! atomic weight
    character(len=255)         :: path    ! path to isotope
    character(len=255)         :: name    ! name of isotope
    logical                    :: thermal ! contains a thermal lib
    ! Added by S. Xu
    logical                    :: doppler !indicte whether doppler braodening is performed

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

    ! display to user
    write(*,*) 'Loading isotope: ',trim(name)

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

print *, this%curr_iso
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
! doppler_broaden_xs
!> @brief routine routine to do the on the fly doppler broadening
!===============================================================================
  subroutine doppler_broaden_xs(this, E)

    use on_the_fly_xs_gen, only: energy_doppler_broadened
    use LinearInterpolation, only: LinInterp
    use constants, only: M_NEUT
    use parameters, only: vmin
    
    ! formal variables
    type(material_type),target :: this ! a material
    real(8)                    :: E    ! neutron energy
    ! local variables
    real(8) :: E_brdn
    integer :: ind ! energy interval
    logical :: stat 
    integer :: i, j  ! loop counter
    integer :: n_sample=2
    real(8) :: v
    real(8) :: v_brdn
    real(8) :: xs_capt_tmp,xs_scat_tmp,xs_fiss_tmp
!    real    :: dE ! for debug

    do i = 1, this%nisotopes
      if (this%isotopes(i)%doppler) then

        v = sqrt(2._8*E/M_NEUT)
        this%isotopes(i)%xs_capt_brdn = 0.0_8
        this%isotopes(i)%xs_scat_brdn = 0.0_8
        this%isotopes(i)%xs_fiss_brdn = 0.0_8

        do j=1,n_sample
        ! sample the relative kinetic energy
          call energy_doppler_broadened(v, this%isotopes(i)%alpha_MB, v_brdn)
          E_brdn = 0.5_8*M_NEUT*v_brdn**2

          ! broaden xs
          call LinInterp(this%isotopes(i)%engy_capt, this%isotopes(i)%xs_capt, &
               & E_brdn, xs_capt_tmp)
          this%isotopes(i)%xs_capt_brdn = this%isotopes(i)%xs_capt_brdn+v_brdn/v*xs_capt_tmp

          call LinInterp(this%isotopes(i)%engy_scat, this%isotopes(i)%xs_scat, &
               & E_brdn, xs_scat_tmp)
          this%isotopes(i)%xs_scat_brdn = this%isotopes(i)%xs_scat_brdn+v_brdn/v*xs_scat_tmp

          call LinInterp(this%isotopes(i)%engy_fiss, this%isotopes(i)%xs_fiss, &
               & E_brdn, xs_fiss_tmp)
          this%isotopes(i)%xs_fiss_brdn = this%isotopes(i)%xs_fiss_brdn+v_brdn/v*xs_fiss_tmp
        end do

        this%isotopes(i)%xs_capt_brdn = this%isotopes(i)%xs_capt_brdn/dble(n_sample)
        this%isotopes(i)%xs_scat_brdn = this%isotopes(i)%xs_scat_brdn/dble(n_sample)
        this%isotopes(i)%xs_fiss_brdn = this%isotopes(i)%xs_fiss_brdn/dble(n_sample)

!        write(997, '(es19.8e3, 3x, es19.8e3)')  E, this%isotopes(i)%xs_capt_brdn
!        write(998, '(es19.8e3, 3x, es19.8e3)')  E, this%isotopes(i)%xs_scat_brdn
!        write(999, '(es19.8e3, 3x, es19.8e3)')  E, this%isotopes(i)%xs_fiss_brdn

      else

        ! find xs
        call LinInterp(this%isotopes(i)%engy_capt, this%isotopes(i)%xs_capt, &
             & E, this%isotopes(i)%xs_capt_brdn)

        call LinInterp(this%isotopes(i)%engy_scat, this%isotopes(i)%xs_scat, &
             & E, this%isotopes(i)%xs_scat_brdn)
  
        call LinInterp(this%isotopes(i)%engy_fiss, this%isotopes(i)%xs_fiss, &
             & E, this%isotopes(i)%xs_fiss_brdn)

      end if

! ! Added by S. Xu for Debug
! if (trim(this%isotopes(i)%name) == 'U-238' .and. E<=10.e-6 .and. E>=1.e-6) then
!       dE = sqrt(2._8*M_NEUT*E/this%isotopes(i)%alpha_MB)
!       write(10, '(2e19.8e3, 2x, f7.3, 3e19.8e3)') E, E_brdn, (E_brdn-E)/dE, &
!    & this%isotopes(i)%xs_capt_brdn, this%isotopes(i)%xs_scat_brdn, this%isotopes(i)%xs_fiss_brdn
! end if

    end do


  end subroutine doppler_broaden_xs

!===============================================================================
! COMPUTE_MACROXS
!> @brief routine to pre-compute macroscopic cross sections 
!===============================================================================

  subroutine compute_macroxs(this)

    ! formal variables
    type(material_type),target :: this ! a material

    ! local variables
    integer                 :: i   ! loop counter
    type(iso_type), pointer :: iso ! pointer to current isotope

! Commented out by S. Xu (Apr. 2012)
!    ! allocate xs arrays
!    if (.not.allocated(this%totalxs))                                          &
!   &                            allocate(this%totalxs(this%npts,this%nisotopes))
!    if (.not.allocated(this%scattxs))                                          &
!   &                            allocate(this%scattxs(this%npts,this%nisotopes))
!    if (.not.allocated(this%absorxs))                                          &
!   &                            allocate(this%absorxs(this%npts,this%nisotopes))
!    if (.not.allocated(this%captuxs))                                          &
!   &                            allocate(this%captuxs(this%npts,this%nisotopes))
!    if (.not.allocated(this%fissixs))                                          &
!   &                            allocate(this%fissixs(this%npts,this%nisotopes))
!    if (.not.allocated(this%transxs))                                          &
!   &                            allocate(this%transxs(this%npts,this%nisotopes))
!
!    ! zero out total xs
!    this%totalxs = 0.0_8
!
!    ! begin loop over isotopes
!    do i = 1,this%nisotopes
!
!      ! set pointer to isotope
!      iso => this%isotopes(i)
!
!      ! multiply microscopic cross section by number density and append
!      this%captuxs(:,i) = iso%N*(iso%xs_capt)
!      this%fissixs(:,i) = iso%N*(iso%xs_fiss)
!      this%scattxs(:,i) = iso%N*(iso%xs_scat)
!      this%absorxs(:,i) = iso%N*(iso%xs_capt + iso%xs_fiss)
!      this%totalxs(:,i) = iso%N*(iso%xs_capt + iso%xs_fiss + iso%xs_scat)
!      this%transxs(:,i) = this%totalxs(:,i) - iso%mubar*this%scattxs(:,i)
!
!    end do

! Added by S. Xu (Apr. 2012)
    ! allocate xs arrays
    if (.not.allocated(this%xs_Total_brdn)) allocate(this%xs_Total_brdn(this%nisotopes))
    if (.not.allocated(this%xs_Capt_brdn)) allocate(this%xs_Capt_brdn(this%nisotopes))
    if (.not.allocated(this%xs_Scat_brdn)) allocate(this%xs_Scat_brdn(this%nisotopes))
    if (.not.allocated(this%xs_Fiss_brdn)) allocate(this%xs_Fiss_brdn(this%nisotopes))
    if (.not.allocated(this%xs_Absb_brdn)) allocate(this%xs_Absb_brdn(this%nisotopes))
    if (.not.allocated(this%xs_Tran_brdn)) allocate(this%xs_Tran_brdn(this%nisotopes))

    ! zero out total xs
    this%xs_Total_brdn = 0.0_8

    ! begin loop over isotopes
    do i = 1,this%nisotopes

      ! set pointer to isotope
      iso => this%isotopes(i)

      ! multiply microscopic cross section by number density and append
      this%xs_Capt_brdn(i) = iso%N*(iso%xs_capt_brdn)
      this%xs_Fiss_brdn(i) = iso%N*(iso%xs_fiss_brdn)
      this%xs_Scat_brdn(i) = iso%N*(iso%xs_scat_brdn)
      this%xs_Absb_brdn(i) = this%xs_Capt_brdn(i)+this%xs_Fiss_brdn(i)
      this%xs_Total_brdn(i) = this%xs_Absb_brdn(i)+this%xs_Scat_brdn(i)
      this%xs_Tran_brdn(i) = this%xs_Total_brdn(i) - iso%mubar*this%xs_Scat_brdn(i)

    end do

  end subroutine compute_macroxs


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
