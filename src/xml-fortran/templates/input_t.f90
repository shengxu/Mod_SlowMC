module xml_data_input_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   save
   integer, private :: lurep_
   logical, private :: strict_

type settings_xml
   integer                                         :: histories
   integer                                         :: seed
   integer                                         :: source_type
   character(len=255)                                :: source_path
   real(kind=kind(1.0d0))                          :: Dancoff
   character(len=255)                                :: res_iso
   real(kind=kind(1.0d0))                          :: radius
   real(kind=kind(1.0d0))                          :: temperature
end type settings_xml

type nuclide_xml
   real(kind=kind(1.0d0))                          :: N
   real(kind=kind(1.0d0))                          :: A
   real(kind=kind(1.0d0))                          :: V
   character(len=255)                                :: path
   logical                                         :: thermal
   character(len=255)                                :: name
   logical                                         :: doppler
end type nuclide_xml

type material_xml
   character(len=255)                                :: type
   real(kind=kind(1.0d0))                          :: V
   type(nuclide_xml), dimension(:), pointer        :: nuclides => null()
end type material_xml

type materials_xml
   type(material_xml), dimension(:), pointer       :: material => null()
end type materials_xml

type tally_xml
   logical                                         :: dv
   real(kind=kind(1.0d0)), dimension(:), pointer   :: Ebins => null()
   character(len=255)                                :: type
   integer                                         :: isotope
   integer                                         :: region
end type tally_xml

type tallies_xml
   type(tally_xml), dimension(:), pointer          :: tally => null()
end type tallies_xml
   type(settings_xml)                              :: settings_
   type(materials_xml)                             :: materials_
   type(tallies_xml)                               :: tallies_
contains
subroutine read_xml_type_settings_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(settings_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(settings_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_settings_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_settings_xml_array

subroutine read_xml_type_settings_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(settings_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_histories
   logical                                         :: has_seed
   logical                                         :: has_source_type
   logical                                         :: has_source_path
   logical                                         :: has_Dancoff
   logical                                         :: has_res_iso
   logical                                         :: has_radius
   logical                                         :: has_temperature
   has_histories                        = .false.
   has_seed                             = .false.
   has_source_type                      = .false.
   has_source_path                      = .false.
   has_Dancoff                          = .false.
   has_res_iso                          = .false.
   has_radius                           = .false.
   has_temperature                      = .false.
   call init_xml_type_settings_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('histories')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%histories, has_histories )
      case('seed')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%seed, has_seed )
      case('source_type')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%source_type, has_source_type )
      case('source_path')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%source_path, has_source_path )
      case('Dancoff')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%Dancoff, has_Dancoff )
      case('res_iso')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%res_iso, has_res_iso )
      case('radius')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%radius, has_radius )
      case('temperature')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%temperature, has_temperature )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_histories ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on histories')
   endif
   if ( .not. has_seed ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on seed')
   endif
   if ( .not. has_source_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on source_type')
   endif
   if ( .not. has_source_path ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on source_path')
   endif
   if ( .not. has_Dancoff ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on Dancoff')
   endif
   if ( .not. has_res_iso ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on res_iso')
   endif
   if ( .not. has_radius ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on radius')
   endif
   if ( .not. has_temperature ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on temperature')
   endif
end subroutine read_xml_type_settings_xml
subroutine init_xml_type_settings_xml_array( dvar )
   type(settings_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_settings_xml_array
subroutine init_xml_type_settings_xml(dvar)
   type(settings_xml) :: dvar
end subroutine init_xml_type_settings_xml
subroutine write_xml_type_settings_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(settings_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_settings_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_settings_xml_array

subroutine write_xml_type_settings_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(settings_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'histories', indent+3, dvar%histories)
   call write_to_xml_integer( info, 'seed', indent+3, dvar%seed)
   call write_to_xml_integer( info, 'source_type', indent+3, dvar%source_type)
   call write_to_xml_word( info, 'source_path', indent+3, dvar%source_path)
   call write_to_xml_double( info, 'Dancoff', indent+3, dvar%Dancoff)
   call write_to_xml_word( info, 'res_iso', indent+3, dvar%res_iso)
   call write_to_xml_double( info, 'radius', indent+3, dvar%radius)
   call write_to_xml_double( info, 'temperature', indent+3, dvar%temperature)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_settings_xml

subroutine read_xml_type_nuclide_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(nuclide_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(nuclide_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_nuclide_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_nuclide_xml_array

subroutine read_xml_type_nuclide_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(nuclide_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_N
   logical                                         :: has_A
   logical                                         :: has_V
   logical                                         :: has_path
   logical                                         :: has_thermal
   logical                                         :: has_name
   logical                                         :: has_doppler
   has_N                                = .false.
   has_A                                = .false.
   has_V                                = .false.
   has_path                             = .false.
   has_thermal                          = .false.
   has_name                             = .false.
   has_doppler                          = .false.
   call init_xml_type_nuclide_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('N')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%N, has_N )
      case('A')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%A, has_A )
      case('V')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%V, has_V )
      case('path')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%path, has_path )
      case('thermal')
         call read_xml_logical( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%thermal, has_thermal )
      case('name')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%name, has_name )
      case('doppler')
         call read_xml_logical( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%doppler, has_doppler )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_N ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on N')
   endif
   if ( .not. has_A ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on A')
   endif
   if ( .not. has_V ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on V')
   endif
   if ( .not. has_path ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on path')
   endif
   if ( .not. has_thermal ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on thermal')
   endif
   if ( .not. has_name ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on name')
   endif
   if ( .not. has_doppler ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on doppler')
   endif
end subroutine read_xml_type_nuclide_xml
subroutine init_xml_type_nuclide_xml_array( dvar )
   type(nuclide_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_nuclide_xml_array
subroutine init_xml_type_nuclide_xml(dvar)
   type(nuclide_xml) :: dvar
end subroutine init_xml_type_nuclide_xml
subroutine write_xml_type_nuclide_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(nuclide_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_nuclide_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_nuclide_xml_array

subroutine write_xml_type_nuclide_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(nuclide_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_double( info, 'N', indent+3, dvar%N)
   call write_to_xml_double( info, 'A', indent+3, dvar%A)
   call write_to_xml_double( info, 'V', indent+3, dvar%V)
   call write_to_xml_word( info, 'path', indent+3, dvar%path)
   call write_to_xml_logical( info, 'thermal', indent+3, dvar%thermal)
   call write_to_xml_word( info, 'name', indent+3, dvar%name)
   call write_to_xml_logical( info, 'doppler', indent+3, dvar%doppler)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_nuclide_xml

subroutine read_xml_type_material_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(material_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(material_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_material_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_material_xml_array

subroutine read_xml_type_material_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(material_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_type
   logical                                         :: has_V
   logical                                         :: has_nuclides
   has_type                             = .false.
   has_V                                = .false.
   has_nuclides                         = .false.
   allocate(dvar%nuclides(0))
   call init_xml_type_material_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('V')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%V, has_V )
      case('nuclide')
         call read_xml_type_nuclide_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nuclides, has_nuclides )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type')
   endif
   if ( .not. has_V ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on V')
   endif
   if ( .not. has_nuclides ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on nuclides')
   endif
end subroutine read_xml_type_material_xml
subroutine init_xml_type_material_xml_array( dvar )
   type(material_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_material_xml_array
subroutine init_xml_type_material_xml(dvar)
   type(material_xml) :: dvar
end subroutine init_xml_type_material_xml
subroutine write_xml_type_material_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(material_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_material_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_material_xml_array

subroutine write_xml_type_material_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(material_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_double( info, 'V', indent+3, dvar%V)
   call write_xml_type_nuclide_xml_array( info, 'nuclide', indent+3, dvar%nuclides)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_material_xml

subroutine read_xml_type_materials_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(materials_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(materials_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_materials_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_materials_xml_array

subroutine read_xml_type_materials_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(materials_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_material
   has_material                         = .false.
   allocate(dvar%material(0))
   call init_xml_type_materials_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('material')
         call read_xml_type_material_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%material, has_material )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_material ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on material')
   endif
end subroutine read_xml_type_materials_xml
subroutine init_xml_type_materials_xml_array( dvar )
   type(materials_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_materials_xml_array
subroutine init_xml_type_materials_xml(dvar)
   type(materials_xml) :: dvar
end subroutine init_xml_type_materials_xml
subroutine write_xml_type_materials_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(materials_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_materials_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_materials_xml_array

subroutine write_xml_type_materials_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(materials_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_xml_type_material_xml_array( info, 'material', indent+3, dvar%material)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_materials_xml

subroutine read_xml_type_tally_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(tally_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_tally_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_tally_xml_array

subroutine read_xml_type_tally_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_dv
   logical                                         :: has_Ebins
   logical                                         :: has_type
   logical                                         :: has_isotope
   logical                                         :: has_region
   has_dv                               = .false.
   has_Ebins                            = .false.
   has_type                             = .false.
   has_isotope                          = .false.
   has_region                           = .false.
   call init_xml_type_tally_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('dv')
         call read_xml_logical( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dv, has_dv )
      case('Ebins')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%Ebins, has_Ebins )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('isotope')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%isotope, has_isotope )
      case('region')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%region, has_region )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_Ebins ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on Ebins')
   endif
   if ( .not. has_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type')
   endif
   if ( .not. has_isotope ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on isotope')
   endif
   if ( .not. has_region ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on region')
   endif
end subroutine read_xml_type_tally_xml
subroutine init_xml_type_tally_xml_array( dvar )
   type(tally_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_tally_xml_array
subroutine init_xml_type_tally_xml(dvar)
   type(tally_xml) :: dvar
   dvar%dv = .false.
end subroutine init_xml_type_tally_xml
subroutine write_xml_type_tally_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_tally_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_tally_xml_array

subroutine write_xml_type_tally_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_logical( info, 'dv', indent+3, dvar%dv)
   call write_to_xml_double_array( info, 'Ebins', indent+3, dvar%Ebins)
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_integer( info, 'isotope', indent+3, dvar%isotope)
   call write_to_xml_integer( info, 'region', indent+3, dvar%region)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_tally_xml

subroutine read_xml_type_tallies_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tallies_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(tallies_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_tallies_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_tallies_xml_array

subroutine read_xml_type_tallies_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tallies_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_tally
   has_tally                            = .false.
   allocate(dvar%tally(0))
   call init_xml_type_tallies_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('tally')
         call read_xml_type_tally_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%tally, has_tally )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_tally ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on tally')
   endif
end subroutine read_xml_type_tallies_xml
subroutine init_xml_type_tallies_xml_array( dvar )
   type(tallies_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_tallies_xml_array
subroutine init_xml_type_tallies_xml(dvar)
   type(tallies_xml) :: dvar
end subroutine init_xml_type_tallies_xml
subroutine write_xml_type_tallies_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tallies_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_tallies_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_tallies_xml_array

subroutine write_xml_type_tallies_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tallies_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_xml_type_tally_xml_array( info, 'tally', indent+3, dvar%tally)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_tallies_xml

subroutine read_xml_file_input_t(fname, lurep, errout)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep
   logical, intent(out), optional         :: errout

   type(XML_PARSE)                        :: info
   logical                                :: error
   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=80), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=100000), dimension(1:1000)  :: data
   integer                                :: nodata
   logical                                         :: has_settings_
   logical                                         :: has_materials_
   logical                                         :: has_tallies_
   has_settings_                        = .false.
   has_materials_                       = .false.
   has_tallies_                         = .false.

   call init_xml_file_input_t
   call xml_open( info, fname, .true. )
   call xml_options( info, report_errors=.false., ignore_whitespace=.true.)
   lurep_ = 0
   if ( present(lurep) ) then
      lurep_ = lurep
      call xml_options( info, report_lun=lurep )
   endif
   do
      call xml_get( info, starttag, endtag, attribs, noattribs, &
         data, nodata)
      if ( starttag /= '!--' ) exit
   enddo
   if ( starttag /= "input" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "input"')
      error = .true.
      call xml_close(info)
      return
   endif
   strict_ = .false.
   error = .false.
   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      if ( xml_error(info) ) then
         write(lurep_,*) 'Error reading input file!'
         error = .true.
         return
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('settings')
         call read_xml_type_settings_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            settings_, has_settings_ )
      case('materials')
         call read_xml_type_materials_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            materials_, has_materials_ )
      case('tallies')
         call read_xml_type_tallies_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            tallies_, has_tallies_ )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_settings_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on settings_')
   endif
   if ( .not. has_materials_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on materials_')
   endif
   if ( .not. has_tallies_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on tallies_')
   endif
   if ( present(errout) ) errout = error
   call xml_close(info)
end subroutine

subroutine write_xml_file_input_t(fname, lurep)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep

   type(XML_PARSE)                        :: info
   integer                                :: indent = 0

   call xml_open( info, fname, .false. )
   call xml_options( info, report_errors=.true.)
   if ( present(lurep) ) then
       call xml_options( info, report_errors=.true.)
   endif
   write(info%lun,'(a)') &
      '<input>'
   call write_xml_type_settings_xml( info, 'settings', indent+3, settings_)
   call write_xml_type_materials_xml( info, 'materials', indent+3, materials_)
   call write_xml_type_tallies_xml( info, 'tallies', indent+3, tallies_)
   write(info%lun,'(a)') '</input>'
   call xml_close(info)
end subroutine

subroutine init_xml_file_input_t

end subroutine

end module
