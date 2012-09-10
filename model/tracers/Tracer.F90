module Tracer_mod
  use Dictionary_mod, only: Dictionary_type
  use TracerSurfaceSource_mod, only: TracerSurfaceSource
  use TracerSource_mod, only: TracerSource3D
  implicit none
  private

  public :: Tracer_type       ! derived type
  public :: Tracer            ! constructor
  public :: clean

  public :: setProperty
  public :: getName
  public :: getProperty
  public :: getProperties
  public :: hasProperty
  public :: merge
  public :: writeUnformatted
  public :: readUnformatted
  public :: readOneTracer

  public :: findSurfaceSources
  public :: addSurfaceSource
  public :: readSurfaceSources

  public :: NTSURFSRCMAX
!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: NTSURFSRCMAX=16
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: NT3DSRCMAX=7

! TODO make Tracer_type EXTEND Dictionary_type when F2003 adequately supported
  type Tracer_type
!!$    private
    type (Dictionary_type) :: properties
    integer :: ntSurfSrc = 0
    type (TracerSurfaceSource) :: surfaceSources(NTSURFSRCMAX)
    type (TracerSource3D) :: sources3D(NT3DSRCMAX)
  end type Tracer_type

  interface Tracer
    module procedure newTracer
    module procedure newTracerName
    module procedure TracerCopy
  end interface

  interface setProperty
    module procedure setProperty_integer, setProperty_integerArr
    module procedure setProperty_real64,  setProperty_real64Arr
    module procedure setProperty_logical, setProperty_logicalArr
    module procedure setProperty_string,  setProperty_stringArr
  end interface

  interface getProperty
    module procedure getProperty_single
  end interface

  interface getProperties
    module procedure getProperties_single
  end interface

  interface hasProperty
    module procedure hasProperty_single
  end interface

  interface merge
    module procedure merge_one
    module procedure merge_multi
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_tracer
  end interface

  interface readUnformatted
    module procedure readUnformatted_tracer
  end interface

  interface clean
    module procedure cleanTracer
  end interface

contains

  function newTracer() result(aTracer)
!@sum Construct empty tracer    
    use Dictionary_mod, only: Dictionary
    type (Tracer_type) :: aTracer
    aTracer%properties = Dictionary()
  end function newTracer

  function newTracerName(name) result(aTracer)
!@sum Construct named tracer
    use Dictionary_mod, only: Dictionary
    character(len=*), intent(in) :: name
    type (Tracer_type) :: aTracer

    aTracer = Tracer()
    call setProperty(aTracer, 'name', name)
  end function newTracerName

  ! Copy constructor
  function TracerCopy(original) result(copy)
    use Dictionary_mod, only: Dictionary
    type (Tracer_type), intent(in) :: original
    type (Tracer_type) :: copy

    copy%properties = Dictionary(original%properties)

  end function TracerCopy

  subroutine setProperty_integer(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_integer

  subroutine setProperty_integerArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_integerArr

  subroutine setProperty_real64(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_real64

  subroutine setProperty_real64Arr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_real64Arr

  subroutine setProperty_logical(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_logical

  subroutine setProperty_logicalArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_logicalArr

  subroutine setProperty_string(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_string

  subroutine setProperty_stringArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_stringArr

  subroutine merge_one(this, pair)
    use KeyValuePair_mod, only: KeyValuePair_type
    use Dictionary_mod, only: merge
    type (Tracer_type), intent(inout) :: this
    type (KeyValuePair_type), intent(in) :: pair
    
    call merge(this%properties, pair)
  end subroutine merge_one

  function getName(this) result (name)
    use GenericType_mod
    use Dictionary_mod, only: lookup
    use KeyValuePair_mod, only: MAX_LEN_KEY
    type (Tracer_type), intent(in) :: this
    character(len=MAX_LEN_KEY) :: name

    name = lookup(this%properties, 'name')
    
  end function getName

  function getProperty_single(this, property) result (propertyValues)
    use GenericType_mod
    use Dictionary_mod, only: lookup
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property
    type (GenericType_type), pointer :: propertyValues(:)

    propertyValues => lookup(this%properties, property)
    
  end function getProperty_single

  function getProperties_single(this) result(properties)
    type (Tracer_type), target, intent(in) :: this
    type (Dictionary_type), pointer :: properties
    properties => this%properties
  end function getProperties_single

  subroutine merge_multi(this, properties)
    use Dictionary_mod, only: merge
    type (Tracer_type), intent(inout) :: this
    type (Dictionary_type), intent(in) :: properties
    
    call merge(this%properties, properties)
  end subroutine merge_multi

  subroutine writeUnformatted_tracer(this, unit)
!@sum Write a tracer to a unit attached to an unformatted sequential file.
    use Dictionary_mod, only: writeUnformatted
    type (Tracer_type), intent(in) :: this
    integer, intent(in) :: unit

    call writeUnformatted(this%properties, unit)
  end subroutine writeUnformatted_tracer

  subroutine readUnformatted_tracer(this, unit)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
    use Dictionary_mod, only: readUnformatted
    type (Tracer_type), intent(out) :: this
    integer, intent(in) :: unit
    call readUnformatted(this%properties, unit)
  end subroutine readUnformatted_tracer

  function readOneTracer(unit, status) result(aTracer)
    use Parser_mod, only: Parser_type
    use Parser_mod, only: setBeginData, setEndData
    use Parser_mod, only: setCommentCharacters, setTokenSeparators
    use Parser_mod, only: parse

    integer, intent(in) :: unit
    integer, intent(out) :: status
    type (Tracer_type) :: aTracer

    type (Parser_type) :: parser

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aTracer = Tracer()
    aTracer%properties = parse(parser, unit, status)
    if (status /= 0) return

  end function readOneTracer

  logical function hasProperty_single(this, property) result(has)
    use Dictionary_mod, only: hasKey
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property
    has = hasKey(this%properties, property)
  end function hasProperty_single

  subroutine cleanTracer(this)
    use Dictionary_mod, only: clean
    type (Tracer_type), intent(inout) :: this
    call clean(this%properties)
  end subroutine cleanTracer

  subroutine findSurfaceSources(tracer, checkname, sect_name) 
!@sum reads headers from emission files to return
!@+ source names and determine the number of sources
!@+ from the number of files in the rundeck of the form:
!@+ trname_##. Then assigns each source to sector(s),
!@+ based on definitions in the rundeck.
!@auth Greg Faluvegi

    use GenericType_mod

    USE SpecialIO_mod, only: write_parallel
    use MpiSupport_mod, only: am_i_root

    implicit none

!@var nsrc number of source to define ntsurfsrc(n)
    type (Tracer_type), intent(inout) :: tracer
    logical, intent(in) :: checkName
    character*10, intent(in):: sect_name(:)

    integer :: n
    character*80 :: fname
    character(len=300) :: out_line
    logical :: fileExists
    integer :: nsrc

    ! loop through potential number of surface sources, checking if
    ! those files exist. If they do, obtain the source name by reading
    ! the header. If not, the number of sources for this tracer has 
    ! been reached.

    nsrc=0
    if (am_i_root()) &
         &  print*,__LINE__,__FILE__,' tracer = ',  &
         &     trim(getName(tracer)), nsrc, ntsurfsrcmax

    loop_n: do n = 1, ntsurfsrcmax

      fname = addIntegerSuffix(getName(tracer), n)
      inquire(file=trim(fname), exist=fileExists)
      if (am_i_root()) print*,'name: ', trim(fname), fileExists

      if (fileExists) then
        nsrc=nsrc+1
        call addSourceFromFile(tracer, fname)
      else
        exit loop_n
      endif
    enddo loop_n

    ! and make sure there isn't a skip:

    n=n+1
    fname = addIntegerSuffix(getName(tracer), n)
    inquire(file=fname,exist=fileExists)

    if (fileExists) then
      write(out_line,*)'problem with num_srf_sources.', &
           &        ' Possibly missing source? n=',n-1
      call write_parallel(trim(out_line))
      call stop_model(trim(out_line),255)
    endif

  contains

    subroutine addSourceFromFile(tracer, fileName)
      use TracerSurfaceSource_mod, only: initSurfaceSource
      type (Tracer_type), intent(inout) :: tracer
      character(len=*), intent(in) :: fileName

      tracer%ntSurfSrc = tracer%ntSurfSrc + 1
      call initSurfaceSource(tracer%surfaceSources(tracer%ntSurfSrc),  &
           &     getName(tracer), fileName, sect_name, checkname)
    end subroutine addSourceFromFile

  end subroutine findSurfaceSources

  ! Use this routine to add a new surface source that
  ! is manipulated by custom logic elsewhere.
  ! Optional sourcename is only used by diagnostics
  subroutine addSurfaceSource(this, sourceName)
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: sourceName

    this%ntSurfSrc = this%ntSurfSrc + 1
    this%surfaceSources(this%ntSurfSrc)%sourceName = sourceName

  end subroutine addSurfaceSource

!TODO - move to string utilities
  function addIntegerSuffix(tracerName, n) result(fullName)
    character(len=*), intent(in) :: tracerName
    character(len=len_trim(tracerName)+3) :: fullName
    integer, intent(in) :: n
    
    character(len=2) :: suffix
      
    write(suffix,'(I2.2)') n
    fullName = trim(tracerName) // '_' // suffix
  end function addIntegerSuffix

  subroutine readSurfaceSources(tracer, n,nsrc,xyear,xday,checkname,itime,itime_tr0,sfc_src)
!@sum reads surface (2D generally non-interactive) sources
!@auth Jean Lerner/Greg Faluvegi
    USE DOMAIN_DECOMP_ATM, only: GRID
    use TracerSurfaceSource_mod, only: readSurfaceSource
    type (Tracer_type), target, intent(inout) :: tracer
    integer, intent(in) :: nsrc,n
    integer, intent(in) :: xyear, xday
    logical, intent(in) :: checkname
    integer, intent(in) :: itime
    integer, intent(in) :: itime_tr0
    real*8, intent(inout) :: sfc_src(grid%i_strt_halo:,grid%j_strt_halo:,:,:)

    integer :: ns

    if (itime < itime_tr0) return
    if (nsrc <= 0) return

    do ns=1,nsrc
      call readSurfaceSource(tracer%surfaceSources(ns), addIntegerSuffix(getName(tracer), ns), checkname, sfc_src(:,:,n,ns), &
           & xyear, xday)
    enddo

    return

  end subroutine readSurfaceSources

end module Tracer_mod
