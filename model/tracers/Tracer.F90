module Tracer_mod
  use Dictionary_mod, only: Dictionary
  use TracerSurfaceSource_mod, only: TracerSurfaceSource
  use TracerSource_mod, only: TracerSource3D
  use AttributeDictionary_mod
  implicit none
  private

  public :: Tracer       ! derived type
  public :: newTracer            ! constructor
  public :: clean

  public :: getName
  public :: writeUnformatted
  public :: readUnformattedTracer
  public :: readOneTracer

  public :: findSurfaceSources
  public :: addSurfaceSource
  public :: readSurfaceSources
  public :: assignment(=)

  public :: NTSURFSRCMAX
  public :: copyInto
!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: NTSURFSRCMAX=16
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: NT3DSRCMAX=7

  type, extends(AttributeDictionary) :: Tracer
!!$    private
    type (Dictionary) :: properties
    integer :: ntSurfSrc = 0
    type (TracerSurfaceSource) :: surfaceSources(NTSURFSRCMAX)
    type (TracerSource3D) :: sources3D(NT3DSRCMAX)
  end type Tracer

  interface newTracer
    module procedure newTracerName
    module procedure newEmptyTracer
    module procedure TracerCopy
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_tracer
  end interface

  interface assignment(=)
    module procedure toTracer
  end interface assignment(=)

  interface clean
    module procedure cleanTracer
  end interface

contains

  ! private constructor
  function newEmptyTracer() result(aTracer)
!@sum Construct empty tracer    
    use Dictionary_mod, only: Dictionary
    type (Tracer), pointer :: aTracer

    allocate(aTracer)
    aTracer%AttributeDictionary = newAttributeDictionary()

  end function newEmptyTracer

  function newTracerName(name) result(aTracer)
!@sum Construct named tracer
    use Dictionary_mod, only: Dictionary
    character(len=*), intent(in) :: name
    type (Tracer), pointer :: aTracer

    aTracer => newEmptyTracer()
    call aTracer%insert('name', trim(name))
    aTracer%ntSurfsrc = 0
    
  end function newTracerName

  ! Copy constructor
  function TracerCopy(original) result(copy)
    use Dictionary_mod, only: Dictionary
    type (Tracer), intent(in) :: original
    type (Tracer) :: copy

    copy%AttributeDictionary = original%AttributeDictionary

  end function TracerCopy

  function getName(this) result (name)
    use AbstractAttribute_mod
    type (Tracer), intent(in) :: this
    character(len=MAX_LEN_ATTRIBUTE_STRING), pointer :: name
    class (AbstractAttribute), pointer :: p

    ! TODO Intel is now struggling with the line below - no idea why.  Worked before other changes.
!!$    name = this%getReference('name')
    p => this%getReference('name')
    name = p
    
  end function getName

  subroutine writeUnformatted_tracer(this, unit)
!@sum Write a tracer to a unit attached to an unformatted sequential file.
    type (Tracer), intent(in) :: this
    integer, intent(in) :: unit

!!$    call this%properties%writeUnformatted(unit)
    call this%writeUnformatted(unit)
    
  end subroutine writeUnformatted_tracer

  subroutine readUnformattedTracer(this, unit)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
!!$    use Dictionary_mod, only: readUnformatted
    type (Tracer), intent(inout) :: this
    integer, intent(in) :: unit
!!$    call readUnformatted(this%properties, unit)
    call this%readUnformatted(unit)
  end subroutine readUnformattedTracer

  function readOneTracer(unit, status) result(aTracer)
    use Parser_mod, only: Parser_type
    use Parser_mod, only: setBeginData, setEndData
    use Parser_mod, only: setCommentCharacters, setTokenSeparators
    use Parser_mod, only: parse

    integer, intent(in) :: unit
    integer, intent(out) :: status
    type (Tracer), pointer :: aTracer

    type (Parser_type) :: parser

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aTracer => newEmptyTracer()
    aTracer%AttributeDictionary = parse(parser, unit, status)

!!$    aTracer%properties = parse(parser, unit, status)
    
    if (status /= 0) return

  end function readOneTracer

  subroutine cleanTracer(this)
    use Dictionary_mod, only: clean
    type (Tracer), intent(inout) :: this
  end subroutine cleanTracer

  subroutine findSurfaceSources(trcer, checkname, sect_name) 
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
    type (Tracer), intent(inout) :: trcer
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
         &     trim(getName(trcer)), nsrc, ntsurfsrcmax

    loop_n: do n = 1, ntsurfsrcmax

      fname = addIntegerSuffix(getName(trcer), n)
      inquire(file=trim(fname), exist=fileExists)
      if (am_i_root()) print*,'name: ', trim(fname), fileExists

      if (fileExists) then
        nsrc=nsrc+1
        call addSourceFromFile(trcer, fname)
      else
        exit loop_n
      endif
    enddo loop_n

    ! and make sure there isn't a skip:

    n=n+1
    fname = addIntegerSuffix(getName(trcer), n)
    inquire(file=fname,exist=fileExists)

    if (fileExists) then
      write(out_line,*)'problem with num_srf_sources.', &
           &        ' Possibly missing source? n=',n-1
      call write_parallel(trim(out_line))
      call stop_model(trim(out_line),255)
    endif

  contains

    subroutine addSourceFromFile(trcer, fileName)
      use TracerSurfaceSource_mod, only: initSurfaceSource
      type (Tracer), intent(inout) :: trcer
      character(len=*), intent(in) :: fileName

      trcer%ntSurfSrc = trcer%ntSurfSrc + 1
      call initSurfaceSource(trcer%surfaceSources(trcer%ntSurfSrc),  &
           &     getName(trcer), fileName, sect_name, checkname)
    end subroutine addSourceFromFile

  end subroutine findSurfaceSources

  ! Use this routine to add a new surface source that
  ! is manipulated by custom logic elsewhere.
  ! Optional sourcename is only used by diagnostics
  subroutine addSurfaceSource(this, sourceName)
    type (Tracer), intent(inout) :: this
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

  subroutine readSurfaceSources(trcer, n,nsrc,xyear,xday,checkname,itime,itime_tr0,sfc_src)
!@sum reads surface (2D generally non-interactive) sources
!@auth Jean Lerner/Greg Faluvegi
    USE DOMAIN_DECOMP_ATM, only: GRID
    use TracerSurfaceSource_mod, only: readSurfaceSource
    type (Tracer), target, intent(inout) :: trcer
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
      call readSurfaceSource(trcer%surfaceSources(ns), addIntegerSuffix(getName(trcer), ns), checkname, sfc_src(:,:,n,ns), &
           & xyear, xday)
    enddo

    return

  end subroutine readSurfaceSources

  subroutine toTracer(pType, pClass)
    type (Tracer), pointer, intent(out) :: pType
    class (Tracer), target, intent(in) :: pClass

    select type (p => pClass)
    type is (Tracer)
      pType => p
    class default
      call stop_model('Illegal conversion in Tracer_mod.',255)
    end select

    
  end subroutine toTracer

  subroutine copyInto(a, b)
    type (Tracer), intent(out) :: a
    type (Tracer), intent(in) :: b

    a%properties = b%properties
    a%AttributeDictionary = b%AttributeDictionary
    a%ntSurfSrc = b%ntSurfSrc
    a%surfaceSources = b%surfaceSources
    a%sources3D = b%sources3D
    
  end subroutine copyInto

end module Tracer_mod
