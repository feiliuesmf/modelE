      module param
!@sum does all the bookkeeping for input parameters
!@auth I. Aleinov
!@ver 1.0
!@use Module param has the following public subroutines:
!@use set_rparam( name, value ) - add parameter <name> to the database
!@use                             and set it to <value>
!@use get_rparam( name, value ) - extract value of <name> from the database
!@use read_param( kunit ) - read parameters from <kunit> unit
!@use write_param( kunit ) - write paramerers to <kunit> unit
!@use 
!@use Here <name> is a character string with length up to 32
!@use      <value> is real*8 or integer variable (recognized implicitly)
!@use Each parameter can be set only once. If an attempt is made to set
!@use a parameter which is already in the database, an error will be 
!@use generated.

      implicit none
      save
      private

      public set_param, get_param, get_pparam, read_param, write_param
      public is_set_param, alloc_param, sync_param

      integer, parameter :: MAX_PARAMS = 256
      integer, parameter :: MAX_RPARAMS = 64
      integer, parameter :: MAX_IPARAMS = 128
      integer, parameter :: MAX_CPARAMS = 64
      integer, parameter :: MAX_NAME_LEN = 32
      integer, parameter :: MAX_CHAR_LEN = 16

      integer, parameter :: MAGIC = 287649201

      type ParamStr
        character(MAX_NAME_LEN) name  ! parameter name
        integer index                 ! storage for its value
        integer dim                   ! number of elements
        character*1 attrib            ! type: real ('r') or int ('i')
      end type ParamStr

      type (ParamStr), target :: Params(MAX_PARAMS)
      real*8, target :: Rdata(MAX_RPARAMS)
      integer, target :: Idata(MAX_IPARAMS)
      character*(MAX_CHAR_LEN), target :: Cdata(MAX_CPARAMS)
      integer :: num_param = 0
      integer :: num_rparam = 0
      integer :: num_iparam = 0
      integer :: num_cparam = 0

      interface set_param
       module procedure set_iparam, set_rparam, set_cparam
       module procedure  set_aiparam, set_arparam, set_acparam
      end interface

      interface get_param
       module procedure get_iparam, get_rparam, get_cparam
       module procedure get_aiparam, get_arparam, get_acparam
      end interface

      interface get_pparam
       module procedure get_piparam, get_prparam, get_pcparam
       module procedure get_paiparam, get_parparam, get_pacparam
      end interface

      interface alloc_param
      module procedure alloc_iparam, alloc_rparam, alloc_cparam
      module procedure alloc_aiparam, alloc_arparam, alloc_acparam
      end interface

      interface sync_param
      module procedure sync_iparam, sync_rparam, sync_cparam
      module procedure sync_aiparam, sync_arparam, sync_acparam
      end interface

      contains

      function is_set_param( name_in )
      implicit none
      character*(*), intent(in) :: name_in
      logical is_set_param
      integer n
      character*(MAX_NAME_LEN) name

      name = name_in
      call lowcase( name )

      do n=1,num_param
        if ( Params(n)%name == name ) exit
      enddo

      if ( n > num_param  ) then
        is_set_param = .false.
      else
        is_set_param = .true.
      endif
      return
      end function is_set_param

      subroutine set_pstr( name_in, dim, attrib, PStr, flag )
      implicit none
      character*(*), intent(in) :: name_in
      integer, intent(in) :: dim
      character*1, intent(in) ::  attrib
      type (ParamStr), pointer :: PStr
      logical, intent(in) :: flag
      integer n
      character*(MAX_NAME_LEN) name

      if ( len(name_in) > MAX_NAME_LEN ) then
        print *, 'PARAM: parameter name too long: ', name_in
        print *, 'PARAM: maximal length allowed: ', MAX_NAME_LEN
        stop 255
      endif

      name = name_in
      call lowcase( name )

      call get_pstr( name, dim, attrib, PStr )

      if ( associated( PStr ) ) then
        if ( .not. flag ) then
          print *, 'PARAM: attempt to set param which is already set'
          print *, 'name: ', name
          stop 255
        else
          return  ! return PStr found by get_pstr
        endif
      endif

      if ( num_param >= MAX_PARAMS ) then
        print *, 'PARAM: Maximal number of parameters exceeded'
        print *, 'PARAM: Please recompile param with bigger MAX_PARAMS'
        stop 255
      endif

      num_param = num_param + 1
      PStr => Params(num_param)
      PStr%name = name
      PStr%attrib = attrib
      PStr%dim = dim

      select case (attrib)
        case ('i')
          if ( num_iparam+dim >= MAX_IPARAMS ) then
            print *, 'PARAM: Maximal number of int parameters exceeded'
            print *, 'PARAM: Recompile param with bigger MAX_IPARAMS'
            stop 255
          endif
          PStr%index = num_iparam + 1
          num_iparam = num_iparam + dim
        case ('r')
          if ( num_rparam+dim >= MAX_RPARAMS ) then
            print *, 'PARAM: Maximal number of real parameters exceeded'
            print *, 'PARAM: Recompile param with bigger MAX_RPARAMS'
            stop 255
          endif
          PStr%index = num_rparam + 1
          num_rparam = num_rparam + dim
        case ('c')
          if ( num_cparam+dim >= MAX_CPARAMS ) then
            print *, 'PARAM: Maximal number of char parameters exceeded'
            print *, 'PARAM: Recompile param with bigger MAX_CPARAMS'
            stop 255
          endif
          PStr%index = num_cparam + 1
          num_cparam = num_cparam + dim
        end select

      return
      end subroutine set_pstr


      subroutine get_pstr( name_in, dim, attrib, PStr )
      implicit none
      character*(*), intent(in) :: name_in
      integer, intent(in) :: dim
      character*1, intent(in) ::  attrib
      type (ParamStr), pointer :: PStr     
      integer n
      character*(MAX_NAME_LEN) name

      name = name_in
      call lowcase( name )

      nullify( PStr )

      do n=1,num_param
        if ( Params(n)%name == name ) exit
      enddo

      if ( n > num_param  ) return   ! not found - return NULL

      if ( Params(n)%attrib /= attrib .or. Params(n)%dim /= dim ) then
        print *, 'PARAM: wrong type or dim of parameter: ', name
        print *, 'ATT: set: ', Params(n)%attrib, ' called: ', attrib
        print *, 'DIM: set: ', Params(n)%dim, ' called: ', dim
        stop 255
      endif
    
      PStr => Params(n)
      return
      end subroutine get_pstr



      !***** integers ******!

      
      subroutine set_iparam( name, value, opt )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: value
      character*1, optional, intent(in) :: opt
      integer v(1)
      v(1) = value
      call set_aiparam( name, v, 1, opt )
      return
      end subroutine set_iparam
 

      subroutine set_aiparam( name, value, np, opt )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: value(np)
      integer, intent(in) :: np
      character*1, optional, intent(in) :: opt
      type (ParamStr), pointer :: PStr
      logical flag

      flag = .false.
      if ( present(opt) .and. opt=='o' ) flag = .true.

      call set_pstr( name, np, 'i', PStr, flag )
      Idata( PStr%index : PStr%index+np-1 ) = value(1:np)
      return
      end subroutine set_aiparam


      subroutine get_iparam( name, value )
      implicit none
      character*(*), intent(in) ::  name
      integer, intent(out) ::  value
      integer v(1)

      call get_aiparam( name, v, 1 )
      value = v(1)
      return
      end subroutine get_iparam


      subroutine get_aiparam( name, value, np )
      implicit none
      character*(*), intent(in) ::  name
      integer, intent(out) ::  value(np)
      integer, intent (in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'i', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      value(1:np) = Idata( PStr%index : PStr%index+np-1 )
      return
      end subroutine get_aiparam


      subroutine get_piparam( name, pvalue )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue
      type (ParamStr), pointer :: PStr

      call get_pstr( name, 1, 'i', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      pvalue => Idata( PStr%index )
      return
      end subroutine get_piparam

      subroutine get_paiparam( name, pvalue, np )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue(:)
      integer, intent(in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'i', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      pvalue => Idata( PStr%index:PStr%index+np-1 )
      return
      end subroutine get_paiparam



      !***** reals ******!

      subroutine set_rparam( name, value, opt )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(in) :: value
      character*1, optional, intent(in) :: opt
      real*8 v(1)
      v(1) = value
      call set_arparam( name, v, 1, opt )
      return
      end subroutine set_rparam
 

      subroutine set_arparam( name, value, np, opt )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(in) :: value(np)
      integer, intent(in) :: np
      character*1, optional, intent(in) :: opt
      type (ParamStr), pointer :: PStr
      logical flag

      flag = .false.
      if ( present(opt) .and. opt=='o' ) flag = .true.

      call set_pstr( name, np, 'r', PStr, flag )
      Rdata( PStr%index : PStr%index+np-1 ) = value(1:np)
      return
      end subroutine set_arparam


      subroutine get_rparam( name, value )
      implicit none
      character*(*), intent(in) ::  name
      real*8, intent(out) ::  value
      real*8 v(1)

      call get_arparam( name, v, 1 )
      value = v(1)
      return
      end subroutine get_rparam


      subroutine get_arparam( name, value, np )
      implicit none
      character*(*), intent(in) ::  name
      real*8, intent(out) ::  value(np)
      integer, intent (in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'r', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      value(1:np) = Rdata( PStr%index : PStr%index+np-1 )
      return
      end subroutine get_arparam


      subroutine get_prparam( name, pvalue )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue
      type (ParamStr), pointer :: PStr

      call get_pstr( name, 1, 'r', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      pvalue => Rdata( PStr%index )
      return
      end subroutine get_prparam

      subroutine get_parparam( name, pvalue, np )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue(:)
      integer, intent(in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'r', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      pvalue => Rdata( PStr%index:PStr%index+np-1 )
      return
      end subroutine get_parparam


      !***** Chars ******!

      subroutine set_cparam( name, value, opt )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(in) :: value
      character*1, optional, intent(in) :: opt
      character*(MAX_CHAR_LEN) v(1)
      if ( len(value) > MAX_CHAR_LEN ) then
        print *, 'PARAM: Char string too long. MAX = ', MAX_CHAR_LEN
        stop 255
      endif
      v(1) = value
      call set_acparam( name, v, 1, opt )
      return
      end subroutine set_cparam

c$$$      subroutine set_cparam( name, value )
c$$$      implicit none
c$$$      character*(*), intent(in) :: name
c$$$      character*(*), intent(in) :: value
c$$$      integer np
c$$$      type (ParamStr), pointer :: PStr
c$$$      np = len (value )
c$$$      call set_pstr( name, np, 'c', PStr )
c$$$      Cdata( PStr%index : PStr%index+np-1 ) = value(1:np)
c$$$      return
c$$$      end subroutine set_cparam
    

      subroutine set_acparam( name, value, np, opt )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(in) :: value(np)
      integer, intent(in) :: np
      character*1, optional, intent(in) :: opt
      type (ParamStr), pointer :: PStr
      integer n
      logical flag

      flag = .false.
      if ( present(opt) .and. opt=='o' ) flag = .true.

      do n=1,np
        if ( len(value(n)) > MAX_CHAR_LEN ) then
          print *, 'PARAM: Char string too long. MAX = ', MAX_CHAR_LEN
          print *, 'You submitted LEN = ', len(value)
          stop 255
        endif
      enddo
      call set_pstr( name, np, 'c', PStr, flag )
      Cdata( PStr%index : PStr%index+np-1 ) = value(1:np)
      return
      end subroutine set_acparam


      subroutine get_cparam( name, value )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), intent(out) ::  value
      character*(MAX_CHAR_LEN) v(1)

      call get_acparam( name, v, 1 )
      value = v(1)
      return
      end subroutine get_cparam


      subroutine get_acparam( name, value, np )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), intent(out) ::  value(np)
      integer, intent (in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'c', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      value(1:np) = Cdata( PStr%index : PStr%index+np-1 )
      return
      end subroutine get_acparam


      subroutine get_pcparam( name, pvalue )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue
      type (ParamStr), pointer :: PStr

      call get_pstr( name, 1, 'c', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      pvalue => Cdata( PStr%index )
      return
      end subroutine get_pcparam

      subroutine get_pacparam( name, pvalue, np )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue(:)
      integer, intent(in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'c', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        stop 255
      endif
      pvalue => Cdata( PStr%index:PStr%index+np-1 )
      return
      end subroutine get_pacparam

      !***** alloc fuctions ******

      subroutine alloc_iparam( name, pvalue, initval )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue
      integer, intent(in) :: initval
      call set_param( name, initval )
      call get_pparam( name, pvalue )
      end subroutine alloc_iparam


      subroutine alloc_aiparam( name, pvalue, initval, dim )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue(:)
      integer, intent(in) :: initval(:)
      integer, intent(in) :: dim
      call set_param( name, initval, dim )
      call get_pparam( name, pvalue, dim )
      end subroutine alloc_aiparam


      subroutine alloc_rparam( name, pvalue, initval )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue
      real*8, intent(in) :: initval
      call set_param( name, initval )
      call get_pparam( name, pvalue )
      end subroutine alloc_rparam


      subroutine alloc_arparam( name, pvalue, initval, dim )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue(:)
      real*8, intent(in) :: initval(:)
      integer, intent(in) :: dim
      call set_param( name, initval, dim )
      call get_pparam( name, pvalue, dim )
      end subroutine alloc_arparam


      subroutine alloc_cparam( name, pvalue, initval )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue
      character*(*), intent(in) :: initval
      call set_param( name, initval )
      call get_pparam( name, pvalue )
      end subroutine alloc_cparam


      subroutine alloc_acparam( name, pvalue, initval, dim )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue(:)
      character*(*), intent(in) :: initval(:)
      integer, intent(in) :: dim
      call set_param( name, initval, dim )
      call get_pparam( name, pvalue, dim )
      end subroutine alloc_acparam

      
      !***** sync functions ******!

      subroutine sync_iparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(inout) :: value
      
      if ( is_set_param( name ) ) then
        call get_param( name, value )
      else
        call set_param( name, value )
      endif
      end subroutine sync_iparam


      subroutine sync_aiparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(inout) :: value(np)
      integer, intent(in) :: np
      
      if ( is_set_param( name ) ) then
        call get_param( name, value, np )
      else
        call set_param( name, value, np )
      endif
      end subroutine sync_aiparam


      subroutine sync_rparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(inout) :: value
      
      if ( is_set_param( name ) ) then
        call get_param( name, value )
      else
        call set_param( name, value )
      endif
      end subroutine sync_rparam


      subroutine sync_arparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(inout) :: value(np)
      integer, intent(in) :: np
      
      if ( is_set_param( name ) ) then
        call get_param( name, value, np )
      else
        call set_param( name, value, np )
      endif
      end subroutine sync_arparam


      subroutine sync_cparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(inout) :: value
      
      if ( is_set_param( name ) ) then
        call get_param( name, value )
      else
        call set_param( name, value )
      endif
      end subroutine sync_cparam


      subroutine sync_acparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(inout) :: value(np)
      integer, intent(in) :: np
      
      if ( is_set_param( name ) ) then
        call get_param( name, value, np )
      else
        call set_param( name, value, np )
      endif
      end subroutine sync_acparam



      !***** input / output ******!

      subroutine read_param( kunit, ovrwrt )
      implicit none
      integer, intent(in) :: kunit
      logical, intent(in) :: ovrwrt
      integer n, np
      type (ParamStr), save :: LParams(MAX_PARAMS)
      real*8, save :: LRdata(MAX_RPARAMS)
      integer, save :: LIdata(MAX_IPARAMS)
      character*(MAX_CHAR_LEN), save :: LCdata(MAX_CPARAMS)
      integer lnum_param, lnum_rparam, lnum_iparam, lnum_cparam
      integer lmagic

      read( kunit, err=10 ) lmagic
      backspace kunit

      if ( lmagic /= MAGIC ) then
        print *, 'WARNING: No parameter data'
        return
      endif

      print *, 'READING PARAMETERS from rsf'
 
      read( kunit, err=10 ) lmagic,
     *     lnum_param, lnum_rparam, lnum_iparam, lnum_cparam,
     *     ( LParams(n), n=1,min(lnum_param,MAX_PARAMS) ),
     *     ( LRdata(n), n=1,min(lnum_rparam,MAX_RPARAMS) ),
     *     ( LIdata(n), n=1,min(lnum_iparam,MAX_IPARAMS) ),
     *     ( LCdata(n), n=1,min(lnum_cparam,MAX_CPARAMS) )

      if (     lnum_param  > MAX_PARAMS 
     *     .or. lnum_rparam > MAX_RPARAMS 
     *     .or. lnum_iparam > MAX_IPARAMS 
     *     .or. lnum_cparam > MAX_CPARAMS 
     *     ) then
        print *, 'PARAM: parameter list in input file too long'
        print *, 'PARAM: please recompile param with bigger MAX_?PARAMS'
        print *, 'PARAM: ',num_param,num_rparam,num_iparam,num_cparam
        stop 255
      endif

      ! now merge the data just read with existing database
      do n=1,lnum_param
        !print *,'read: ',LParams(n)
        if ( (.not. is_set_param(LParams(n)%name)) .or. ovrwrt ) then
          select case( LParams(n)%attrib )
          case ('i')
            call set_aiparam( LParams(n)%name, LIdata(LParams(n)%index),
     *           LParams(n)%dim, 'o' )
          case ('r')
            call set_arparam( LParams(n)%name, LRdata(LParams(n)%index),
     *           LParams(n)%dim, 'o' )
          case ('c')
            call set_acparam( LParams(n)%name, LCdata(LParams(n)%index),
     *           LParams(n)%dim, 'o' )
          end select
        endif
      enddo
      return
 10   print *, 'PARAM: Error reading, unit = ', kunit
      stop 255
      end subroutine read_param
  

      subroutine write_param( kunit )
      implicit none
      integer, intent(in) :: kunit
      integer n
      
      write( kunit, err=10 ) MAGIC,
     *     num_param, num_rparam, num_iparam, num_cparam,
     *     ( Params(n), n=1,min(num_param,MAX_PARAMS) ),
     *     ( Rdata(n), n=1,min(num_rparam,MAX_RPARAMS) ),
     *     ( Idata(n), n=1,min(num_iparam,MAX_IPARAMS) ),
     *     ( Cdata(n), n=1,min(num_cparam,MAX_CPARAMS) )
      return
 10   print *, 'PARAM: Error writing, unit = ', kunit
      stop 255
      end subroutine write_param


      subroutine lowcase( str )
      implicit none
      character*(*) str
      integer n, i
      integer A, Z, shift, c
 
      A = iachar( 'A' )
      Z = iachar( 'Z' )
      shift = iachar( 'a' ) - iachar( 'A' )
      n = len_trim(str)
      do i=1,n
         c = iachar( str(i:i) )
        if ( c>=A .and. c<=Z ) str(i:i) = achar( c + shift )
      enddo
      end subroutine lowcase


      end module param


#ifdef TEST_MODULE_PARAM
ccc this part is for testing - will be removed soon

  ! trying to add some code 

      program foo
      use param
      implicit none
      real*8 a, b, w(10)
      integer i, ai(3)
      integer, pointer :: pi(:)
      character*32 c
      character*16 cc(4)
      integer, pointer :: ppi(:), ppi1
      real*8, pointer :: ppr(:), ppr1
      character*4, pointer :: ppc(:), ppc1
      integer :: d6(2,4) 
      data d6 /63,17, 17,34, 37,27, 13,23/
      integer dd(8)

      a = 1234.77d0
      b = 122

      w(1) = 9.d0
      w(3) = 11.d0
      w(10) = 13.d0

      ai(1) = 11
      ai(2) = 12
      ai(3) = 13+100
      
      cc(1) = 'fd'
      cc(2) = '12'
      cc(3) = 'a'
      
      !call read_param( 99 )

      dd(1:8) =  d6(1:8,1)
      call set_param( "d6", d6(1:8,1), 8 )

      d6(:,:) = 0
      print *, dd, d6

      call get_param( "d6", d6(1:8,1), 8 )
      print *, dd, d6

      call alloc_param( "ppi", ppi, (/1,2,3/), 3)
      call alloc_param( "ppi1", ppi1, 7)

      print *, 'ppi = ', (ppi(i), i=1,3), ' :: ',ppi1

      ppi(2) = 22

      print *, 'ppi = ', (ppi(i), i=1,3), ' :: ',ppi1

      call get_param( 'ppi', ai, 3 )
      
      print *, 'ppi = ', ai



      call alloc_param( "ppr", ppr, (/11.d0,12.d0,13.d0/), 3)
      call alloc_param( "ppr1", ppr1, 77.d0)

      print *, 'ppr = ', (ppr(i), i=1,3), ' :: ',ppr1

      ppr(2) = 222.

      print *, 'ppr = ', (ppr(i), i=1,3), ' :: ',ppr1

      call get_param( 'ppr', w , 3 )
      
      print *, 'ppr = ', w


      call alloc_param( "ppc", ppc, (/'asdf','ghjk','qwer','tyui'/), 4)
      call alloc_param( "ppc1", ppc1, 'abc')
      print *, 'ppc = ', (ppc(i), i=1,4), ' :: ',ppc1

      ppc(3) = 'JAN'

      print *, 'ppc = ', (ppc(i), i=1,4), ' :: ',ppc1

      call get_param( 'ppc', cc, 4 )
      
      print *, 'ppc = ', cc


      a =  999999

      !call set_param( 'par_a', a )
      call set_param( 'par_c', 'test23_oOOO', 'o' )
      call set_param( 'parbb_b', b )
      !call set_param( 'par_i', ai, 3, 'o' )
      call set_param( 'par_cc', cc, 3 )

      call read_param( 99, .false. )

      print *, is_set_param('par_c'), is_set_param('parbb_b'),
     *     is_set_param('par_a'), is_set_param('abc'),
     $     is_set_param('par_c')

      call get_param( 'par_c', c )
      print *, 'c = ', c

      print *, 'a = ', a
      call sync_param( 'par_a', a )
      print *, 'a = ', a

      call get_param( 'par_a', a )
      print *, 'a = ', a

      call get_param( 'par_cc', cc, 3 )
      print *, 'cc = ', cc
      

      !call get_pparam( 'par_i', pi, 3 )

      call get_param( 'parbb_b', a )
      print *, 'b = ', a

      call set_param( 'par_d', 12301, 'o' )
     

      call get_param( 'par_c', c )
      print *, 'c = ', c

      !call get_param( 'par_i', ai, 3 )
      print *, 'i = ', ai

      !pi(1) = 2001
      !pi(2) = 2002

      ai(1) = 1024
      ai(2) = 2048
      print *, 'ai = ', ai
      call sync_param( 'par_i', ai, 3 )
      print *, 'ai = ', ai
      call get_param( 'par_i', ai, 3 )
      print *, 'ai = ', ai


      call get_param( 'par_d', i )
      print *, 'par_d = ', i

      call write_param( 98 )

      end program foo
#endif




