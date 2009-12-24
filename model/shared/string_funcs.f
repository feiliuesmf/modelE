!@sum string manipulation functions not provided by fortran

      module string_funcs
      implicit none

      private

      public :: lowercase

      integer, parameter :: iA=iachar('A'),iZ=iachar('Z'),
     &     jump_upper_to_lower=iachar('a')-iA
      contains

      function lowercase(str)
      implicit none
      character(len=*) :: str
      character(len=len_trim(str)) :: lowercase
      integer :: p,ic
      character :: c
      do p=1,len_trim(str)
        c = str(p:p)
        ic = iachar(c)
        if(ic >= iA .and. ic <= iZ) c=achar(ic+jump_upper_to_lower) 
        lowercase(p:p) = c
      enddo
      end function lowercase

      end module string_funcs

c      program testit
c      use string_funcs
c      implicit none
c      character(len=10) :: str
c      str='Abc'
c      write(6,*) lowercase(str)
c      end program testit
