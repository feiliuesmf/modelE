module TimeFormatUtilities_mod
   implicit none
   private

   public :: formatHMS
   public :: formatSeconds

   integer, parameter :: r64 =selected_real_kind(14)

contains

   function formatHMS(time) result(hms)
      character(len=15) :: hms
      real (kind=r64), intent(in) :: time

      integer, parameter :: PRECISION = 2

      integer :: hours 
      integer :: minutes
      integer :: seconds
      integer :: fraction

      integer :: totalSeconds
      integer :: totalFraction
      real (kind=r64) :: roundTime
      integer :: remainder

      integer, parameter :: SECONDS_PER_HOUR   = 3600
      integer, parameter :: SECONDS_PER_MINUTE = 60

      call splitTime(time, PRECISION, totalSeconds, fraction)

      hours = totalSeconds / SECONDS_PER_HOUR
      remainder = totalSeconds - hours * SECONDS_PER_HOUR
      minutes = remainder / SECONDS_PER_MINUTE
      seconds = remainder - minutes * SECONDS_PER_MINUTE

      write(hms,"(I3.1,':',I2.2,':',i2.2,'.',i2.2)") &
           & hours, minutes, seconds, fraction

   contains
      
      subroutine splitTime(time, precision, seconds, fraction)
         real (kind=r64), intent(in) :: time
         integer, intent(in)  :: precision
         integer, intent(out) :: seconds
         integer, intent(out) :: fraction

         integer :: denominator
         integer :: numTicks

         denominator = 10**precision
         numTicks = nint(time * denominator) 
         seconds = numTicks / denominator
         fraction = numTicks - denominator*seconds
         
      end subroutine splitTime

   end function formatHMS

   function formatSeconds(x,decimalsBeforePoint,decimalsAfterPoint) result(str)
      real (kind=r64), intent(in) :: x
      integer, intent(in), optional :: decimalsBeforePoint
      integer, intent(in), optional :: decimalsAfterPoint
      character(len=20) :: str
      integer :: whole
      integer :: fraction

      integer :: before, after
      character(len=65) :: fmt

      before = 5
      after  = 2

      if (present(decimalsBeforePoint)) before = decimalsBeforePoint
      if (present(decimalsAfterPoint)) after = decimalsAfterPoint

      write(fmt,'(a,i1,a,i1,a,i1,a)') "(i", before, ".1,'.',i",after,".",after,")"
 
      whole = floor(x)
      fraction = nint(10**after * (x - whole))

      write(str,fmt) whole, fraction
      
   end function formatSeconds

end module TimeFormatUtilities_mod
