module SpecialFunctions_mod
   implicit none
   private

   public :: erf 

   interface erf
     module procedure erf_single
     module procedure erf_double
   end interface erf

   integer, parameter :: DP = selected_real_kind(14)

contains
  
  ! Based upon Abromowitz and Stegun 7.1.26
  ! http://en.wikipedia.org/wiki/Abramowitz_and_Stegun
  ! Original source:
  ! C. Hastings, Jr.  Approximations for Digital Computers. Princeton University Press, 1955.
 
  real(kind=DP) function erf_double(x)
    real(kind=DP), intent(in) :: x

    real(kind=DP), parameter ::  p = +0.3275911
    real(kind=DP), parameter :: a1 = +0.254829592
    real(kind=DP), parameter :: a2 = -0.284496736
    real(kind=DP), parameter :: a3 = +1.421413741
    real(kind=DP), parameter :: a4 = -1.453152027
    real(kind=DP), parameter :: a5 = +1.061405429

    real(kind=DP) :: t

    t = 1/(1 + p *abs(x))
    erf_double = 1 - (a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5)*exp(-x**2)

    if (x < 0) erf_double = -erf_double

  end function erf_double


  real function erf_single(x)
    real, intent(in) :: x

    erf_single = erf(real(x,kind=DP))

  end function erf_single


end module SpecialFunctions_mod
