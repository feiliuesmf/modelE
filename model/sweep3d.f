
      subroutine sweep3diag( a, b, c, f, x, n)
ccc solves 3-diagonal system:
ccc
ccc a x_{i-1} + c x_i + b_{i+1} = f_i
ccc
ccc !!! caution: some input arrays are changed inside the program !!!
ccc
      implicit none
      REAL*8 EPS
      parameter( EPS=1.d-30 )
ccc input:
      REAL*8 a(0:1), b(0:1), c(0:1), f(0:1)
      integer n
ccc output:
      REAL*8 x(0:1)
ccc working variables:
      integer k
      REAL*8 denom

ccc corresponding variables in Samarskii book:
ccc  x(0) = beta(1)   f(0) = alpha(1)

      if ( abs(c(0)) .lt. EPS )
     &   call error_abort('sweep3diag: zero denominator  ',30)
      x(0) = f(0)/c(0)
      if ( n .lt. 2 ) return
      f(0) = -b(0)/c(0)

      do k=1,n-2
        denom = c(k) + a(k)*f(k-1)
        if ( abs(denom) .lt. EPS )
     &      call error_abort('sweep3diag: zero denominator  ',30)
        x(k) = ( f(k) - a(k)*x(k-1) )/denom
        f(k) = -b(k)/denom
        enddo

      denom = c(n-1) + a(n-1)*f(n-2)
      if ( abs(denom) .lt. EPS )
     &      call error_abort('sweep3diag: zero denominator  ',30)
      x(n-1) = ( f(n-1) - a(n-1)*x(n-2) )/denom

      do k=n-2,0,-1
        x(k) = f(k)*x(k+1) + x(k)
        enddo

      return
      end




      subroutine error_abort( str , len)
      implicit none
      integer len
      character str(len)
      print *, 'ERROR: ', str
      stop 66
      end

