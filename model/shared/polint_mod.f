      module polint_mod
      implicit none
      contains

      SUBROUTINE polint2dlin(x1a,x2a,ya,m,n,x1,x2,y,dy)

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: m,n
      REAL*8,INTENT(IN) :: x1a(m),x2a(n),ya(m,n),x1,x2
      REAL*8,INTENT(OUT) :: dy,y
      INTEGER, PARAMETER :: nmax=2
      INTEGER j,k,jjj,iii,xx,yy
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax)

      call locate(x1a,m,x1,xx)
      call locate(x2a,n,x2,yy)

      do k=1,nmax
        iii=k
        do j=1,nmax
          jjj=j
          yntmp(j)=ya(xx+jjj-1,yy+iii-1)
          x11(j)=x1a(xx+jjj-1)
        enddo
        if (yntmp(1).eq.-1000.) then
          ymtmp(k)=-1000.
          x22(k)=x2a(yy+iii-1)
        else
          call polint(x11,yntmp,nmax,x1,ymtmp(k),dy)
          x22(k)=x2a(yy+iii-1)
        endif
      enddo
      if (ymtmp(1).eq.-1000.) then
        y=-1000.
      else
        call polint(x22,ymtmp,nmax,x2,y,dy)
      endif

      return
      END SUBROUTINE polint2dlin
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint2dcub(x1a,x2a,ya,m,n,x1,x2,y,dy)

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: m,n
      REAL*8,INTENT(IN) :: x1a(m),x2a(n),ya(m,n),x1,x2
      REAL*8,INTENT(OUT) :: dy,y
      INTEGER, PARAMETER :: nmax=4
      INTEGER j,k,jjj,iii,xx,yy
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax)

      call locate(x1a,m,x1,xx)
      call locate(x2a,n,x2,yy)

      do k=1,nmax
        iii=k-1
        do j=1,nmax
          jjj=j-1
          yntmp(j)=ya(xx+jjj-1,yy+iii-1)
          x11(j)=x1a(xx+jjj-1)
        enddo
        if (yntmp(1).eq.-1000..or.yntmp(2).eq.-1000.) then
          ymtmp(k)=-1000.
          x22(k)=x2a(yy+iii-1)
        else
          call polint(x11,yntmp,nmax,x1,ymtmp(k),dy)
          x22(k)=x2a(yy+iii-1)
        endif
      enddo
      if (ymtmp(1).eq.-1000..or.ymtmp(2).eq.-1000.) then
        y=-1000.
      else
        call polint(x22,ymtmp,nmax,x2,y,dy)
      endif

      return
      END SUBROUTINE polint2dcub
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint3dlin(x1a,x2a,x3a,ya,m,n,lkm,x1,x2,x3,y,dy) 
 
      implicit none 
      INTEGER, INTENT(IN) :: m,n,lkm 
      REAL*8, INTENT(IN) :: x1,x2,x3,x1a(m),x2a(n),x3a(lkm),ya(m,n,lkm) 
      REAL*8, INTENT(OUT):: y,dy 
      INTEGER, PARAMETER :: nmax=2
      INTEGER i,j,k,jjj,iii,xx,yy,zz,kkk 
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax),x33(nmax) 
      real*8 yotmp(nmax) 
 
      call locate(x1a,m,x1,xx) 
      call locate(x2a,n,x2,yy) 
      call locate(x3a,lkm,x3,zz) 

      do i=1,nmax
         kkk=i
         x33(i)=x3a(zz+kkk-1)
         do k=1,nmax
            iii=k
            x22(k)=x2a(yy+iii-1)  
            do j=1,nmax
               jjj=j
               x11(j)=x1a(xx+jjj-1)
               yntmp(j)=ya(xx+jjj-1,yy+iii-1,zz+kkk-1)
            enddo
            if (yntmp(1).eq.-1000) then
               ymtmp(k)=-1000.
            else
               call polint(x11,yntmp,nmax,x1,ymtmp(k),dy)
            endif
         enddo
         if (ymtmp(1).eq.-1000) then
            yotmp(i)=-1000.
         else
            call polint(x22,ymtmp,nmax,x2,yotmp(i),dy)
         endif
      enddo
      if (yotmp(2).eq.-1000)  then
         y=-1000.
      else
         call polint(x33,yotmp,nmax,x3,y,dy)
      endif

      return
      END SUBROUTINE POLINT3DLIN
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint3dcub(x1a,x2a,x3a,ya,m,n,lkm,x1,x2,x3,y,dy)

      implicit none
      INTEGER, INTENT(IN) :: m,n,lkm
      REAL*8, INTENT(IN) :: x1,x2,x3,x1a(m),x2a(n),x3a(lkm),ya(m,n,lkm)
      REAL*8, INTENT(OUT):: y,dy
      INTEGER, PARAMETER :: nmax=4
      INTEGER i,j,k,jjj,iii,xx,yy,zz,kkk
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax),x33(nmax)
      real*8 yotmp(nmax)

      call locate(x1a,m,x1,xx)
      call locate(x2a,n,x2,yy)
      call locate(x3a,lkm,x3,zz)

      do i=1,nmax
         kkk=i-1
         x33(i)=x3a(zz+kkk-1)
         do k=1,nmax
            iii=k-1
            x22(k)=x2a(yy+iii-1)
            do j=1,nmax
               jjj=j-1
               x11(j)=x1a(xx+jjj-1)
               yntmp(j)=ya(xx+jjj-1,yy+iii-1,zz+kkk-1)
            enddo
            if (yntmp(1).eq.-1000.or.yntmp(2).eq.-1000.) then
               ymtmp(k)=-1000.
            else
               call polint(x11,yntmp,nmax,x1,ymtmp(k),dy)
            endif
         enddo
         if (ymtmp(1).eq.-1000.or.ymtmp(2).eq.-1000.) then
            yotmp(i)=-1000.
         else
            call polint(x22,ymtmp,nmax,x2,yotmp(i),dy)
         endif
      enddo
      if (yotmp(3).eq.-1000.or.yotmp(4).eq.-1000.)  then
         y=-1000.
      else
         call polint(x33,yotmp,nmax,x3,y,dy)
      endif

      return
      END SUBROUTINE POLINT3DCUB
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint(xa,ya,n,x,y,dy)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(OUT) :: y,dy
      REAL*8, INTENT(IN) :: x,xa(n),ya(n)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(n),d(n)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0) CALL stop_model('failure in polint',255)
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue

      return
      END SUBROUTINE POLINT
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      subroutine locate(xx,n,x,j)
!@sum locates parameters of integration in lookup table

      implicit none
      INTEGER, INTENT(IN):: n
      INTEGER, INTENT(OUT):: j
      REAL*8, INTENT(IN):: x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl

      return
      end subroutine locate

      end module polint_mod
