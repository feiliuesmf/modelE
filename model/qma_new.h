c
c Routines:
c read_qma:       reads stratospheric WV input file,
c                 and vertically regrids to model resolution
c lat_interp_qma: horizontally regrids to model latitudes
c 

      subroutine read_qma (iu,plb)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
!@ver  1.0
      use domain_decomp, only : write_parallel
      use rad_com, only : dH2O,jma=>jm_dh2o,lat_dh2o
      use model_com, only : lm
      use constant, only : radian
      implicit none
      integer, parameter:: lma=24
      integer m,iu,j,l,ll,ldn(lm),lup(lm)
      real*8 :: plb(lm+1)
      real*4 pb(0:lma+1),h2o(jma,0:lma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,dh,fracl
      character(len=300) :: out_line

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (lat_dh2o(j),j=1,jma)
      lat_dh2o(:) = lat_dh2o(:)*radian

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) vertically
        do j=1,jma
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+h2o(j,ll)*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/365. !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine read_qma

      subroutine lat_interp_qma (rlat,lev,mon,dh2o_interp)
!@sum  interpolate CH4->H2O production rates in latitude
!@auth R. Ruedy
!@ver  1.0
      use rad_com, only : jma=>jm_dh2o,xlat=>lat_dh2o,dh2o
      implicit none
      real*8 :: rlat ! input latitude (radians)
      integer :: lev,mon ! input level, month
      real*8 :: dh2o_interp  ! output
      real*8 w1,w2
      integer :: j1,j2

C**** Interpolate (extrapolate) horizontally
      j2 = 2+(jma-1)*(rlat-xlat(1))/(xlat(jma)-xlat(1)) ! first guess
      j2 = min(max(2,j2),jma)
      j1 = j2-1
      if(rlat.gt.xlat(j2)) then ! j guess was too low
         do while (j2.lt.jma .and. rlat.gt.xlat(j2))
          j2 = j2+1
         end do
         j1 = j2-1
      elseif(rlat.lt.xlat(j1)) then ! j guess was too high
         do while (j1.gt.1 .and. rlat.lt.xlat(j1))
          j1 = j1-1
         end do
         j2 = j1+1
      endif
C**** coeff. for latitudinal linear inter/extrapolation
      w1 = (xlat(j2)-rlat)/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
      if(w1.gt.1.) w1=.5+.5*w1
      if(w1.lt.0.) w1=.5*w1
      w2 = 1.-w1
      dh2o_interp = w1*dh2o(j1,lev,mon)+w2*dh2o(j2,lev,mon)
      return
      end subroutine lat_interp_qma
