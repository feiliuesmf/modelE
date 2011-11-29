#include "rundeck_opts.h"

! A WIP version of OCN_TRACER appropriate for RUNTIME_NTM_OCEAN and
! tracers whose ICs are defined in a netcdf file containing arrays
! whose names match those specified in ocean_trname in the rundeck.

      subroutine tracer_ic_ocean(atmocn)
      use model_com, only: itime,itimei
      use ocn_tracer_com, only : ntm, trname
      use ocean, only : im,jm,lmo
      use ocean, only : dxypo,mo
      use ocean, only : trmo
      use ocean, only : txmo,tymo,tzmo
      use ocean, only : nbyzm,i1yzm,i2yzm
      use domain_decomp_1d, only : get, am_i_root
      use oceanr_dim, only : grid=>ogrid
      use exchange_types, only : atmocn_xchng_vars
      use dictionary_mod
      use pario, only : read_dist_data,par_open,par_close
      implicit none
      type(atmocn_xchng_vars) :: atmocn
c
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     tr_ic
      integer n,i,j,l,nt,fid
      integer :: j_0s, j_1s, j_0, j_1

      if(itime.ne.itimei) return

      call get(grid, j_strt_skp = j_0s, j_stop_skp = j_1s,
     &     j_strt = j_0, j_stop = j_1)

      fid = par_open(grid,'OCN_TRACER_IC','read')

! Loop over tracers, read the IC for each, convert to extensive units (kg).
! straits IC not an option yet.
      do nt=1,ntm
        call read_dist_data(grid,fid,trim(trname(nt)),tr_ic)
        do l=1,lmo
        do j=j_0,j_1
        trmo(:,j,l,nt) = 0.
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          trmo(i,j,l,nt) = tr_ic(i,j,l)*mo(i,j,l)*dxypo(j)
        enddo
        enddo
        enddo
        enddo
      enddo

      call par_close(grid,fid)

      return
      end subroutine tracer_ic_ocean

      subroutine oc_tdecay(dts)
      implicit none
      real*8, intent(in) :: dts
      return
      end subroutine oc_tdecay

      subroutine diagtco (m,nt0,atmocn)
      use exchange_types, only : atmocn_xchng_vars
      implicit none
      integer, intent(in) :: m
      integer, intent(in) :: nt0
      type(atmocn_xchng_vars) :: atmocn
      return
      end subroutine diagtco
