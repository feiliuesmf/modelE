#include "rundeck_opts.h"

      subroutine diag_trac_prep
      implicit none
      call gather_zonal_trdiag
      call diagjlt_prep
      call diagijt_prep
      call diagijlt_prep
      call diagtcp_prep
      return
      end subroutine diag_trac_prep

      SUBROUTINE DIAGJLT_prep
! comments to be added
      use constant, only : teeny,grav
      USE MODEL_COM, only: lm,idacc,fim
      USE TRACER_COM, only : ntm,ntm_power,trname
     &     ,n_Water,n_CH4,n_O3
#ifdef TRACERS_WATER
     &     ,trw0,dowetdep
#endif
      USE DIAG_COM, only: jm=>jm_budg,dxyp=>dxyp_budg,
     &     ia_dga,ia_src,ajl
     &     ,jl_dpa,jl_dpasrc,jl_dwasrc
      USE TRDIAG_COM, only : tajln, tajls, lname_jln, sname_jln,
     *     units_jln,  scale_jln, lname_jls, sname_jls, units_jls,
     *     scale_jls, jls_power, jls_ltop, ia_jls, jwt_jls, jgrid_jls,
     *     jls_3Dsource, jlnt_conc, jlnt_mass, jlnt_nt_tot, jlnt_nt_mm,
     *     jlnt_lscond,  jlnt_turb,  jlnt_vt_tot, jlnt_vt_mm, jlnt_mc,
     *     jgrid_jlq, ia_jlq, scale_jlq, jlq_power, ktajls, jls_source
#ifdef TRACERS_WATER
     *     ,jlnt_cldh2o
#endif
     &     ,tajl=>tajl_out, ktajl_, ktajl_out, cdl_tajl
     &     ,denom_tajl,ia_tajl,sname_tajl,lname_tajl,ltop_tajl
     &     ,units_tajl,scale_tajl,pow_tajl,jgrid_tajl,lgrid_tajl
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      USE TRDIAG_COM, only : to_per_mil
#endif
#if !defined(CUBED_SPHERE) && !defined(CUBE_GRID)
      USE BDJLT
#endif
      use domain_decomp_atm, only : am_i_root
      IMPLICIT NONE
      real*8 :: bydxyp(jm),byapo(jm),onespo(jm)
      INTEGER :: J,L,N,K,KK,KKK,jtpow,n1,n2,k_dpa,k_dwa,k_vap,k_cnd
      REAL*8 :: dD, d18O, d17O, byiacc
      real*8, dimension(:,:,:), allocatable :: tajl_tmp
      character(len=10) :: zstr

      if(.not. am_i_root()) return

#if !defined(CUBED_SPHERE) && !defined(CUBE_GRID)
      call JLt_TITLEX ! needed for some extra titles
#endif

      onespo = 1.
      onespo(1) = fim
      onespo(jm) = fim
      bydxyp = 1d0/dxyp
      byapo = bydxyp*onespo/fim

      do k=1,ktajl_
        denom_tajl(k) = 0
        ia_tajl(k) = ia_src
        sname_tajl(k) = 'unused'
        lname_tajl(k) = 'unused'
        units_tajl(k) = 'unused'
        scale_tajl(k) = 1.
        pow_tajl(k) = 0
        jgrid_tajl(k) = 1
        lgrid_tajl(k) = 1
        ltop_tajl(k) = lm
      enddo

      k = 0

c
      k = k + 1
      k_dpa = k
      tajl(:,:,k) = ajl(:,:,jl_dpasrc)
c
      k = k + 1
      k_dwa = k
      tajl(:,:,k) = ajl(:,:,jl_dwasrc)

#ifdef TRACERS_WATER
c
      k = k + 1
      k_vap = k
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,jlnt_mass,n_water)*bydxyp(:)
      enddo
c
      k = k + 1
      k_cnd = k
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,jlnt_cldh2o,n_water)*bydxyp(:)
      enddo

#endif


C****
C**** LOOP OVER TRACERS
C****

      DO N=1,NTM
C****
C**** TRACER CONCENTRATION
C****
      k = k + 1
      kk = jlnt_conc
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,jlnt_mass,n)*bydxyp(:)
      enddo

#ifdef TRACERS_WATER
      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
        denom_tajl(k) = k_vap
        tajl(:,:,k) = 1d3*(tajl(:,:,k)/trw0(n)-tajl(:,:,k_vap))
      else
#endif

      denom_tajl(k) = k_dpa
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jln(n)*scale_jlq(kk)*10.**(-jtpow)

#ifdef TRACERS_WATER
      end if
#endif

#ifdef TRACERS_COSMO
      if(n.eq.n_Be7) k_Be7 = k
      if(n.eq.n_Pb210) k_Pb210 = k
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_OM_SP)
C****
C**** Mass diagnostic (this is saved for everyone, but only output
C**** for Dorothy and Drew for the time being)
C****
      k = k + 1
      kk=jlnt_mass

      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      tajl(:,:,k) = tajln(:,:,kk,n)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP)
      jtpow = jtpow+13
#else
      denom_tajl(k) = k_dpa
      byiacc = 1d0/(idacc(ia_tajl(k))+teeny)
      do l=1,lm
        tajl(:,l,k) = tajl(:,l,k)*byapo(:)*byiacc*tajl(:,l,k_dpa)
      enddo
#endif
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)

#endif

#ifdef TRACERS_WATER
C****
C**** TRACER CLOUD WATER CONCENTRATION
C****
      if (dowetdep(n)) then
      k = k + 1
      kk = jlnt_cldh2o

      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,kk,n)*bydxyp(:)
      enddo

      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
        denom_tajl(k) = k_cnd
        tajl(:,:,k) = 1d3*(tajl(:,:,k)/trw0(n)-tajl(:,:,k_cnd))
      else
        denom_tajl(k) = k_dwa
        jtpow = ntm_power(n)+jlq_power(kk)
        scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      end if
      end if
#endif
C****
C**** NORTHWARD TRANSPORTS: Total and eddies
C****
#if !defined(CUBED_SPHERE) && !defined(CUBE_GRID)
      k = k + 1
      kk = jlnt_nt_tot
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      jgrid_tajl(k) = 2
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(2:jm,:,k) = tajln(1:jm-1,:,kk,n)
c
      k = k + 1
      kk = jlnt_nt_eddy
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      jgrid_tajl(k) = 2
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(2:jm,:,k) = tajln(1:jm-1,:,jlnt_nt_tot,n)
     &                -tajln(1:jm-1,:,jlnt_nt_mm ,n)
#endif
C****
C**** VERTICAL TRANSPORTS: Total and eddies
C****
#if !defined(CUBED_SPHERE) && !defined(CUBE_GRID)
      k = k + 1
      kk = jlnt_vt_tot
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      lgrid_tajl(k) = 2
      ltop_tajl(k) = lm-1
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(:,:,k) = tajln(:,:,kk,n)
c
      k = k + 1
      kk = jlnt_vt_eddy
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      lgrid_tajl(k) = 2
      ltop_tajl(k) = lm-1
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(:,:,k) = tajln(:,:,jlnt_vt_tot,n)-tajln(:,:,jlnt_vt_mm,n)
#endif
c
c tendencies from various processes
c
      do kkk=1,3
        k = k + 1
        if(kkk.eq.1) kk = jlnt_mc
        if(kkk.eq.2) kk = jlnt_lscond
        if(kkk.eq.3) kk = jlnt_turb
        sname_tajl(k) = sname_jln(kk,n)
        lname_tajl(k) = lname_jln(kk,n)
        units_tajl(k) = units_jln(kk,n)
        ia_tajl(k) = ia_jlq(kk)
        jtpow = ntm_power(n)+jlq_power(kk)
        scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
        do l=1,lm
          tajl(:,l,k) = tajln(:,l,kk,n)*onespo(:)
        enddo
      enddo

      enddo ! end loop over tracers

C****
C**** JL Specials (incl. Sources and sinks)
C**** Partial move towards correct units (kg/(mb m^2 s)).
C**** Plot depends on jwt_jls.
C**** Note that only jwt_jls=3 is resolution independent.
C****
      do kk=1,ktajls
        if (sname_jls(kk).eq."daylight" .or. sname_jls(kk).eq."H2O_mr"
     *       .or. lname_jls(kk).eq."unused") cycle

        k = k+1
        sname_tajl(k) = sname_jls(kk)
        lname_tajl(k) = lname_jls(kk)
        units_tajl(k) = units_jls(kk)
        ia_tajl(k) = ia_jls(kk)
        ltop_tajl(k) = jls_ltop(kk)
        jtpow = jls_power(kk)
        scale_tajl(k) = scale_jls(kk)*10.**(-jtpow)
        tajl(:,:,k) = tajls(:,:,kk)
        byiacc = 1d0/(idacc(ia_tajl(k))+teeny)
        if(jwt_jls(kk).eq.1) then
          do l=1,lm
            tajl(:,l,k) = tajl(:,l,k)*onespo(:)
          enddo
        elseif(jwt_jls(kk).eq.2) then
          denom_tajl(k) = k_dpa
          do l=1,lm
            tajl(:,l,k) = tajl(:,l,k)*bydxyp(:)*byiacc*tajl(:,l,k_dpa)
          enddo
        elseif(jwt_jls(kk).eq.3) then
          denom_tajl(k) = k_dpa
          do l=1,lm
            tajl(:,l,k) = tajl(:,l,k)*bydxyp(:)*100./grav
          enddo
        endif

c        select case (jwt_jls(kk))
c        case (1)   !  simple sum (like kg/s),
c          CALL JLMAP_t (lname_jls(kk),sname_jls(kk),units_jls(kk),plm,
c     *         tajls(1,1,kk),scalet,onespo,ones,jls_ltop(kk),jwt_jls(kk)
c     *         ,jgrid_jls(kk))
c        case (2)   !  area weighting (like kg/m^2 s)
c          CALL JLMAP_t (lname_jls(kk),sname_jls(kk),units_jls(kk),plm
c     *         ,tajls(1,1,kk),scalet,bydxyp,ones,jls_ltop(kk),jwt_jls(kk)
c     *         ,jgrid_jls(kk))
c        case (3)   !  area + pressure weighting (like kg/mb m^2 s)
c          CALL JLMAP_t (lname_jls(kk),sname_jls(kk),units_jls(kk),plm
c     *         ,tajls(1,1,kk),scalet,byapo,ones,jls_ltop(kk),jwt_jls(kk)
c     *         ,jgrid_jls(kk))
c        end select

        end do

#ifdef TRACERS_SPECIAL_Lerner
C**** some special combination diagnostics

C**** total chemical change for CH4
      if (n_CH4.gt.0) then
        k = k + 1
        kk=jls_3Dsource(1,n_CH4)
        sname_tajl(k) = 'Total_Chem_change'//trname(n_CH4)
        lname_tajl(k) = 'TOTAL CHANGE OF CH4 BY CHEMISTRY'
        units_tajl(k) = units_jls(kk)
        ia_tajl(k) = ia_jls(kk)
        ltop_tajl(k) = jls_ltop(kk)
        scale_tajl(k) = scale_jls(kk)*10.**(-jls_power(kk))
        tajl(:,:,k) = tajls(:,:,jls_3Dsource(1,n_CH4))
     &              + tajls(:,:,jls_3Dsource(2,n_CH4))
        do l=1,lm
          tajl(:,l,k) = tajl(:,l,k)*onespo(:)
        enddo
      end if
C**** total chemical change for O3
      if (n_O3.gt.0) then
        k = k + 1
        kk=jls_3Dsource(1,n_O3)
        sname_tajl(k) = 'Total_change_chem+depo'//trname(n_O3)
        lname_tajl(k) = 'Total Change of O3 by Chemistry and deposition'
        units_tajl(k) = units_jls(kk)
        ia_tajl(k) = ia_jls(kk)
        ltop_tajl(k) = jls_ltop(kk)
        scale_tajl(k) = scale_jls(kk)*10.**(-jls_power(kk))
        tajl(:,:,k) = tajls(:,:,jls_3Dsource(1,n_O3))
     &              + tajls(:,:,jls_3Dsource(2,n_O3))
     &              + tajls(:,:,jls_3Dsource(3,n_O3))
     &              + tajls(:,:,jls_source(1,n_O3))
        do l=1,lm
          tajl(:,l,k) = tajl(:,l,k)*onespo(:)
        enddo
      end if
#endif

#ifdef TRACERS_COSMO
C**** ratios : Be7/Pb210 and Be10/Be7
      if (n_Be7.gt.0 .and. n_Be10.gt.0 .and. n_Pb210.gt.0) then
C*** ratio Be10/Be7
        k = k + 1
        lname_tajl(k)="Be10 to Be7 ratio"
        sname_tajl(k)="be10be7"
        units_tajl(k)=" "
        denom_tajl(k)=k_Be7
        tajl(:,:,k) = tajln(:,:,jlnt_mass,n_Be10)

C*** ratio Be7/Pb210
        k = k + 1

        lname_tajl(k)="Be7 to Pb210 ratio"
        sname_tajl(k)="be7pb210"
        units_tajl(k)="mBq/mBq"        !be sure this is 1/scalet
        denom_tajl(k)=k_Pb210
C*** scale by (Be7decay/mm_Be7)/(Pb210decay/mm_Pb210) to convert to mBq
        scale_tajl(k)=
     &       trdecay(n_Be7)*tr_mm(n_Pb210)/trdecay(n_Pb210)/tr_mm(n_Be7)
        tajl(:,:,k) = tajln(:,:,jlnt_mass,n_Be7)

      end if

#endif

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-8*d18O)
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** Concentration in water vapour
        k = k + 1
        kk=jlnt_mass
        lname_tajl(k)="Deuterium excess"
        sname_tajl(k)="dexcess"
        units_tajl(k)="per mil"
        denom_tajl(k)=k_vap
        do l=1,lm
          do j=1,jm
            d18O=tajln(j,l,kk,n1)/trw0(n1)-tajln(j,l,kk,n_water)
            dD=tajln(j,l,kk,n2)/trw0(n2)-tajln(j,l,kk,n_water)
            tajl(j,l,k)=1d3*(dD-8.*d18O)*bydxyp(j)
          enddo
        enddo

C**** Concentration in cloud water
        k = k + 1
        kk=jlnt_cldh2o
        lname_tajl(k)="Deuterium excess in cloud water"
        sname_tajl(k)="dexcess_cldh2o"
        units_tajl(k)="per mil"
        denom_tajl(k)=k_cnd
        do l=1,lm
          do j=1,jm
            d18O=tajln(j,l,kk,n1)/trw0(n1)-tajln(j,l,kk,n_water)
            dD=tajln(j,l,kk,n2)/trw0(n2)-tajln(j,l,kk,n_water)
            tajl(j,l,k)=1d3*(dD-8.*d18O)*bydxyp(j)
          enddo
        enddo

      end if

C****
C**** Calculations of 17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17

C**** Concentration in water vapour
        k = k + 1
        kk=jlnt_mass

        lname_tajl(k)="D17O excess"
        sname_tajl(k)="D17O_excess"
        units_tajl(k)="per meg"
        denom_tajl(k)=k_vap

        do l=1,lm
          do j=1,jm
            if (tajln(j,l,kk,n_water).gt.0) then
              d18O=tajln(j,l,kk,n1)/(trw0(n1)*tajln(j,l,kk,n_water))
              d17O=tajln(j,l,kk,n2)/(trw0(n2)*tajln(j,l,kk,n_water))
              tajl(j,l,k)=1d6*(log(d17O)-0.529d0*log(d18O))*
     &             tajl(j,l,k_vap)
            else
              tajl(j,l,k)=0.
            end if
          end do
        end do

C**** Concentration in cloud water
        k = k + 1
        kk=jlnt_cldh2o
        lname_tajl(k)="D17O excess in cloud water"
        sname_tajl(k)="D17O_excess_cldh2o"
        units_tajl(k)="per meg"
        denom_tajl(k)=k_cnd

        do l=1,lm
          do j=1,jm
            if (tajln(j,l,kk,n_water).gt.0) then
              d18O=tajln(j,l,kk,n1)/(trw0(n1)*tajln(j,l,kk,n_water))
              d17O=tajln(j,l,kk,n2)/(trw0(n2)*tajln(j,l,kk,n_water))
              tajl(j,l,k)=1d6*(log(d17O)-0.529d0*log(d18O))*
     &             tajl(j,l,k_cnd)
            else
              tajl(j,l,k)=0.
            end if
          end do
        end do

      end if
#endif

      ktajl_out = k

c
c if necessary, reallocate tajl to be the right size
c
      if(ktajl_out.lt.size(tajl,3)) then
        allocate(tajl_tmp(jm,lm,ktajl_out))
        tajl_tmp(:,:,1:ktajl_out) = tajl(:,:,1:ktajl_out)
        deallocate(tajl)
        allocate(tajl(jm,lm,ktajl_out))
        tajl = tajl_tmp
        deallocate(tajl_tmp)
      elseif(ktajl_out.gt.size(tajl,3)) then
        call stop_model('error: ktajl_out > size(tajl,3)',255)
      endif

#ifdef TRACERS_SPECIAL_Lerner
c adjust sname elements beginning with 14C
      do k=1,ktajl_out
        if(sname_tajl(k)(1:3).eq.'14C') then
          sname_tajl(k) = 'a'//trim(sname_tajl(k))
        endif
      enddo
#endif

c
c Declare the dimensions and metadata of TAJL output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      if(.not.allocated(cdl_tajl)) then
      allocate(cdl_tajl(6*ktajl_out))
      cdl_tajl = ''
      cdl_tajl(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_tajl(3),'(a,i3,a)') '   lat_budg = ',jm,' ;'
      write(cdl_tajl(4),'(a,i3,a)') '   lat_budg_plus3 = ',jm+3,' ;'
      write(cdl_tajl(5),'(a,i3,a)') '   plm = ',lm,' ;'
      write(cdl_tajl(6),'(a,i3,a)') '   ple = ',lm,' ;'
      write(cdl_tajl(7),'(a,i3,a)') '   shnhgm = 3 ;'
      cdl_tajl(8:14)(:) = (/
     &     'variables:                           ',
     &     'float lat_budg(lat_budg) ;           ',
     &     '   lat_budg:units = "degrees_north" ;',
     &     'float plm(plm) ;                     ',
     &     '   plm:units = "mb" ;                ',
     &     'float ple(ple) ;                     ',
     &     '   ple:units = "mb" ;                '
     &     /)
      kk = count(len_trim(cdl_tajl).gt.0)
      do k=1,ktajl_out
        if(trim(units_tajl(k)).eq.'unused') cycle
        if(lgrid_tajl(k).eq.1) then
          zstr='plm'
        else
          zstr='ple'
        endif
        kk = kk + 1
        cdl_tajl(kk) = 'float '//trim(sname_tajl(k))//'('//
     &       trim(zstr)//',lat_budg) ;'
        kk = kk + 1
        cdl_tajl(kk) = '   '//trim(sname_tajl(k))//':long_name = "'//
     &       trim(lname_tajl(k))//'" ;'
        kk = kk + 1
        cdl_tajl(kk) = '   '//trim(sname_tajl(k))//':units = "'//
     &       trim(units_tajl(k))//'" ;'
        if(pow_tajl(k).ne.0) then
          kk = kk + 1
          write(cdl_tajl(kk),'(a,i3,a)')
     &         '   '//trim(sname_tajl(k))//':prtpow = ',pow_tajl(k),' ;'
        endif
        kk = kk + 1
        cdl_tajl(kk) = 'float '//trim(sname_tajl(k))//'_hemis('//
     &       trim(zstr)//',shnhgm) ;'
        if(denom_tajl(k).gt.0) then
          kk = kk + 1
          cdl_tajl(kk) = 'float '//trim(sname_tajl(k))//
     &         '_vmean(lat_budg_plus3) ;'
        endif
        if(ltop_tajl(k).ne.lm) then
          kk = kk + 1
          write(cdl_tajl(kk),'(a,i3,a)')
     &         trim(sname_tajl(k))//':ltop = ',ltop_tajl(k),' ;'
        endif
      enddo
      kk = kk + 1
      cdl_tajl(kk) = '}'
      endif

      RETURN
      END SUBROUTINE DIAGJLT_prep

      subroutine diagijt_prep
! comments to be added
      use model_com, only: im,jm,lm,idacc,focean
      use tracer_com
      use diag_com
      use trdiag_com, only : taijn=>taijn_loc, taijs=>taijs_loc,
     &     ktaij_,ktaij_out,taij=>taij_out,scale_taij,cdl_taij,
     &     ir_taij,ia_taij,denom_taij,lname_taij,sname_taij,units_taij,
     &     sname_tij, lname_tij,
     &     units_tij, scale_tij, tij_mass, lname_ijts,  sname_ijts,
     &     units_ijts,  scale_ijts,  ia_ijts, ktaij, ktaijs, 
     &     tij_drydep, tij_gsdep, tij_surf, tij_grnd, tij_prec, 
     &     tij_uflx, tij_vflx, tij_kw, ijs_NO2_1330c, ijs_NO2_1030c, 
     &     ijs_NO2_1330, ijs_NO2_1030, tij_alpha
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
     &     ,to_per_mil
#endif
      use constant, only : teeny
      use domain_decomp_atm, only : grid,am_i_root
      use geom, only : byaxyp
      implicit none
      integer ::  i,j,k,kk,kx,n,n1,n2
      integer :: k_water(ktaij),k_Be7,k_Pb210,k_clr,
     & k_no2_1030=0,k_no2_1330=0,k_no2_1030c=0,k_no2_1330c=0
      logical :: div_by_area
      integer :: i_0,i_1,j_0,j_1, i_0h,i_1h,j_0h,j_1h
      real*8, dimension(:,:,:), allocatable :: taij_tmp
      character(len=16) :: ijstr

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

C**** Fill in the undefined pole box duplicates
      if(grid%have_south_pole) then
        do i=2,im
          taijn(i,1,:,:) = taijn(1,1,:,:)
          taijs(i,1,:) = taijs(1,1,:)
        enddo
      endif
      if(grid%have_north_pole) then
        do i=2,im
          taijn(i,jm,:,:) = taijn(1,jm,:,:)
          taijs(i,jm,:) = taijs(1,jm,:)
        enddo
      endif

      do k=1,ktaij_
        denom_taij(k) = 0
        ia_taij(k) = ia_src
        sname_taij(k) = 'unused'
        lname_taij(k) = 'unused'
        units_taij(k) = 'unused'
        scale_taij(k) = 1.
      enddo

      k = 0

      k = k + 1
      k_clr = k
      do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k)=real(idacc(ia_rad))-aij_loc(i,j,ij_cldcv)
        enddo
      enddo

c
c Tracer sums/means and ground conc
c
      do n=1,ntm
      do kx=1,ktaij
        if (index(lname_tij(kx,n),'unused').gt.0) cycle
        k = k+1

        taij(i_0:i_1,j_0:j_1,k) = taijn(i_0:i_1,j_0:j_1,kx,n)

        sname_taij(k) = sname_tij(kx,n)
        lname_taij(k) = lname_tij(kx,n)
        units_taij(k) = units_tij(kx,n)
        scale_taij(k) = scale_tij(kx,n)
        ir_taij(k) = ir_log2

#ifdef TRACERS_COSMO
        if(kx.eq.tij_surf .and. n.eq.n_Be7) k_Be7 = k
        if(kx.eq.tij_surf .and. n.eq.n_Pb210) k_Pb210 = k
#endif

#ifdef TRACERS_WATER
        if(n.eq.n_water) then
          k_water(kx) = k  ! water must go before any per_mil tracers
        endif
        if (to_per_mil(n).gt.0) then
          if((kx-tij_mass)*(kx-tij_uflx)*(kx-tij_vflx).ne.0) then
            do j=j_0,j_1; do i=i_0,i_1
              taij(i,j,k) =
     &             1d3*(taij(i,j,k)/trw0(n)-taijn(i,j,kx,n_water))
            enddo       ; enddo
            denom_taij(k) = k_water(kx)
            if(n.eq.n_water) then ! save a non-per-mil denom
              k = k+1
              k_water(kx) = k
              denom_taij(k-1) = k_water(kx)
              do j=j_0,j_1; do i=i_0,i_1
                taij(i,j,k) =  taijn(i,j,kx,n)
              enddo       ; enddo
              ia_taij(k) = ia_taij(k-1)
            endif
          endif
        endif
#endif

      enddo ! end loop over k

#if (defined TRACERS_WATER) && (defined TRACERS_DRYDEP)
c
c dry/(dry+wet) deposition fraction
c
      if (dodrydep(n).and.dowetdep(n)) then
        k=k+1
        sname_taij(k) = "pc_dry_dep_"//trim(trname(n))
        lname_taij(k) = trim(trname(n))//" Percent Dry Deposition"
        units_taij(k) = "%"
        ir_taij(k) = ir_pct
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = taijn(i,j,tij_drydep,n)+taijn(i,j,tij_gsdep,n)
          taij(i,j,k+1) = taij(i,j,k) + taijn(i,j,tij_prec,n)
c ijt_mapk does not apply the scale factor to ratios.  workaround.
          taij(i,j,k) = taij(i,j,k)*100.
        enddo
        enddo
        denom_taij(k) = k+1
        k=k+1
      endif
#endif

      end do ! end loop over tracers

#ifdef TRACERS_COSMO
      if (n_Be7.gt.0 .and. n_Be10.gt.0 .and. n_Pb210.gt.0) then
C**** Be10/Be7
        k=k+1
        sname_taij(k) = "be10be7_ij"
        lname_taij(k) = "surface ratio Be10 to Be7"
        units_taij(k) = " "
        ir_taij(k) = ir_0_18
        ia_taij(k) = ia_srf
        denom_taij(k) = k_Be7
        taij(i_0:i_1,j_0:j_1,k) = taijn(i_0:i_1,j_0:j_1,tij_surf,n_Be10)
C**** Be7/Pb210
        k=k+1
        sname_taij(k) = "be7pb210_ij"
        lname_taij(k) = "surface ratio Be7 to Pb210"
        units_taij(k) = "mBq/mBq "
        ir_taij(k) = ir_0_180
        ia_taij(k) = ia_srf
        denom_taij(k) = k_Pb210
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = taij(i,j,tij_surf,n_Be7)
C*** scale by (Be7decay/mm_Be7)/(Pb210decay/mm_Pb210) to convert to mBq
          taij(i,j,k)=taij(i,j,k)*trdecay(n_Be7)*tr_mm(n_Pb210)
     &             /trdecay(n_Pb210)/tr_mm(n_Be7)
        enddo
        enddo
      endif
#endif

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-d18O)
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** precipitation
        k=k+1
        denom_taij(k) = k_water(tij_prec)
        sname_taij(k) = "prec_ij_dex"
        lname_taij(k) = "Deuterium excess in precip"
        units_taij(k) = "per mil"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = 1d3*(taijn(i,j,tij_prec,n2)/trw0(n2)-
     &                  8.*taijn(i,j,tij_prec,n1)/trw0(n1)+
     &                  7.*taijn(i,j,tij_prec,n_water))
        enddo
        enddo
C**** ground concentration
        k=k+1
        denom_taij(k) = k_water(tij_grnd)
        sname_taij(k) = "grnd_ij_dex"
        lname_taij(k) = "Deuterium excess at Ground"
        units_taij(k) = "per mil"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = 1d3*(taijn(i,j,tij_grnd,n2)/trw0(n2)-
     &                  8.*taijn(i,j,tij_grnd,n1)/trw0(n1)+
     &                  7.*taijn(i,j,tij_grnd,n_water))
        enddo
        enddo
      end if

C****
C**** Calculations of D17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
C**** precipitation
        k=k+1
        denom_taij(k) = k_water(tij_prec)
        sname_taij(k) = "prec_ij_D17O"
        lname_taij(k) = "D17O excess in precip"
        units_taij(k) = "per meg"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
          do i=i_0,i_1
            if (taijn(i,j,tij_prec,n_water).gt.0) then
              taij(i,j,k) = 1d6*taijn(i,j,tij_prec,n_water)*
     &             (log(taijn(i,j,tij_prec,n2)/trw0(n2))-
     &             0.529d0*log(taijn(i,j,tij_prec,n1)/trw0(n1))-
     &             0.471d0*log(taijn(i,j,tij_prec,n_water)))
            else
              taij(i,j,k)=0.
            end if
          end do
        end do
C**** ground concentration
        k=k+1
        denom_taij(k) = k_water(tij_grnd)
        sname_taij(k) = "grnd_ij_D17O"
        lname_taij(k) = "D17O excess at Ground"
        units_taij(k) = "per meg"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
          do i=i_0,i_1
            if (taijn(i,j,tij_grnd,n_water).gt.0) then
              taij(i,j,k) = 1d6*taijn(i,j,tij_grnd,n_water)*
     &             (log(taijn(i,j,tij_grnd,n2)/trw0(n2))-
     &             0.529d0*log(taijn(i,j,tij_grnd,n1)/trw0(n1))-
     &             0.471d0*log(taijn(i,j,tij_grnd,n_water)))
            else
              taij(i,j,k)=0.
            end if
          end do
        end do
      end if
#endif

c
c tracer sources and sinks
c
      do kx=1,ktaijs
        if (index(lname_ijts(kx),'unused').gt.0) cycle

        k = k+1

        sname_taij(k) = sname_ijts(kx)
        lname_taij(k) = lname_ijts(kx)
        units_taij(k) = units_ijts(kx)
        ir_taij(k) = ir_log2    ! should be the correct default
        ia_taij(k) = ia_ijts(kx)
        scale_taij(k) = scale_ijts(kx)

        taij(i_0:i_1,j_0:j_1,k) = taijs(i_0:i_1,j_0:j_1,kx)

        div_by_area = .true.

        if(sname_taij(k)(1:3).eq.'tau') div_by_area = .false.
        if(sname_taij(k)(1:3).eq.'swf') div_by_area = .false.
        if(sname_taij(k)(1:3).eq.'lwf') div_by_area = .false.
        if(sname_taij(k)(1:3).eq.'no_') div_by_area = .false.
        if(sname_taij(k)(1:5).eq.'wtrsh') div_by_area = .false.
        if(sname_taij(k)(1:8).eq.'ext_band') div_by_area = .false.
        if(sname_taij(k)(1:8).eq.'sct_band') div_by_area = .false.
        if(sname_taij(k)(1:8).eq.'asf_band') div_by_area = .false.
        if(sname_taij(k)(1:4).eq.'DIAM') div_by_area = .false.
        if(sname_taij(k)(1:8).eq.'DMS_con_') div_by_area = .false.
        if(sname_taij(k)(1:8).eq.'SO2_con_') div_by_area = .false.
        if(sname_taij(k)(1:8).eq.'SO4_con_') div_by_area = .false.
        select case (trim(sname_taij(k)))
        case('Ox_loss','Ox_prod','OH_vmr','OH_con','NO3_con','HO2_con',
     &  'J_H2O2','NO2_1030c','NO2_1330c','NO2_1030','NO2_1330',
     &  'COprod','COdest','Oxprod','Oxdest','CH4dest','OxpHO2',
     &  'OxpCH3O2','OxpRO2','OxlOH','OxlHO2','OxlALK','phO1d','pO1d',
     &  'pOH','NO_vmr','NO2_vmr')
          div_by_area = .false.
        end select

        if(div_by_area) then
          do j=j_0,j_1
          do i=i_0,i_1
            taij(i,j,k) = taij(i,j,k)*byaxyp(i,j)
          enddo
          enddo
        endif

        if(sname_taij(k)(5:6).eq.'CS') then
          do j=j_0,j_1
          do i=i_0,i_1
c ijt_mapk does not apply the scale factor to ratios.  workaround.
            taij(i,j,k) = taij(i,j,k)*scale_taij(k)
          enddo
          enddo
          scale_taij(k) = 1.
          denom_taij(k) = k_clr
        endif

        if(sname_taij(k)=='NO2_1030c') k_no2_1030c = k
        if(sname_taij(k)=='NO2_1030') k_no2_1030 = k
        if(sname_taij(k)=='NO2_1330c') k_no2_1330c = k
        if(sname_taij(k)=='NO2_1330') k_no2_1330 = k

      end do

      ktaij_out = k

c set this denominator by introducing denom_ijts to tracer code?
      if(k_no2_1030.gt.0) then
        k = k_no2_1030
        denom_taij(k) = k_no2_1030c
        do j=j_0,j_1
        do i=i_0,i_1
c ijt_mapk does not apply scale factor to ratios.  workaround.
          taij(i,j,k) = taij(i,j,k)*scale_taij(k)
        enddo
        enddo
        scale_taij(k) = 1.
      endif

      if(k_no2_1330.gt.0) then
        k = k_no2_1330
        denom_taij(k) = k_no2_1330c
        do j=j_0,j_1
        do i=i_0,i_1
c ijt_mapk does not apply scale factor to ratios.  workaround.
          taij(i,j,k) = taij(i,j,k)*scale_taij(k)
        enddo
        enddo
        scale_taij(k) = 1.
      endif

c
c if necessary, reallocate taij to be the right size
c
      if(ktaij_out.lt.size(taij,3)) then
        allocate(taij_tmp(i_0:i_1,j_0:j_1,ktaij_out))
        taij_tmp(i_0:i_1,j_0:j_1,1:ktaij_out) =
     &      taij(i_0:i_1,j_0:j_1,1:ktaij_out)
        deallocate(taij)
        allocate(taij(i_0h:i_1h,j_0h:j_1h,ktaij_out))
        taij(i_0:i_1,j_0:j_1,1:ktaij_out) =
     &       taij_tmp(i_0:i_1,j_0:j_1,1:ktaij_out)
        deallocate(taij_tmp)
      elseif(ktaij_out.gt.size(taij,3)) then
        call stop_model('error: ktaij_out > size(taij,3)',255)
      endif

#ifdef TRACERS_SPECIAL_Lerner
c adjust sname elements beginning with 14C
      do k=1,ktaij_out
        if(sname_taij(k)(1:3).eq.'14C') then
          sname_taij(k) = 'a'//trim(sname_taij(k))
        endif
      enddo
#endif

c
c Declare the dimensions and metadata of TAIJ output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      if(.not.allocated(cdl_taij)) then
      allocate(cdl_taij(6*ktaij_out))
      cdl_taij = ''
      cdl_taij(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      cdl_taij(3) = '   shnhgm = 3 ;'
#if defined(CUBED_SPHERE) || defined(CUBE_GRID)
      ijstr='(tile,y,x) ;'
      write(cdl_taij(4),'(a,i3,a)') '   x = ',im,' ;'
      write(cdl_taij(5),'(a,i3,a)') '   y = ',im,' ;'
      write(cdl_taij(6),'(a,i3,a)') '   tile = 6 ;'
      do k=4,6
        cdl_taij(k)=trim(cdl_taij(k))//' // remove_from_latlon'
      enddo
      cdl_taij(7) = '// add_to_latlon    lon = xxx ;'
      cdl_taij(8) = '// add_to_latlon    lat = xxx ;'
      cdl_taij(9) = 'variables:'
      cdl_taij(10:17)(:) = (/
     &     'float x(x) ;                                        ',
     &     '   x:long_name = "nondimensional cube coordinate" ; ',
     &     'float y(y) ;                                        ',
     &     '   y:long_name = "nondimensional cube coordinate" ; ',
     &     'float lon(tile,y,x) ;                               ',
     &     '   lon:units = "degrees_east" ;                     ',
     &     'float lat(tile,y,x) ;                               ',
     &     '   lat:units = "degrees_north" ;                    '
     &     /)
      do k=10,17
        cdl_taij(k)=trim(cdl_taij(k))//' // remove_from_latlon'
      enddo
      cdl_taij(18:21)(:) = (/
     &     '// add_to_latlon float lon(lon) ;               ',
     &     '// add_to_latlon   lon:units = "degrees_east" ; ',
     &     '// add_to_latlon float lat(lat) ;               ',
     &     '// add_to_latlon   lat:units = "degrees_north" ;'
     &     /)
#else
      ijstr='(lat,lon) ;'
      write(cdl_taij(4),'(a,i3,a)') '   lon = ',im,' ;'
      write(cdl_taij(5),'(a,i3,a)') '   lat = ',jm,' ;'
      cdl_taij(6:10)(:) = (/
     &     'variables:                      ',
     &     'float lon(lon) ;                ',
     &     '   lon:units = "degrees_east" ; ',
     &     'float lat(lat) ;                ',
     &     '   lat:units = "degrees_north" ;'
     &     /)
#endif
      kk = count(len_trim(cdl_taij).gt.0)
      do k=1,ktaij_out
        if(trim(sname_taij(k)).eq.'unused') cycle
        kk = kk + 1
        cdl_taij(kk) = 'float '//trim(sname_taij(k))//trim(ijstr)
#if defined(CUBED_SPHERE) || defined(CUBE_GRID)
        cdl_taij(kk)=trim(cdl_taij(kk))//' // remove_from_latlon'
        kk = kk + 1
        cdl_taij(kk) = '// add_to_latlon float '//
     &       trim(sname_taij(k))//'(lat,lon);'
#endif
        kk = kk + 1
        cdl_taij(kk) = '   '//trim(sname_taij(k))//':long_name = "'//
     &       trim(lname_taij(k))//'" ;'
        kk = kk + 1
        cdl_taij(kk) = '   '//trim(sname_taij(k))//':units = "'//
     &       trim(units_taij(k))//'" ;'
        kk = kk + 1
        cdl_taij(kk) = 'float '//trim(sname_taij(k))//
     &       '_hemis(shnhgm) ;'
      enddo
      kk = kk + 1
      cdl_taij(kk) = '}'
      endif

      return
      end subroutine diagijt_prep

      SUBROUTINE DIAGIJLt_prep
! comments to be added
      use model_com, only: im,jm,lm,idacc
      use tracer_com
      use diag_com
      use trdiag_com, only : taijln=>taijln_loc, taijls=>taijls_loc,
     &     ktaijl_,ktaijl_out,taijl=>taijl_out,scale_taijl,ir_taijl,
     &     ia_taijl,denom_taijl,lname_taijl,sname_taijl,units_taijl,
     &     cdl_taijl, sname_ijlt, lname_ijlt,
     &     units_ijlt, sname_ijt, lname_ijt, units_ijt, scale_ijt,
     &     ir_ijlt, ia_ijlt, scale_ijlt, ktaijl
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
     &     ,to_per_mil
#endif
      use constant, only : teeny
      use domain_decomp_atm, only : grid,am_i_root
      use geom, only : byaxyp
      implicit none
      integer i,j,l,k,kx,kk,n,n1,n2
      real*8 :: r1,r2
      integer :: k_water
      integer :: i_0,i_1,j_0,j_1, i_0h,i_1h,j_0h,j_1h
      real*8, dimension(:,:,:,:), allocatable :: taijl_tmp
      character(len=16) :: zstr,hstr,tstr

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

C**** Fill in the undefined pole box duplicates
      if(grid%have_south_pole) then
        do i=2,im
          taijln(i,1,:,:) = taijln(1,1,:,:)
          taijls(i,1,:,:) = taijls(1,1,:,:)
        enddo
      endif
      if(grid%have_north_pole) then
        do i=2,im
          taijln(i,jm,:,:) = taijln(1,jm,:,:)
          taijls(i,jm,:,:) = taijls(1,jm,:,:)
        enddo
      endif

      do k=1,ktaijl_
        denom_taijl(k) = 0
        ia_taijl(k) = ia_src
        sname_taijl(k) = 'unused'
        lname_taijl(k) = 'unused'
        units_taijl(k) = 'unused'
        scale_taijl(k) = 1.
      enddo

      k = 0

C**** Tracer concentrations
      do n=1,ntm
        k = k+1
        sname_taijl(k) = sname_ijt(n)
        lname_taijl(k) = lname_ijt(n)
        units_taijl(k) = units_ijt(n)
        ir_taijl(k) = ir_log2
        ia_taijl(k) = ia_src
        scale_taijl(k) = scale_ijt(n)
        do l=1,lm
        do j=j_0,j_1; do i=i_0,i_1
          taijl(i,j,l,k) = taijln(i,j,l,n)*byaxyp(i,j)
        enddo       ; enddo
        enddo
#ifdef TRACERS_WATER
        if(n.eq.n_water) k_water = k
        if(to_per_mil(n).gt.0) then
          do l=1,lm
          do j=j_0,j_1; do i=i_0,i_1
            taijl(i,j,l,k) = 1d3*(taijl(i,j,l,k)/trw0(n)
     &           -byaxyp(i,j)*taijln(i,j,l,n_water))
          enddo       ; enddo
          enddo
          denom_taijl(k) = k_water
          if(n.eq.n_water) then ! save a non-per-mil denom
            k = k+1
            k_water = k
            denom_taijl(k-1) = k_water
            do l=1,lm
            do j=j_0,j_1; do i=i_0,i_1
              taijl(i,j,l,k) =  taijln(i,j,l,n)*byaxyp(i,j)
            enddo       ; enddo
            enddo
            ia_taijl(k) = ia_taijl(k-1)
          endif
        endif
#endif
      enddo

C**** Tracer specials 
      do kx=1,ktaijl
        if (index(lname_ijlt(kx),'unused').gt.0) cycle
        k = k+1
        sname_taijl(k) = sname_ijlt(kx)
        lname_taijl(k) = lname_ijlt(kx)
        units_taijl(k) = units_ijlt(kx)
        ir_taijl(k) = ir_ijlt(kx)
        ia_taijl(k) = ia_ijlt(kx)
        scale_taijl(k) = scale_ijlt(kx)
        do l=1,lm
        do j=j_0,j_1; do i=i_0,i_1
          taijl(i,j,l,k) = taijls(i,j,l,kx)
        enddo       ; enddo
        enddo
      enddo

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-d18O)
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** water vapour
        k=k+1
        denom_taijl(k) = k_water
        lname_taijl(k) = "Deuterium excess water vapour"
        sname_taijl(k) = "wvap_ij_dex"
        units_taijl(k) = "per mil"
        ir_taijl(k) = ir_m45_130
        ia_taijl(k) = ia_src
        do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          taijl(i,j,l,k) = 1d3*(taijln(i,j,l,n2)/trw0(n2)-
     &         8.*taijln(i,j,l,n1)/trw0(n1)+
     &         7.*taijln(i,j,l,n_water))*byaxyp(i,j)
        enddo
        enddo
        enddo
      end if

C****
C**** Calculations of D17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
C**** water vapour
        k=k+1
        denom_taijl(k) = k_water
        lname_taijl(k) = "D17O excess water vapour"
        sname_taijl(k) = "wvap_ij_D17O"
        units_taijl(k) = "per meg"
        ir_taijl(k) = ir_m45_130
        ia_taijl(k) = ia_src
        do l=1,lm
          do j=j_0,j_1
            do i=i_0,i_1
              if (taijln(i,j,l,n_water).gt.0) then
                r1 = taijln(i,j,l,n1)/(trw0(n1)*taijln(i,j,l,n_water))
                r2 = taijln(i,j,l,n2)/(trw0(n2)*taijln(i,j,l,n_water))
                taijl(i,j,l,k) = 1d6*taijl(i,j,l,k_water)*
     &               (log(r2)-0.529d0*log(r1))
              else
                taijl(i,j,l,k)=0.
              end if
            end do
          end do             
        end do
      end if
#endif

      ktaijl_out = k

c
c if necessary, reallocate taijl to be the right size
c
      if(ktaijl_out.lt.size(taijl,4)) then
        allocate(taijl_tmp(i_0:i_1,j_0:j_1,lm,ktaijl_out))
        taijl_tmp(i_0:i_1,j_0:j_1,:,1:ktaijl_out) =
     &      taijl(i_0:i_1,j_0:j_1,:,1:ktaijl_out)
        deallocate(taijl)
        allocate(taijl(i_0h:i_1h,j_0h:j_1h,lm,ktaijl_out))
        taijl(i_0:i_1,j_0:j_1,:,1:ktaijl_out) =
     &       taijl_tmp(i_0:i_1,j_0:j_1,:,1:ktaijl_out)
        deallocate(taijl_tmp)
      elseif(ktaijl_out.gt.size(taijl,4)) then
        call stop_model('error: ktaijl_out > size(taijl,4)',255)
      endif

#ifdef TRACERS_SPECIAL_Lerner
c adjust sname elements beginning with 14C
      do k=1,ktaijl_out
        if(sname_taijl(k)(1:3).eq.'14C') then
          sname_taijl(k) = 'a'//trim(sname_taijl(k))
        endif
      enddo
#endif

c
c Declare the dimensions and metadata of TAIJL output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      if(.not.allocated(cdl_taijl)) then
      allocate(cdl_taijl(50+6*ktaijl_out))
      cdl_taijl = ''
      cdl_taijl(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_taijl(3),'(a,i3,a)') '   plm = ',lm,' ;'
#if defined(CUBED_SPHERE) || defined(CUBE_GRID)
      tstr='(tile,'
      hstr=',y,x) ;'
      write(cdl_taijl(4),'(a,i3,a)') '   x = ',im,' ;'
      write(cdl_taijl(5),'(a,i3,a)') '   y = ',im,' ;'
      write(cdl_taijl(6),'(a,i3,a)') '   tile = 6 ;'
      do k=4,6
        cdl_taijl(k)=trim(cdl_taijl(k))//' // remove_from_latlon'
      enddo
      cdl_taijl(7)  = '// add_to_latlon    lon = xxx ;'
      cdl_taijl(8)  = '// add_to_latlon    lat = xxx ;'
      cdl_taijl(9) = 'variables:'
      cdl_taijl(10:17)(:) = (/
     &     'float x(x) ;                                        ',
     &     '   x:long_name = "nondimensional cube coordinate" ; ',
     &     'float y(y) ;                                        ',
     &     '   y:long_name = "nondimensional cube coordinate" ; ',
     &     'float lon(tile,y,x) ;                               ',
     &     '   lon:units = "degrees_east" ;                     ',
     &     'float lat(tile,y,x) ;                               ',
     &     '   lat:units = "degrees_north" ;                    '
     &     /)
      do k=10,17
        cdl_taijl(k)=trim(cdl_taijl(k))//' // remove_from_latlon'
      enddo
      cdl_taijl(18:21)(:) = (/
     &     '// add_to_latlon float lon(lon) ;               ',
     &     '// add_to_latlon   lon:units = "degrees_east" ; ',
     &     '// add_to_latlon float lat(lat) ;               ',
     &     '// add_to_latlon   lat:units = "degrees_north" ;'
     &     /)
#else
      tstr='('
      hstr=',lat,lon) ;'
      write(cdl_taijl(4),'(a,i3,a)') '   lon = ',im,' ;'
      write(cdl_taijl(5),'(a,i3,a)') '   lat = ',jm,' ;'
      cdl_taijl(6:10)(:) = (/
     &     'variables:                      ',
     &     'float lon(lon) ;                ',
     &     '   lon:units = "degrees_east" ; ',
     &     'float lat(lat) ;                ',
     &     '   lat:units = "degrees_north" ;'
     &     /)
#endif
      kk = 1+count(len_trim(cdl_taijl).gt.0)
      cdl_taijl(kk) = '// vertical_coords: plm'
      kk = kk + 1
      cdl_taijl(kk:kk+1)(:) = (/
     &     'float plm(plm) ;                ',
     &     '   plm:units = "mb" ;           '
     &     /)
      kk = count(len_trim(cdl_taijl).gt.0)
      do k=1,ktaijl_out
        if(trim(sname_taijl(k)).eq.'unused') cycle
        zstr='plm'
        kk = kk + 1
        cdl_taijl(kk) = 'float '//trim(sname_taijl(k))//
     &       trim(tstr)//trim(zstr)//trim(hstr)
#if defined(CUBED_SPHERE) || defined(CUBE_GRID)
        cdl_taijl(kk)=trim(cdl_taijl(kk))//' // remove_from_latlon'
        kk = kk + 1
        cdl_taijl(kk) = '// add_to_latlon float '//
     &       trim(sname_taijl(k))//'('//trim(zstr)//',lat,lon);'
#endif
        kk = kk + 1
        cdl_taijl(kk) = '   '//trim(sname_taijl(k))//':long_name = "'//
     &       trim(lname_taijl(k))//'" ;'
        kk = kk + 1
        cdl_taijl(kk) = '   '//trim(sname_taijl(k))//':units = "'//
     &       trim(units_taijl(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_taijl(kk) = '}'
      endif

      return
      end subroutine diagijlt_prep

      SUBROUTINE DIAGTCP_prep
! comments to be added
      USE CONSTANT, only: teeny
      USE MODEL_COM, only : fim,idacc
      USE GEOM, only: areag
      USE DIAG_COM, only: jm=>jm_budg,dxyp=>dxyp_budg,lat_dg=>lat_budg
      USE TRACER_COM, only: ntm
      USE TRDIAG_COM, only : tconsrv,
     &     ktcon,ktcon_out,ntmxcon,nsum_tcon,scale_tcon,
     &     ia_tcon,title_tcon,hemis_tconsrv,name_tconsrv,tconsrv_out,
     &     ia_tcon_out,scale_tcon_out,sname_tconsrv_out,cdl_tconsrv
      use domain_decomp_atm, only : am_i_root
      IMPLICIT NONE
      INTEGER :: j,n,k,kk,KTCON_max,j1,j2
      real*8 :: hemfac
      character*80 :: sname,titles(ktcon*ntmxcon)

      if(.not.am_i_root()) return

      if(.not.allocated(tconsrv_out)) then
c
c determine how many actual output fields there are
c
        ktcon_out = 0
        do n=1,ntm
          do k=ktcon,1,-1
            if(nsum_tcon(k,n).gt.0) ktcon_max = nsum_tcon(k,n)
          enddo
          do k=1,ktcon_max
            ktcon_out = ktcon_out + 1
          enddo
        enddo
c
c allocate space for actual number of outputs
c
        allocate(tconsrv_out(jm,ktcon_out))
        allocate(hemis_tconsrv(3,ktcon_out))
        allocate(ia_tcon_out(ktcon_out))
        allocate(scale_tcon_out(ktcon_out))
        allocate(sname_tconsrv_out(ktcon_out))
        allocate(cdl_tconsrv(5*ktcon_out))
c
c copy metadata
c
        kk = 0
        do n=1,ntm
          do k=ktcon,1,-1
            if(nsum_tcon(k,n).gt.0) ktcon_max = nsum_tcon(k,n)
          enddo
          do k=1,ktcon_max
            kk = kk + 1
            ia_tcon_out(kk) = ia_tcon(k,n)
            scale_tcon_out(kk) = scale_tcon(k,n)
            sname_tconsrv_out(kk) = name_tconsrv(k,n)
            titles(kk) = title_tcon(k,n)
          enddo
        enddo

c
c Declare the dimensions and metadata of TCONSRV output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
        cdl_tconsrv = ''
        cdl_tconsrv(1:2)(:) = (/
     &       'netcdf xxx { ', 'dimensions:  ' /)
        write(cdl_tconsrv(3),'(a,i3,a)') '   lat_budg = ',jm,' ;'
        write(cdl_tconsrv(4),'(a,i3,a)') '   shnhgm = 3 ;'
        cdl_tconsrv(5:9)(:) = (/
     &       'variables:                           ',
     &       'float lat_budg(lat_budg) ;           ',
     &       '   lat_budg:units = "degrees_north" ;',
     &       'float area_budg(lat_budg) ;          ',
     &       '   area_budg:units = "m^2" ;         '
     &       /)
        kk = count(len_trim(cdl_tconsrv).gt.0)
        do k=1,ktcon_out
          sname = trim(sname_tconsrv_out(k))
          kk = kk + 1
          cdl_tconsrv(kk) = 'float '//trim(sname)//'(lat_budg) ;'
          kk = kk + 1
          cdl_tconsrv(kk) = '   '//trim(sname)//':long_name = "'//
     &         trim(titles(k))//'" ;'
          kk = kk + 1
          cdl_tconsrv(kk) = 'float '//trim(sname)//'_hemis(shnhgm) ;'
        enddo
        kk = kk + 1
        cdl_tconsrv(kk) = '}'

      endif ! memory allocation and setup

c
c copy the nonzero contents of tconsrv into tconsrv_out
c also calculate sums of changes, hemispheric/global avgs
c
      hemfac = 2./(fim*sum(dxyp))
      kk = 0
      do n=1,ntm
        do k=ktcon,1,-1
C**** LOOP BACKWARDS SO THAT INITIALIZATION IS DONE BEFORE SUMMATION!
          if(nsum_tcon(k,n).eq.0) then
            tconsrv(:,k,n)=0.
          elseif(nsum_tcon(k,n).gt.0) then
            tconsrv(:,nsum_tcon(k,n),n)=
     &      tconsrv(:,nsum_tcon(k,n),n)+tconsrv(:,k,n)
     &           *scale_tcon(k,n)*idacc(12)/(idacc(ia_tcon(k,n))+teeny)
            ktcon_max = nsum_tcon(k,n)
          endif
        enddo
        do k=1,ktcon_max
          kk = kk + 1
          tconsrv_out(:,kk) = tconsrv(:,k,n)
          j1 = 1; j2 = jm/2
          hemis_tconsrv(1,kk) = hemfac*sum(tconsrv_out(j1:j2,kk))
          j1 = jm/2+1; j2 = jm
          hemis_tconsrv(2,kk) = hemfac*sum(tconsrv_out(j1:j2,kk))
          hemis_tconsrv(3,kk) = .5*sum(hemis_tconsrv(1:2,kk))
          tconsrv_out(:,kk) = tconsrv_out(:,kk)/(fim*dxyp)
        enddo
      enddo

      RETURN
      END SUBROUTINE DIAGTCP_prep
