#include "rundeck_opts.h"

      module icedyn_com_dummy
c dummy module just to test compilation
      implicit none
      real*8, dimension(:,:), allocatable :: rsix,rsiy,rsisave
      end module icedyn_com_dummy

      SUBROUTINE ADVSIcs ! change name to test compilation
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
!@auth rewrite for cubed sphere: M. Kelley
c NOTE: CURRENTLY ASSUMING THAT THERE IS NO TRANSPORT OF ICE TO/FROM
c EQUATORIAL CUBE FACES.  WILL UPGRADE AS NEEDED.
      USE CONSTANT, only : byshi,lhm,grav
      USE MODEL_COM, only : im,focean,p,ptop,kocean
      USE DOMAIN_DECOMP_ATM, only : grid, GET, HALO_UPDATE
      USE GEOM, only : axyp,byaxyp,
     &     dlxsina,dlysina, ull2ucs,vll2ucs, ull2vcs,vll2vcs
      USE ICEDYN_COM_dummy, only : rsix,rsiy,rsisave
c     &     ,icij,ij_musi,ij_mvsi,ij_husi,ij_hvsi,ij_susi,ij_svsi
c#ifdef TRACERS_WATER
c     *     ,ticij,ticij_tusi,ticij_tvsi
c#endif
c      USE ICEDYN, only : grid_MIC ! needed for latlon usi,vsi?
      USE SEAICE, only : ace1i,xsi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,lmi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,apress,msicnv,fwsim
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DIAG_COM, only : oa
      IMPLICIT NONE
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=2+2*LMI
#else
      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
      INTEGER ITR
      REAL*8 TRSNOW(NTM), TRICE(NTM)
#endif
      INTEGER I,J,L,K
      REAL*8 DMHSI,ASI,YRSI,XRSI,FRSI,SICE,COUR,FAO,CNEW,
     &     ullavg,vllavg
C****
C**** USIDT_ll, VSIDT_ll  latlon-oriented U,V components of
C****                     time integrated sea ice velocity (m)
C****
C**** FAW    flux of surface water area (m^2) = USIDT*DYP
C**** FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C**** FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)
C****
      REAL*8, DIMENSION(ntrice,grid%I_STRT_HALO:grid%I_STOP_HALO) ::
     &     FMSI, FMSI_jm1
      REAL*8, DIMENSION(ntrice) :: AMSI, FMSI_im1
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO) ::
     &     FASI, FXSI, FYSI, FAW,
     &     FASI_jm1, FXSI_jm1, FYSI_jm1, FAW_jm1
      REAL*8 ::
     &     FASI_im1, FXSI_im1, FYSI_im1, FAW_im1

!@var MHS mass/heat/salt content of sea ice
      REAL*8 MHS(NTRICE,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO)

      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     FOA,BYFOA, USIDT_ll,VSIDT_ll, UDYDT,VDXDT

      INTEGER I_0,I_1,J_0,J_1, I_0Y,I_1Y

      call stop_model('testing',255)

C**** Get grid parameters
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

      i_0y = max(1 ,i_0-1)
      i_1y = min(im,i_1+1)

C**** Reduce ice concentration gradients if ice amounts decreased
      DO J=J_0,J_1
      DO I=I_0,I_1
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
C**** update RSISAVE for diagnostics
        RSISAVE(I,J)=RSI(I,J)
C**** define area arrays
        FOA(I,J)=AXYP(I,J)*FOCEAN(I,J)
        IF(FOCEAN(I,J).gt.0) THEN
          BYFOA(I,J)=BYAXYP(I,J)/FOCEAN(I,J)
        ELSE
          BYFOA(I,J)=0.
        END IF
C**** set up local MHS array to contain all advected quantities
C**** MHS(1:2) = MASS, MHS(3:2+LMI) = HEAT, MHS(3+LMI:2+2*LMI)=SALT
C**** Currently this is on atmospheric grid
        MHS(1,I,J) = ACE1I + SNOWI(I,J)
        MHS(2,I,J) = MSI(I,J)
        DO L=1,LMI
          MHS(L+2,I,J) = HSI(L,I,J)
          MHS(L+2+LMI,I,J) = SSI(L,I,J)
        END DO
#ifdef TRACERS_WATER
C**** add tracers to advected arrays
        DO ITR=1,NTM
          IF (SNOWI(I,J)*XSI(2).gt.XSI(1)*ACE1I) THEN ! layer 1:all snow
            SICE=SSI(1,I,J)+SSI(2,I,J)
            TRSNOW(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J)*MAX(1.
     *           -(ACE1I-SICE)/(XSI(2)*(ACE1I+SNOWI(I,J))-SICE),0d0)
          ELSE                  ! first layer is snow and some ice
            TRSNOW(ITR) = TRSI(ITR,1,I,J)*MIN(SNOWI(I,J)/(XSI(1)*(ACE1I
     *           +SNOWI(I,J))-SSI(1,I,J)),1d0)
          END IF
          TRICE(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J) - TRSNOW(ITR)
          MHS(1+2+(1+ITR)*LMI,I,J)=TRSNOW(ITR)
          MHS(2+2+(1+ITR)*LMI,I,J)=TRICE(ITR)
          DO L=3,LMI
            MHS(L+2+(1+ITR)*LMI,I,J)=TRSI(ITR,L,I,J)
          END DO
        END DO
#endif
      ENDDO ! i
      ENDDO ! j


C****
C**** Interpolate to obtain latlon-oriented ice velocities at
C**** cell centers.  Or get collocated latlon u,v at cell edges
C**** and skip the A-to-C step.
C****
c      call interp_routine(uice,vice,usidt_ll,vsidt_ll)

C****
C**** A-to-C velocity average and transformation to CS orientation.
C**** Move later into the transport loops.
C****
      CALL HALO_UPDATE(grid, USIDT_ll) ! not necessary if interp does it
      CALL HALO_UPDATE(grid, VSIDT_ll) ! not necessary if interp does it
      do j=j_0-1,j_1
      do i=i_0y,i_1y
        ullavg = .5*(usidt_ll(i,j)+usidt_ll(i,j+1))
        vllavg = .5*(vsidt_ll(i,j)+vsidt_ll(i,j+1))
        vdxdt(i,j) = dlxsina(i,j+1)* ! note offset
     &       (ullavg*ull2vcs(i,j+1)+vllavg*vll2vcs(i,j+1))
      enddo
      enddo
      do j=j_0,j_1
      do i=i_0-1,i_1
        ullavg = .5*(usidt_ll(i,j)+usidt_ll(i+1,j))
        vllavg = .5*(vsidt_ll(i,j)+vsidt_ll(i+1,j))
        udydt(i,j) = dlysina(i+1,j)* ! note offset
     &       (ullavg*ull2ucs(i+1,j)+vllavg*vll2ucs(i+1,j))
      enddo
      enddo

C****
C**** Update halos of transported quantities
C****
      CALL HALO_UPDATE(grid, FOCEAN) ! not needed every time.
      CALL HALO_UPDATE(grid, RSI)
      CALL HALO_UPDATE(grid, RSIX)
      CALL HALO_UPDATE(grid, RSIY)
      CALL HALO_UPDATE(grid, MHS, jdim=3)


C****
C**** Transport in the Y direction
C****
      j = j_0-1
      do i=i_0y,i_1y ! updating i-halo cells for subsequent x-sweep
        faw(i) = vdxdt(i,j)  ! should compute vdxdt here
        if(faw(i).le.0.) then
c**** sea ice velocity is southward at grid box edge
          cour = faw(i)*byaxyp(i,j+1)
          fao = faw(i)*focean(i,j+1)
          fasi(i)=fao*(rsi(i,j+1)-(1d0+cour)*rsiy(i,j+1))
          fxsi(i)=fao*rsix(i,j+1)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j+1)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j+1)
        else
c**** sea ice velocity is northward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsiy(i,j))
          fxsi(i)=fao*rsix(i,j)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        end if
        faw_jm1(i)  = faw(i)
        fasi_jm1(i) = fasi(i)
        fxsi_jm1(i) = fxsi(i)
        fysi_jm1(i) = fysi(i)
        fmsi_jm1(1:ntrice,i) = fmsi(1:ntrice,i)
      enddo

      do j=j_0,j_1
      do i=i_0y,i_1y ! updating i-halo cells for subsequent x-sweep

c _jm1 qtys are already zero. no need to re-zero.
        if(vdxdt(i,j-1).eq.0. .and. vdxdt(i,j).eq.0.) cycle

        faw(i) = vdxdt(i,j)  ! should compute vdxdt here
        if(faw(i).le.0.) then
c**** sea ice velocity is southward at grid box edge
          cour = faw(i)*byaxyp(i,j+1)
          fao = faw(i)*focean(i,j+1)
          fasi(i)=fao*(rsi(i,j+1)-(1d0+cour)*rsiy(i,j+1))
          fxsi(i)=fao*rsix(i,j+1)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j+1)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j+1)
        else
c**** sea ice velocity is northward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsiy(i,j))
          fxsi(i)=fao*rsix(i,j)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        end if
c put these into atm-grid acc arrays
c        icij(i,j,ij_mvsi)=icij(i,j,ij_mvsi)+sum(fmsi(1:2,i))
c        icij(i,j,ij_hvsi)=icij(i,j,ij_hvsi)+sum(fmsi(3:2+lmi,i))
c        icij(i,j,ij_svsi)=icij(i,j,ij_svsi)+sum(fmsi(3+lmi:2+2*lmi,i))
c#ifdef TRACERS_WATER
c        do itr=1,ntm
c          ticij(i,j,ticij_tvsi,itr)=ticij(i,j,ticij_tvsi,itr)+
c     &         sum(fmsi(3+(1+itr)*lmi:2+(2+itr)*lmi,i))
c        end do
c#endif

        if(faw_jm1(i).gt.0. .or. faw(i).lt.0.) then
! when there is inflow, use general-case formulas
          asi = rsi(i,j)*foa(i,j)
          do k=1,ntrice
            amsi(k) = asi*mhs(k,i,j) + (fmsi_jm1(k,i)-fmsi(k,i))
          enddo
          asi = asi + (fasi_jm1(i)-fasi(i))
          if(asi.le.foa(i,j)) then
            yrsi = (rsiy(i,j)*axyp(i,j)*foa(i,j)+
     &           (fysi_jm1(i)-fysi(i)) + 3d0*((faw_jm1(i)+
     &           faw(i))*asi-axyp(i,j)*(fasi_jm1(i)+fasi(i))))
     &           / (axyp(i,j) + (faw_jm1(i)-faw(i)))
            rsi(i,j)  = asi*byfoa(i,j)
            rsiy(i,j) = yrsi*byfoa(i,j)
            rsix(i,j) = rsix(i,j) + (fxsi_jm1(i)-fxsi(i))*byfoa(i,j)
            if(asi.gt.0) mhs(1:ntrice,i,j) = amsi(1:ntrice)/asi
          else
c**** sea ice crunches into itself and completely covers grid box
            rsi(i,j)   = 1d0
            rsix(i,j)  = 0.
            rsiy(i,j)  = 0.
            mhs(1,i,j) = amsi(1)/asi
            mhs(2,i,j) =(amsi(1)+amsi(2))*byfoa(i,j) - mhs(1,i,j)
            do k=1,(ntrice-2)/lmi
              mhs(3+lmi*(k-1),i,j) = amsi(3+lmi*(k-1)) / asi
              mhs(4+lmi*(k-1),i,j) = amsi(4+lmi*(k-1)) / asi
              dmhsi = (amsi(3+lmi*(k-1))+amsi(4+lmi*(k-1))
     &             +amsi(5+lmi*(k-1))
     &             +amsi(6+lmi*(k-1)))*(byfoa(i,j) -1d0 / asi )
              mhs(5+lmi*(k-1),i,j) = amsi(5+lmi*(k-1)) / asi +
     &             xsi(3)*dmhsi
              mhs(6+lmi*(k-1),i,j) = amsi(6+lmi*(k-1)) / asi +
     &             xsi(4)*dmhsi
            end do

          endif
        else
! when there is only outflow, use simpler formulas.
! why is mhs not updated here.
          rsi(i,j)  =  rsi(i,j) + (fasi_jm1(i)-fasi(i))*byfoa(i,j)
          cnew = 1d0+(faw_jm1(i)*focean(i,j-1)
     &               -faw(i)*focean(i,j))*byfoa(i,j)
          rsix(i,j) = rsix(i,j)*cnew
          rsiy(i,j) = rsiy(i,j)*cnew**2
        endif

        faw_jm1(i)  = faw(i)
        fasi_jm1(i) = fasi(i)
        fxsi_jm1(i) = fxsi(i)
        fysi_jm1(i) = fysi(i)
        fmsi_jm1(1:ntrice,i) = fmsi(1:ntrice,i)

      enddo ! i
      enddo ! j

C****
C**** Transport in the X direction
C****

      do j=j_0,j_1

      i=i_0-1
      faw(i) = udydt(i,j)  ! should compute udydt here
      if(faw(i).le.0.) then
c**** sea ice velocity is westward at grid box edge
c put limit on rsix here.
        cour = faw(i)*byaxyp(i+1,j)
        fao = faw(i)*focean(i+1,j)
        fasi(i)=fao*(rsi(i+1,j)-(1d0+cour)*rsix(i+1,j))
        fxsi(i)=faw(i)*(cour*fao*rsix(i+1,j)-3d0*fasi(i))
        fysi(i)=fao*rsiy(i+1,j)
        fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i+1,j)
      else
c**** sea ice velocity is eastward at grid box edge
c put limit on rsix here.
        cour = faw(i)*byaxyp(i,j)
        fao = faw(i)*focean(i,j)
        fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsix(i,j))
        fxsi(i)=faw(i)*(cour*fao*rsix(i,j)-3d0*fasi(i))
        fysi(i)=fao*rsiy(i,j)
        fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
      endif
      faw_im1 = faw(i)
      fasi_im1 = fasi(i)
      fxsi_im1 = fxsi(i)
      fysi_im1 = fysi(i)
      fmsi_im1(1:ntrice) = fmsi(1:ntrice,i)

      do i=i_0,i_1

c _im1 qtys are already zero. no need to re-zero.
        if(udydt(i-1,j).eq.0. .and. udydt(i,j).eq.0.) cycle

        faw(i) = udydt(i,j) ! should compute udydt here
        if(faw(i).le.0.) then
c**** sea ice velocity is westward at grid box edge
c put limit on rsix here.
          cour = faw(i)*byaxyp(i+1,j)
          fao = faw(i)*focean(i+1,j)
          fasi(i)=fao*(rsi(i+1,j)-(1d0+cour)*rsix(i+1,j))
          fxsi(i)=faw(i)*(cour*fao*rsix(i+1,j)-3d0*fasi(i))
          fysi(i)=fao*rsiy(i+1,j)
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i+1,j)
        else
c**** sea ice velocity is eastward at grid box edge
c put limit on rsix here.
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsix(i,j))
          fxsi(i)=faw(i)*(cour*fao*rsix(i,j)-3d0*fasi(i))
          fysi(i)=fao*rsiy(i,j)
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        endif
c put these into atm-grid acc arrays
c        icij(i,j,ij_musi)=icij(i,j,ij_musi)+sum(fmsi(1:2,i))
c        icij(i,j,ij_husi)=icij(i,j,ij_husi)+sum(fmsi(3:2+lmi,i))
c        icij(i,j,ij_susi)=icij(i,j,ij_susi)+sum(fmsi(3+lmi:2+2*lmi,i))
c#ifdef TRACERS_WATER
c        do itr=1,ntm
c          ticij(i,j,ticij_tusi,itr)=ticij(i,j,ticij_tusi,itr)+
c     &         sum(fmsi(3+(1+itr)*lmi:2+(2+itr)*lmi,i))
c        enddo

        if(faw_im1.gt.0. .or. faw(i).lt.0.) then
! when there is inflow, use general-case formulas
          asi = rsi(i,j)*foa(i,j)
          do k=1,ntrice
            amsi(k) = asi*mhs(k,i,j) + (fmsi_im1(k)-fmsi(k,i))
          enddo
          asi = asi + (fasi_im1-fasi(i))
          if(asi.le.foa(i,j)) then
            xrsi = (rsix(i,j)*axyp(i,j)*foa(i,j)+
     &        (fxsi_im1-fxsi(i)) + 3d0*((faw_im1+faw(i))*asi-
     &           axyp(i,j)*(fasi_im1+fasi(i))))
     &           / (axyp(i,j) + (faw_im1-faw(i)))
            rsi(i,j)  = asi*byfoa(i,j)
            rsix(i,j) = xrsi*byfoa(i,j)
            rsiy(i,j) = rsiy(i,j) + (fysi_im1-fysi(i))*byfoa(i,j)
            if (asi.gt.0) mhs(1:ntrice,i,j) = amsi(1:ntrice)/asi
          else
c**** sea ice crunches into itself and completely covers grid box
            rsi(i,j)   = 1d0
            rsix(i,j)  = 0.
            rsiy(i,j)  = 0.
            mhs(1,i,j) = amsi(1)/asi
            mhs(2,i,j) =(amsi(1)+amsi(2))*byfoa(i,j) - mhs(1,i,j)
            do k=1,(ntrice-2)/lmi
              mhs(3+lmi*(k-1),i,j) = amsi(3+lmi*(k-1)) / asi
              mhs(4+lmi*(k-1),i,j) = amsi(4+lmi*(k-1)) / asi
              dmhsi = (amsi(3+lmi*(k-1))+amsi(4+lmi*(k-1))
     &             +amsi(5+lmi*(k-1))
     &             +amsi(6+lmi*(k-1)))*(byfoa(i,j) -1d0/ asi)
              mhs(5+lmi*(k-1),i,j) = amsi(5+lmi*(k-1)) / asi +
     &             xsi(3)*dmhsi
              mhs(6+lmi*(k-1),i,j) = amsi(6+lmi*(k-1)) / asi +
     &             xsi(4)*dmhsi
            enddo
          endif
        else
! when there is only outflow, use simpler formulas
! why is mhs not updated here.
          rsi(i,j)  =  rsi(i,j) + (fasi_im1-fasi(i))*byfoa(i,j)
          cnew = 1d0+(faw_im1*focean(i-1,j)
     &               -faw(i )*focean(i  ,j))*byfoa(i,j)
          rsix(i,j) = rsix(i,j)*cnew**2
          rsiy(i,j) = rsiy(i,j)*cnew
        endif
        faw_im1  = faw(i)
        fasi_im1 = fasi(i)
        fxsi_im1 = fxsi(i)
        fysi_im1 = fysi(i)
        fmsi_im1(1:ntrice) = fmsi(1:ntrice,i)
        
C**** Limit RSIX and RSIY so that sea ice is positive at the edges.
c why is it necessary to do this after the advection?
        rsi(i,j) = max(0d0,rsi(i,j))
        if(rsi(i,j)-rsix(i,j).lt.0.)  rsix(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsix(i,j).lt.0.)  rsix(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsix(i,j).gt.1d0) rsix(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsix(i,j).gt.1d0) rsix(i,j) =1d0-rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).lt.0.)  rsiy(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsiy(i,j).lt.0.)  rsiy(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).gt.1d0) rsiy(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsiy(i,j).gt.1d0) rsiy(i,j) =1d0-rsi(i,j)

      enddo ! i
      enddo ! j

      IF (KOCEAN.ge.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
C**** Fresh water sea ice mass convergence (needed for qflux model)
            MSICNV(I,J) = RSI(I,J)*(MHS(1,I,J)+MHS(2,I,J)-SUM(MHS(3
     *           +LMI:2*LMI+2,I,J))) - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     *           +MSI(I,J)-SUM(SSI(1:LMI,I,J)))
C**** sea ice prognostic variables
            SNOWI(I,J)= MAX(0d0,MHS(1,I,J) - ACE1I)
            MSI(I,J)  = MHS(2,I,J)
            DO L=1,LMI
              HSI(L,I,J) = MHS(L+2,I,J)
            END DO
C**** ensure that salinity is only associated with ice
            SICE=MHS(1+2+LMI,I,J)+MHS(2+2+LMI,I,J)
            IF (SNOWI(I,J).gt.XSI(2)*(ACE1I+SNOWI(I,J))) THEN
              SSI(1,I,J)=0.
            ELSE
              SSI(1,I,J)=SICE*(XSI(1)*ACE1I-XSI(2)*SNOWI(I,J))/ACE1I
            END IF
            SSI(2,I,J)=SICE-SSI(1,I,J)
C**** correction of heat energy to compensate for salinity fix
            HSI(1,I,J)=HSI(1,I,J)-(MHS(1+2+LMI,I,J)-SSI(1,I,J))*LHM
            HSI(2,I,J)=HSI(2,I,J)+(MHS(1+2+LMI,I,J)-SSI(1,I,J))*LHM
            DO L=3,LMI
               SSI(L,I,J) = MHS(L+2+LMI,I,J)
            END DO
#ifdef TRACERS_WATER
C**** reconstruct tracer arrays
            DO ITR=1,NTM
              IF (ACE1I.gt.XSI(2)*(SNOWI(I,J)+ACE1I)) THEN
                TRSI(ITR,1,I,J)= MHS(2+2+(1+ITR)*LMI,I,J) *(ACE1I
     *               -XSI(2)*(SNOWI(I,J)+ACE1I))/ACE1I +MHS(1+2+(1+ITR)
     *               *LMI,I,J)
              ELSE
                TRSI(ITR,1,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)*XSI(1)*(ACE1I
     *               +SNOWI(I,J))/SNOWI(I,J)
              END IF
              TRSI(ITR,2,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)+MHS(2+2+(1+ITR)
     *             *LMI,I,J)-TRSI(ITR,1,I,J)
              DO L=3,LMI
                TRSI(ITR,L,I,J)=MHS(L+2+(1+ITR)*LMI,I,J)
              END DO
            END DO
#endif
            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     *           SUM(SSI(1:LMI,I,J)))
            END IF
          END DO
        END DO
C**** Set atmospheric arrays
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
C**** set total atmopsheric pressure anomaly in case needed by ocean
              APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+RSI(I,J)
     *             *(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV
              GTEMP(1:2,2,I,J)=((HSI(1:2,I,J)-SSI(1:2,I,J)*LHM)/
     *             (XSI(1:2)*(SNOWI(I,J)+ACE1I))+LHM)*BYSHI
#ifdef TRACERS_WATER
              GTRACER(:,2,I,J)=TRSI(:,1,I,J)/(XSI(1)*MHS(1,I,J)
     *             -SSI(1,I,J))
#endif
            END IF
          END DO
        END DO
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
              OA(I,J,13)=OA(I,J,13)+(RSI(I,J)*SUM(MHS(3:2+LMI,I,J))
     *             -RSISAVE(I,J)*SUM(HSI(1:LMI,I,J)))
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****
      RETURN
      END SUBROUTINE ADVSIcs
