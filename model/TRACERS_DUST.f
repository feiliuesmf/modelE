#include "rundeck_opts.h"
      SUBROUTINE tracers_dust_old
!@sum soil dust sources and sinks
!auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)

#ifndef TRACERS_WATER
      CALL dust_wet
#endif
#endif

      RETURN
      END SUBROUTINE tracers_dust_old

      SUBROUTINE dust_wet
!@sum  Computes dust wet deposition
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE constant,ONLY : Grav
      USE resolution,ONLY : Im,Jm,Lm
      USE geom,ONLY : bydxyp
      USE model_com,ONLY : Dtsrc,zatmo,jhour
      USE fluxes,ONLY : prec,trsrfflx,tr3Dsource
      USE diag_com,ONLY : ndiupt,ijdd,adiurn,adiurn_dust
#ifdef TRACERS_DUST
     &     ,idd_wet
#endif
      USE tracer_com,ONLY : n_clay,Ntm_dust,trm,trmom
      USE trdiag_com,ONLY : ijts_source,jls_3Dsource,nDustWetij,
     &     nDustWet3Djl
      USE trdiag_com,ONLY : taijs=>taijs_loc
      USE trdiag_com,ONLY : tajls=>tajls_loc
      USE tracers_dust,ONLY : prelay
      USE CLOUDS_COM, ONLY : cldmc,cldss
      USE DYNAMICS, ONLY : gz

      IMPLICIT NONE

      REAL*8,PARAMETER :: Z=700.
      INTEGER :: i,j,ih,kr,l,n,n1,naij,najl,layer
      INTEGER,DIMENSION(Jm) :: lwdep
      REAL*8,DIMENSION(jm) :: h
      REAL*8 :: y
      REAL*8 :: cloudlay,sum1,cloudfrac(im,jm),height,cldmc1(im,jm,lm)

      ih=jhour+1

#ifdef WET_DEPO_Ina
      DO J=1,6
         LWDEP(J) = 3
         H(J) = 2800
         LWDEP(JM+1-J) = 3
         H(JM+1-J) = 2800
      ENDDO
      DO J=7,12
         LWDEP(J) = 4
         H(J) = 4900
         LWDEP(JM+1-J) = 4
         H(JM+1-J) = 4900
      ENDDO
      DO J=13,16
         LWDEP(J) = 5
         H(J) = 7400
         LWDEP(JM+1-J) = 5
         H(JM+1-J) = 7400
      ENDDO
      DO J=17,23
         LWDEP(J) = 6
         H(J) = 10300
         LWDEP(JM+1-J) = 6
         H(JM+1-J) = 10300
      ENDDO
      
c**** Wet Deposition
      DO n=1,Ntm_dust
        n1=n_clay+n-1
        tr3Dsource(:,:,:,nDustWet3Djl,n1)=0D0
        trsrfflx(:,:,n1)=0D0
        DO j=1,Jm
          DO l=1,lwdep(j)
            DO i=1,Im
              y = Z*(prec(i,j)/h(j))
              IF (y > 1.) y=1.
              tr3Dsource(i,j,l,nDustWet3Djl,n1)=y*trm(i,j,l,n1)/Dtsrc
              trsrfflx(i,j,n1)=trsrfflx(i,j,n1)+
     &             tr3Dsource(i,j,l,nDustWet3Djl,n1)
              trm(i,j,l,n1)=trm(i,j,l,n1)*(1.-y)
              trmom(:,i,j,l,n1)=trmom(:,i,j,l,n1)*(1.-y)
            END DO
          END DO
        END DO
      END DO

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        naij=ijts_source(nDustWetij,n1)
        najl=jls_3Dsource(nDustWet3Djl,n1)
        taijs(:,:,naij)=taijs(:,:,naij)+trsrfflx(:,:,n1)*Dtsrc
        DO j=1,Jm
          DO l=1,lwdep(j)
            tajls(j,l,najl)=tajls(j,l,najl)+
     &           SUM(tr3Dsource(:,j,l,nDustWet3Djl,n1))*Dtsrc
          END DO
        END DO
#ifdef TRACERS_DUST
        IF (adiurn_dust == 1) THEN
          DO kr=1,Ndiupt
            IF(i == ijdd(1,kr) .AND. j == ijdd(2,kr)) THEN
              adiurn(ih,idd_wet,kr)=adiurn(ih,idd_wet,kr)
     &             +trsrfflx(i,j,n1)*bydxyp(j)
            END IF
          END DO
        END IF
#endif
      END DO
#else

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        tr3Dsource(:,:,:,nDustWet3Djl,n1)=0D0
        trsrfflx(:,:,n1)=0D0
        DO j=1,Jm
          DO i=1,im

            layer=0
            do l=lm,1,-1
              if (prelay(i,j,l).ne.0.) then 
                layer=l
                exit
              endif
            enddo

            if (layer.eq.1) then
              height=(gz(i,j,layer)-zatmo(i,j))/Grav
              y = Z*(prec(i,j)/height)
              IF (y > 1.) y=1.
              DO l=1,layer
                tr3Dsource(i,j,l,nDustWet3Djl,n1)=y*trm(i,j,l,n1)/Dtsrc
                trsrfflx(i,j,n1)=trsrfflx(i,j,n1)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)
                trm(i,j,l,n1)=trm(i,j,l,n1)*(1.-y)
                trmom(:,i,j,l,n1)=trmom(:,i,j,l,n1)*(1.-y)
              END DO

              naij=ijts_source(nDustWetij,n1)   
              najl=jls_3Dsource(nDustWet3Djl,n1)
              taijs(i,j,naij)=taijs(i,j,naij)+trsrfflx(i,j,n1)*Dtsrc
              DO l=1,Lm
                tajls(j,l,najl)=tajls(j,l,najl)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)*Dtsrc
              END DO
#ifdef TRACERS_DUST
              IF (adiurn_dust == 1) THEN
                DO kr=1,Ndiupt
                  IF(i == ijdd(1,kr) .AND. j == ijdd(2,kr)) THEN
                    adiurn(ih,idd_wet,kr)=adiurn(ih,idd_wet,kr)
     &                   +trsrfflx(i,j,n1)*bydxyp(j)
                  END IF
                END DO
              END IF
#endif

            else if (layer.ne.0) then
              height=(gz(i,j,layer)-gz(i,j,1))/Grav
              y = Z*(prec(i,j)/height)
              IF (y > 1.) y=1.
              DO l=1,layer
                tr3Dsource(i,j,l,nDustWet3Djl,n1)=y*trm(i,j,l,n1)/Dtsrc
                trsrfflx(i,j,n1)=trsrfflx(i,j,n1)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)
                trm(i,j,l,n1)=trm(i,j,l,n1)*(1.-y)
                trmom(:,i,j,l,n1)=trmom(:,i,j,l,n1)*(1.-y)
              END DO
              naij=ijts_source(nDustWetij,n1)   
              najl=jls_3Dsource(nDustWet3Djl,n1)
              taijs(i,j,naij)=taijs(i,j,naij)+trsrfflx(i,j,n1)*Dtsrc
              DO l=1,Lm
                tajls(j,l,najl)=tajls(j,l,najl)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)*Dtsrc
              END DO
#ifdef TRACERS_DUST
              IF (adiurn_dust == 1) THEN
                DO kr=1,Ndiupt
                  IF(i == ijdd(1,kr) .AND. j == ijdd(2,kr)) THEN
                    adiurn(ih,idd_wet,kr)=adiurn(ih,idd_wet,kr)
     &                   +trsrfflx(i,j,n1)*bydxyp(j)
                  END IF
                END DO
              END IF
#endif
            endif

          END DO
        END DO
      END DO

#endif

#endif

      RETURN
      END SUBROUTINE dust_wet
