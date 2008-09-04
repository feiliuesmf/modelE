C****   
C**** OGEOM2.f    Spherical Geometry Used by Ocean    2006/10/10
C****
      Subroutine GEOMO
!@sum  GEOMO Calculates the spherical geometry for the C grid
!@auth Gary Russell
!@ver  2.0
      Use CONSTANT, Only: TWOPI,RADIUS,OMEGA

      Use GEOM,  Only: DXYPA=>DXYP, zDXYPA=>BYDXYP
      Use OCEAN, Only: IM,JM, DLON,DLAT,FJEQ, 
     *                 DXYP=>DXYPO, DXYVO, DXYS=>DXYSO, DXYN=>DXYNO,
     *                 DXP=>DXPO, DYP=>DYPO, DXV=>DXVO, DYV=>DYVO,
     *                 RLAT, COSV=>COSVO, SINP=>SINPO, COSP=>COSPO,
     *                 RAMVS,RAMVN, zDXYP=>BYDXYPO, RATOC,ROCAT,
     *                 COSM,COSQ, SINxY,TANxY, DXPGF,DYPGF,
     *                 SINI=>SINIC, COSI=>COSIC, SINU,COSU,
     *                 J1O, JMPF=>J40S, IMAXJ
      Implicit None
      Integer*4 I,J
      Real*8 LATS, !  LATitude in radians at South edge of primary cell
     *       LATN, !  LATitude in radians at North edge of primary cell
     *       SINS, !  SINe of LATS
     *       SINN, !  SINe of LATN
     *       PCOS, !  = s[cos(LAT)^2 dLAT] / S[cos(LAT) dLAT]
     *   PLAT(JM)  !  = s[LAT cos(LAT) dLAT] / S[cos(LAT) dLAT]
C****
C**** Calculate geometric parameters defined at V latitudes
C****
      Do 10 J=1,JM-1
      COSV(J)  = Cos (DLAT*(J+.5-FJEQ))
      DXV(J)   = RADIUS*DLON*COSV(J)
   10 DYV(J)   = RADIUS*DLAT
      DXV(JM)  = 0
      DYV(JM)  = 0
      COSV(0)  = 0
      COSV(JM) = 0
C****
C**** Calculate geometric parameters defined at primary latitudes
C****
      Do 20 J=1,JM
      LATN = DLAT*(J+.5-FJEQ)  ;  If(J==JM) LATN =  TWOPI/4
      LATS = DLAT*(J-.5-FJEQ)  ;  If(J==1 ) LATS = -TWOPI/4
      SINN = Sin (LATN)
      SINS = Sin (LATS)
      RLAT(J)  = DLAT*(J-FJEQ)
      SINP(J)  = Sin (RLAT(J))
      COSP(J)  = Cos (RLAT(J))
      DXYP(J)  = RADIUS*RADIUS*DLON*(SINN-SINS)
      zDXYP(J) = 1 / DXYP(J)
      DXP(J)   = .5*RADIUS*DLON*(COSV(J-1)+COSV(J))
      DYP(J)   = RADIUS*(LATN-LATS)
      COSM(J)  = .5*(COSV(J-1)+COSV(J))
      COSQ(J)  = .5*(COSV(J-1)**2 + COSV(J)**2)
      DXYS(J)  = .5*DXYP(J)
      DXYN(J)  = .5*DXYP(J)
      SINxY(J) = COSV(J-1) - COSV(J)
      TANxY(J) = SINxY(J)  / COSM(J)
      PCOS     = .5*(LATN-LATS+.5*(SIN(2*LATN)-SIN(2*LATS)))/(SINN-SINS)
      DYPGF(J) = DXYP(J) / (RADIUS*DLON*PCOS)
   20 PLAT(J)  = (COS(LATN)+LATN*SINN-COS(LATS)-LATS*SINS) / (SINN-SINS)
      RLAT(1)  = -TWOPI/4
      SINP(1)  = -1
      PLAT(1)  = -TWOPI/4
      RLAT(JM) =  TWOPI/4
      SINP(JM) =  1
      PLAT(JM) =  TWOPI/4
C****
      Do 30 J=1,JM-1
   30 DXYVO(J) = DXYN(J) + DXYS(J+1)
C****
C**** Calculate area ratios for applying source terms to V winds
C****
      Do 40 J=1,JM-1
      RAMVN(J  ) = DXYN(J  ) / (DXYN(J)+DXYS(J+1))
   40 RAMVS(J+1) = DXYS(J+1) / (DXYN(J)+DXYS(J+1))
      RAMVS(1  ) = 0
      RAMVN(JM ) = 0
C****
C**** Calculate DXPGF = DXYV/DY used by V component of PGF
C****
      DO 50 J=1,JM-1
   50 DXPGF(J)  = (DXYN(J)+DXYS(J+1)) / (RADIUS*(PLAT(J+1)-PLAT(J)))
      DXPGF(0)  =  DXYN(1 )           / (RADIUS*(PLAT(2  )-PLAT(1)))
      DXPGF(JM) =  DXYS(JM)           / (RADIUS*(PLAT(JM )-PLAT(JM-1)))
C****
C**** Calculate sines and cosines of longitude
C****
      DO 60 I=1,IM
      SINI(I)  = Sin ((I-.5)*TWOPI/IM)
      COSI(I)  = Cos ((I-.5)*TWOPI/IM)
      SINU(I)  = Sin (I*TWOPI/IM)
   60 COSU(I)  = Cos (I*TWOPI/IM)
      SINU(IM) = 0
      COSU(IM) = 1
C**** Calculate area ratios for converting atmospheric grid to ocean
      Do J=1,JM
        RATOC(J) = DXYPA(J)*zDXYP(J)
        ROCAT(J) = DXYP(J)*zDXYPA(J)  ;  EndDo
C**** Calculate JMPF = greatest J in SH where polar filter is applied
      Do J=1,JM/2
        If (RLAT(J) > -40*TWOPI/360)  Then  !  find first J north of 40S
          JMPF = J-1  ;  Exit  ;  EndIf  ;  EndDo
C****
C**** Conditions at the poles
      DO J=1,JM,JM-1
        IMAXJ(J)=1
      END DO
C**** Conditions at non-polar points
      DO J=2,JM-1
        IMAXJ(J)=IM
      END DO
      
      Return
      End
