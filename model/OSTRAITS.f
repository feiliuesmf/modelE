      SUBROUTINE STPGF
!@sum  STPGF calculates pressure gradient force through strait
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : grav,rrt12=>byrt12
      USE OCEAN, only : lmo,lmm,mo,opress,g0m,gzmo,s0m,szmo,dxypo,hocean
     *     ,dts
      USE STRAITS, only : nmst,lmst,ist,jst,wist,must,distpg
      IMPLICIT NONE
      INTEGER I,J,K,L,N
      REAL*8 PE,PHIE,GUP,GDN,SUP,SDN,DP,PUP,PDN,VUP,VDN
      REAL*8 VOLGSP
      REAL*8, DIMENSION(LMO,2) :: MEND,PEND,PHI,DH
C****
      DO N=1,NMST
C****
C**** Accelerate the channel speed due to the pressure gradient
C**** force, reduce the channel speed due to friction
C****
      DO K=1,2
      I=IST(N,K)
      J=JST(N,K)
C**** Calculate pressure by integrating from the top down
      PE = OPRESS(I,J)
      DO L=1,LMM(I,J)
        MEND(L,K) = MO(I,J,L)
        PEND(L,K) = PE + GRAV*MO(I,J,L)*5d-1
        PE        = PE + GRAV*MO(I,J,L)
      END DO
C**** Calculate geopotential by integrating from the bottom up,
C**** also calculate the specific volume
      PHIE = -HOCEAN(I,J)*GRAV
      DO L=LMM(I,J),1,-1
      GUP = (G0M(I,J,L)-2d0*RRT12*GZMO(I,J,L))/(MO(I,J,L)*DXYPO(J))
      GDN = (G0M(I,J,L)+2d0*RRT12*GZMO(I,J,L))/(MO(I,J,L)*DXYPO(J))
      SUP = (S0M(I,J,L)-2d0*RRT12*SZMO(I,J,L))/(MO(I,J,L)*DXYPO(J))
      SDN = (S0M(I,J,L)+2d0*RRT12*SZMO(I,J,L))/(MO(I,J,L)*DXYPO(J))
      DP  = GRAV*MEND(L,K)
      PUP = PEND(L,K) - RRT12*DP
      PDN = PEND(L,K) + RRT12*DP
      VUP = VOLGSP(GUP,SUP,PUP)
      VDN = VOLGSP(GDN,SDN,PDN)
      PHI(L,K) = PHIE + (VUP*(5d-1-RRT12)+VDN*(5d-1+RRT12))*5d-1*DP
       DH(L,K) = MEND(L,K)*(VUP+VDN)*5d-1
C**** Calculate PHI at top of current layer (or bottom of next layer)
       PHIE = PHIE + (VDN+VUP)*5d-1*DP
      END DO
      END DO
C**** Subtract Pressure Gradient Force from the mass flux
      DO L=1,LMST(N)
        MUST(L,N) = MUST(L,N) - 5d-1*DTS*WIST(N)*
     *       (( DH(L,2)+ DH(L,1))*(PEND(L,2)-PEND(L,1)) +
     +       (PHI(L,2)-PHI(L,1))*(MEND(L,2)+MEND(L,1)))/DISTPG(N)
      END DO
      END DO
C****      
      RETURN
      END

      SUBROUTINE STADV
!@sum  STADV advects tracers and water mass through the straits
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : dts,mo,dxypo,bydxypo,g0m,gxmo,gymo,gzmo,s0m,sxmo
     *     ,symo,szmo
      USE STRAITS, only : nmst,lmst,ist,jst,wist,must,distpg,g0mst,s0mst
     *     ,gxmst,sxmst,gzmst,szmst
      USE ODIAG, only : olnst,ln_mflx,ln_gflx,ln_sflx
      IMPLICIT NONE
      INTEGER I1,J1,I2,J2,N,L

      REAL*8 MM1,MM2,AM
C****
      DO N=1,NMST
      I1=IST(N,1)
      J1=JST(N,1)
      I2=IST(N,2)
      J2=JST(N,2)
      DO L=1,LMST(N)
      AM = DTS*MUST(L,N)
      MM1 = MO(I1,J1,L)*DXYPO(J1)
      MM2 = MO(I2,J2,L)*DXYPO(J2)
      CALL STADVT (N,L,AM,MM1,MM2,G0MST(L,N),GXMST(L,N),GZMST(L,N),
     *     G0M,GXMO,GYMO,GZMO,OLNST(L,N,LN_GFLX))
      CALL STADVT (N,L,AM,MM1,MM2,S0MST(L,N),SXMST(L,N),SZMST(L,N),
     *     S0M,SXMO,SYMO,SZMO,OLNST(L,N,LN_SFLX))
      MO(I1,J1,L) = MO(I1,J1,L) - AM*BYDXYPO(J1)
      MO(I2,J2,L) = MO(I2,J2,L) + AM*BYDXYPO(J2)
        OLNST(L,N,LN_MFLX) = OLNST(L,N,LN_MFLX) + AM
      END DO
      END DO
      RETURN
      END

      SUBROUTINE STADVT(N,L,AM,MM1,MM2,RMST,RXST,RZST,RM,RX,RY,RZ,OLN)
!@sum  STADVT advects tracers through the straits
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst,lmst,ist,jst,xst,yst,mmst
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: AM,MM1,MM2
      REAL*8, INTENT(INOUT) :: RMST,RXST,RZST
      REAL*8, DIMENSION(IM,JM,LMO), INTENT(INOUT) :: RX,RY,RZ,RM
      REAL*8, INTENT(INOUT) :: OLN
      REAl*8 A1,A2,FM1,FZ1,FM2,FZ2
      INTEGER I1,I2,J1,J2,N,L
C****
      I1=IST(N,1)
      J1=JST(N,1)
      I2=IST(N,2)
      J2=JST(N,2)
      IF(AM.LT.0.)  GO TO 200
C****
C**** Water flux is moving from grid box 1 to grid box 2  (AM > 0)
C****
       A1 = AM/MM1
      FM1 = A1*(RM(I1,J1,L) +
     *      (1d0-A1)*(XST(N,1)*RX(I1,J1,L)+YST(N,1)*RY(I1,J1,L)))
      FZ1 = A1*RZ(I1,J1,L)
       A2 = AM/MMST(L,N)
      FM2 = A2*(RMST + (1d0-A2)*RXST)
      FZ2 = A2*RZST
C**** Calculate first moments of tracer mass for grid boxes
      RX(I1,J1,L) = RX(I1,J1,L)*(1d0-A1)*(1d0-A1*XST(N,1)*XST(N,1))
      RY(I1,J1,L) = RY(I1,J1,L)*(1d0-A1)*(1d0-A1*YST(N,1)*YST(N,1))
      RXST        = RXST*(1d0-2d0*A2) - FM1 + FM2
      RX(I2,J2,L) = RX(I2,J2,L) +
     *  XST(N,2)*(FM2 - (RM(I2,J2,L)-XST(N,2)*RX(I2,J2,L))*AM/MM2)
      RY(I2,J2,L) = RY(I2,J2,L) +
     *  YST(N,2)*(FM2 - (RM(I2,J2,L)-YST(N,2)*RY(I2,J2,L))*AM/MM2)
      GO TO 300
C****
C**** Water flux is moving from grid box 2 to grid box 1  (AM < 0)
C****
  200  A1 = AM/MMST(L,N)
      FM1 = A1*(RMST - (1d0+A1)*RXST)
      FZ1 = A1*RZST
       A2 = AM/MM2
      FM2 = A2*(RM(I2,J2,L) +
     *      (1d0+A2)*(XST(N,2)*RX(I2,J2,L)+YST(N,2)*RY(I2,J2,L)))
      FZ2 = A2*RZ(I2,J2,L)
C**** Calculate first moments of tracer mass for grid boxes
      RX(I1,J1,L) = RX(I1,J1,L) -
     *  XST(N,1)*(FM1 - (RM(I1,J1,L)-XST(N,1)*RX(I1,J1,L))*AM/MM1)
      RY(I1,J1,L) = RY(I1,J1,L) -
     *  YST(N,1)*(FM1 - (RM(I1,J1,L)-YST(N,1)*RY(I1,J1,L))*AM/MM1)
      RXST        = RXST*(1d0+2d0*A1) + FM1 - FM2
      RX(I2,J2,L) = RX(I2,J2,L)*(1d0+A2)*(1d0+A2*XST(N,2)*XST(N,2))
      RY(I2,J2,L) = RY(I2,J2,L)*(1d0+A2)*(1d0+A2*YST(N,2)*YST(N,2))
C****
C**** Calculate new tracer mass, vertical moment of tracer mass,
C**** and horizontal moment of tracer mass for the straits
C****
  300 RM(I1,J1,L) = RM(I1,J1,L) -  FM1
      RZ(I1,J1,L) = RZ(I1,J1,L) -  FZ1
      RM(I2,J2,L) = RM(I2,J2,L) +  FM2
      RZ(I2,J2,L) = RZ(I2,J2,L) +  FZ2
      RMST        = RMST   + (FM1-FM2)
      RZST        = RZST   + (FZ1-FZ2)
       OLN        =  OLN   + (FM1+FM2)
      RETURN
      END

      SUBROUTINE STBDRA
!@sum  STBDRA exerts a bottom drag on the lowest layer and a side drag
!@sum  on each layer of each strait.
!@auth Gary Russell
!@ver  1.0
C**** MUST = MUST*(1-x)  is approximated by  MUST = MUST/(1+x)
C****
C**** STRACB  MMST  Mass of sea water in strait (kg)
C****         MUST  Mass flow through strait (kg/s)
C****         WIST  Width of strait (m)
C****         DIST  Length of strait (m)
C****
      USE OCEAN, only : lmo,dts
      USE STRAITS, only : nmst,lmst,must,mmst,wist,dist
      IMPLICIT NONE
      REAL*8, PARAMETER :: BDRAGX=1d0, SDRAGX=1d-1
      INTEGER N,L

      DO N=1,NMST
C**** Apply bottom drag
      L=LMST(N)
      MUST(L,N) = MUST(L,N) * MMST(L,N)**2 /
     *  (MMST(L,N)**2 + DTS*BDRAGX*ABS(MUST(L,N))*WIST(N)*DIST(N)**2)
C**** Apply side drag
      DO L=1,LMST(N)
        MUST(L,N) = MUST(L,N) * MMST(L,N)*WIST(N) /
     *       (MMST(L,N)*WIST(N) + DTS*SDRAGX*ABS(MUST(L,N))*DIST(N))
      END DO
      END DO
C****
      RETURN
      END SUBROUTINE STBDRA

      SUBROUTINE init_STRAITS(iniOCEAN)
!@sum  init_STRAITS initializes strait variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: twopi,radius
      USE OCEAN, only : im,jm,lmo,ze,mo,s0m,sxmo,symo,szmo,g0m,gxmo,gymo
     *     ,gzmo,dxypo
      USE STRAITS, only : dist,wist,nmst,distpg,rsist,rsixst,msist,hsist
     *     ,mmst,g0mst,s0mst,gxmst,sxmst,gzmst,szmst,must,lmst,ist,jst
     *     ,xst,yst
      IMPLICIT NONE
      INTEGER I,J,L,N,I1,J1,I2,J2
      REAL*8 DLON,DLAT,FLAT,DFLON,DFLAT,SLAT,DSLON,DSLAT,TLAT,DTLON
     *     ,DTLAT,G01,GZ1,S01,SZ1,G02,GZ2,S02,SZ2,FJEQ
      LOGICAL, INTENT(IN) :: iniOCEAN
C****
C**** Calculate distance of strait, distance between centers of
C**** ocean grid boxes through strait for pressure gradient force,
C**** and mass of water in the strait
C****
      DLON = TWOPI/IM
      DLAT = NINT(180d0/(JM-1))*TWOPI/360d0
      FJEQ = (JM+1d0)/2d0
      DO N=1,NMST
      DSLON = (IST(N,2)+5d-1*XST(N,2)-IST(N,1)-5d-1*XST(N,1))*DLON
      DSLAT = (JST(N,2)+5d-1*YST(N,2)-JST(N,1)-5d-1*YST(N,1))*DLAT
       SLAT = (JST(N,2)+5d-1*YST(N,2)+JST(N,1)+5d-1*YST(N,1)-2d0*FJEQ)
     *     *DLAT/2d0
      DIST(N) = RADIUS*SQRT((DSLON*COS(SLAT))**2 + DSLAT*DSLAT)
      DFLON = 5d-1*XST(N,1)*DLON
      DFLAT = 5d-1*YST(N,1)*DLAT
       FLAT = (JST(N,1)+2.5d-1*YST(N,1)-FJEQ)*DLAT
      DTLON = 5d-1*XST(N,2)*DLON
      DTLAT = 5d-1*YST(N,2)*DLAT
       TLAT = (JST(N,2)+2.5d-1*YST(N,2)-FJEQ)*DLAT
      DISTPG(N) = DIST(N) +
     +   RADIUS*SQRT((DFLON*COS(FLAT))**2 + DFLAT*DFLAT) +
     +   RADIUS*SQRT((DTLON*COS(TLAT))**2 + DTLAT*DTLAT)
      DO L=1,LMST(N)
        MMST(L,N) = DIST(N)*WIST(N)*(ZE(L)-ZE(L-1))/9.7d-4
      END DO
      DO L=LMST(N)+1,LMO
        MMST(L,N) = 0.
      END DO
      END DO
C****
      IF(iniOCEAN) THEN
C****
C**** Initialize the mass flux, potential heat, and salinity in
C**** each strait
C****
      DO N=1,NMST
      I1=IST(N,1)
      J1=JST(N,1)
      I2=IST(N,2)
      J2=JST(N,2)
      DO L=1,LMST(N)
      MUST(L,N) = 0.
      G01= (G0M(I1,J1,L)+XST(N,1)*GXMO(I1,J1,L)+YST(N,1)*GYMO(I1,J1,L))/
     /       (MO(I1,J1,L)*DXYPO(J1))
      GZ1=  GZMO(I1,J1,L)/(MO(I1,J1,L)*DXYPO(J1))
      G02= (G0M(I2,J2,L)+XST(N,2)*GXMO(I2,J2,L)+YST(N,2)*GYMO(I2,J2,L))/
     /       (MO(I2,J2,L)*DXYPO(J2))
      GZ2=  GZMO(I2,J2,L)/(MO(I2,J2,L)*DXYPO(J2))
      G0MST(L,N) = (G01+G02)*MMST(L,N)*5d-1
      GXMST(L,N) = (G02-G01)*MMST(L,N)*5d-1
      GZMST(L,N) = (GZ1+GZ2)*MMST(L,N)*5d-1
      S01= (S0M(I1,J1,L)+XST(N,1)*SXMO(I1,J1,L)+YST(N,1)*SYMO(I1,J1,L))/
     /       (MO(I1,J1,L)*DXYPO(J1))
      SZ1=  SZMO(I1,J1,L)/(MO(I1,J1,L)*DXYPO(J1))
      S02= (S0M(I2,J2,L)+XST(N,2)*SXMO(I2,J2,L)+YST(N,2)*SYMO(I2,J2,L))/
     /       (MO(I2,J2,L)*DXYPO(J2))
      SZ2=  SZMO(I2,J2,L)/(MO(I2,J2,L)*DXYPO(J2))
      S0MST(L,N) = (S01+S02)*MMST(L,N)*5d-1
      SXMST(L,N) = (S02-S01)*MMST(L,N)*5d-1
      SZMST(L,N) = (SZ1+SZ2)*MMST(L,N)*5d-1
      END DO
      DO L=LMST(N)+1,LMO
       MUST(L,N) = 0.
      G0MST(L,N) = 0.
      GXMST(L,N) = 0.
      GZMST(L,N) = 0.
      S0MST(L,N) = 0.
      SXMST(L,N) = 0.
      SZMST(L,N) = 0.
      END DO
      END DO
C**** Initialize sea ice in straits
      RSIST=0
      RSIXST=0
      MSIST=0
      HSIST=0
      END IF
C****
      RETURN
      END

      SUBROUTINE STADVI                                                 
!@sum  STADVI advects sea ice through the straits                        
!@auth Gary Russell/Gavin Schmidt                                       
!@ver  1.0                                                              
      USE OCEAN, only : im,jm,dts,dxypo,bydxypo
      USE SEAICE_COM, only : lmi,rsi,hsi,msi2=>msi,snowi
      USE SEAICE, only : xsi,ace1i
      USE STRAITS
      USE ICEDYN, only : rsix,rsiy
      USE ODIAG, only : olnst,ln_icfl
      IMPLICIT NONE

      REAL*8, SAVE ::  WDST(NMST),BYWDST(NMST)
      REAL*8 HMSI(LMI), MSI(2,IM,JM)                                        
      REAL*8 MMSI1,MMSI2,USTDT,DRSIST,ASI,ASIST,DHSI,ALPHA,RMEAN
      INTEGER I,J,K,L,N
      INTEGER, SAVE :: IFIRST=1
C****
      IF(IFIRST.EQ.1) THEN
        IFIRST = 0
        DO N=1,NMST
          WDST(N)    = WIST(N)*DIST(N)
          BYWDST(N)  = 1d0 / WDST(N)
        END DO
      END IF
C**** SET UP LOCAL ICE ARRAY
      MSI(1,:,:) = ACE1I + SNOWI
      MSI(2,:,:) = MSI2
C**** 
C**** Loop over all straits                                             
C****                                                                   
      DO 900 N=1,NMST                                                   
      IF(MUST(1,N).GT.0.) GO TO 500                                     
      IF(RSIST(N).LE.0.)  GO TO 300                                     
C****                                                                   
C**** MUST < 0: sea ice may be leaving western end of strait            
C****                                                                   
      USTDT = DTS*MUST(1,N)*DIST(N)/MMST(1,N)                           
      IF(RSIXST(N).LE.0.)  GO TO 120                                    
      IF(USTDT+2d0*RSIXST(N) .LT. 0.)  GO TO 110                        
C**** RSIXST > 0 and no sea ice reaches western end of strait           
      RSIXST(N) = RSIXST(N) + USTDT                                     
      GO TO 300                                                         
C**** RSIXST > 0 but some sea ice flows into western ocean box          
  110 DRSIST = -RSIST(N)*(USTDT+2d0*RSIXST(N))/(DIST(N)-2d0*RSIXST(N))  
      RSIST(N) =  RSIST(N) - DRSIST                                     
      RSIXST(N) = 5d-1*USTDT                                            
      GO TO 200                                                         
  120 IF(USTDT+2d0*RSIXST(N)+DIST(N) .LE. 0.)  GO TO 130                
C**** RSIXST < 0 and some sea ice flows into western ocean box          
      DRSIST   = -RSIST(N)*USTDT / (DIST(N)+2d0*RSIXST(N))              
      RSIST(N) =  RSIST(N) - DRSIST                                     
      RSIXST(N) =  RSIXST(N) + 5d-1*USTDT                               
      GO TO 200                                                         
C**** RSIXST < 0 and all sea ice in strait flows into western ocean box 
  130 DRSIST   = RSIST(N)                                               
      RSIST(N) = 0.                                                     
      RSIXST(N) = 0.                                                    
C****                                                                   
C**** Western ocean box receives sea ice from strait                    
C****                                                                   
  200 I = IST(N,1)                                                      
      J = JST(N,1)                                                      
      ASI   = RSI(I,J)*DXYPO(J)                                         
      ASIST = DRSIST*WDST(N)                                            
      MMSI1 = ASI*MSI(1,I,J) + ASIST*MSIST(1,N)                         
      MMSI2 = ASI*MSI(2,I,J) + ASIST*MSIST(2,N)                         
      DO K=1,LMI                                                       
        HMSI(K) = ASI*HSI(K,I,J) + ASIST*HSIST(K,N)                     
      END DO                                                            
C     USI(I,J) = USI(I,J)*ASI / (ASI+ASIST)                             
C     VSI(I,J) = VSI(I,J)*ASI / (ASI+ASIST)                             
        OLNST(2,N,LN_ICFL) = OLNST(2,N,LN_ICFL) - ASIST*(MSIST(1,N)
     *     +MSIST(2,N))    
      IF(ASI+ASIST.GT.DXYPO(J))  GO TO 210                              
      RSI(I,J)   = RSI(I,J)  + ASIST*BYDXYPO(J)                         
      RSIY(I,J)  = RSIY(I,J) + ASIST*BYDXYPO(J)*2d0*YST(N,1)            
      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)               
      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)               
      IF(RSI(I,J)-RSIY(I,J).GT.1d0) RSIY(I,J) =  RSI(I,J) - 1d0         
      IF(RSI(I,J)+RSIY(I,J).GT.1d0) RSIY(I,J) =  1d0 - RSI(I,J)         
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =  RSI(I,J)-1d0           
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =  1d0-RSI(I,J)           
      MSI(1,I,J) = MMSI1 / (ASI+ASIST)                                  
      MSI(2,I,J) = MMSI2 / (ASI+ASIST)                                  
      DO K=1,LMI                                                       
        HSI(K,I,J) = HMSI(K) / (ASI+ASIST)                              
      END DO                                                            
      GO TO 300                                                         
C**** Sea ice crunches into itself and completely covers grid box       
  210 RSI(I,J)  = 1d0                                                   
      RSIY(I,J) = 0.                                                    
      MSI(1,I,J) = MMSI1   / (ASI+ASIST)                                
      HSI(1,I,J) = HMSI(1) / (ASI+ASIST)                                
      HSI(2,I,J) = HMSI(2) / (ASI+ASIST)                                
      MSI(2,I,J) = (MMSI1+MMSI2)*BYDXYPO(J) - MSI(1,I,J)                
      DHSI = (HMSI(1)+HMSI(2)+HMSI(3)+HMSI(4))*(BYDXYPO(J) -            
     *     1d0/(ASI+ASIST))                                             
      HSI(3,I,J) = HMSI(3) / (ASI+ASIST) + XSI(3)*DHSI                    
      HSI(4,I,J) = HMSI(4) / (ASI+ASIST) + XSI(4)*DHSI                    
C****                                                                   
C**** Eastern ocean box may send sea ice into strait                    
C****                                                                   
  300 I = IST(N,2)                                                      
      J = JST(N,2)                                                      
      IF(RSI(I,J).LE.0.)  GO TO 900                                     
      USTDT = DTS*MUST(1,N)*DIST(N)/MMST(1,N)                           
      ALPHA = -USTDT*WIST(N)*BYDXYPO(J)                                 
      RMEAN = RSI(I,J) + YST(N,2)*RSIY(I,J)*(1d0-ALPHA)                 
      ASI   = ALPHA*RMEAN*DXYPO(J)                                      
      RSI(I,J)  = RSI(I,J) - ALPHA*RMEAN                                
      RSIY(I,J) = RSIY(I,J)*(1d0-ALPHA)*(1d0-ALPHA)                     
      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)               
      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)               
      IF(RSI(I,J)-RSIY(I,J).GT.1d0) RSIY(I,J) =  RSI(I,J) - 1d0         
      IF(RSI(I,J)+RSIY(I,J).GT.1d0) RSIY(I,J) =  1d0 - RSI(I,J)         
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =  RSI(I,J)               
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) = -RSI(I,J)               
C****                                                                   
C**** Sea ice is entering strait from eastern ocean box                 
C****                                                                   
  400 DRSIST   = ASI*BYWDST(N)                                          
      RSIXST(N) = (RSIST(N)*(MSIST(1,N)+MSIST(2,N))*RSIXST(N) +         
     *            DRSIST*(MSI(1,I,J)+MSI(2,I,J))*5d-1*(DIST(N)+USTDT))/ 
     *           (RSIST(N)*(MSIST(1,N)+MSIST(2,N)) +                    
     *             DRSIST*(MSI(1,I,J)+MSI(2,I,J)))                      
      MSIST(1,N) = (RSIST(N)*MSIST(1,N) + DRSIST*MSI(1,I,J)) /          
     *             (RSIST(N) + DRSIST)                                  
      MSIST(2,N) = (RSIST(N)*MSIST(2,N) + DRSIST*MSI(2,I,J)) /          
     *             (RSIST(N) + DRSIST)                                  
      DO K=1,LMI                                                       
        HSIST(K,N) = (RSIST(N)*HSIST(K,N) + DRSIST*HSI(K,I,J)) /        
     *       (RSIST(N) + DRSIST)                                        
      END DO                                                            
      RSIST(N) = RSIST(N) + DRSIST                                      
         OLNST(1,N,LN_ICFL) = OLNST(2,N,LN_ICFL) -
     *     ASI*(MSI(1,I,J)+MSI(2,I,J))      
      GO TO 900                                                         
C****                                                                   
  500 IF(RSIST(N).LE.0.)  GO TO 700                                     
C****                                                                   
C**** MUST > 0: sea ice may be leaving eastern end of strait            
C****                                                                   
      USTDT = DTS*MUST(1,N)*DIST(N)/MMST(1,N)                           
      IF(RSIXST(N).GE.0.)  GO TO 520                                    
      IF(USTDT+2d0*RSIXST(N) .GT. 0.)  GO TO 510                        
C**** RSIXST < 0 and no sea ice reaches eastern end of strait           
      RSIXST(N) = RSIXST(N) + USTDT                                     
      GO TO 700                                                         
C**** RSIXST < 0 but some sea ice flows into eastern ocean box          
  510 DRSIST = RSIST(N)*(USTDT+2d0*RSIXST(N)) / (DIST(N)+2d0*RSIXST(N)) 
      RSIST(N) = RSIST(N) - DRSIST                                      
      RSIXST(N) = 5d-1*USTDT                                            
      GO TO 600                                                         
  520 IF(USTDT+2d0*RSIXST(N)-DIST(N) .GE. 0.)  GO TO 530                
C**** RSIXST > 0 and some sea ice flows into eastern ocean box          
      DRSIST   = RSIST(N)*USTDT / (DIST(N)-2d0*RSIXST(N))               
      RSIST(N) = RSIST(N) - DRSIST                                      
      RSIXST(N) = RSIXST(N) + 5d-1*USTDT                                
      GO TO 600                                                         
C**** RSIXST > 0 and all sea ice in strait flows into eastern ocean box 
  530 DRSIST   = RSIST(N)                                               
      RSIST(N) = 0.                                                     
      RSIXST(N) = 0.                                                    
C****                                                                   
C**** Eastern ocean box receives sea ice from strait                    
C****                                                                   
  600 I = IST(N,2)                                                      
      J = JST(N,2)                                                      
      ASI   = RSI(I,J)*DXYPO(J)                                         
      ASIST = DRSIST*WDST(N)                                            
      MMSI1 = ASI*MSI(1,I,J) + ASIST*MSIST(1,N)                         
      MMSI2 = ASI*MSI(2,I,J) + ASIST*MSIST(2,N)                         
      DO K=1,LMI                                                       
        HMSI(K) = ASI*HSI(K,I,J) + ASIST*HSIST(K,N)                     
      END DO                                                            
C     USI(I,J) = USI(I,J)*ASI / (ASI+ASIST)                             
C     VSI(I,J) = VSI(I,J)*ASI / (ASI+ASIST)                             
         OLNST(2,N,LN_ICFL) = OLNST(2,N,LN_ICFL) +
     *     ASIST*(MSIST(1,N)+MSIST(2,N))    
      IF(ASI+ASIST.GT.DXYPO(J))  GO TO 610                              
      RSI(I,J)   = RSI(I,J)  + ASIST*BYDXYPO(J)                         
      RSIY(I,J)  = RSIY(I,J) + ASIST*BYDXYPO(J)*2d0*YST(N,2)            
      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)               
      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)               
      IF(RSI(I,J)-RSIY(I,J).GT.1d0) RSIY(I,J) =  RSI(I,J) - 1d0         
      IF(RSI(I,J)+RSIY(I,J).GT.1d0) RSIY(I,J) =  1d0 - RSI(I,J)         
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =  RSI(I,J) - 1d0         
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =  1d0 - RSI(I,J)         
      MSI(1,I,J) = MMSI1 / (ASI+ASIST)                                  
      MSI(2,I,J) = MMSI2 / (ASI+ASIST)                                  
      DO K=1,LMI                                                       
        HSI(K,I,J) = HMSI(K) / (ASI+ASIST)                              
      END DO                                                            
      GO TO 700                                                         
C**** Sea ice crunches into itself and completely covers grid box       
  610 RSI(I,J)  = 1d0                                                   
      RSIY(I,J) = 0.                                                    
      MSI(1,I,J) = MMSI1 / (ASI+ASIST)                                  
      HSI(1,I,J) = HMSI(1) / (ASI+ASIST)                                
      HSI(2,I,J) = HMSI(2) / (ASI+ASIST)                                
      MSI(2,I,J) = (MMSI1+MMSI2)*BYDXYPO(J) - MSI(1,I,J)                
      DHSI = (HMSI(1)+HMSI(2)+HMSI(3)+HMSI(4))*(BYDXYPO(J) -            
     *     1d0/(ASI+ASIST))                                             
      HSI(3,I,J) = HMSI(3) / (ASI+ASIST) + XSI(3)*DHSI                    
      HSI(4,I,J) = HMSI(4) / (ASI+ASIST) + XSI(4)*DHSI                    
C****                                                                   
C**** Western ocean box may send sea ice into strait                    
C****                                                                   
  700 I = IST(N,1)                                                      
      J = JST(N,1)                                                      
      IF(RSI(I,J).LE.0.)  GO TO 900                                     
      USTDT = DTS*MUST(1,N)*DIST(N)/MMST(1,N)                           
      ALPHA = USTDT*WIST(N)*BYDXYPO(J)                                  
      RMEAN = RSI(I,J) + YST(N,1)*RSIY(I,J)*(1d0-ALPHA)                 
      ASI   = ALPHA*RMEAN*DXYPO(J)                                      
      RSI(I,J)  = RSI(I,J) - ALPHA*RMEAN                                
      RSIY(I,J) = RSIY(I,J)*(1d0-ALPHA)*(1d0-ALPHA)                     
      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)               
      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)               
      IF(RSI(I,J)-RSIY(I,J).GT.1d0) RSIY(I,J) =  RSI(I,J) - 1d0         
      IF(RSI(I,J)+RSIY(I,J).GT.1d0) RSIY(I,J) =  1d0 - RSI(I,J)         
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =  RSI(I,J)               
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) = -RSI(I,J)               
C****                                                                   
C**** Sea ice is entering strait from western ocean box                 
C****                                                                   
  800 DRSIST   = ASI*BYWDST(N)                                          
      RSIXST(N) = (RSIST(N)*(MSIST(1,N)+MSIST(2,N))*RSIXST(N) +         
     *            DRSIST*(MSI(1,I,J)+MSI(2,I,J))*5d-1*(USTDT-DIST(N)))/ 
     *           (RSIST(N)*(MSIST(1,N)+MSIST(2,N)) +                    
     *              DRSIST*(MSI(1,I,J)+MSI(2,I,J)))                     
      MSIST(1,N) = (RSIST(N)*MSIST(1,N) + DRSIST*MSI(1,I,J)) /          
     *             (RSIST(N) + DRSIST)                                  
      MSIST(2,N) = (RSIST(N)*MSIST(2,N) + DRSIST*MSI(2,I,J)) /          
     *             (RSIST(N) + DRSIST)                                  
      DO K=1,LMI                                                       
        HSIST(K,N) = (RSIST(N)*HSIST(K,N) + DRSIST*HSI(K,I,J)) /        
     *             (RSIST(N) + DRSIST)                                  
      END DO                                                            
      RSIST(N)   =  RSIST(N) + DRSIST                                   
         OLNST(1,N,LN_ICFL) = OLNST(2,N,LN_ICFL) +
     *     ASI*(MSI(1,I,J)+MSI(2,I,J))      
C****                                                                   
  900 CONTINUE                                                          
C**** UPDATE GLOBAL ICE ARRAY
      SNOWI = MAX(0d0,MSI(1,:,:) - ACE1I)
      MSI2  = MSI(2,:,:)
C****
      RETURN                                                            
      END                                                               

