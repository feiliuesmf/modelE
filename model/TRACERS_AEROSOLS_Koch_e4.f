
      MODULE AEROSOL_SOURCES
!@sum repository for Koch aerosol sources, features, etc.
!@auth Dorothy Koch
      USE TRACER_COM
      INTEGER, PARAMETER :: ndmssrc  = 1
!@var DMS_src           DMS ocean source (kg/s)
      real*8 DMS_src(im,jm,ndmssrc),DMSinput(im,jm,12)
      INTEGER, PARAMETER :: nso2src  = 2
!@var SO2_src    SO2 industry, biomass: surface sources (kg/s)
      real*8 SO2_src(im,jm,nso2src)
      INTEGER, PARAMETER :: nso2src_3d  = 2
!@var SO2_src_3d SO2 volcanic and aircraft sources (kg/s)
      real*8 SO2_src_3d(im,jm,lm,nso2src_3d)
!@var PBLH boundary layer height
!@var MDF is the mass of the downdraft flux
      real*8, DIMENSION(IM,JM):: PBLH = 0,shdtt = 0.   ! ,MDF
c note that tno3,tno3r had dimension (im,jm,lm,12) for Wang source
      real*8, DIMENSION(IM,JM,LM):: ohr,dho2r,perjr,
     *   tno3r,oh,dho2,perj,tno3
      END MODULE AEROSOL_SOURCES

      subroutine read_SO2_source(nt,iact)
!@sum reads in industrial and biomass SO2 sources
!@auth Koch
c want kg SO2/m2/s
      USE CONSTANT, only: sday
      USE MODEL_COM, only: im,jm,jmon,jday,dtsrc,zatmo,t
      USE GEOM, only: dxyp
      USE TRACER_COM
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
      USE FILEMANAGER, only: openunit,closeunit
      USE AEROSOL_SOURCES, only: SO2_src,SO2_src_3d,nso2src,
     * nso2src_3d
      USE DYNAMICS, only: pmid,pk
      implicit none
c biomass burning parameters:
Cewg Above 56N, allow biomass burning during the months of May(1.5%),    
Cewg June(36.9%), July(44.8%), August(16.2%), and September(0.6%).  
Cewg Divide monthly totals by number of time steps in the month.  
      real*8, parameter :: B60N(12)=(/0.0, 0.0, 0.0, 0.0, 0.015, 
     *  0.369, 0.448, 0.162, 0.006, 0.0, 0.0, 0.0/)
Cewg Between 40N and 56N, allow biomass burning during May(25.8%),   
Cewg June(36.6%), July(23.5%), August(12.5%), and September(1.7%).
      real*8, parameter :: B40N(12)=(/0.0, 0.0, 0.0, 0.0, 0.258, 
     *  0.366, 0.235, 0.125, 0.017, 0.0, 0.0, 0.0/)
Cewg Between 24N and 40N, allow biomass burning from September to April.
      real*8, parameter :: B24N(12)=(/0.125,0.125,0.125,0.125,
     * 0.0, 0.0,  0.0,0.0, 0.125, 0.125, 0.125, 0.125/)

      logical :: ifirst=.true.
      REAL*8 tx,txx,ty,tyy,txy,tmas,
     *  so2_ind_input(im,jm,6),cfac,
     *  tb,tbx,tbxx,tby,tbyy,tbxy,sdday,endday,fdry,addtc
     *  ,vemis,zg,vht,zh,pl1,pl2,te,hight,vtot
      real*4 craft(im,jm,lm)
      character*56  titleg
      integer i,j,iuc,ij,iact,nt,
     *  iuc2,ib,jb,iburn,
     *  iv,jv,ivht,ir,l,najl
      save ifirst,so2_ind_input

      if (ifirst) then
      call openunit('SO2_IND',iuc,.false.)
      DO 10 ij = 1,9999
      read(iuc,902)   I,J,TMAS,TX,TXX,TY,TYY,TXY      
 902  FORMAT(3X,2(I4),E11.3,5F9.2)  
      if (i.eq.0) goto 12     
      so2_ind_input(i,j,1)=tmas
      so2_ind_input(i,j,2)=tx  
      so2_ind_input(i,j,3)=ty  
      so2_ind_input(i,j,4)=txx 
      so2_ind_input(i,j,5)=tyy 
      so2_ind_input(i,j,6)=txy   
 10    continue  
 12   continue  
      call closeunit(iuc)
      endif
c I think units are: Kg S/box/yr
c We need kg SO2/m2/s
      do i=1,im
      do j=1,jm
      cfac=tr_mm(nt)/32.d0/365.d0/sday  !*dtsrc
      so2_src(i,j,1)=so2_ind_input(i,j,1)*cfac
      end do
      end do

c??      DTT=REAL(NDYN)*DT/(86400.*30.) (ADDTC=TB*B60N*DTT)

C  Read in emissions from biomass burning 
      call openunit('SO2_BIOMASS',iuc2,.false.)

      DO 154 ij=1,9999   
Cewg  SDDAY -- first day (as day of the year) of the 90 driest days 
      READ (iuc2,9051) IB,JB,TB,TBX,TBXX,TBY,TBYY,TBXY,SDDAY   
 9051  FORMAT(4X,2I5,6E10.3,F5.0)   
      IF (IB.EQ.0) GO TO 155     
C Flux is tons SO2 / 4x5 box / year. Convert to kg SO2 / 4x5 box / year 6369.1  
C Flux is tons SO2 / 4x5 box / year. Convert to kg SO2 / 4x5 box /sec
c Must be tons SO2/4x5 box/month. Convert to kg SO2/4x5 box/sec
       TB = TB * 1.e3/30.4d0/24.d0/3600.d0  
        IF ((SDDAY+89.) .LE. 365.) THEN     
          ENDDAY = SDDAY + 89.    
        ELSE       
          ENDDAY = (SDDAY + 89.) - 365.   
        ENDIF     
 158  CONTINUE    
Cewg Allow burning above 40N only in May, June, July, Aug, and Sept.  
      IBURN = 0           
      IF (JMON .GE. 5  .AND.  JMON .LE. 9) IBURN = 1    
        IF (JB .GE. 34  .AND.  IBURN .EQ. 0) GO TO 154    
        IF (JB .GE. 39  .AND.  IBURN .EQ. 1) THEN   
          ADDTC = TB*B60N(JMON)   
           GO TO 165 
        ENDIF       
        IF ((JB.GE.34 .AND. JB.LE.38)  .AND.  IBURN .EQ. 1) THEN   
          ADDTC = TB*B40N(JMON) 
           GO TO 165    
        ENDIF  
Cewg Allow burning between 32N & 40N from September through April  
        IF (JB .GE. 30  .AND.  JB .LE. 33) THEN     
          IF (JMON .LE. 4  .OR.  JMON .GE. 9) THEN    
             ADDTC = TB*B24N(JMON)
             GO TO 165    
          ELSE    
            GO TO 154  
          ENDIF   
        ENDIF    
Cewg Allow burning south of 32N on the 90 driest days of the year 
        FDRY = 1./3.     
        IF (ENDDAY .LT. SDDAY) THEN   
          IF (REAL(JDAY).GT.ENDDAY.AND.REAL(JDAY).LT.SDDAY) GOTO 154  
        ELSE      
          IF (REAL(JDAY).LT.SDDAY.OR.REAL(JDAY).GT.ENDDAY) GOTO 154   
        ENDIF     
        ADDTC = TB*FDRY    
 165    so2_src(IB,JB,2) =  ADDTC   
 400  CONTINUE   
 154  CONTINUE   
 155  call closeunit(iuc2)

c continuously erupting volcanic emissions 
c     from GEIA (Andres and Kasgnoc, 1998)
      if (ifirst) then
      so2_src_3d(:,:,:,1)=0.
      call openunit('SO2_VOLCANO',iuc2,.false.)
      do 116 ir=1,49                        
      read(iuc2,*) iv,jv,ivht,vemis    
C Convert emissions from Mg SO2/day to Kg/sec  
       vemis = vemis * 1.e3 /sday   
C Find layer to put emissions in    
C ZG is height of ground above sea level, in meters  
       zg = zatmo(iv,jv)  
C VHT is height in meters of volcanic plume above sea level  
C Assume all emissions occur at top of plume.    
       vht = real(ivht)  
       zh = zg      
       do 21 l=1,lm-1
        pl1=pmid(l,iv,jv)
        pl2=pmid(l+1,iv,jv)
        te=pk(l,iv,jv)*t(iv,jv,l)
        HIGHT = 2.9271E+01*TE*LOG(PL1/PL2)   
C ZH is height in meters of top of layer above sea level   
        ZH = ZH + HIGHT     
        IF (VHT .LT. ZH) GO TO 24    
 21    CONTINUE       
       GO TO 100     
 24    CONTINUE   
       so2_src_3d(iv,jv,l,1)=so2_src_3d(iv,jv,l,1)+vemis
c        vtot=vemis*dtsrc+vtot
c        write(6,*) 'volcsrc',l,iv,jv,vemis,zg,sday,vtot
 100   CONTINUE    
 116  CONTINUE                                                         
      call closeunit(iuc2) 
      endif
c Aircraft emissions
      if (ifirst) then
      so2_src_3d(:,:,:,2)=0.
      call openunit('AIRCRAFT',iuc2,.true.)
      DO L=1,LM           
      READ(iuc2) titleg,((craft(i,j,l),i=1,im),j=1,jm) 
      END DO            
      call closeunit(iuc2)
c craft is Kg fuel/day. Convert to Kg SO2/s. 2.3 factor
c adjusts 2015 source to 1990 source.
c 4.d0d-4 converts to kg S
c (for BC multiply this by 0.1)
      so2_src_3d(:,:,:,2)=craft(:,:,:)*4.0d-4*tr_mm(n_SO2)/32.d0
     *  /2.3d0/sday
      endif


      end subroutine read_SO2_source



cg      subroutine apply_SO2_3Dsrc
cg      USE MODEL_COM, only: im,jm,dtsrc
cg      USE TRACER_COM
cg      USE TRACER_DIAG_COM, only: jls_3Dsource,tajls,itcon_3Dsrc
cg      USE AEROSOL_SOURCES, only: SO2_src_3d
cg      integer najl,i,j,l
cg      real*8 vtot
cg
cg        najl = jls_3Dsource(1,n_SO2)
cg        do l=1,lm
cg        do j=1,jm
cg        do i=1,imaxj(j)
cg        tajls(j,l,najl)=tajls(j,l,najl)+so2_src_3d(i,j,l,1)*dtsrc
cg        trm(i,j,l,n_so2)=trm(i,j,l,n_so2)
cg     *       +so2_src_3d(i,j,l,1)*dtsrc
cg        end do
cg        end do
cg        end do
cg      call DIAGTCA(itcon_3Dsrc(1,n_SO2),n_SO2)
cg
cg        najl = jls_3Dsource(2,n_SO2)
cg        do l=1,lm
cg        do j=1,jm
cg        do i=1,imaxj(j)
cg        tajls(j,l,najl)=tajls(j,l,najl)+so2_src_3d(i,j,l,2)*dtsrc
cg        trm(i,j,l,n_so2)=trm(i,j,l,n_so2)
cg     *       +so2_src_3d(i,j,l,2)*dtsrc
cg        end do
cg        end do
cg        end do
cg      call DIAGTCA(itcon_3Dsrc(2,n_SO2),n_SO2)
cg
cg      end subroutine apply_SO2_3Dsrc

      subroutine read_DMS_sources(nt,iact)
!@sum reads in DMS ocean source
!@auth Koch
c Monthly DMS ocean concentration sources are read in and combined
c  with wind and ocean temperature functions to get DMS air surface
c  concentrations
c want kg DMS/m2/s
      USE CONSTANT, only: sday,grav,rgas,teeny
      USE MODEL_COM, only: im,jm,jmon,focean,t,dtsrc
      USE GEOM, only: bydxyp,imaxj
      USE TRACER_COM
      USE TRACER_DIAG_COM, only: tajls,jls_source
      USE FILEMANAGER, only: openunit,closeunit
      USE SEAICE_COM, only:rsi
      USE PBLCOM, only: wsavg,eabl,tsavg
      USE CLOUDS_COM, only : airx
      USE FLUXES, only : gtemp
      USE AEROSOL_SOURCES, only: DMSinput,DMS_src,PBLH,SHDTT ! ,MDF
      USE DYNAMICS, only: BYAM,pmid,pk
      implicit none
      logical :: ifirst=.true.
      integer, parameter :: nanns=0, nmons=1
      integer najl
      REAL*8 swind,ot,schm,scrat,akw,erate,conc,steppd,
     *  foc,f_ice_free,grdi,dms1(jm)
      real*8 wt,aa,wtke,wd,arho,cp,wm,sig,wdf,
     * x,bess,dydx,yy,erate2,w1,w2,besf,bessi0,exx

      integer i,j,jj,nt,iact,iu,k,m,mont,ii

       cp=1005.20d0 !J/kg/K

      steppd = 1./sday
      dms1(:)=0.d0
c
      do j=1,jm
      do i=1,im
      if (dmsinput(i,j,jmon).gt.teeny) then
      dms1(j)=dms1(j)+1
      endif
      end do
      end do

      do j=1,jm
      do i=1,imaxj(j)
       DMS_src(i,j,1)=0
       erate2=0.d0
       erate=0.d0
         swind=wsavg(i,j)  !m/s
         ot=gtemp(1,1,i,j)   !mixed layer temp, C
         foc=focean(i,j) !fraction of gridbox with water
         f_ice_free=1.-rsi(i,j)
c Liss and Merlivat         
         schm = 2674.d0-147.12d0*ot+3.726d0*ot**2-0.038d0*ot**3d0   
         scrat = 600.d0/schm 
         IF (swind.LE.3.6d0) THEN                  
          akw = 0.041d0*swind*(scrat**0.667d0)      
         ELSE IF (swind.GT.3.6d0.AND.swind.LE.13.d0) THEN   
          akw = (0.68d0*swind - 2.31d0)*SQRT(scrat)     
         ELSE                           
          akw = (1.42d0*swind - 11.8d0)*SQRT(scrat) 
         ENDIF         
c My function  
c     AKW=0.0  
c     IF (SWIND.GT.2.) AKW=0.86*(SWIND-2.)*(1.-ot)*SQRT(SCRAT)
c Wanninkov                           
c     SCRAT = 660./SCHM     
c     AKW = 0.074*SWIND*SWIND*SQRT(SCRAT)*(1.-ot) 
      erate=akw*dmsinput(i,j,jmon)*1.d-9*tr_mm(nt)*
     *  steppd*foc*f_ice_free  !*dtsrc
c use this for Tans et al. source     
c     AKW = 0.0                     
c     IF (SWIND.GT.3.0)          
c    *  AKW = 0.016*(SWIND -3.)  
c     ERATE = AKW*CONC*1.E-6*TCMASS(5)*NDYN*DT/(3600.*24.*365.)  
c    *      *(1.-ODATA(I,J,2))*(ODATA(I,J,1)+273.15)*.082 
cccccccccccccccccccccc   
c Subgrid surface wind parameterization, using Liss and Merlivat 
         IF (swind.LE.3.6d0) THEN
         wt=0d0
         aa=0.041d0*(scrat**0.667d0)
         w1=0d0
         w2=3.6d0
         ELSE IF (swind.GT.3.6d0.AND.swind.LE.13.d0) THEN
         wt=3.40d0
         aa=0.68d0*SQRT(scrat)
         w1=3.6d0
         w2=13.d0
         ELSE 
         wt=8.31d0
         aa=1.42d0*SQRT(scrat)
         w1=13.d0
         w2=50.d0
         ENDIF
c TKE contribution
c        wtke=sqrt(2d0/3d0*eabl(1,i,j,1)*byam(1,i,j)*bydxyp(j))
         wtke=sqrt(2d0/3d0*eabl(1,i,j,1))
c dry convection contribution
         arho=1d2*PMID(1,I,J)/(RGAS*T(I,J,1)*PK(1,I,J))   !kg/m3
c I'm not sure about the sign here
         if (shdtt(i,j).lt.0d0) then
         wdf=(-shdtt(i,j)*1000.d0*grav*PBLH(I,J)
     *              /arho/cp/tsavg(I,J))
        wd=wdf**(0.33333333333d0)
        else
        wd=0.d0
        endif
c moist convection contribution
c        wm=200.d0*MDF(i,j)/dtsrc/arho*100.d0/grav
C**** use already saved AIRX array - are you sure this is what you want????
        wm=200.d0*(AIRX(i,j)/3d0)/dtsrc/arho*100.d0/grav
c sigma
        sig=wm+wd+wtke
c integrate
         yy=0d0
         do ii=1,100
         x=w1+(w2-w1)*(ii-1)/99.d0
         if (sig.eq.0) go to 22   ! added by gavin to prevent NaN
         besf=x*swind/(sig*sig)
         if (besf.gt.700.d0) go to 22
         exx=dexp(-(x*x+swind*swind)/(2.d0*sig*sig))
         bess=BESSI0(besf,i,j,exx)
         dydx=x*(x-wt)/(sig*sig)*bess

         yy=yy+dydx*(w2-w1)/99.d0
         if (bess.lt.teeny) go to 23
         enddo
 23    continue
       erate2=aa*yy*dmsinput(i,j,jmon)*1.d-9*tr_mm(nt)*
     *  steppd*foc*f_ice_free
c we can augment this by a factor of order 1
       erate2=erate2*4.d0
c       if (erate2.gt.erate) write(6,*) 'SubDMS!',
c    * erate2,erate
c       if (dmsinput(i,j,1,jmon).gt.15.) then
c       write(6,*) 'DMSdiag',erate,erate2,sig,
c    * wm,wd,wtke,mdf(i,j)
c       endif
c       if (dmsinput(i,j,1,jmon).gt.5.and.mdf(i,j).gt.teeny) then
c       write(6,*) 'DMS WD',sig,
c    * wm,wd,wtke,mdf(i,j),arho
c       endif
c       write(6,*) 'DMSwind',swind,yy,sig,
c    *  dmsinput(i,j,1,jmon),i,j,wtke,wd,wm
c       write(6,*) aa,tr_mm(nt),foc,f_ice_free
c       write(6,*) mdf(i,j),arho,
c    *  shdtt(i,j),pblh(i,j),tsavg(i,j),wdf
c
        if (sig.gt.teeny.and.dmsinput(i,j,jmon).gt.teeny) then 
        najl = jls_source(2,nt)
        tajls(j,1,najl) = tajls(j,1,najl)+wtke/sig*100.d0/dms1(j)
        najl = jls_source(3,nt)
        tajls(j,1,najl) = tajls(j,1,najl)+wm/sig*100.d0/dms1(j)
        najl = jls_source(4,nt)
        tajls(j,1,najl) = tajls(j,1,najl)+wd/sig*100.d0/dms1(j)
        endif
        najl = jls_source(5,nt)
        if (erate2.gt.erate.and.dmsinput(i,j,jmon).gt.teeny) then
        tajls(j,1,najl)=tajls(j,1,najl)
     *      +(erate2*1.D12-erate*1.D12)/(1.D12*erate)*100.d0/dms1(j)
        erate=erate2
        endif

 22    continue
ccccccccccccccccccccccccc
          DMS_src(i,j,1)=erate                 
c source -> return
      end do
      end do

      return
      end subroutine read_DMS_sources

      REAL*8 FUNCTION BESSI0(X,i,j,exx)
      IMPLICIT NONE
      real*8 ax,y
      real*8, parameter:: P1=1.0D0,P2=3.5156229D0,
     * P3=3.0899424D0,
     * P4=1.2067492D0,
     * P5=0.2659732D0,P6=0.360768D-1,P7=0.45813D-2
      real*8, parameter:: Q1=0.39894228D0,Q2=0.1328592D-1,
     *   Q3=0.225319D-2,Q4=-0.157565D-2,Q5=0.916281D-2,
     *   Q6=-0.2057706D-1,
     *    Q7=0.2635537D-1,Q8=-0.1647633D-1,Q9=0.392377D-2
      real*8, intent(in)::x,exx
      integer i,j
      IF (ABS(X).LT.3.75d0) THEN
        Y=(X/3.75d0)**2.d0
        BESSI0=exx*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=dABS(X)
        Y=3.75d0/AX
        BESSI0=(dEXP(AX)/dSQRT(AX))*exx*
     *     (Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
c     if (i.eq.1.and.j.eq.4) write(6,*) 'DMSbug1',
c    * Y,BESSI0,X,AX
      end function bessi0



      subroutine aerosol_gas_chem
!@sum aerosol gas phase chemistry
!@auth Dorothy Koch
      USE TRACER_COM
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
     *     ,jls_OHcon,jls_HO2con,jls_NO3,jls_phot
      USE MODEL_COM, only: im,jm,jmon,lm,jhour,dtsrc,t,q,jday
      USE DYNAMICS, only: pmid,am,pk
      USE GEOM, only: dxyp,imaxj
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit
      USE AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,
     & dho2,perj,
     * tno3
c Aerosol chemistry
      implicit none
      logical :: ifirst=.true.
      real*8 ohx(36,24,lm),dho2x(36,24,lm),perjx(36,24,lm) 
cg      real*8 told(im,jm,lm,ntm),ttemp(im,jm,lm)
      real*8 ppres,te,tt,mm,dmm,ohmc,r1,d1,r2,d2,ttno3,r3,d3,
     * ddno3,dddms,ddno3a,fmom,dtt
      real*8 rk4,ek4,r4,d4
      real*8 r6,d6,ek9,ek9t,ch2o,eh2o,dho2mc,dho2kg,eeee,xk9,
     * r5,d5,dmssink
      integer i,j,l,n,iuc,iun,itau,ixx1,ixx2,ichemi,itopen,itt,
     * ittime,isp,iix,jjx,llx,ii,jj,ll,iuc2,it,nm,najl
      save ifirst,itopen,iuc

c Do I need special treatment at the poles (like before?)

cg       told(:,:,:,:) = trm(:,:,:,:)
ccccccccccccccccc
c Use these for Wang NO3 and Spivakovsky for others
c open and read in NO3 file
c       if(ifirst) then
c       call openunit('AER_NO3',iuc2,.true.)
c       read(iuc2) tno3r
c       call closeunit(iuc2)
c       endif
c 204   format(8(E9.3,x))
c read in chemistry;Leave OH, HO2 in units of molecules/cm^3
c
c       itau=(jday-1)*24+jhour 
c 3     IF (itau.LT.13200) THEN  
c       itau = itau + 8760    
c       GO TO 3     
c       ENDIF  
c       ixx1=MOD(itau,8760)    
c       ixx2=INT(REAL(ixx1)/120.)   
c       ichemi = 120*ixx2+8760    
c 4     IF (ichemi.LT.13200) THEN     
c       ichemi = ichemi + 8760   
c       GO TO 4    
c       ENDIF                       
c 5     IF (ichemi.GT.21840) THEN    
c       ichemi = ichemi - 8760   
c       GO TO 5           
c       ENDIF     
c      write(6,*) 'hourly S',jhour,itau,ixx1,ixx2,
c    * ichemi,itopen
c       IF (itopen.EQ.ichemi) GO TO 27  
c       write(6,*) 'S opening',jhour,ichemi,itau
c      if(ifirst) call openunit('AER_CHEM',iuc,.true.)
c        call openunit('AER_CHEM',iuc,.true.)
c       ifirst = .false.
c       DO I=1,im*jm*lm  
c        ohr(I,1,1)=0.0      
c        perjr(I,1,1)=0.0    
c        dho2r(I,1,1)=0.0   
c       END DO   
c       DO 334 itt=1,73   
c       READ(iuc) ittime  
c      write(6,*) ittime,ichemi
c       IF (ittime.NE.ichemi) THEN 
c        READ(iuc) (ohx(i,1,1),i=1,36*24*lm)
c        READ(iuc) (dho2x(i,1,1),i=1,36*24*lm)
c        READ(iuc) (perjx(i,1,1),i=1,36*24*lm)
c        GO TO 334    
c       ELSE    
c        READ(iuc) (ohx(i,1,1),i=1,36*24*lm)
c        READ(iuc) (dho2x(i,1,1),i=1,36*24*lm)
c        READ(iuc) (perjx(i,1,1),i=1,36*24*lm)
c        do ll=1,lm
c        do ii=1,36
c        do jj=1,24
c        IF (jj.EQ.1) THEN   
c         ohr(2*ii,1,ll)=ohx(ii,jj,ll)    
c         ohr(2*ii-1,1,ll)=ohx(ii,jj,ll)
c         dho2r(2*ii,1,ll)=dho2x(ii,jj,ll)
c         dho2r(2*ii-1,1,ll)=dho2x(ii,jj,ll)
c         perjr(2*ii,1,ll)=perjx(ii,jj,ll)
c         perjr(2*ii-1,1,ll)=perjx(ii,jj,ll)
c        else if (jj.eq.24) then   
c         ohr(2*ii,46,ll)=ohx(ii,jj,ll)
c         ohr(2*ii-1,46,ll)=ohx(ii,jj,ll)    
c         dho2r(2*ii,46,ll)=dho2x(ii,jj,ll)    
c         dho2r(2*ii-1,46,ll)=dho2x(ii,jj,ll)   
c         perjr(2*ii,46,ll)=perjx(ii,jj,ll)  
c         perjr(2*ii-1,46,ll)=perjx(ii,jj,ll)   
c        else           
c         ohr(2*ii,2*jj-1,ll)=ohx(ii,jj,ll) 
c         ohr(2*ii-1,2*jj-1,ll)=ohx(ii,jj,ll)   
c         ohr(2*ii,2*jj-2,ll)=ohx(ii,jj,ll)
c         ohr(2*ii-1,2*jj-2,ll)=ohx(ii,jj,ll)   
c         dho2r(2*ii,2*jj-1,ll)=dho2x(ii,jj,ll)    
c         dho2r(2*ii-1,2*jj-1,ll)=dho2x(ii,jj,ll)    
c         dho2r(2*ii,2*jj-2,ll)=dho2x(ii,jj,ll)  
c         dho2r(2*ii-1,2*jj-2,ll)=dho2x(ii,jj,ll)   
c         perjr(2*ii,2*jj-1,ll)=perjx(ii,jj,ll) 
c         perjr(2*ii-1,2*jj-1,ll)=perjx(ii,jj,ll)    
c         perjr(2*ii,2*jj-2,ll)=perjx(ii,jj,ll) 
c         perjr(2*ii-1,2*jj-2,ll)=perjx(ii,jj,ll)   
c        endif       
c        end do    
c        end do
c        end do
c       endif  
c 334   continue    
c 203    FORMAT(X,3I3,X,3(E9.3,1X))    
c 999   call closeunit(iuc)
c       itopen=ichemi 
cccccccccccccc
c Use this for chem inputs from B4360C0M23, from Drew
        call openunit('AER_CHEM',iuc,.true.)
        do ii=1,jmon
        read(iuc) ichemi
        read(iuc) ohr
        read(iuc) dho2r
        read(iuc) perjr
        read(iuc) tno3r
        end do
        call closeunit(iuc)

c need to scale TNO3, OH and PERJ using cosine of zenith angle   
 27   CALL SCALERAD    

       dtt=dtsrc           
      do 20 l=1,lm       
      do 21 j=1,jm   
      do 22 i=1,imaxj(j)    

C**** initialise source arrays
       tr3Dsource(i,j,l,1,n_DMS)=0.  ! DMS chem sink
       tr3Dsource(i,j,l,1,n_MSA)=0.  ! MSA chem sink
       tr3Dsource(i,j,l,3,n_SO2)=0.  ! SO2 chem source
       tr3Dsource(i,j,l,4,n_SO2)=0.  ! SO2 chem sink
       tr3Dsource(i,j,l,1,n_SO4)=0.  ! SO4 chem source
       tr3Dsource(i,j,l,1,n_H2O2_s)=0. ! H2O2 chem source      
       tr3Dsource(i,j,l,2,n_H2O2_s)=0. ! H2O2 chem sink

c ptop,psf(surface),psfmpt,sige,sig
c I used to have to treat these differently above the tropopause??
c pmid=plij*sig(l)+ptop ;plij=p or psf-ptop 
c pk=pmid**kapa 
      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm=am(l,i,j)*dxyp(j)
      tt = 1.d0/te         
c here are the old ones, in case I have a problem
c      IF (L.GE.LS1) THEN               
c      PPRES = (SIG(L)*(PSF-PTOP)+PTOP)*9.869d-4 !in atm   
c      TE = T(I,J,L)*(SIG(L)*(PSF-PTOP)+PTOP)**KAPA 
c      MM = (PSF-PTOP)*DSIG(L)*100./GRAV*DXYP(J) 
c      ELSE                      
c      PPRES = (SIG(L)*P(I,J) + PTOP)*9.869E-4 !in atm   
c      TE = EXPBYK(SIG(L)*P(I,J) + PTOP)*T(I,J,L)  
c      MM = P(I,J)*DSIG(L)*100./GRAV*DXYP(J)   
c      ENDIF               
c DMM is number density of air in molecules/cm3  
      dmm=ppres/(.082d0*te)*6.02d20          
       ohmc = oh(i,j,l)  !oh is alread in units of molecules/cm3 

       do 23 n=1,ntm

       select case (trname(n))

        case ('DMS')
C***1.DMS + OH -> 0.75SO2 + 0.25MSA        
C***2.DMS + OH -> SO2         
C***3.DMS + NO3 -> HNO3 + SO2    
       r1 = 1.7d-22*dmm*0.21d0*1.d-20*exp(7810.d0*tt)/    
     *(1.d0+5.5d-20*exp(7460.d0*tt)*dmm*0.21d0*1.d-11)*ohmc     
       d1 = exp(-r1*dtsrc)                   
       r2 = 9.6d-12 * exp(-234.d0*tt)*ohmc   
       d2 = exp(-r2*dtsrc)   
c NO3 is in mixing ratio: convert to molecules/cm3  
c - not necessary for Shindell source  
       if (l.gt.8) then               
       ttno3=0.d0       
       go to 87   
       endif   
       ttno3 = tno3(i,j,l) !*6.02d20*ppres/(.082056d0*te) 
 87    r3 = ttno3*1.9d-13*exp(520.d0*tt)        
       d3= exp(-r3*dtsrc) 
       ddno3=r3*trm(i,j,l,n)/tr_mm(n)*1000.d0*dtsrc     
       dddms=trm(i,j,l,n)/tr_mm(n)*1000.d0    
       if (ddno3.gt.dddms) ddno3=dddms    
c      if (l.gt.7) then      
c      ddno3a=0.0    
c      go to 89   
c      endif   
c       ddno3a=tno3(i,j,l,jmon)*mm/.2897d0   
c 89    if (ddno3.gt.ddno3a) ddno3=ddno3a   
       ddno3=ddno3*0.9      
C DMS losses: eqns 1, 2 ,3  

cg       trm(i,j,l,n) = trm(i,j,l,n)*d1*d2   
       tr3Dsource(i,j,l,1,n) = trm(i,j,l,n)*(d1*d2-1.)/dtsrc

       dmssink=ddno3*tr_mm(n)/1000.d0
cg       if (dmssink.gt.trm(i,j,l,n)) dmssink=trm(i,j,l,n)
cg       trm(i,j,l,n)=trm(i,j,l,n)-dmssink  
cg       if (trm(i,j,l,n).lt.0.) write(6,*)'sssdmschem',i,j,l,
cg     * d1,d2,dmssink,trm(i,j,l,n),ddno3

       if (dmssink.gt.trm(i,j,l,n)+tr3Dsource(i,j,l,1,n)*dtsrc)
     *      dmssink=trm(i,j,l,n)+tr3Dsource(i,j,l,1,n)*dtsrc
       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) - dmssink/dtsrc

c  31   TSUM(8) = TSUM(8) + TrM(I,J,L,n) - TCO(n)   
 
cg this is now done automatically in apply_tracer_3Dsource
cg        najl = jls_3Dsource(1,n)
cg        tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-told(i,j,l,n))

        case ('MSA')
C MSA gain: eqn 1                
cg       TrM(I,J,L,n) = TrM(I,J,L,n) +      
cg     *0.25d0*Tr_mm(n)/Tr_mm(n_dms)*told(i,j,l,n_dms)*(1.d0 -D1)*SQRT(D2)   

          tr3Dsource(i,j,l,1,n) = 0.25d0*Tr_mm(n)/Tr_mm(n_dms)*trm(i,j
     *         ,l,n_dms)*(1.d0 -D1)*SQRT(D2)/dtsrc

cg this is now done automatically in apply_tracer_3Dsource
cg        najl = jls_3Dsource(1,n)
cg        tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-told(i,j,l,n))

       case ('SO2')
c SO2 production from DMS
cg       trm(i,j,l,n) = trm(i,j,l,n)     
cg     * +0.75*tr_mm(n)/tr_mm(n_dms)*told(i,j,l,n_dms)
cg     * *(1.d0 - d1)*sqrt(d2)    
cg     * + tr_mm(n)/tr_mm(n_dms)*told(i,j,l,n_dms)*(1.d0 - d2)*sqrt(d1)   
cg       trm(i,j,l,n)=trm(i,j,l,n)+dmssink*tr_mm(n)/tr_mm(n_dms)    

         tr3Dsource(i,j,l,3,n) = (0.75*tr_mm(n)/tr_mm(n_dms)*trm(i,j,l
     *        ,n_dms)*(1.d0 - d1)*sqrt(d2)+ tr_mm(n)/tr_mm(n_dms)*trm(i
     *        ,j,l,n_dms)*(1.d0 - d2)*sqrt(d1)+dmssink*tr_mm(n)
     *        /tr_mm(n_dms))/dtsrc

cg this is now done automatically in apply_tracer_3Dsource
cg       najl = jls_3Dsource(3,n)
cg       tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-told(i,j,l,n))
   
cg this stays since it is a simple local diagnostic
         najl = jls_NO3
         tajls(j,l,najl) = tajls(j,l,najl)+ttno3
       end select

 23    CONTINUE
 22    CONTINUE                 
 21    CONTINUE            
 20    CONTINUE       

cg this is now done automatically in apply_tracer_3Dsource
cg       call DIAGTCA(itcon_3Dsrc(3,n_SO2),n_SO2)


      do 30 l=1,lm       
      do 31 j=1,jm   
      do 32 i=1,imaxj(j)

      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm=am(l,i,j)*dxyp(j)
      tt = 1.d0/te         
      dmm=ppres/(.082d0*te)*6.02d20          
       ohmc = oh(i,j,l)  !oh is alread in units of molecules/cm3 

       do 33 n=1,ntm

       select case (trname(n))


       case ('SO2')
c oxidation of SO2 to make SO4: SO2 + OH -> H2SO4
       rk4 = 4.0d-20 *((tt*300.d0)**(3.3d0))*dmm*1.d-11  
       ek4 = 1.d0/(1.d0 + ((log10(rk4/2.0d-12))**2.d0))  
       r4 = ohmc * (rk4/(1.d0 + rk4/2.0d-12))*(0.45d0**ek4)   
       d4 = exp(-r4*dtsrc)     
c      IF (I.EQ.30.AND.J.EQ.30.and.L.EQ.2) WRITE(6,*)'msulf',TE,DMM, 
c    *  PPRES,RK4,EK4,R4,D4,ohmc
       IF (d4.GE.1.) d4=0.99999d0   
cg       trm(i,j,l,n) = trm(i,j,l,n)-told(i,j,l,n)*(1.d0-d4) 

       tr3Dsource(i,j,l,4,n) = -trm(i,j,l,n)*(1.d0-d4)/dtsrc 

c diagnostics to save oxidant fields
        najl = jls_OHcon
        tajls(j,l,najl) = tajls(j,l,najl)+oh(i,j,l)
        najl = jls_HO2con
        tajls(j,l,najl) = tajls(j,l,najl)+dho2(i,j,l)

       case('SO4')
C SO4 production   
cg       trm(i,j,l,n) = trm(i,j,l,n) + tr_mm(n)/tr_mm(n_so2)               
cg     *   *told(i,j,l,n_so2)*(1.d0 -d4)    
cg        najl = jls_3Dsource(1,n)
cg        tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-told(i,j,l,n))

         tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n)+tr_mm(n)
     *        /tr_mm(n_so2)*trm(i,j,l,n_so2)*(1.d0 -d4)/dtsrc

c      if (i.eq.30.and.j.eq.40.and.l.eq.2) 
c    *   write(6,*)'mkso4',te,dmm,ppres,rk4,ek4,r4,d4,ohmc,
c    *   told(i,j,l,n_so2)
       case('H2O2_s')
c hydrogen peroxide formation and destruction:
C***5.H2O2 +hv -> 2OH                        
C***6.H2O2 + OH -> H2O + HO2     
C***7.H2O2 + SO2 -> H2O + SO3 (in-cloud, in CB)   
C***9.HO2 + HO2 ->H2O2 + O2   
C     HO2 + HO2 + M ->        
C     HO2 + HO2 + H2O ->    
C     HO2 + HO2 + H2O + M ->    

       r6 = 2.9d-12 * exp(-160.d0*tt)*ohmc     
       d6 = exp(-r6*dtsrc)     
       ek9 = 2.2d-13*exp(600.d0*tt)  
       ek9t = 1.9d-20*dmm*0.78d0*exp(980.d0*tt)*1.d-13   
       ch2o = q(i,j,l)*6.02d20*28.97d0/18.d0*ppres/(.082d0*te) 
       eh2o = 1.+1.4d-21*exp(2200.d0*tt)*ch2o  
       dho2mc = dho2(i,j,l)      !/mm*1.292/.033*6.02e17  
       dho2kg = dho2(i,j,l)*mm*te*.082056d0/(ppres*28.97d0*6.02d20)
       eeee = eh2o*(ek9+ek9t)*dtt*dho2mc  
       xk9 = dho2kg*eeee       
c H2O2 production: eqn 9       
cg       trm(i,j,l,n) = trm(i,j,l,n) + tr_mm(n)*xk9       

       tr3Dsource(i,j,l,1,n) = tr_mm(n)*xk9/dtsrc

cg       ttemp(i,j,l)=trm(i,j,l,n)
cg        najl = jls_3Dsource(1,n)
cg        tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-told(i,j,l,n))
c H2O2 losses:5 and 6        
       r5 = perj(i,j,l)        
       d5 = exp(-r5*dtsrc)  
cg       trm(i,j,l,n) = trm(i,j,l,n)*d5*d6  
cg        najl = jls_3Dsource(2,n)
cg        tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-ttemp(i,j,l))

c      tr3Dsource(i,j,l,2,n)=(trm(i,j,l,n)+tr_mm(n)*xk9)*(d5*d6-1.d0)
c    *      /dtsrc
       
       tr3Dsource(i,j,l,2,n)=(trm(i,j,l,n))*(d5*d6-1.d0)
     *      /dtsrc

c      if (i.eq.30.and.j.eq.35.and.l.eq.2) write(6,*) 'hchemn',
c    * d5,d6,r5,r6,xk9,dho2kg,eeee,eh2o,ek9,ek9t,dho2mc,
c    * dho2(i,j,l),mm,te,ppres,q(i,j,l)
c      if (i.eq.30.and.j.eq.10.and.l.eq.2) write(6,*) 'hchems',
c    * d5,d6,r5,r6,xk9,dho2kg,eeee,eh2o,ek9,ek9t,dho2mc,
c    * dho2(i,j,l),mm,te,ppres,q(i,j,l)

        najl = jls_phot
        tajls(j,l,najl) = tajls(j,l,najl)+perj(i,j,l)

       end select

cgc adjust moments  DONE IN APPLY_TRACER_3DSOURCE
cg       if (trm(i,j,l,n).lt.0.0) trm(i,j,l,n)=0.0  
cg       if (told(i,j,l,n).gt.trm(i,j,l,n).and.
cg     *      told(i,j,l,n).gt.1.d-10) then  
cg       fmom = trm(i,j,l,n)/told(i,j,l,n)      
cg       else         
cg       fmom=1.d0      
cg       endif     
cg       do nm=1,nmom
cg        trmom(nm,i,j,l,n)=trmom(nm,i,j,l,n)*fmom
cg      end do
    
 33    CONTINUE
 32    CONTINUE                 
 31    CONTINUE            
 30    CONTINUE       

cg       call DIAGTCA(itcon_3Dsrc(1,n_DMS),n_DMS)
cg       call DIAGTCA(itcon_3Dsrc(1,n_MSA),n_MSA)
cg       call DIAGTCA(itcon_3Dsrc(4,n_SO2),n_SO2)
cg       call DIAGTCA(itcon_3Dsrc(1,n_SO4),n_SO4)
cg       call DIAGTCA(itcon_3Dsrc(1,n_H2O2_s),n_H2O2_s)

c BC: insoluble -> soluble                                              6643.4  
c      DBC=10.*DSO4*TCMASS(6)/TCMASS(3)                                 6643.5  
c      DBC=T0M(I,J,L,6)*(1.-DEXP(-9.9E-6*DT*NDYN))                      6643.55 
c      DBC=T0M(I,J,L,6)*9.736E-6*DT*NDYN  !this used last               6643.56 
c      IF (DBC.GT.T0M(I,J,L,6)) DBC=T0M(I,J,L,6)*.985                   6643.6  
c      T0M(I,J,L,6)=T0M(I,J,L,6)-DBC  !0.95*DBC                         6643.7  
c      T0M(I,J,L,7)=T0M(I,J,L,7)+DBC  !0.95*DBC                         6643.8  

   
       RETURN          
       END subroutine aerosol_gas_chem    


      SUBROUTINE SCALERAD     
      use MODEL_COM, only: im,jm,lm,jday,jhour,jmon
      use AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,dho2,perj,tno3
      use CONSTANT, only: radian
      implicit none
      real*8 ang1,xnair,vlon,vlat,ctime,timec,p1,p2,p3,fact,rad,
     *  rad1,rad2,rad3
      real*8 suncos(im,jm),tczen(im,jm)
      integer i,j,hrstrt,jdstrt,ihr,l
      integer nradn(im,jm)

      ang1 = 90.d0/91.3125d0 
      xnair = 6.022d23 * 1.d3 / 28.9644d0
      jdstrt=0
      hrstrt=0  
C*     
      DO 100 j = 1,jm   
      DO 100 i = 1,im    
C*** Calculate cos theta (RAD) at the beginning of the time step (RAD)
      vlon = 180.d0 - (i-1)*5.d0       
      vlat=-90.d0+(j-1)*4.d0     
      if (j.eq.1) vlat=-88.d0   
      if (j.eq.46) vlat=88.d0    
      ctime = (jday*24.d0) + jhour 
      timec = ctime*3600.d0    
      p1 = 15.d0*(timec/3600.d0 + hrstrt - vlon/15.d0 - 12.d0)    
      fact = (jdstrt + timec/86400.d0 - 81.1875d0)*ang1   
      p2 = 23.5d0*sin(fact*radian)      
      p3 = vlat       
C*       
      rad = (SIN(p3*radian)*SIN(p2*radian)) +    
     1            (COS(p1*radian)*COS(p2*radian)*COS(p3*radian))  
      if (rad.lt.0.d0) rad = 0.d0     
      suncos(I,J) = rad    
c     if (i.eq.30.and.j.eq.40) write(6,*) 'ohrad',ctime,jday,jhour,
c    * rad,p1,p2,p3
 100  CONTINUE     
c Scale OH and PERJ depending on time of day    
c   [OH] is approximately proportional to cos(zenith angle), 
c   therefore OH is calculated as   
c **                                                               
c **    OH4HR(i,j,l) = OH5DY(i,j,l) * 6 * SUNCOS(i,j) / TCZEN(i,j)    
c **                    
c ** where TCZEN(i,j) is the total cos(zenith angle) calculated at the   
c ** beginning of the day.                                 
c etc for NO3, PERJ. The 6 is because TCZEN is average of 6 times    
c    NO3 is different since we want to have it nonzero only during dark  
C----------------------------------------------------------------------  
        do j = 1, jm             
           do i = 1, im       
             tczen(i,j) = 0.      
             nradn(i,j)=0  
             vlon = 180.d0 - (i-1)*5.d0      
             vlat=-90.d0+(j-1)*4.d0     
             if (j.eq.1) vlat=-88.d0    
             if (j.eq.46) vlat=88.d0    
              do ihr = 0,20,4               
                 ctime = (jday*24.d0) + ihr     
                 p1 = 15.d0*(ctime - vlon/15.d0 - 12.d0)     
                 fact = (ctime/24.d0 - 81.1875d0)*ang1    
                 p2 = 23.5d0 * sin(fact*radian)     
                 p3 = vlat       
                 rad1 = sin(p3*radian) * sin(p2*radian) +    
     m             cos(p3*radian) * cos(p2*radian) * cos(p1*radian)    
                 if (rad1.lt.0.d0) rad1 = 0.d0      
                 if (rad1.eq.0.0) nradn(i,j)=nradn(i,j)+1    
c                            
                 ctime = (jday*24.d0) + ihr+2      
                 p1 = 15.d0*(ctime - vlon/15.d0 - 12.d0)    
                 fact = (ctime/24.d0 - 81.1875d0)*ang1   
                 p2 = 23.5d0 * sin(fact*radian)    
                 p3 = vlat                 
                 rad2 = sin(p3*radian) * sin(p2*radian) +    
     m             cos(p3*radian) * cos(p2*radian) * cos(p1*radian)
                 if (rad2.lt.0.d0) rad2 = 0.d0      
                 if (rad2.eq.0.0) nradn(i,j)=nradn(i,j)+2 
c                   
                 ctime = (jday*24) + ihr+4    
                 p1 = 15.d0*(ctime - vlon/15.d0 - 12.d0)   
                 fact = (ctime/24.d0 - 81.1875d0)*ang1      
                 p2 = 23.5d0 * sin(fact*radian)      
                 p3 = vlat          
                 rad3 = sin(p3*radian) * sin(p2*radian) +    
     m              cos(p3*radian) * cos(p2*radian) * cos(p1*radian)   
                 if (rad3.lt.0.d0) rad3 = 0.d0    
                 if (rad3.eq.0.0) nradn(i,j)=nradn(i,j)+1  
c                          
                 rad = (rad1 + 2.d0*rad2 + rad3)/4.d0   
                 tczen(i,j) = tczen(i,j) + rad     
                 if(tczen(i,j).eq.0.) tczen(i,j) = 1.d-32     
              end do             
           end do        
        end do      
c**************************
      do l = 1,lm      
          do j = 1,jm      
            do i = 1,im    
c Get NO3 only if dark, weighted by number of dark hours      
c        if (I.EQ.1.AND.L.EQ.1) write(6,*)'NO3R',TAU,J,NRADN(I,J)  
            if (suncos(i,j).eq.0.and.nradn(i,j).gt.0) then    
            tno3(i,j,l)=tno3r(i,j,l)*24.d0/real(nradn(i,j))    
            else        
            tno3(i,j,l)=0.d0     
            endif      
  88          if (tczen(i,j).eq.1.d-32) then       
                 oh(i,j,l) = 0.d0        
                 perj(i,j,l)=0.d0      
                 dho2(i,j,l)=0.d0   
              else    
                 oh(i,j,l)=ohr(i,j,l)*6.d0*suncos(i,j)/tczen(i,j)     
                 perj(i,j,l)=perjr(i,j,l)*6.d0*suncos(i,j)/tczen(i,j)  
c     if (i.eq.30.and.j.eq.40.and.l.eq.2) write(6,*) 'ohrad2',
c    *  suncos(i,j),tczen(i,j)
                 dho2(i,j,l)=dho2r(i,j,l)*6.d0*suncos(i,j)/tczen(i,j) 
              end if         
           end do      
         end do   
      end do    
      RETURN    
      END subroutine SCALERAD


      subroutine simple_dry_dep
!@sum simple dry deposition for aerosols
!@auth Dorothy Koch
      USE TRACER_COM
      USE TRACER_DIAG_COM, only : tajls,jls_source,itcon_dd
      USE MODEL_COM, only: im,jm,dtsrc,fland,flice,t,p
      USE DYNAMICS, only: pmid,pk
      USE GEOM, only: dxyp
c Constant dep velocity:
      implicit none
      integer i,j,nm,n
      real*8 dvz,p1,p2,te,thik,dryloss,tarea,
     * dvz_l,dvz_o,dvz_i

      real*8, parameter :: dvz_p_l=0.002
      real*8, parameter :: dvz_p_o=0.001
      real*8, parameter :: dvz_p_i=0.0003

      real*8, parameter :: dvz_so2_l=0.003
      real*8, parameter :: dvz_so2_o=0.01
      real*8, parameter :: dvz_so2_i=0.003

      real*8, parameter :: dvz_h2o2_l=0.002
      real*8, parameter :: dvz_h2o2_o=0.002
      real*8, parameter :: dvz_h2o2_i=0.002

      DO 21 J=1,JM  
      DO 20 I=1,IM  
c dvz is the dep vel in m/s
      do 22 n=1,ntm
        dvz_l=0.
        dvz_o=0.
        dvz_i=0.


      SELECT CASE(tr_wd_TYPE(N))
        CASE(nPart) 
        dvz_l=dvz_p_l
        dvz_o=dvz_p_o
        dvz_i=dvz_p_i
      END SELECT
       select case (trname(n))
       case('SO2')
        dvz_l=dvz_so2_l
        dvz_o=dvz_so2_o
        dvz_i=dvz_so2_i

       case('H2O2_s')
        dvz_l=dvz_h2o2_l
        dvz_o=dvz_h2o2_o
        dvz_i=dvz_h2o2_i

       end select

       dvz=fland(i,j)*(1.-flice(i,j))*dvz_l  !over ice-free land
     *     +(1.-fland(i,j))*dvz_o            ! over ocean
     *     +fland(i,j)*flice(i,j)*dvz_i      !over icy land
       tarea=fland(i,j)*(1.-flice(i,j))+(1.-fland(i,j))
     *  +fland(i,j)*flice(i,j)
      if (tarea.ne.1.) write(6,*)'lfrac',tarea,fland(i,j),flice(i,j)
      p1=pmid(1,i,j)
      p2=pmid(2,i,j)
      te=pk(1,i,j)*t(i,j,1)
      THIK = 2.9271E+01*TE*LOG(P1/P2)     
       DRYLOSS = dtsrc/THIK*DVZ*2.  
       trm(i,j,1,n)=trm(i,j,1,n)*(1.-dryloss)
       do nm=1,nmom
        trmom(nm,i,j,1,n)=trmom(nm,i,j,1,n)*(1.-dryloss)
       end do
 22   CONTINUE 
 20   CONTINUE  
 21   CONTINUE  

       call DIAGTCA(itcon_dd(n_MSA),n_MSA)
       call DIAGTCA(itcon_dd(n_SO2),n_SO2)
       call DIAGTCA(itcon_dd(n_SO4),n_SO4)
       call DIAGTCA(itcon_dd(n_h2o2_s),n_h2o2_s)

      RETURN
      END subroutine simple_dry_dep

      SUBROUTINE GET_SULFATE(L,temp,fcloud,
     *  tr_conv,wa_vol,wmxtr,sulfin,sulfinc,sulfout,tr_left,
     *  tm,tmcl,airm,LHX)
!@sum  GET_SULFATE calculates formation of sulfate from SO2 and H2O2
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch 
!@ver  1.0 (based on CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,mb2kg,gasc,LHE
      USE TRACER_COM, only: tr_RKD,tr_DHD,n_H2O2_s,n_SO2
     *     ,trname,ntm,tr_mm,lm,n_SO4
      USE CLOUDS, only: PL,NTIX,NTX,DXYPJ,DT_SULF_MC,DT_SULF_SS
      USE MODEL_COM, only: dtsrc
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var L index for altitude loop
!@var N index for tracer number loop
!@var RKD dummy variable (= tr_RKD*EXP[ ])
!@var WA_VOL is the cloud water volume in L
!@var WMX_INC is the change in cloud water ratio
!@var CLWC is the cloud liquid water content: L water/L air
!@var sulfin is the amount of SO2 and H2O2 used to make sulfate
!@var sulfout is the amount of sulfate generated
!@var tr_left is the amount of SO2 and H2O2 left after sulfate is made
!@+    and is now available to condense
!@var PPH is the partial pressure of the gas in M/kg 
!@var amass is airmass in kg airm(l)*mb2kg
!@var trd is dissolved portion in moles/L
!@var sulfinc is the change in dissolved SO2 and H2O2 as we use those
!@+  to form sulfate (in addition to what is dissolved)
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD, Henry_const(ntm)
      real*8 clwc,rk1f,rkdm(ntm),amass,trd(ntm),trdr(ntm),
     * dso4g,dso4d,pph(ntm),trmol(ntm),trdmol(ntm),dso4gt,dso4dt
      integer n,ih,is,is4
      real*8, parameter :: rk1=1.3E-2 !M
      real*8, parameter :: dh1=-1.6736E4 !J/mol
      real*8, parameter :: rk=6.357E14    !1/(M*M*s)
      real*8, parameter :: ea=3.95E4 !J/mol
      REAL*8,  INTENT(IN) :: fcloud,temp,wa_vol,wmxtr,LHX
      real*8 tm(lm,ntm), tmcl(ntm,lm), airm(lm)
c     REAL*8,  INTENT(OUT):: 
      real*8 sulfin(ntm),sulfout(ntm),tr_left(ntm)
     *  ,sulfinc(ntm)
      INTEGER, INTENT(IN) :: L
      LOGICAL tr_conv
      do n=1,ntx
        sulfin(N)=0.
        sulfinc(N)=0.
        sulfout(N)=0.
        tr_left(N)=1. 
      end do
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      if (LHX.NE.LHE.or.fcloud.lt.teeny) go to 333
c First allow for formation of sulfate from SO2 and H2O2. Then remaining
c  gases may be allowed to dissolve (amount given by tr_left)
C H2O2 + SO2 -> H2O + SO3 -> H2SO4 
      amass=airm(l)*mb2kg*DXYPJ
      Ppas = PL(L)*1.D2
      tfac = (1./temp - by298k)*bygasc  !mol/J
c  cloud liquid water content
      clwc=wmxtr*mair*ppas/temp*bygasc/1.D6/fcloud
      rk1f=rk1*dexp(-dh1*tfac)
 
      do n=1,ntx
       select case (trname(ntix(n)))
       case('SO2')
       is=ntix(n)
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(is)=tr_rkd(is)*(1.+ rk1f/3.2e-5)  
c mole of tracer, used to limit so4 production
      trmol(is)=1000.*tm(l,is)/tr_mm(is)*fcloud
c partial pressure of gas x henry's law coefficient                 
      pph(is)=mair*1.d-3*ppas/tr_mm(is)/amass*           
     *   tr_rkd(is)*dexp(-tr_dhd(is)*tfac)    
c the following is from Phil:                                  
c      reduction in partial pressure as species dissolves     
      henry_const(is)=rkdm(is)*dexp(-tr_dhd(is)*tfac)   
      pph(is)=pph(is)/(1+(henry_const(is)*clwc*gasc*temp))
c again all except tmcl(n,l)
      trdr(is)=mair*ppas/tr_mm(is)/amass*bygasc
     *    /temp*1.D-3  !M/kg
c dissolved moles
      trdmol(is)=trdr(is)*1000./tr_mm(is)
       case('H2O2_s')
       ih=ntix(n)
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(ih)=tr_rkd(ih)
c mole of tracer, used to limit so4 production
      trmol(ih)=1000.*tm(l,ih)/tr_mm(ih)*fcloud
c partial pressure of gas x henry's law coefficient                 
      pph(ih)=mair*1.D-3*ppas/tr_mm(ih)/amass*           
     *   tr_rkd(ih)*dexp(-tr_dhd(ih)*tfac)    
c the following is from Phil:                                  
c      reduction in partial pressure as species dissolves     
      henry_const(ih)=rkdm(ih)*dexp(-tr_dhd(ih)*tfac)   
      pph(ih)=pph(ih)/(1+(henry_const(ih)*clwc*gasc*temp))   
c all except tmcl(n,l)
      trdr(ih)=mair*ppas/tr_mm(ih)
     *    /amass*bygasc/temp*1.D-3  !M/kg
c dissolved moles
      trdmol(ih)=trdr(ih)*1000./tr_mm(ih)
      end select
      end do
      if (tm(l,ih).lt.teeny.or.tm(l,is).lt.teeny) then
      dso4g=0.
      go to 21
      endif
c this part from gas phase:moles/kg/kg
      dso4g=rk*dexp(-ea/(gasc*temp))*rk1f
     *    *pph(ih)*pph(is)*dtsrc*wa_vol
c check to make sure no overreaction: moles of production:
      dso4gt=dso4g*tm(l,ih)*tm(l,is)
c can't be more than moles going in:
      if (dso4gt.gt.trmol(is)) then
        dso4g=trmol(is)/(tm(l,ih)*tm(l,is))
      endif
      dso4gt=dso4g*tm(l,ih)*tm(l,is)
      if (dso4gt.gt.trmol(ih)) then
        dso4g=trmol(ih)/(tm(l,ih)*tm(l,is))
      endif
c this part from dissolved gases
 21    dso4d=rk*dexp(-ea/(gasc*temp))*rk1f
     *    *trdr(ih)*trdr(is)*dtsrc*wa_vol

      if (trdr(ih).lt.teeny.or.trdr(is).lt.teeny) then
      dso4d=0.
      go to 22
      endif

c check to make sure no overreaction: moles of production:
      dso4dt=dso4d*trdr(ih)*trdr(is)
c can't be more than moles going in:
      if (dso4dt.gt.trdmol(is)) then
        dso4d=trdmol(is)/(trdr(ih)*trdr(is))
      endif
      dso4dt=dso4d*trdr(ih)*trdr(is)
      if (dso4dt.gt.trdmol(ih)) then
        dso4d=trdmol(ih)/(trdr(ih)*trdr(is))
      endif

 22   continue
      do n=1,ntx
       select case (trname(ntix(n)))
       case('SO4')
       is4=ntix(n)
       sulfout(is4)=tr_mm(is4)/1000.*(dso4g*tm(l,is)*tm(l,ih)
     *  +dso4d*tmcl(is,l)*tmcl(ih,l)) !kg
       if (tr_conv) then
        dt_sulf_mc(l,is4)=dt_sulf_mc(l,is4)+sulfout(is4)
       else
        dt_sulf_ss(l,is4)=dt_sulf_ss(l,is4)+sulfout(is4)
       endif
       case('SO2')
       is=ntix(n)
       sulfin(is)=-dso4g*tm(l,ih)*tr_mm(is)/1000. !dimnless
       sulfinc(is)=-dso4d*tmcl(ih,l)*tr_mm(is)/1000.
       sulfinc(is)=max(-1d0,sulfinc(is))
       sulfin(is)=max(-1d0,sulfin(is))
       tr_left(is)=0.
       if (fcloud.gt.dabs(sulfin(is))) then
       tr_left(is)=(fcloud+sulfin(is))
       endif
       if (tr_conv) then
        dt_sulf_mc(l,is)=dt_sulf_mc(l,is)+sulfin(is)*tm(l,is)
     *   +sulfinc(is)*trdr(is)
       else
        dt_sulf_ss(l,is)=dt_sulf_ss(l,is)+sulfin(is)*tm(l,is)
     *   +sulfinc(is)*trdr(is)
       endif

       case('H2O2_s')
       ih=ntix(n)
       sulfin(ih)=-dso4g*tm(l,is)*tr_mm(ih)/1000.
       sulfinc(ih)=-dso4d*tmcl(is,l)*tr_mm(ih)/1000.
       sulfinc(ih)=max(-1d0,sulfinc(ih))
       sulfin(ih)=max(-1d0,sulfin(ih))
       tr_left(ih)=0.
       if (fcloud.gt.dabs(sulfin(ih))) then
       tr_left(ih)=fcloud+sulfin(ih)
       endif
       if (tr_conv) then
        dt_sulf_mc(l,ih)=dt_sulf_mc(l,ih)+sulfin(ih)*tm(l,ih)
     *    +sulfinc(ih)*trdr(ih)
       else
        dt_sulf_ss(l,ih)=dt_sulf_ss(l,ih)+sulfin(ih)*tm(l,ih)
     *  +sulfinc(ih)*trdr(ih)
       endif

      end select
      END DO    
   
 333  RETURN
      END SUBROUTINE GET_SULFATE
      SUBROUTINE HETER
!@sum heterogeneous production of sulfate
!@auth Dorothy Koch
      USE CONSTANT, only: BYGASC, MAIR,teeny,mb2kg,gasc,pi
      USE TRACER_COM, only:tr_RKD,tr_DHD,n_H2O2_s,n_SO2
     *     ,trname,ntm,tr_mm,n_SO4,trm,trmom
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
      USE MODEL_COM, only: im,jm,lm,dtsrc,t,p
      USE CLOUDS_COM, only:rhsav
      USE DYNAMICS, only: pmid,am,pk
      USE GEOM, only: dxyp,imaxj
      USE FLUXES, only : tr3Dsource

      IMPLICIT NONE
      integer i,j,l,n,najl
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      real*8, parameter :: rk1=1.3E-2 !M
      real*8, parameter :: dh1=-1.6736E4 !J/mol
      real*8, parameter :: rk=6.357E14    !1/(M*M*s)
      real*8, parameter :: ea=3.95E4 !J/mol
      real*8 te,amass,ppas,rk1f,tfac,ssfac,RKD,Henry_const(ntm),
     * pph(ntm),r1,A,B,pp,rr,aa,bb,BA,B2,xx,y,pn,tv,avol,
     * dso4g,dso4gt,tso2,th2o2,sulfout,sulfin,wv,ss,clwc,
     * rkdm(ntm),trmol(ntm),tt1,tt2,tt3,ptr

      ptr=0.3d0  !um
      DO 19 L=1,LM
      DO 21 J=1,JM  
      DO 20 I=1,IMAXJ(J)

C**** initialise source arrays
        tr3Dsource(i,j,l,2,n_SO4)=0.
        tr3Dsource(i,j,l,5,n_SO2)=0.
        tr3Dsource(i,j,l,3,n_H2O2_s)=0.

      amass=am(l,i,j)*DXYP(j)   !kg
      ppas=pmid(l,i,j)*100.    !Pa
      te=pk(l,i,j)*t(i,j,l)
      tfac = (1./te - by298k)*bygasc  !mol/J
c  water from deliquescence; assume complete dissociation of (NH)2SO4
c   Assume average radius of 0.1 um
      ss=rhsav(l,i,j)
      if (ss.gt.0.99d0) go to 88
         R1=log(SS)
         A=3.2E-7/te*1E6
c        B=1.09E-4*7.4 !7.4 is particle mass for 0.1 um particle
         B=1.09E-4*4.d0*pi/3.d0*1760.d0*ptr*ptr*ptr
         pp=-A/R1
         rr=B/R1
         aa=-pp*pp/3.
         bb=(2.*pp*pp*pp+27.*rr)/27.
         BA=(-bb/2.+sqrt(bb*bb/4.+aa*aa*aa/27.))**(1./3.)
         B2=(-(bb/2.+sqrt(bb*bb/4.+aa*aa*aa/27.)))**(1./3.)
         xx=BA+B2
         y=xx-pp/3.
c y is wet diameter. wv is volume of water per particle
       wv=4./3.*pi*(y-ptr)**3.d0*1000.d0/1.E18
c pn is number of particles
       pn=trm(i,j,l,n_so4)/1760.d0*3.d0/(4.d0*pi)*1.E18/(ptr)**3
c tv is total volume
       tv=wv*pn
c clwc is vol water/volume air
c      clwc=wmxtr*mair*ppas/temp*bygasc/1.D6/fcloud
c avol is air volume
      avol=amass/mair*1000.d0*gasc*te/ppas*1000.d0   !L
      clwc=tv/avol
c
      rk1f=rk1*dexp(-dh1*tfac)

      if (trm(i,j,l,n_so2).lt.teeny.or.
     *  trm(i,j,l,n_h2o2_s).lt.teeny.or.tv.lt.teeny) then 
      go to 88
      endif

      DO N=1,NTM
       select case (trname(n))
       case('SO2')
      tso2=trm(i,j,l,n)

c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(n)=tr_rkd(n)*(1.+ rk1f/3.2d-5)  
c mole of tracer, used to limit so4 production
      trmol(n)=1000.*trm(i,j,l,n)/tr_mm(n)
c partial pressure of gas x henry's law coefficient                 
      pph(n)=mair*1.d-3*ppas/tr_mm(n)/amass*           
     *   tr_rkd(n)*dexp(-tr_dhd(n)*tfac)    
c the following is from Phil:                                  
c      reduction in partial pressure as species dissolves     
      henry_const(n)=rkdm(n)*dexp(-tr_dhd(n)*tfac)   
      pph(n)=pph(n)/(1+(henry_const(n)*clwc*gasc*te))

       case('H2O2_s')
      th2o2=trm(i,j,l,n)
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(n)=tr_rkd(n)*(1.+ rk1f/3.2d-5)  
c mole of tracer, used to limit so4 production
      trmol(n)=1000.*trm(i,j,l,n)/tr_mm(n)
c partial pressure of gas x henry's law coefficient                 
      pph(n)=mair*1.d-3*ppas/tr_mm(n)/amass*           
     *   tr_rkd(n)*dexp(-tr_dhd(n)*tfac)    
c the following is from Phil:                                  
c      reduction in partial pressure as species dissolves     
      henry_const(n)=rkdm(n)*dexp(-tr_dhd(n)*tfac)   
      pph(n)=pph(n)/(1+(henry_const(n)*clwc*gasc*te))
      end select
      END DO 

c this part from gas phase:moles/kg/kg
      dso4g=rk*dexp(-ea/(gasc*te))*rk1f
     *    *pph(n_h2o2_s)*pph(n_so2)*dtsrc*tv
c check to make sure no overreaction: moles of production:
      dso4gt=dso4g*tso2*th2o2
c can't be more than moles going in:
      if (dso4gt.gt.trmol(n_so2)) then
        dso4g=trmol(n_so2)/(tso2*th2o2)
      endif
      dso4gt=dso4g*tso2*th2o2
      if (dso4gt.gt.trmol(n_h2o2_s)) then
        dso4g=trmol(n_h2o2_s)/(tso2*th2o2)
      endif


      do n=1,ntm
       select case (trname(n))
       case('SO4')
       sulfin=0.
       sulfout=tr_mm(n)/1000.*(dso4g*tso2*th2o2) !kg
c diagnostic
cg        najl = jls_3Dsource(4,n)
cg        tajls(j,l,najl)=tajls(j,l,najl)+sulfout
cg        trm(i,j,l,n)=trm(i,j,l,n)+sulfout

       tr3Dsource(i,j,l,2,n)=sulfout/dtsrc

c      tt1=tt1+sulfout 
c      if (l.eq.2.and.j.eq.34) write(6,*)'Het',i,sulfout,
c    * tv,pn,wv,y,ss
       sulfin=0.
       case('SO2')
       sulfin=-dso4g*th2o2*tr_mm(n)/1000. !dimnless
       sulfin=max(-1d0,sulfin)
cg       trm(i,j,l,n)=trm(i,j,l,n)*(1.d0+sulfin)
cg       trmom(:,i,j,l,n)=trmom(:,i,j,l,n)*trm(i,j,l,n)/tso2
       if (sulfin.eq.-1d0) sulfin = -.998d0
       tr3Dsource(i,j,l,5,n)=trm(i,j,l,n)*sulfin/dtsrc
c       tt2=tt2+trm(i,j,l,n)-tso2 
       sulfin=0.
       case('H2O2_s')
       sulfin=-dso4g*tso2*tr_mm(n)/1000. !dimnless
       sulfin=max(-1d0,sulfin)
cg       trm(i,j,l,n)=trm(i,j,l,n)*(1.d0+sulfin)
cg       trmom(:,i,j,l,n)=trmom(:,i,j,l,n)*trm(i,j,l,n)/th2o2

       tr3Dsource(i,j,l,3,n)=trm(i,j,l,n)*sulfin/dtsrc

c      if (l.eq.2.and.j.eq.34) write(6,*)'H2O2',i,sulfin,th2o2,
c    * trm(i,j,l,n)
c      tt3=tt3+trm(i,j,l,n)-th2o2 
      end select
      END DO    
   
 
  88  continue

  20  CONTINUE
  21  CONTINUE
  19  CONTINUE
c     write(6,*)'Hethet',tt1,tt2,tt3
cg      call DIAGTCA(itcon_3Dsrc(4,n_SO4),n_SO4)
cg      call DIAGTCA(itcon_3Dsrc(6,n_SO2),n_SO2)
cg      call DIAGTCA(itcon_3Dsrc(5,n_H2O2_s),n_H2O2_s)
      return
      END subroutine HETER
