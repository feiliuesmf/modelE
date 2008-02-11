      subroutine obio_ptend(vrbos,kmax,i,j)

c  Computes tendencies of biological particles (units/hr)
c  tracer
c  P(1) = nitrate (uM)
c  P(2) = ammonium (uM)
c  P(3) = silica (uM)
c  P(4) = iron (nM)
c  P(5) = diatoms (mg chl m-3)
c  P(6) = chlorophytes (mg chl m-3)
c  P(7) = cyanobacteria (mg chl m-3)
c  P(8) = coccolithophores (mg chl m-3)
c  P(9) = herbivores (mg chl m-3)
 
      USE obio_dim
      USE obio_incom,only: cnratio,cfratio,remin,obio_wsh
     .                    ,wsdeth,rkn,rks,rkf,Rm
      USE obio_forc, only: tirrq
      USE obio_com, only : dp1d,obio_P,obio_ws,P_tend,D_tend
     .                    ,gro,rlamz,dratez1,dratez2,bn
     .                    ,greff,pnoice,drate,tfac,regen,Fescav
     .                    ,wshc,rikd,rmuplsr,det
     .                    ,gcmax1d,covice_ij,atmFe_ij
     .                    ,temp1d,wsdet,tzoo,p1d
     .                    ,pnoice2,rhs

      implicit none

#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"

      integer nt,kmax

      real rmu4(nchl)     !growth on ammonium, NH4
      real rmu3(nchl)     !growth on nitrate, NO4
      real rmu5(nchl)     !growth on silica
      real rmuf(nchl)     !growth on iron
      real zoo(ntyp)      !herbivores (zooplankton)
      real dphy(ntyp)     !death rate of phytoplankton
      real viscfac(kdm)



      real bs,cchlratio,bf,ptot
      real Pzoo,Gzoo,Dzoo1,Dzoo2,exc
      real tirrqice,upn,upa,upf,ups

      real rnut2,tnit,tmp,rnut1,rmmn,framm,rmml,rmmlice,rmms
      real rmmf,rlim,rlimice,grate,rlimnfix,rlimrkn,rfix,gratenfix1
      real graterkn,gratenlim,gratenfix,gron,gronfix
      real dphyt

      logical vrbos

      real term
      character bio_var*7(12)
      data (bio_var(nt),nt=1,12) /
     .  ' NO3   ',' NH4   ',' SiO2  ',' Iron  ','Diatom ','Chlphy ',
     .  'Cyanob ','Coccol ','Herbiv ','N/Cdet ','Silica ','Fe_det '/


       rhs=0.0
       obio_ws = 0.0
       P_tend = 0.0
       D_tend = 0.0
       wsdet = 0.0
       rmu4 = 0.0
       rmu3 = 0.0
       rmu5 = 0.0
       rmuf = 0.0
       zoo  = 0.0
       dphy = 0.0
       gro = 0.0
       viscfac = 0.0



!define no ice points based on covice (here: covice_ij)
      pnoice=1.-covice_ij
      pnoice2=pnoice
      !pnoice=1.  !Watson prefers that we do not change pnoice


!tendency terms are computed in the mid-level (m)

c  Start Model Space Loop
!      m = indext2     !index of "past" (t-1)

!Iron + atm iron: disperse in layer and convert to nM
       k = 1
       term = atmFe_ij*1.d-3/max(p1d(k+1),1.e-3)
       rhs(k,4,4) = term
       P_tend(k,4) = term

       do k=2,kmax
          rhs(k,4,4) = 0.0
          P_tend(k,4) = 0.0
       enddo

!Grazing/regeneration of ammonium
       do k = 1,kmax

         bs = 2.0*bn(k)
         cchlratio = bn(k)*cnratio
         bf = cchlratio/cfratio  !Fe:chl ratio (nM/ugl)
 
         ptot = 0.0
         do nt = nnut+1,ntyp-nzoo
          ptot = ptot + obio_P(k,nt)
         enddo
         ptot = max(ptot,1.0E-36)

         Pzoo = obio_P(k,ntyp)
         gzoo = tzoo*Rm*(1.0-exp(-rlamz*ptot))*Pzoo
         dzoo1 = dratez1*Pzoo
         dzoo2 = dratez2*Pzoo*Pzoo

         term = ((1.0-greff)*gzoo-dzoo1-dzoo2) * pnoice
         rhs(k,ntyp,ntyp) = term    !!!!!!!this also takes into accout P(5,6,7,8)
         P_tend(k,ntyp) = term

cdiag    if (vrbos)
cdiag.   write(*,*)'ptend1: ',
cdiag.        nstep,i,j,k,ntyp,obio_P(k,ntyp),tzoo,Rm,rlamz,
cdiag.        ptot,Pzoo,dratez1,dratez2,greff,P_tend(k,ntyp)

         do nt = nnut+1,ntyp-nzoo
          !fraction of grazing for this group
          zoo(nt) = gzoo*obio_P(k,nt)/ptot
          dphy(nt) = drate*obio_P(k,nt)

          term = -zoo(nt) * pnoice
          rhs(k,nt,ntyp) = term
          P_tend(k,nt) = term      !!nonlin term takes into accountP(ntyp),P(5:8)

cdiag     if (vrbos) write(*,*)
cdiag.     'obio_ptend2: ',nstep,nt,i,j,k,gzoo,obio_P(k,nt),ptot,
cdiag.      zoo(nt),P_tend(k,nt)

          term = -dphy(nt) * pnoice           !death of phytoplankton
          rhs(k,nt,nt) = term
          P_tend(k,nt) = P_tend(k,nt) + term

cdiag     if (vrbos) write(*,*)
cdiag.     'obio_ptend3: ',nstep,nt,i,j,k,drate,dphy(nt),
cdiag.     obio_P(k,nt),P_tend(k,nt)

cdiag      if (vrbos) then
cdiag       if(k.eq.1 .and. nt.eq.nnut+1) write(101,'(a,a)')
cdiag.      'nstep    k   nt   dp1d   zoo(nt) dphy(nt)   obio_P(k,nt)'
cdiag.     ,'gzoo   ptot    drate  P_tend(k,nt)'
cdiag       write(101,'(3i5,8(1x,es8.2))')nstep,k,nt,dp1d(k),
cdiag.        zoo(nt),dphy(nt),obio_P(k,nt),gzoo,ptot,drate,
cdiag.        P_tend(k,nt)
cdiag      endif

         enddo

         exc = greff*gzoo

         term = tfac(k)*remin(1)*det(k,1)/cnratio * pnoice
         rhs(k,1,ntyp+1) = term    !put this in diff column
         P_tend(k,1) = term

         term = bn(k)*(exc + regen*dzoo2) * pnoice
         rhs(k,2,ntyp) = term
         P_tend(k,2) = term

         term = tfac(k)*remin(2)*det(k,2) * pnoice
         rhs(k,3,ntyp+2) = term    !put this in diff column
         P_tend(k,3) = term

         term = bf*(exc + regen*dzoo2) * pnoice    
         rhs(k,4,ntyp) = term
         P_tend(k,4) = P_tend(k,4) + term

         term = tfac(k)*remin(3)*det(k,3) * pnoice
         rhs(k,4,ntyp+3) = term                        !put this in diff column
         P_tend(k,4) = P_tend(k,4) + term

         term = -Fescav(k)
         rhs(k,4,13) = term             !this should actually be rhs(4,4) but we
                                        !have already defined it so let it be rhs(4,13)
         P_tend(k,4) = P_tend(k,4) + term
 
cdiag    if (vrbos) then
cdiag     if(k.eq.1) write(102,'(a,a)')
cdiag.    'nstep    k  dp1d(k) tfac(k)   det(k,1) bn(k)    exc   ',
cdiag.    ' det(3)   P_t(1)   P_t(2)   P_t(3) P_t(4) '
cdiag       write(102,'(2i5,10(1x,es8.2))')nstep,k,
cdiag.         dp1d(k),tfac(k),det(k,1),bn(k),exc,det(k,3),
cdiag.         P_tend(k,1),P_tend(k,2),P_tend(k,3),P_tend(k,4)
cdiag    endif
cdiag    if (vrbos) write(*,*)'ptend4: ',
cdiag.     nstep,i,j,k,dp1d(k),tfac(k),det(k,1),bn(k),exc,det(k,3),
cdiag.     P_tend(k,1),P_tend(k,2),P_tend(k,3),P_tend(k,4)


         dphyt = 0.0
         do nt = nnut+1,ntyp-nzoo
          dphyt = dphyt + dphy(nt)
         enddo
 
!1st detrital fraction is carbon
         term = dphyt*cchlratio * pnoice
         rhs(k,ntyp+1,nnut+1) = term
         D_tend(k,1) = term

         term = dzoo1*cchlratio * pnoice
     .        + (1.0-regen)*dzoo2*cchlratio * pnoice
         rhs(k,ntyp+1,ntyp) = term
         D_tend(k,1) = D_tend(k,1) + term

         term = -tfac(k)*remin(1)*det(k,1) * pnoice
         rhs(k,ntyp+1,ntyp+1) = term
         D_tend(k,1) = D_tend(k,1) + term

!2nd detrital fraction is silica
         term = bs*dphy(nnut+1) * pnoice
         rhs(k,ntyp+2,nnut+1) = term
         D_tend(k,2) = term

         term = bs*zoo(nnut+1) * pnoice
         rhs(k,ntyp+2,ntyp) = term
         D_tend(k,2) = D_tend(k,2) + term

         term = -tfac(k)*remin(2)*det(k,2) * pnoice
         rhs(k,ntyp+2,ntyp+2) = term
         D_tend(k,2) = D_tend(k,2) + term

!3rd detrital fraction is iron
         term = bf*dphyt * pnoice
         rhs(k,ntyp+3,nnut+1) = term
         D_tend(k,3) = term

         term = bf*dzoo1 * pnoice
     .        + bf*(1.0-regen)*dzoo2 * pnoice
         rhs(k,ntyp+3,ntyp) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = -tfac(k)*remin(3)*det(k,3) * pnoice
         rhs(k,ntyp+3,ntyp+3) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = Fescav(k)
         rhs(k,ntyp+3,4) = term
         D_tend(k,3) = D_tend(k,3) + term

cdiag    if (vrbos) then
cdiag     if(k.eq.1) write(103,'(a,a)')
cdiag.    'nstep    k  dp1d(k)   dphyt   dzoo1     dzoo2   dphy(5) '
cdiag.   ,'zoo(5)   tfac(k)  D_t(1)     D_t(2) D_t(3)'
cdiag    write(103,'(2i5,10(1x,es8.2))')nstep,k,dp1d(k),
cdiag.   dphyt,dzoo1,dzoo2,dphy(nnut+1),zoo(nnut+1),tfac(k),
cdiag.   D_tend(k,1),D_tend(k,2),D_tend(k,3)
cdiag    endif
cdiag    if (vrbos) write(*,*)'ptend5 :',
cdiag.   nstep,i,j,k,
cdiag.   dp1d(k),dphyt,dzoo1,dzoo2,dphy(nnut+1),zoo(nnut+1),tfac(k),
cdiag.   D_tend(k,1),D_tend(k,2),D_tend(k,3)

       enddo !kmax

c Day: Grow
      do k = 1,kmax

cdiag if(vrbos)
cdiag.write(*,*)'obio_ptend6: ',nstep,i,j,k,tirrq(k)

      if (tirrq(k) .gt. 0.0)then
        bs = 2.0*bn(k)
        cchlratio = bn(k)*cnratio
        tirrqice = tirrq(k)*0.01  !reduce light in ice by half
        bf = cchlratio/cfratio  !Fe:chl ratio (nM/ugl)

c Light-regulated growth
!!#if NCHL_DEFINED > 0
      if (nchl > 0) then

        ! Diatoms
        nt = 1

cdiag   if (vrbos) then
cdiag     if (k.eq.1)write(104,'(a,3i5)')
cdiag.                'ptend: diatoms growth, nstep,k,nt=',nstep,k,nt
cdiag     if (k.eq.1)write(104,'(5x,a,12x,a)')
cdiag.      'dp    tirrq    P(1)  P(2)  P(3)  P(4)','rkn    rks    rkf'
cdiag     write(104,'(9(1x,es8.2))')dp1d(k),tirrq(k),
cdiag.                       obio_P(k,1),obio_P(k,2),obio_P(k,3),
cdiag.                       obio_P(k,4),rkn(nt),rks(nt),rkf(nt)
cdiag   endif
        if (vrbos)
     .    write(*,*)'ptend7 :', nstep,i,j,k,dp1d(k),tirrq(k),
     .                       obio_P(k,1),obio_P(k,2),obio_P(k,3),
     .                       obio_P(k,4),rkn(nt),rks(nt),rkf(nt)

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
         tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))     !nitrate
          tmp = 1.0 - rnut2
        ! Enforce preferential utilization of ammonium
          rnut1 = min(tmp,tnit)
           rmmn = rnut1 + rnut2
          framm = rnut2/rmmn
           rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmlice = tirrqice/(tirrqice+0.5*rikd(k,nt))


        ! Silicate
        rmms = obio_P(k,3)/(rks(nt)+obio_P(k,3)) 
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
           !rhs(k,nt+nnut,1) = rmml
           !rhs(k,nt+nnut,2) = rmmn
           !rhs(k,nt+nnut,3) = rmms
           !rhs(k,nt+nnut,4) = rmmf
        rlim = min(rmml,rmmn,rmms,rmmf)
        rlimice = min(rmmlice,rmmn,rmms,rmmf)
        grate = rmuplsr(k,nt)*rlim*pnoice2
     .          + rmuplsr(k,nt)*rlimice*(1.0-pnoice2)
        rmu4(nt) = grate*framm
        rmu3(nt) = grate*(1.0-framm)
        rmu5(nt) = grate*rmms
        rmuf(nt) = grate*rmmf
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

      endif

!!#endif


!!#if NCHL_DEFINED > 1
      if (nchl > 1) then

! Chlorophytes
        nt = 2

cdiag   if (vrbos) then
cdiag     if (k.eq.1)write(105,'(a,3i5)')
cdiag.                'ptend: chlorop growth, nstep,k,nt=',nstep,k,nt
cdiag     if (k.eq.1)write(105,'(5x,a,12x,a)')
cdiag.      'dp    tirrq    P(1)  P(2)  P(3)  P(4)','rkn    rks    rkf'
cdiag     write(105,'(9(1x,es8.2))')dp1d(k),tirrq(k),
cdiag.               obio_P(k,1),obio_P(k,2),obio_P(k,4),rkn(nt),
cdiag.                       rkf(nt), gro(k,nt),obio_P(k,nt+nnut)
cdiag   endif
        if (vrbos)
     .    write(*,*)'ptend8: ',nstep,i,j,k,dp1d(k),tirrq(k),
     .               obio_P(k,1),obio_P(k,2),obio_P(k,4),rkn(nt),
     .                       rkf(nt), gro(k,nt),obio_P(k,nt+nnut)


        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit  = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp   = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2                         
        framm = rnut2/rmmn                           
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmlice = tirrqice/(tirrqice+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
           !rhs(k,nt+nnut,1) = rmml
           !rhs(k,nt+nnut,2) = rmmn
           !rhs(k,nt+nnut,4) = rmmf
        rlim = min(rmml,rmmn,rmmf)
        rlimice = min(rmmlice,rmmn,rmmf)
        grate = rmuplsr(k,nt)*rlim * pnoice2
     .        + rmuplsr(k,nt)*rlimice * (1.0-pnoice2)
        rmu4(nt) = grate*framm
        rmu3(nt) = grate*(1.0-framm)
        rmuf(nt) = grate*rmmf
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term  
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term
      endif
!!#endif


!!#if NCHL_DEFINED > 2
      if (nchl > 2) then
! Cyanobacteria
        nt = 3
        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        framm = rnut2/rmmn
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
           !rhs(k,nt+nnut,1) = rmml
           !rhs(k,nt+nnut,2) = rmmn
           !rhs(k,nt+nnut,4) = rmmf
        rlim = min(rmml,rmmn,rmmf)
        rlimnfix = min(rmml,rmmf)         !limitation for N2 fixation
        rlimrkn = min(rmml,rkn(nt),rmmf)   !limitation at kn

        grate = rmuplsr(k,nt)*rlim
        rmu4(nt) = grate*framm*pnoice2
        rmu3(nt) = grate*(1.0-framm)*pnoice2
        rfix = 0.25*exp(-(75.0*obio_P(k,nt+nnut)))
        rfix = max(rfix,0.0)
c        rfix = min(rfix,0.2)

        gratenfix1 = rmuplsr(k,nt)*rlimnfix*rfix  !N fix
        graterkn = rmuplsr(k,nt)*rlimrkn
        gratenlim = graterkn - grate
        gratenfix = min(gratenlim,gratenfix1)  !N fix cannot exceed kn
        gratenfix = max(gratenfix,0.0)

        rmuf(nt) = (grate+gratenfix)*rmmf * pnoice
        gron = grate*obio_P(k,nt+nnut)
        gronfix = gratenfix*obio_P(k,nt+nnut)
        gro(k,nt) = gron + gronfix

        term = gro(k,nt) * pnoice
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

cdiag   if (vrbos) write(*,*)
cdiag.   'obio_ptend9: ',nt,i,j,k,rmuplsr(k,nt),rlim,rlimnfix,rlimrkn,
cdiag.    rfix,gratenfix,gronfix,obio_P(k,nt+nnut),gron,gro(k,nt),
cdiag.    P_tend(k,nt+nnut)

cdiag   if (vrbos) write(*,*)
cdiag.  'obio_ptend10: ',nt,i,j,k,tirrq(k),rikd(k,nt),rmml,obio_P(k,4),
cdiag.   rkf(nt),rmmf,rkn(nt),rnut1,rnut2,tmp,tnit,rmmn,obio_P(k,1),
cdiag.   obio_P(k,2)

      endif
!!#endif


!!#if NCHL_DEFINED > 3
      if (nchl > 3) then
! Coccolithophores
        nt = 4

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        framm = rnut2/rmmn
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
           !rhs(k,nt+nnut,1) = rmml
           !rhs(k,nt+nnut,2) = rmmn
           !rhs(k,nt+nnut,4) = rmmf
        rlim = min(rmml,rmmn,rmmf)
        grate = rmuplsr(k,nt)*rlim
        rmu4(nt) = grate*framm * pnoice2
        rmu3(nt) = grate*(1.0-framm) * pnoice2
        rmuf(nt) = grate*rmmf * pnoice2
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt) * pnoice
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term
        gcmax1d(k) = max(gcmax1d(k),grate)
      endif
!!#endif


!!#if NCHL_DEFINED > 4
      if (nchl > 4) then
! Dinoflagellates
        nt = 5

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        framm = rnut2/rmmn
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
           !rhs(k,nt+nnut,1) = rmml
           !rhs(k,nt+nnut,2) = rmmn
           !rhs(k,nt+nnut,4) = rmmf
        rlim = min(rmml,rmmn,rmmf)
        grate = rmuplsr(k,nt)*rlim
        rmu4(nt) = grate*framm * pnoice2
        rmu3(nt) = grate*(1.0-framm) * pnoice2
        rmuf(nt) = grate*rmmf * pnoice2
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt) * pnoice2
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term
      endif
!!#endif

cdiag   if (vrbos) then
cdiag   do nt = 1,nchl
cdiag     write(801,'(i7,2e12.4)')
cdiag.          nstep,gro(1,nt),obio_P(1,nt+nnut)
cdiag   enddo
cdiag   endif

cdiag  if(vrbos)write(*,*)'ptend11: ',
cdiag.  nstep,i,j,k,tirrq(k),P_tend(k,5),P_tend(k,6),P_tend(k,7) 
cdiag.  ,P_tend(k,8),P_tend(k,9)


! Nutrient uptake
        upn = 0.0
        upa = 0.0
        upf = 0.0
        do nt = 1,nchl
         term = bn(k)*(rmu3(nt)*obio_P(k,nnut+nt))
         rhs(k,1,nnut+nt) = -term   !use - sign here because that is how it goes in P_tend
         upn = upn + term

         term = bn(k)*(rmu4(nt)*obio_P(k,nnut+nt))
         rhs(k,2,nnut+nt) = -term   !use - sign here because that is how it goes in P_tend
         upa = upa + term

         term = bf*(rmuf(nt)*obio_P(k,nnut+nt))
         rhs(k,4,nnut+nt) = -term   !use - sign here because that is how it goes in P_tend
         upf = upf + term
        enddo

        term = bs*(rmu5(1)*obio_P(k,nnut+1))    
        rhs(k,3,nnut+1) = -term   !use - sign here because that is how it goes in P_tend
        ups = term

cdiag  if(vrbos)then
cdiag     if(k.eq.1) write(106,'(a,a)')
cdiag.'nstep    k  dp1d(k)   bn(k)   rmu3(1)  rmu3(2)   rmu3(3)',
cdiag.'rmu3(4)  P(5)   P(6)   P(7)   P(8)     P_t(1)bef P_t(1) '
cdiag   write(106,'(2i5,12(1x,es8.2))')
cdiag. nstep,k,dp1d(k),bn(k),rmu3(1),rmu3(2),rmu3(3),
cdiag. rmu3(4),obio_P(k,nnut+1),obio_P(k,nnut+2),obio_P(k,nnut+3),
cdiag. obio_P(k,nnut+4),P_tend(k,1),P_tend(k,1) - upn
cdiag  endif

cdiag  if(vrbos)write(*,*)'obio_ptend12: ',
cdiag. nstep,i,j,k,dp1d(k),bn(k),rmu3(1),rmu3(2),rmu3(3),rmu3(4),
cdiag. obio_P(k,5),obio_P(k,6),obio_P(k,7),obio_P(k,8),
cdiag. P_tend(k,1),P_tend(k,2),P_tend(k,3),P_tend(k,4),
cdiag. upn,upa,upf,ups

        P_tend(k,1) = P_tend(k,1) - upn
        P_tend(k,2) = P_tend(k,2) - upa
        P_tend(k,3) = P_tend(k,3) - ups
        P_tend(k,4) = P_tend(k,4) - upf

       endif !tirrq(k) .gt. 0.0
      enddo  !kmax
  
cdiag if(vrbos) then
cdiag do k=1,kdm
cdiag  write(*,*)'obio_ptend13: ',
cdiag. nstep,i,j,k,tirrq(k),gro(k,1),gro(k,2),
cdiag. gro(k,3),gro(k,4),P_tend(k,8),P_tend(k,9)
cdiag enddo
cdiag endif

      call obio_carbon(vrbos,kmax,i,j)

cdiag if(vrbos) then
cdiag do k=1,kdm
cdiag  write(*,*)'obio_ptend14: ',
cdiag.  nstep,i,j,k,P_tend(k,8),P_tend(k,9)
cdiag enddo
cdiag endif

 
cdiag if (vrbos) then
cdiag  do k=1,1
cdiag   write (*,107) k,(bio_var(l),l=1,9),
cdiag.   (bio_var(nt),obio_P(k,nt),P_tend(k,nt),
cdiag.    (rhs(k,nt,l),l=1,16),nt=1,9),
cdiag.   (bio_var(nt),det(k,nt-9),D_tend(k,nt-9),
cdiag.    (rhs(k,nt,l),l=1,6),nt=10,12)
cdiag
cdiag   write (*,107) k,(bio_var(l),l=7,12),
cdiag.   (bio_var(nt),obio_P(k,nt),P_tend(k,nt),
cdiag.    (rhs(k,nt,l),l=9,16),nt=1,9),
cdiag.   (bio_var(nt),det(k,nt-9),D_tend(k,nt-9),
cdiag.    (rhs(k,nt,l),l=7,12),nt=10,12)
cdiag  end do
cdiag end if
 107  format (/'lyr',i3,4x,'amount   tndcy   ',
     .    9(2x,a7)/(a7,2es9.1,2x,6es9.1))

!compute here rates but all detritus update done inside the update routine

c Sinking rate temperature (viscosity) dependence (also convert to /hr)
      do k = 1,kmax
        viscfac(k) = 0.451 + 0.0178*temp1d(k)
      enddo
      do nt = 1,nchl
       do k = 1,kmax
         obio_ws(k,nt) = obio_wsh(nt)*viscfac(k)*pnoice
       enddo
      enddo
      nt = 4
      do k = 1,kmax
        obio_ws(k,nt) = wshc(k)*viscfac(k)*pnoice
      enddo
      do nt = 1,nchl
        obio_ws(kmax+1,nt)=obio_ws(kmax,nt)
      enddo
      do nt = 1,ndet
       do k = 1,kmax
         wsdet(k,nt) = wsdeth(nt)*viscfac(k)*pnoice
       enddo
      enddo

cdiag  if(vrbos .and. diagno)then
cdiag   do k=1,kmax
cdiag    write(905,'(2i5,11e12.4)')
cdiag.     nstep,k,temp1d(k),viscfac(k),pnoice,
cdiag.             (obio_wsh(nt),nt=1,nchl),(obio_ws(k,nt),nt=1,nchl)
cdiag      write(906,'(2i5,4e12.4)')
cdiag.       nstep,k,wshc(k),viscfac(k),pnoice,obio_ws(k,4)
cdiag      write(*,'(a,2i5,8e12.4)')'obio_ptend, wsdet ',
cdiag.       nstep,k,(wsdeth(nt),nt=1,ndet),viscfac(k),pnoice,
cdiag.              (wsdet(k,nt),nt=1,ndet)
cdiag   enddo
cdiag  endif

c Save method for hard boundary condition (no flux)
c      srate = 0.0 - obio_wsh(n)*tracer(i,k-1,m,n)

      return
      end
