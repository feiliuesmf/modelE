c compilation:
c f90 -c -64 -mips4 -O2 acc2nc.f
c f90 -64 -mips4 -O2 MODEL_COM.o DAGCOM.o RADNCB.o NCACC.o DEFACC.o
c             GEOM_B.o acc2nc.o /usr/local/netcdf-3.4/lib64/libnetcdf.a
c f90 -64 -mips4 -O2 acc2nc.o ../model/gcmlib.a
c             /usr/local/netcdf-3.4/lib64/libnetcdf.a -o acc2nc
      program acc2nc
      use MODEL_COM
      use DAGCOM
      use GEOM
      implicit none
      character(len=80) :: accfile,outfile,cmdstr,stopstr
      INTEGER :: IARGC
      INTEGER :: IUNIT
      INTEGER :: K,L
      LOGICAL :: EX
      INTEGER :: IT
      REAL :: DUMMY4
      CHARACTER(len=8) :: MODULE_HEADER
      integer :: k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,
     &           k11,k12,k13,k14,k15,k16,k17,k18,k19   
      real ::
     &     aj_odd,areg_odd,apj_odd,ajl_odd,asjl_odd,aij_odd,
     &     ail_odd,aijg_odd,energy_odd,consrv_odd,speca_odd,
     &     atpe_odd,adaily_odd,wave_odd,ajk_odd,aijk_odd,
     &     aijl_odd,ajlsp_odd,tsfrez_odd

      call getarg(0,cmdstr)
      if(iargc().ne.2) then
         stopstr=trim(cmdstr)//' infile outfile'
         stop trim(stopstr)
      endif
      call getarg(1,accfile)
      inquire(file=accfile,exist=ex)
      if(.not.ex) then
         stopstr='file does not exist: '//trim(accfile)
         stop trim(stopstr)
      endif
      call getarg(2,outfile)
C****
C**** READ THE SINGLE-PRECISION ACC FILE
C**** the ugly code is necessary bc no ioread_single option,
C**** and the acc arrays no longer stored adjacent in a common block


      k01 = size(aj)
      k02 = size(areg)
      k03 = size(apj)
      k04 = size(ajl)
      k05 = size(asjl)
      k06 = size(aij)
      k07 = size(ail)
      k08 = size(aijg)
      k09 = size(energy)
      k10 = size(consrv)
      k11 = size(speca)
      k12 = size(atpe)
      k13 = size(adaily)
      k14 = size(wave)
      k15 = size(ajk)
      k16 = size(aijk)
      k17 = size(aijl)
      k18 = size(ajlsp)
      k19 = size(tsfrez)

      iunit=30
      open(iunit,file=accfile,form='unformatted')
      read(iunit,err=840,end=850) it,jc,
     &     xlabel,namd6,amon,amon0,rc
      read(iunit) ! temporary blank record
      read(iunit,err=840,end=850) module_header,keynr,
     &     (tsfrez(k,1,1) ,k=1,k19/2) ,(tsfrez_odd ,k=1,mod(k19,2)),
     &     (aj(k,1,1)     ,k=1,k01/2) ,(aj_odd     ,k=1,mod(k01,2)),
     &     (areg(k,1)     ,k=1,k02/2) ,(areg_odd   ,k=1,mod(k02,2)),
     &     (apj(k,1)      ,k=1,k03/2) ,(apj_odd    ,k=1,mod(k03,2)),
     &     (ajl(k,1,1)    ,k=1,k04/2) ,(ajl_odd    ,k=1,mod(k04,2)),
     &     (asjl(k,1,1)   ,k=1,k05/2) ,(asjl_odd   ,k=1,mod(k05,2)),
     &     (aij(k,1,1)    ,k=1,k06/2) ,(aij_odd    ,k=1,mod(k06,2)),
     &     (ail(k,1,1)    ,k=1,k07/2) ,(ail_odd    ,k=1,mod(k07,2)),
     &     (aijg(k,1,1)   ,k=1,k08/2) ,(aijg_odd   ,k=1,mod(k08,2)),
     &     (energy(k,1)   ,k=1,k09/2) ,(energy_odd ,k=1,mod(k09,2)),
     &     (consrv(k,1)   ,k=1,k10/2) ,(consrv_odd ,k=1,mod(k10,2)),
     &     (speca(k,1,1)  ,k=1,k11/2) ,(speca_odd  ,k=1,mod(k11,2)),
     &     (atpe(k,1)     ,k=1,k12/2) ,(atpe_odd   ,k=1,mod(k12,2)),
     &     (adaily(k,1,1) ,k=1,k13/2) ,(adaily_odd ,k=1,mod(k13,2)),
     &     (wave(k,1,1,1) ,k=1,k14/2) ,(wave_odd   ,k=1,mod(k14,2)),
     &     (ajk(k,1,1)    ,k=1,k15/2) ,(ajk_odd    ,k=1,mod(k15,2)),
     &     (aijk(k,1,1,1) ,k=1,k16/2) ,(aijk_odd   ,k=1,mod(k16,2)),
     &     (aijl(k,1,1,1) ,k=1,k17/2) ,(aijl_odd   ,k=1,mod(k17,2)),
     &     (ajlsp(k,1,1,1),k=1,k18/2) ,(ajlsp_odd  ,k=1,mod(k18,2)),
     &     it
      close(iunit)

C**** Convert ACC single to double precision
      call CNV428(aj      ,aj_odd      ,k01)
      call CNV428(areg    ,areg_odd    ,k02)
      call CNV428(apj     ,apj_odd     ,k03)
      call CNV428(ajl     ,ajl_odd     ,k04)
      call CNV428(asjl    ,asjl_odd    ,k05)
      call CNV428(aij     ,aij_odd     ,k06)
      call CNV428(ail     ,ail_odd     ,k07)
      call CNV428(aijg    ,aijg_odd    ,k08)
      call CNV428(energy  ,energy_odd  ,k09)
      call CNV428(consrv  ,consrv_odd  ,k10)
      call CNV428(speca   ,speca_odd   ,k11)
      call CNV428(atpe    ,atpe_odd    ,k12)
      call CNV428(adaily  ,adaily_odd  ,k13)
      call CNV428(wave    ,wave_odd    ,k14)
      call CNV428(ajk     ,ajk_odd     ,k15)
      call CNV428(aijk    ,aijk_odd    ,k16)
      call CNV428(aijl    ,aijl_odd    ,k17)
      call CNV428(ajlsp   ,ajlsp_odd   ,k18)
      CALL CNV428(tsfrez  ,tsfrez_odd  ,k19)

      call geom_b
      call init_DIAG
      call write_nc_acc(outfile)

      stop
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP OR PARAMETER CONFLICTS
C****
  840 WRITE(6,*) 'READ ERROR ON RESTART FILE TAPE'
      IF(IM.NE.JC(1).OR.JM.NE.JC(2).OR.LM.NE.JC(3)) then
         WRITE(6,*) 'IM,JM,LM mismatch: ',IM,JM,LM,' vs ',JC(1:3)
      ENDIF
      IF(KACC.NE.JC(5)) WRITE(6,*) 'KACC mismatch: ',KACC,' vs ',JC(5)
      STOP 3
  850 WRITE(6,'('' END OF TAPE REACHED...LAST TAU ON TAPE '',I6)')IT
      STOP 3

      end program acc2nc

      subroutine cnv428(arr,arr_odd,n)
      implicit none
      integer :: n
      double precision, dimension(n) :: arr
      real*4 :: arr_odd
      integer :: i
      if(mod(n,2).eq.1) arr(n)=dble(arr_odd)
      do i=n/2,1,-1
         arr(2*i-1:2*i)=dble(transfer(arr(i),arr_odd,size=2))
      enddo
      return
      end subroutine cnv428
