#include "rundeck_opts.h"

      module CMP
      implicit none
      private

      public store, guess_dims, do_compare

      interface store
        module procedure store_r, store_i
      end interface

      integer, parameter :: TYPE_DOUBLE=0, TYPE_INT=1
      integer, parameter :: max_num=170
      real*8, public :: EPS=1.d-36
      integer, public :: max_err=10

      type db_str
        character*32 name
        integer type
        real*8, dimension(:), pointer :: buf1, buf2
        integer, dimension(:), pointer :: ibuf1, ibuf2
        integer im,jm,km,lm,n
      end type db_str

      integer, save :: num=0
      type (db_str) db(max_num)

      contains


      subroutine store_r(ind, name, a, n)
      implicit none
      integer ind,  n
      character*(*) name
      real*8 a(:)
      ! local
      integer i

      if (ind==2) then
        do i=1,num 
          if ( db(i)%name == name ) exit
        enddo
        if (i>num) call stop_model("name not present in db",255)
        db(i)%buf2(1:n) = a(1:n)
        return
      endif

      do i=1,num 
        if ( db(i)%name == name ) then
          print *,"ERROR: You are trying to check array twice: ",name
          call stop_model("Do not check the same array more than once!",255)
        endif
      enddo

      num=num+1
      if (num>max_num) call stop_model("too many arrays",255)
      db(num)%name=name
      db(num)%type=TYPE_DOUBLE
      db(num)%im=1
      db(num)%jm=1
      db(num)%km=1
      db(num)%lm=1
      db(num)%n=n
      allocate( db(num)%buf1(n) )
      allocate( db(num)%buf2(n) )
      db(num)%buf1(1:n) = a(1:n)

      return
      end subroutine store_r

      subroutine store_i(ind, name, a, n)
      implicit none
      integer ind,  n
      character*(*) name
      integer a(:)
      ! local
      integer i

      if (ind==2) then
        do i=1,num 
          if ( db(i)%name == name ) exit
        enddo
        if (i>num) call stop_model("name not present in db",255)
        db(i)%ibuf2(1:n) = a(1:n)
        return
      endif

      num=num+1
      if (num>max_num) call stop_model("too many arrays",255)
      db(num)%name=name
      db(num)%type=TYPE_INT
      db(num)%im=1
      db(num)%jm=1
      db(num)%km=1
      db(num)%lm=1
      db(num)%n=n
      allocate( db(num)%ibuf1(n) )
      allocate( db(num)%ibuf2(n) )
      db(num)%ibuf1(1:n) = a(1:n)

      return
      end subroutine store_i


      subroutine guess_dims(ind, name, dims )
      implicit none
      character*(*) name
      integer ind, dims(:)
      ! local
      integer n,i
!!!     write(0,*) 'checking ',name
      do i=1,num 
        if ( db(i)%name == name ) exit
      enddo
      if (i>num) call stop_model("rs: name not present in db",255)
      n = db(i)%n
      db(i)%im = dims(1)
      n = n/dims(1)
      if ( n<=1 ) return
      db(i)%jm = dims(2)
      n = n/dims(2)
      if ( n<=1 ) return
      db(i)%km = dims(3)
      n = n/dims(3)
      if ( n<=1 ) return
      db(i)%lm = n ! all other dimensions combined

      return
      end subroutine guess_dims

!The detection of NaN''s is based on the following info from
!Compaq Fortran documentation:
!(http://h18009.www1.hp.com/
!     fortran/docs/unix-um/dfumdatarep.htm#sec_fltng_pt_exc)
!
!Infinity(+)Z'7FF0000000000000'
!Infinity(--)Z'FFF0000000000000'
!Zero(+0)Z'0000000000000000'
!Zero(-0)Z'8000000000000000'
!QuietNaN(+)FromZ'7FF8000000000000'toZ'7FFFFFFFFFFFFFFF'
!QuietNaN(--)FromZ'FFF8000000000000'toZ'FFFFFFFFFFFFFFFF'
!SignalingNaN(+)FromZ'7FF0000000000001'toZ'7FF7FFFFFFFFFFFF'
!SignalingNaN(--)FromZ'FFF0000000000001'toZ'FFF7FFFFFFFFFFFF'

      function isNaN(a)
      real*8 a
      logical isNaN
      real*8 b
      integer*4 c(2),x
      equivalence (c(1),b)
      integer*4, parameter :: mask=Z'FFF00000'
      integer*4, parameter :: nan_1=Z'FFF00000', nan_2=Z'7FF00000'

      isNaN = .false.
      b = a
#ifdef CONVERT_BIGENDIAN
      ! for little-endian machines
      x = iand( c(2), mask )
      if ( x==nan_1 .or. x==nan_2 ) isNaN = .true.
#else
      ! for big-endian machines
      x = iand( c(1), mask )
      if ( x==nan_1 .or. x==nan_2 ) isNaN = .true.
#endif
      return
      end function isNaN


      subroutine do_compare
!     real*8, parameter :: EPS=1.d-36
      integer m, nerr, n
      real*8 err, rel_err, abs_val
      integer i,j,k,l,nn
      real*8 v1,v2
      integer*4 iv1(2), iv2(2)
      equivalence (v1,iv1(1)), (v2,iv2(1))


      do m=1,num 
        nerr = 0
        print '(" dims=  ",i6,3i4,"  name=  ",a16,"     tot. points=",i8)', &
           db(m)%im, db(m)%jm, db(m)%km,  db(m)%lm, trim(db(m)%name), &
           db(m)%n
        do n=1,db(m)%n
          select case( db(m)%type )
          case (TYPE_DOUBLE)
            v1 = db(m)%buf2(n)
            v2 = db(m)%buf1(n)
          case (TYPE_INT)
            v1 = db(m)%ibuf2(n)
            v2 = db(m)%ibuf1(n)
          case default
            call stop_model("wrong type in DB",255)
          end select
          !if values are binary identical skip the rest of the test
          !this should work also for NaN's
          if ( iv1(1) == iv2(1) .and. iv1(2) == iv2(2) ) cycle
          if ( isNaN(v1) .or. isNaN(v2) ) then
            !NaN
            err = 0 ; rel_err = 1.d30
            if ( isNaN(v1) ) err = err + 10
            if ( isNaN(v2) ) err = err +  2
          else
            !normal number
            err = v2 -v1
            abs_val = .5*( abs(v1) + abs(v2) )
            rel_err = 0.d0
            if ( abs_val > 0.d0 ) rel_err = abs(err)/abs_val
          endif
          if (rel_err>EPS) then
            nn = n-1
            l = nn/(db(m)%im*db(m)%jm*db(m)%km)
            nn = mod( nn, db(m)%im*db(m)%jm*db(m)%km )
            k = nn/(db(m)%im*db(m)%jm)
            nn = mod( nn, db(m)%im*db(m)%jm )
            j = nn/(db(m)%im)
            i = mod( nn, db(m)%im )
            print '(i6,": ",i6,3i4,"    ",4e24.16)', &
               n, i+1,j+1,k+1,l+1, v1,v2,err,rel_err
            nerr = nerr+1
            if ( max_err > 0 .and. nerr > max_err ) exit
          endif
        enddo
      enddo

      end subroutine do_compare

      end module CMP

#ifdef USE_ENT
      module ent_data_for_cmp

      private
      public ent_data, get_ent_data_for_cmp
      real*8, pointer :: ent_data(:)

      contains

      subroutine get_ent_data_for_cmp
      use ent_com, only : entcells
      use ent_mod

      call ent_cell_pack( ent_data, entcells )

      end subroutine get_ent_data_for_cmp

      end module ent_data_for_cmp
#endif

#define check(y,x) call store(i,y,pack(x,tt),size(x)); \
                   call guess_dims(i,y,shape(x))

#define checks(y,x) call store(i,y,(/x/),1); \
                   call guess_dims(i,y,(/1/))

      program compare
      use CMP
      use filemanager
!ccc module to allocate dynamic arrays
!AOO use statements added for domain_decomp and dynamics to pull in
!AOO dynamically allocated arrays
      use domain_decomp, only : init_app, grid, finish_app
      use model_com, only : ioread,ioread_nodiag
      use model_com, only : im,jm,lm
!ccc  modules with data to compare
      use model_com, only : u,v,t,q,p
#ifdef CHECK_OCEAN
      use ocean, only : mo,uo,vo,g0m,gxmo,gymo,gzmo,s0m,sxmo,symo,szmo &
           ,ogeoz,ogeoz_sv
      use straits, only : must,g0mst,gxmst,gzmst,s0mst,sxmst,szmst &
           ,rsist,rsixst,msist,hsist,ssist
#endif
      use lakes_com, only : mldlk,mwl,tlake,gml
      use seaice_com, only : rsi,hsi,snowi,msi,ssi,pond_melt,flag_dsws
      use ghy_com, only : snowe,tearth,wearth,aiearth,snoage &
           ,evap_max_ij,fr_sat_ij,qg_ij
      use ghy_com, only : w_ij,ht_ij,snowbv, &
        nsn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij
#ifndef USE_ENT
      use veg_com, only : Cint,Qfol,cnc_ij
#endif
      use landice_com, only : snowli,tlandi,MdwnImp,EdwnImp
      use landice, only : accpda,accpdg, eaccpda,eaccpdg
      use pblcom, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg, &
           ustar_pbl,egcm,w2gcm,tgvavg,qgavg
      use pblcom, only : uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs,ipbl
      use clouds_com, only : ttold,qtold,svlhx,rhsav,cldsav,airx,lmc
      use somtq_com, only : tmom,qmom
      use rad_com, only : tchg,rqt,kliq,  s0,srhr,trhr,fsf, &
           fsrdir,srvissurf,srdn,cfrac,rcld,salb
#ifdef TRACERS_SPECIAL_Shindell
      use rad_com, only : O3_trac=>O3_tracer_save, rad_to_chem
#endif
#ifdef TRACERS_DUST
      use rad_com, only : srnflb_save,trnflb_save,ttausv_save,ttausv_cs_save
#endif
      use icedyn_com, only : rsix,rsiy,usi,vsi,icij
      use icedyn, only : imic

#ifndef SKIP_DIAG
      use diag_com, only : keynr,tsfrez,tdiurn,oa
      use diag_com, only : aj,aregj,apj,ajl,asjl,aij,ail,energy,consrv &
           ,speca,atpe,adiurn,wave,ajk,aijk,aisccp
#ifndef NO_HDIURN
      use diag_com, only : hdiurn
#endif
      use model_com, only : idacc

#ifdef CHECK_OCEAN
      use odiag, only: oij,oijl,ol,olnst
#endif
#endif

!ccc  include tracers data here
#ifdef TRACERS_ON
      use pblcom, only : trabl
      use tracer_com, only : trm,trmom
#ifndef SKIP_DIAG
      use trdiag_com, only: taijln,taijn,taijs,tajln,tajls,tconsrv
#endif

#  ifdef TRACERS_WATER
      use lakes_com, only : trlake
      use seaice_com, only : trsi
      use ghy_com, only : tr_w_ij, tr_wsn_ij, trsnowbv0
      use landice_com, only : trsnowli,trlndi
      use tracer_com, only : trwm
      use icedyn_com, only : ticij
#  endif

#ifdef TRACERS_SPECIAL_Shindell
      use TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3 &
       ,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,ss,corrOx,JPPJ
#ifdef SHINDELL_STRAT_CHEM
      use TRCHEM_Shindell_COM, only: &
       SF3,SF2,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      use TRACER_SOURCES, only: day_ncep,DRA_ch4,sum_ncep,PRS_ch4, &
       HRA_ch4,iday_ncep,i0_ncep,iHch4,iDch4,i0ch4,first_ncep,first_mod &
       ,max_days,nra_ncep,nra_ch4,maxHR_ch4
#endif
#endif

#  ifdef CHECK_OCEAN
#    ifdef TRACERS_WATER
       use straits, only : trsist
#      ifdef TRACERS_OCEAN
         use ocean, only : trmo,txmo,tymo,tzmo
         use straits, only : trmst,txmst,tzmst
#      endif
#    endif
#  endif

#endif
#ifdef USE_ENT
      use ent_data_for_cmp, only : ent_data, get_ent_data_for_cmp
#endif

      implicit none
      logical, parameter :: tt=.true.
      integer, external :: iargc
      integer Itime(2)
      integer ioerr
      character*120 file_name(2)
      integer fd,i

      IF(IARGC().le.1) then
         print *,"CMPE002 compares data in two restart files"
         print *,"Usage: CMPE002 file_1 file_2 tolerance(1.d-36) #_shown(10,0=all)"
         call stop_model("Incorrect arguments",255)
      endif

!AOO added calls ti init routines for dynamically allocated arrays.
      call init_app(grid,im,jm,lm)
      call alloc_drv()

!C****
!C**** Read ReStartFiles
!C****
      if(IARGC().ge.3) then
        call getarg ( 3, file_name(1) )
        read(file_name(1),*) EPS
      end if
      if(IARGC().ge.4) then
        call getarg ( 4, file_name(1) )
        read(file_name(1),*) max_err
      end if
      write(*,*) 'listing differences >',eps
      call getarg ( 1, file_name(1) )
      call getarg ( 2, file_name(2) )

      do i=1,2
        call openunit( file_name(i), fd, .true., .true. )
#ifndef SKIP_DIAG
        call io_rsf( fd, Itime(i), ioread, ioerr )
#else
        call io_rsf( fd, Itime(i), ioread_nodiag, ioerr )
#endif
        call closeunit( fd )
        if ( ioerr == 1 ) then
           print *, 'There was an error while reading input file.'
           print *, 'You are probably using incompatible version'
           print *, 'of CMPE002. Try to recompile it.'
           ! stop
        endif
        print *,"read file: ", trim(file_name(i)),"  time= ",Itime(i)
#ifdef USE_ENT
        ! get Ent data
        call get_ent_data_for_cmp
#endif
        ! data from model_com
        check("u",u)
        check("v",v)
        check("t",t)
        check("q",q)
        check("p",p)
        ! strat
        check("airx",airx)
        check("lmc",lmc)
        ! ocean
#ifdef CHECK_OCEAN
        check("mo",mo)
        check("uo",uo)
        check("vo",vo)
        check("g0m",g0m)
        check("gxmo",gxmo)
        check("gymo",gymo)
        check("gzmo",gzmo)
        check("s0m",s0m)
        check("sxmo",sxmo)
        check("symo",symo)
        check("szmo",szmo)
        check("ogeoz",ogeoz)
        check("ogeoz_sv",ogeoz_sv)
        ! straits
        check("must",must)
        check("g0mst",g0mst)
        check("gxmst",gxmst)
        check("gzmst",gzmst)
        check("s0mst",s0mst)
        check("sxmst",sxmst)
        check("szmst",szmst)
        check("rsist",rsist)
        check("rsixst",rsixst)
        check("msist",msist)
        check("hsist",hsist)
        check("ssist",ssist)
#endif
        ! lakes
        check("mldlk",mldlk)
        check("mwl",mwl)
        check("tlake",tlake)
        check("gml",gml)
        ! sea ice
        check("rsi",rsi)
        check("hsi",hsi)
        check("snowi",snowi)
        check("msi",msi)
        check("ssi",ssi)
        check("pond_melt",pond_melt)
        !check("flag_dsws",flag_dsws)
        ! earth
        check("snowe",snowe)
        check("tearth",tearth)
        check("wearth",wearth)
        check("aiearth",aiearth)
        check("snoage",snoage)
        check("evap_max_ij",evap_max_ij)
        check("fr_sat_ij",fr_sat_ij)
        check("qg_ij",qg_ij)
        ! ground hydrology data from ghy_com
        check("w_ij",w_ij)
        check("ht_ij",ht_ij)
        check("snowbv",snowbv)
#ifdef USE_ENT
        check("ent_data",ent_data)
#else
        check("Cint",Cint)
        check("Qfol",Qfol)
        check("cnc_ij",cnc_ij)
#endif
        check("nsn_ij",nsn_ij)
        check("dzsn_ij",dzsn_ij)
        check("wsn_ij",wsn_ij)
        check("hsn_ij",hsn_ij)
        check("fr_snow_ij",fr_snow_ij)
        ! land ice
        check("snowli",snowli)
        check("tlandi",tlandi)
        check("MdwnImp",MdwnImp)
        check("EdwnImp",EdwnImp)
        !check("accpda",accpda)
        !check("accpdg",accpdg)
        !check("Eaccpda",Eaccpda)
        !check("Eaccpdg",Eaccpdg)
        ! bldat
        check("wsavg",wsavg)
        check("tsavg",tsavg)
        check("qsavg",qsavg)
        check("dclev",dclev)
        check("usavg",usavg)
        check("vsavg",vsavg)
        check("tauavg",tauavg)
        check("ustar_pbl",ustar_pbl)
        check("egcm",egcm)
        check("w2gcm",w2gcm)
        check("tgvavg",tgvavg)
        check("qgavg",qgavg)
        ! pbl data from pblcom
        check("uabl",uabl)
        check("vabl",vabl)
        check("tabl",tabl)
        check("qabl",qabl)
        check("eabl",eabl)
        check("cmgs",cmgs)
        check("chgs",chgs)
        check("cqgs",cqgs)
        check("ipbl",ipbl)
        ! clouds
        check("ttold",ttold)
        check("qtold",qtold)
        check("svlhx",svlhx)
        check("rhsav",rhsav)
        check("cldsav",cldsav)
        ! somtq
        check("tmom",tmom)
        check("qmom",qmom)
        ! radiation
        check("Tchg",Tchg)
        check("rqt",rqt)
        check("kliq",kliq)
        !check("s0",s0)
        check("srhr",srhr)
        check("trhr",trhr)
        check("fsf",fsf)
        check("fsrdir",fsrdir)
        check("srvissurf",srvissurf)
        check("salb",salb)
        check("srdn",srdn)
        check("cfrac",cfrac)
        check("rcld",rcld)
#ifdef TRACERS_SPECIAL_Shindell
        check("rad_to_chem",rad_to_chem)
        check("O3_trac",O3_trac)
#endif
#ifdef TRACERS_DUST
        check("srnflb_save",srnflb_save)
        check("trnflb_save",trnflb_save)
        check("ttausv_save",ttausv_save)
        check("ttausv_cs_save",ttausv_cs_save)
#endif
        ! icedyn
        if(imic.gt.0) then
          check("RSIX",RSIX)
          check("RSIY",RSIY)
          check("USI",USI)
          check("VSI",VSI)
        end if

#ifndef SKIP_DIAG
        ! diagnostics from diag_com
        check("aj",aj)
        check("aregj",aregj)
        check("apj",apj)
        check("ajl",ajl)
        check("asjl",asjl)
        check("aij",aij)
        check("ail",ail)
        check("energy",energy)
        check("consrv",consrv)
        check("speca",speca)
        check("atpe",atpe)
        check("adiurn",adiurn)
        check("wave",wave)
        check("ajk",ajk)
        check("aijk",aijk)
        check("aisccp",aisccp)
#ifndef NO_HDIURN
        check("hdiurn",hdiurn)
#endif
        ! diags
        check("KEYNR",KEYNR)
        check("TSFREZ",TSFREZ)
        check("idacc",idacc)
        check("TDIURN",TDIURN)
        check("OA",OA)
        !check("it",it)
        ! icdiag
        if(imic.gt.0) then
          check("ICIJ",ICIJ)
        end if

#ifdef CHECK_OCEAN
        check("oij",oij)
        check("oijl",oijl)
        check("ol",ol)
        check("olnst",olnst)
#endif
#endif

!ccc    compare tracers data here
#ifdef TRACERS_ON
        check("TRM",TRM)
        check("TRmom",TRmom)
#ifndef SKIP_DIAG
        check("taijln",taijln)
        check("taijn",taijn)
        check("taijs",taijs)
        check("tajln",tajln)
        check("tajls",tajls)
        check("tconsrv",tconsrv)
#endif

#  ifdef TRACERS_WATER
        check("trlake",trlake)
        check("trsi",trsi)
        check("tr_w_ij",tr_w_ij)
        check("tr_wsn_ij",tr_wsn_ij)
        check("trsnowbv0",trsnowbv0)
        check("trsnowli",trsnowli)
        check("trlndi",trlndi)
        check("trabl",trabl)
        check("trwm",trwm)
        check("ticij",ticij)
#  endif

#  ifdef TRACERS_SPECIAL_Shindell
        check("yNO3",yNO3)
        check("pHOx",pHOx)
        check("pNOx",pNOx)
        check("pOx",pOx)
        check("yCH3O2",yCH3O2)
        check("yC2O3",yC2O3)
        check("yROR",yROR)
        check("yXO2",yXO2)
        check("yAldehyde",yAldehyde)
        check("yXO2N",yXO2N)
        check("yRXPAR",yRXPAR)
        check("corrOx",corrOx)
        check("ss",ss)
#  ifdef SHINDELL_STRAT_CHEM
        check("SF3",SF3)
        check("SF2",SF2)
        check("pClOx",pClOx)
        check("pClx",pClx)
        check("pOClOx",pOClOx)
        check("pBrOx",pBrOx)
        check("yCl2",yCl2)
        check("yCl2O2",yCl2O2)
#  endif
#ifdef INTERACTIVE_WETLANDS_CH4
        check("day_ncep",day_ncep)
        check("DRA_ch4",DRA_ch4)
        check("sum_ncep",sum_ncep)
        check("PRS_ch4",PRS_ch4)
        check("HRA_ch4",HRA_ch4)
        check("iday_ncep",iday_ncep)
        check("i0_ncep",i0_ncep)
        check("iHch4",iHch4)
        check("iDch4",iDch4)
        check("i0ch4",i0ch4)
        check("first_ncep",first_ncep)
        check("first_mod",first_mod)
#  endif
#  endif

#  ifdef CHECK_OCEAN
#    ifdef TRACERS_WATER
       check("trsist",trsist)
#      ifdef TRACERS_OCEAN
         trmo,txmo,tymo,tzmo
         trmst,txmst,tzmst
         check("trmo",trmo)
         check("txmo",txmo)
         check("tymo",tymo)
         check("tzmo",tzmo)
         check("trmst",trmst)
         check("txmst",txmst)
         check("tymo",tymo)
         check("tzmo",tzmo)
#      endif
#    endif
#  endif
#endif

      enddo

      print *," ------     Comparing data     -----"
      call do_compare

!** not sure if this is needed, but just in case ...
      call finish_app()

      end

