#include "rundeck_opts.h"

      module CMP
      implicit none
      private

      public store, guess_dims, do_compare

      interface store
        module procedure store_r, store_i
      end interface

      integer, parameter :: TYPE_DOUBLE=0, TYPE_INT=1
      integer, parameter :: max_num=140


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
      db(i)%lm = dims(4)

      return
      end subroutine guess_dims


      subroutine do_compare
      real*8, parameter :: EPS=1.d-36
      integer m, nerr, n
      real*8 err, rel_err, abs_val, v1,v2
      integer i,j,k,l,nn

      do m=1,num
        nerr = 0
        print '(" dims=  ",i6,3i4,"  name=  ",a16,"     tot. points= ",i7)', &
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
          err = v2 -v1
          abs_val = .5*( abs(v1) + abs(v2) )
          rel_err = 0.d0
          if ( abs_val > 0.d0 ) rel_err = abs(err)/abs_val
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
            if ( nerr > 10 ) exit
          endif
        enddo
      enddo

      end subroutine do_compare

      end module CMP


#define check(y,x) call store(i,y,pack(x,tt),size(x)); \
                   call guess_dims(i,y,shape(x))

      program compare
      use CMP
      use filemanager
!ccc module to allocate dynamic arrays
!AOO use statements added for domain_decomp and dynamics to pull in
!AOO dynamically allocated arrays
      use domain_decomp, only : init_app, grid, finish_app
      use model_com, only : ioread
      use model_com, only : im,jm
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
      use ghy_com, only : wbare,wvege,htbare,htvege,snowbv, &
        nsn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij
      use veg_com, only : Cint,Qfol,cnc_ij
      use landice_com, only : snowli,tlandi
      use pblcom, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg, &
           ustar_pbl,egcm,w2gcm,tgvavg,qgavg
      use pblcom, only : uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs,ipbl
      use clouds_com, only : ttold,qtold,svlhx,rhsav,cldsav,airx,lmc
      use somtq_com, only : tmom,qmom
      use rad_com, only : tchg,rqt,kliq,  s0,srhr,trhr,fsf, &
           fsrdir,srvissurf,srdn,cfrac,rcld,salb,O3_rad_save, &
           O3_trac=>O3_tracer_save
      use icedyn_com, only : rsix,rsiy,usi,vsi,icij
      use icedyn, only : imic
      use diag_com, only : keynr,tsfrez,tdiurn,oa
      use diag_com, only : aj,areg,apj,ajl,asjl,aij,ail,energy,consrv &
           ,speca,atpe,adiurn,wave,ajk,aijk,aisccp,hdiurn
      use model_com, only : idacc

!ccc  include tracers data here
#ifdef TRACERS_ON
      use pblcom, only : trabl
      use tracer_com, only : trm,trmom
      use trdiag_com, only: taijln,taijn,taijs,tajln,tajls,tconsrv

#  ifdef TRACERS_WATER
      use lakes_com, only : trlake
      use seaice_com, only : trsi
      use ghy_com, only : tr_wbare, tr_wvege, tr_wsn_ij
      use landice_com, only : trsnowli,trlndi
      use tracer_com, only : trwm
      use icedyn_com, only : ticij
#  endif

#  ifdef TRACERS_SPECIAL_Shindell
      use TRCHEM_Shindell_COM, only : &
           yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2 &
           ,yAldehyde,yXO2N,yRXPAR,corrOx,ss &
#  ifdef SHINDELL_STRAT_CHEM
           ,SF3,pClOx,pClx,pOClOx,pBrOx
#  endif
#  endif

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

      implicit none
      logical, parameter :: tt=.true.
      integer, external :: iargc
      integer Itime(2)
      integer ioerr
      character*120 file_name(2)
      integer fd,i

      IF(IARGC().NE.2) then
         print *,"CMPE002 compares data in two restart files"
         print *,"Usage: CMPE002 file_1 file_2"
         call stop_model("Incorrect arguments",255)
      endif

!AOO added calls ti init routines for dynamically allocated arrays.
      call init_app(grid,im,jm)
      call alloc_drv()

!C****
!C**** Read ReStartFiles
!C****
      call getarg ( 1, file_name(1) )
      call getarg ( 2, file_name(2) )

      do i=1,2
        call openunit( file_name(i), fd, .true., .true. )
        call io_rsf( fd, Itime(i), ioread, ioerr )
        call closeunit( fd )
        if ( ioerr == 1 ) then
           print *, 'There was an error while reading input file.'
           print *, 'You are probably using incompatible version'
           print *, 'of CMPE002. Try to recompile it.'
           ! stop
        endif
        print *,"read file: ", trim(file_name(i)),"  time= ",Itime(i)
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
        check("wbare",wbare)
        check("wvege",wvege)
        check("htbare",htbare)
        check("htvege",htvege)
        check("snowbv",snowbv)
        check("Cint",Cint)
        check("Qfol",Qfol)
        check("cnc_ij",cnc_ij)
        check("nsn_ij",nsn_ij)
        check("dzsn_ij",dzsn_ij)
        check("wsn_ij",wsn_ij)
        check("hsn_ij",hsn_ij)
        check("fr_snow_ij",fr_snow_ij)
        ! land ice
        check("snowli",snowli)
        check("tlandi",tlandi)
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
        check("O3_rad_save",O3_rad_save)
        check("O3_trac",O3_trac)
        ! icedyn
        if(imic.gt.0) then
          check("RSIX",RSIX)
          check("RSIY",RSIY)
          check("USI",USI)
          check("VSI",VSI)
        end if

        ! diagnostics from diag_com
        check("aj",aj)
        check("areg",areg)
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
        check("hdiurn",hdiurn)

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

!ccc    compare tracers data here
#ifdef TRACERS_ON
        check("TRM",TRM)
        check("TRmom",TRmom)
        check("taijln",taijln)
        check("taijn",taijn)
        check("taijs",taijs)
        check("tajln",tajln)
        check("tajls",tajls)
        check("tconsrv",tconsrv)

#  ifdef TRACERS_WATER
        check("trlake",trlake)
        check("trsi",trsi)
        check("tr_wbare",tr_wbare)
        check("tr_wvege",tr_wvege)
        check("tr_wsn_ij",tr_wsn_ij)
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
        check("pClOx",pClOx)
        check("pClx",pClx)
        check("pOClOx",pOClOx)
        check("pBrOx",pBrOx)
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

