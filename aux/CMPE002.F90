#include "rundeck_opts.h"

      module CMP
      implicit none
      private

      public store, guess_dims, do_compare

      interface store
        module procedure store_r, store_i
      end interface

      integer, parameter :: TYPE_DOUBLE=0, TYPE_INT=1
      integer, parameter :: max_num=128


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
        print '(" dims=  ",4i4,"  name=  ",a16,"     tot. points= ",i6)', &
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
            print '(i6,": ",4i4,"    ",3e24.16)', &
               n, i+1,j+1,k+1,l+1, v1,v2,err
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
      use model_com, only : ioread
!ccc  modules with data to compare
      use ghycom, only : wbare,wvege,htbare,htvege,snowbv,Cint,Qfol, &
        nsn_ij,isn_ij,dzsn_ij,hsn_ij,fr_snow_ij
      use model_com, only : u,v,t,q,p
      use pblcom, only : uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs
      use dagcom, only : aj, areg, apj, ajl, asjl, aij, ail
!ccc  include tracers data here
#ifdef TRACERS_ON
      use pblcom, only : trabl
      use ghycom, only : tr_wbare, tr_wvege, tr_wsn_ij
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
!C****
!C**** Read ReStartFiles
!C****
      call getarg ( 1, file_name(1) )
      call getarg ( 2, file_name(2) )

      do i=1,2
        call openunit( file_name(i), fd, .true., .true. )
        call io_rsf( fd, Itime(i), ioread, ioerr )
        call closeunit( fd )
        print *,"read file: ", trim(file_name(i)),"  time= ",Itime(i)
        ! data from model_com
        check("u",u)
        check("v",v)
        check("t",t)
        check("q",q)
        check("p",p)
         ! ground hydrology data from ghycom
        check("wbare",wbare)
        check("wvege",wvege)
        check("htbare",htbare)
        check("htvege",htvege)
        check("snowbv",snowbv)
        check("Cint",Cint)
        check("Qfol",Qfol)
        check("nsn_ij",nsn_ij)
        check("isn_ij",isn_ij)
        check("dzsn_ij",dzsn_ij)
        check("hsn_ij",hsn_ij)
        check("fr_snow_ij",fr_snow_ij)
        ! diagnostics from dagcom
        check("aj",aj)
        check("areg",areg)
        check("apj",apj)
        check("ajl",ajl)
        check("asjl",asjl)
        check("aij",aij)
        check("ail",ail)
        ! pbl data from pblcom
        check("uabl",uabl)
        check("vabl",vabl)
        check("tabl",tabl)
        check("qabl",qabl)
        check("eabl",eabl)
        check("cmgs",cmgs)
        check("chgs",chgs)
        check("cqgs",cqgs)

!ccc    compare tracers data here
#ifdef TRACERS_ON
        check("trabl",trabl)
        check("tr_wbare",tr_wbare)
        check("tr_wvege",tr_wvege)
        check("tr_wsn_ij",tr_wsn_ij)
#endif

     enddo

      print *," ------     Comparing data     -----"
      call do_compare

      end
 
