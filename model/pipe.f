      subroutine pipe_init(master)
c
c --- this set of routines facilitates output comparison from two micom
c --- versions running side by side. one model, the 'slave', writes its
c --- output into a named pipe. the other model, the 'master', reads from
c --- the pipe and compares. differences are recorded in 'base.out'.
c
c --- call 'pipe_init' initially from both code versions undergoing testing.
c --- one version must set master=.true., the other must set master=.false.
c
      logical master,slave
c
      common/cmp_pipe/iunit,lpunit,slave
      iunit=39
      lpunit=38
      slave=.not.master
c
c --- open the pipe and some output files
c
      open (unit=iunit,file='cmp_pipe',status='unknown',
     .   form='unformatted')
      if (master) then
        open (unit=lpunit,file='base.out',status='unknown')
      else
        open (unit=lpunit,file='test.out',status='unknown')
      end if
c
      return
      end
c
c
      subroutine compare(field,mask,what)
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
c --- call this routine from anywhere in the code (from both versions, of
c --- course) to check whether data stored in 'field' are identical
c
      real field(idm,jdm),field1(idm,jdm)
      integer mask(idm,jdm)
      character*20 what,which
      logical slave,fail
      common/cmp_pipe/iunit,lpunit,slave
c
      if (nstep.le.450) return		!  don't start right away
c
      if (slave) then
      write (lpunit,'(2a)') 'writing for comparison: ',what
      write (lp    ,'(2a)') 'writing for comparison: ',what
      write (iunit) what,field
c
      else				!  slave = .false.
c
      read (iunit) which,field1
      write (lpunit,'(2a)') 'reading for comparison: ',which
      write (lp    ,'(2a)') 'reading for comparison: ',which
      if (what.ne.which) then
        write (lpunit,'(4a)') 'out of sync -- trying to compare ',what,
     .     '  to  ',which
        stop
      end if
c
      fail=.false.
      do 1 j=1,jdm
      do 1 i=1,idm
      if (mask(i,j).gt.0 .and. field(i,j).ne.field1(i,j)) then
        write (lpunit,'(a,2i5,1p,3(a,e14.7))') 'i,j=',i,j,
     .    '  master:',field(i,j),'  slave:',field1(i,j),
     .    '  diff:',field(i,j)-field1(i,j)
        fail=.true.
        if (fail) return		!  optional
      end if
 1    continue
      if (fail) stop			!  optional
c
      end if
      return
      end
c
c
      subroutine comparall(m,n,mm,nn,info)
c
c --- write out a standard menu of arrays for testing
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      common/cmp_pipe/iunit,lpunit,slave
      character text*20,info*(*)
c
      if (nstep.le.450) return		!  don't start right away
c
      write (lpunit,'(2a)') 'comparall called: ',info
      do 1 k=1,kk
      km=k+mm
      kn=k+nn
 100  format (a9,i3,a8)
      write (text,100) 'u(kn)  k=',k,info(13:20)
      call compare(u(1,1,kn),iu,text)
      write (text,100) 'v(kn)  k=',k,info(13:20)
      call compare(v(1,1,kn),iv,text)
      write (text,100) 'dp(kn) k=',k,info(13:20)
      call compare(dp(1,1,kn),ip,text)
      write (text,100) 'temp(kn) ',k,info(13:20)
      call compare(temp(1,1,kn),ip,text)
      write (text,100) 'saln(kn) ',k,info(13:20)
      call compare(saln(1,1,kn),ip,text)
      write (text,100) 'th3d(kn) ',k,info(13:20)
      call compare(th3d(1,1,kn),ip,text)
 1    continue
c
      return
      end
