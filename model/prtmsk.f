      subroutine prtmsk(mask,array,work,idm,ii,jj,offset,scale,title)
c
c --- Delete 'array' elements outside 'mask'. Then
c --- break 'array' into sections, each 'nchar' characters wide, for printing.
c
      common/linepr/lp
c
      character title*(*)
      real array(idm,1),work(idm,1)
      integer mask(idm,1)
      data nchar/76/
ccc   data nchar/80/
ccc   data nchar/132/
c
      cvmgp(a,b,c)=a*(.5+sign(.5,c))+b*(.5-sign(.5,c))
      cvmgz(a,b,ic)=cvmgp(a,b,-1.*iabs(ic))
c
      ncols=nchar/4
      do 1 n=1,jj/ncols+1
      j1=ncols*(n-1)+1
      j2=min0(ncols*n,jj)
      if (j1.gt.j2) go to 1
      write (lp,'(/'' Sec.'',i2,'' (cols'',i4,'' -'',i4,'') -- '',a)')
     .   n,j1,j2,title
ccc      if (j2.lt.j1+5) then
ccc      write (lp,'('' (Not printed. Too few columns. Save paper.)'')')
ccc      go to 1
ccc      end if
      do 2 i=1,ii
      do 3 j=j1,j2
 3    work(i,j)=cvmgz(0.,array(i,j),mask(i,j))
      do 4 j=j1,j2
 4    work(i,j)=cvmgz(0.,(work(i,j)-offset)*scale,mask(i,j))
ccc 2    write (lp,'(32i4)') (int(work(i,j)),j=j1,j2)
 2    write (lp,'(32i4)') i,(int(work(i,j)),j=j1,j2)
 1    continue
      return
      end
