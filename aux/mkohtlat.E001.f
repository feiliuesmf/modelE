      parameter (im=72,jm=46)
      real*8 COT(im,jm)
      dimension dxyp(jm),plk(im,jm),po(im,jm),onht(jm)
      character*80 title,runID*10

      call getarg(1,runID)

      read(25)
      read(25)
      read(25) dxyp
      rewind 25

      read(26) title,po
      read(26) title,plk

c     do 10 j=1,jm
c     do 10 i=1,im
c  10 po(i,j)=po(i,j)+plk(i,j)  ! ocn+lakes

      read (31) (COT,i=1,9)

      areag=dxyp(1)+dxyp(jm)
      onht(jm)=cot(1,jm)*im*dxyp(jm)*po(1,jm)
         write(*,*) cot(1,jm),1.e-15*onht(jm),dxyp(jm)
      do 20 j=jm-1,2,-1
      areag=areag+dxyp(j)
      onht(j)=onht(j+1)
      do 20 i=1,im
   20 onht(j)=onht(j)+COT(i,j)*dxyp(j)*po(i,j)
      write(*,*) 'areag=',im*areag

      write(99,*) 'Global Northward Ocean Heat Transport '
      write(99,*) 'latitude'
      write(99,*) '10**15 W'
      write(99,500) 'lat',RunID 
      do 30 j=2,jm
   30 write(99,*) -92+(j-1)*4,1.e-15*onht(j)
      write(99,*) ' '
 500  format (1x,a,2x,a)
      stop
      end












