      subroutine linout(value,char,index)
c
      common/linepr/lp
c
      parameter (length=77)
      character*1 char,line(length)
      common/lnout/line
c
c --- replace n-th element of array 'line' by character 'char', where
c --- n = 'value' modulo 'length'
c --- index < 0  -- initialize 'line' by blanks before adding 'char'
c --- index > 0  -- output 'line' after adding 'char'
c
      if (index.lt.0) then
        do 1 l=1,length
 1      line(l)=' '
      end if
      if (value.gt.0.) then
        n=int(mod(value,float(length)))
        line(n)=char
      end if
      if (index.gt.0) write (lp,'(''=+='',80a1)') (line(l),l=1,length)
      return
      end
