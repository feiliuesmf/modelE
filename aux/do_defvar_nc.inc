      integer :: fid
      integer, intent(in) :: npx,npy,ntiles
      character(len=*), intent(in), optional :: defby
      character(len=*) :: varinfo
      integer :: nd,shp(7),rc,vid
      nd = size(shape(arr))
      shp(1:nd) = shape(arr)
      call define_var(fid,dtype,varinfo,nd,shp,
     &     npx,npy,ntiles,rc,vid)
      if(present(defby)) then
        rc = nf_put_att_text(fid,vid,'defined_by',len_trim(defby),
     &       defby)
      endif
