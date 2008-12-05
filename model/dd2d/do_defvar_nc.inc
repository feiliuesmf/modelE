      type(dd2d_grid) :: grid
      integer :: fid
      character(len=*) :: varinfo
      integer :: nd,shp(7),rc
      nd = size(shape(arr))
      shp(1:nd) = shape(arr)
      if(grid%am_i_globalroot) then
      call define_var(fid,dtype,varinfo,nd,shp,
     &     grid%npx,grid%npy,grid%ntiles,rc)
      endif
      call stoprc(rc,0)
