      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: varname
      integer :: vid,rc
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(grid%am_i_globalroot .and. rc.ne.nf_noerr) then
        write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      rc = nfmpi_begin_indep_data(fid)
      if(grid%am_i_globalroot) then
        rc = nfmpi_put_var_double(fid,vid,arr)
      endif
      rc = nfmpi_end_indep_data(fid)
