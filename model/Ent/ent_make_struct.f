      module ent_make_struct

      use ent_types
      use ent_const
      use cohorts
      use patches
      use entcells

      public ent_struct_readcsv
      
      contains

!************************************************************************
      subroutine skipstar (iu_entstruct)
      !* Skip a comment line in an Ent structure csv file.
      integer :: iu_entstruct
      !---
      character :: check

      read(iu_entstruct,*) check
      write(*,*) check
      do while (check.eq.'*') 
        read(iu_entstruct,*) check
      end do
      backspace(iu_entstruct)
      end subroutine skipstar

!************************************************************************


      subroutine read_entcell_struct ( ecp, iu_entstruct )
      type(entcelltype) :: ecp
      integer :: iu_entstruct
      !---
      !real*8 :: stext1,stext2,stext3,stext4,stext5 !soil_texture
      
      !read(iu_entstruct,*) stext1,stext2,stext3,stext4,stext5 !soil_texture
      !ecp%soil_texture(1) = stext1
      !ecp%soil_texture(2) = stext2
      !ecp%soil_texture(3) = stext3
      !ecp%soil_texture(4) = stext4
      !ecp%soil_texture(5) = stext5
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) ecp%soil_texture 
      end subroutine read_entcell_struct

!************************************************************************
      

      subroutine read_patch_struct ( pp,iu_entstruct )
      type(patch),pointer :: pp
      integer :: iu_entstruct
      !---
      integer :: layer
!      real*8 :: age, area
!      real*8 :: tc1,tc2,tc3,tc4,tc5,tc6,tc7,tc8,tc9,tc10,tc11,tc12 !TpoolC
!      real*8 :: tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10,tn11,tn12 !TpoolN

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) pp%age, pp%area, pp%soil_type
      do layer = 1,N_CASA_LAYERS
        call skipstar(iu_entstruct)
        read(iu_entstruct,*) pp%Tpool(CARBON,:,layer)
        read(iu_entstruct,*) pp%Tpool(NITROGEN,:,layer)
      end do
      end subroutine read_patch_struct

!************************************************************************

      subroutine read_cohort_struct ( cop, iu_entstruct )
      type(cohort) :: cop
      integer :: iu_entstruct
      !---
!      real*8 :: nm,Ntot,LAI,LMA,h,crown_dx,crown_dy,dbh,root_d,clump
!      real*8 :: fracroot !skip
!      real*8 :: C_fol,N_fol,C_sw,N_sw,C_hw,N_hw,C_lab,N_lab,
!     &     C_froot,N_froot,C_croot,N_croot

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%pft
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%n,cop%h,cop%crown_dx,cop%crown_dy,
     &     cop%dbh,cop%root_d,cop%clump
      write(*,*) cop%n,cop%h,cop%crown_dx,cop%crown_dy,
     &     cop%dbh,cop%root_d,cop%clump
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%C_fol,cop%N_fol,cop%C_sw,cop%N_sw,
     &     cop%C_hw,cop%N_hw,cop%C_lab,cop%N_lab,cop%
     &     C_froot,cop%N_froot,cop%C_croot,cop%N_croot

      end subroutine read_cohort_struct
!************************************************************************

! left the old program for a reference - should be removed
cddd      subroutine ent_struct_readcsv (IM,JM,I0,I1,J0,J1)
cddd!@sum Read ent vegetation structure from ASCII CSV file
cddd!@sum Order must be:
cddd!@sum First line:  N_CASA_LAYERS
cddd!@sum e, entcells in any order spanning dimension IM,JM
cddd!@sum p, patches in order from oldest to youngest
cddd!@sum c, cohorts in order from tallest to shortest
cddd!@sum
cddd!@sum This ordering is not critical, because the insert routines will
cddd!@sum sort the patches and cohorts, given age and height.
cddd
cddd      use FILEMANAGER
cddd      implicit none
cddd      integer :: IM,JM,I0,I1,J0,J1
cddd      type(entcelltype) :: cells(I0:I1,J0:J1)
cddd      type(entcelltype) :: cellsbuf(I0:I1,J0:J1)
cddd      !----
cddd      integer :: iu_entstruct
cddd      character*80, parameter :: ENTSTRUCT="ENTSTRUCT"
cddd      integer :: N_CASAlayers    !N_CASA_layers
cddd      integer :: i,j
cddd      type(entcelltype),pointer :: ecp
cddd      type(patch),pointer :: pp
cddd      type(cohort),pointer :: cop
cddd      real*8 :: fracroot(N_DEPTH)
cddd      logical :: MORE
cddd      character :: next
cddd
cddd      fracroot(1:N_DEPTH) = 0.d0
cddd      MORE = .true.  !initialize
cddd
cddd!      open(iu_entstruct,file="ent_struct_init.csv",status='unknown')
cddd      call openunit(trim(ENTSTRUCT),iu_entstruct,.false.,.true.)
cddd      read(iu_entstruct,*) N_CASAlayers
cddd      write(*,*) 'N_CASAlayers',N_CASAlayers
cddd      do while(MORE)
cddd        read(iu_entstruct,*) next
cddd        if (next.eq.'$') then !End of data
cddd          MORE = .false.
cddd        else if (next.eq.'*') then !skip comment
cddd        else if (next.eq.'e') then !new entcell
cddd          write(*,*) 'e'
cddd          read(iu_entstruct,*) i,j
cddd          write(*,*) i,j
cddd          call check_ij_bounds(i,j,I0,I1,J0,J1)
cddd          call entcell_construct(ecp)
cddd          call zero_entcell( ecp )
cddd          call read_entcell_struct( ecp,iu_entstruct )
cddd          cells(i,j) = ecp
cddd        else if (next.eq.'p') then !new patch
cddd          write(*,*) 'p'
cddd          call insert_patch(ecp,0.d0,0)
cddd          pp => ecp%youngest
cddd          call read_patch_struct(pp,iu_entstruct)
cddd          write(*,*) 'patch area age soil',pp%area,pp%age,pp%soil_type
cddd        else if (next.eq.'c') then !new cohort
cddd          write(*,*) 'c'
cddd          call insert_cohort(pp,0,
cddd     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
cddd     &         fracroot,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
cddd     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
cddd     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
cddd          cop => pp%shortest
cddd          call read_cohort_struct(cop,iu_entstruct)
cddd          write(*,*) 'cohort pft height',cop%pft,cop%h
cddd        end if
cddd        cells(i,j) = ecp
cddd      end do
cddd
cddd!      write(*,*) 'Got here 0.'
cddd!      do i = I0,I1
cddd!        do j = J0,J1
cddd!          call entcell_print(6,cells(i,j))
cddd!        end do
cddd!      end do
cddd
cddd      write(*,*) 'Summarizing entcells...'
cddd      do i = I0,I1
cddd        do j = J0,J1
cddd          call summarize_entcell(cells(i,j))
cddd          write(*,*) '...after summarize_entcell'
cddd          call entcell_print(6,cells(i,j))
cddd          cellsbuf(i,j) = cells(i,j)
cddd          write(*,*) '...after assign cellsbuf(i,j)'
cddd        end do
cddd      end do
cddd
cddd      write(*,*) 'Printing cellsbuf...'
cddd      do i = I0,I1
cddd         do j = J0,J1
cddd            call entcell_print(6, cellsbuf(i,j))
cddd         end do
cddd      end do 
cddd      
cddd      call closeunit(iu_entstruct)
cddd      write(*,*) 'Done with creating Ent vegetation structure.'
cddd      end subroutine ent_struct_readcsv


      subroutine ent_struct_readcsv (ecp, iu_entstruct)
!@sum Read ent vegetation structure from ASCII CSV file
!@sum Called by do_ent_struct:
!@sum    - do_ent_struct loops through entcells and checks ij bounds
!@sum    - ent_struct_read_csv allocates an entcell and loops 
!@sum         through patches and cohorts
!@sum Order in CSV file must be:
!@sum First line:  N_CASA_LAYERS
!@sum e, entcells in any order spanning dimension IM,JM
!@sum p, patches in order from oldest to youngest
!@sum c, cohorts in order from tallest to shortest
!@sum *, comment line
!@sum $, end of entcell
!@sum #, end of file
!@sum The ordering is not critical, because the CSV file must give the
!@sum ij indices of the entcell, and the insert routines will
!@sum sort the patches and cohorts, given age and height.

      implicit none
      type(entcelltype),pointer :: ecp
      integer :: iu_entstruct      
      !----
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      real*8 :: fracroot(N_DEPTH)
      character :: next
      integer :: counter

      fracroot(1:N_DEPTH) = 0.d0

      call entcell_construct(ecp)
      call zero_entcell( ecp )
      call read_entcell_struct( ecp,iu_entstruct )
      
      counter = 0
      do
        counter = counter + 1
        if ( counter >= 1000) 
     &       call stop_model("ent_readcsv: infinite loop",255)
        read(iu_entstruct,*) next
        if (next.eq.'$') then
           write(*,*) 'Ent of entcell.'
           exit
        else if (next.eq.'*') then !skip comment
        else if (next.eq.'p') then !new patch
          write(*,*) 'p'
          call insert_patch(ecp,0.d0,0)
          pp => ecp%youngest
          call read_patch_struct(pp,iu_entstruct)
          write(*,*) 'patch area age soil',pp%area,pp%age,pp%soil_type
        else if (next.eq.'c') then !new cohort
          write(*,*) 'c'
          call insert_cohort(pp,0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         fracroot,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
          cop => pp%shortest
          call read_cohort_struct(cop,iu_entstruct)
          write(*,*) 'cohort pft height',cop%pft,cop%h
        end if
      end do

      call summarize_entcell(ecp)

      end subroutine ent_struct_readcsv


      end module ent_make_struct
