      program convSOIL
c reads soil textures from  
c C.A. Reynolds, T. J. Jackson, and W.J. Rawls. 1999. 
c Estimated Available Water Content from the FAO Soil Map of 
c the World, Global Soil Profile Databases, and Pedo-transfer Function
c http://www.ngdc.noaa.gov/ecosys/cdroms/reynolds/reynolds/reynolds.htm
c and replaces missing data by data from 
c  * GSFC GLDAS soil database at 
c    http://ldas.gsfc.nasa.gov/GLDAS/SOILS/GLDASsoils.shtml
c  * for bedrock and peat, we use data from Cynthia Rosenzweig original
c    modelE file S360X180_0098M.rep
c the slope sl and layer thickness originate from S360X180_0098M.rep
c
c     
      integer, parameter :: imf=4320, jmf=2160,nx=1440,ny=600
      real*4 :: dz(360,180,6),
     &     ftext(360,180,5,6),      ! coarse grid, 5 textures, 6 layers
     &     ftextk(360,180,5,6),     ! coarse grid, 5 textures, 6 layers
     &     sl(360,180),sum(360,180,6),sum_rock_peat(360,180,6),
     &     ftext_fine(imf,jmf,3,2), ! fine grid, 3 textures, 2 layers
     &     sumfine(imf,jmf),landclss(imf,jmf),
     &     sand_top(imf,jmf),silt_top(imf,jmf),clay_top(imf,jmf),
     &     sand_bot(imf,jmf),silt_bot(imf,jmf),clay_bot(imf,jmf),
     &     sand_top_coarse(720,360),silt_top_coarse(720,360),
     &     clay_top_coarse(720,360),sumfixed(720,360),
     &     sand_bot_coarse(720,360),silt_bot_coarse(720,360),
     &     clay_bot_coarse(720,360),
     &     sand_lay(720,360,6),
     &     silt_lay(720,360,6),
     &     clay_lay(720,360,6),
     &     peat_lay(720,360,6),
     &     rock_lay(720,360,6),
     &     sl_targ(720,360),dz_targ(720,360,6),
     &     ftext_targ(720,360,5,6)
      real*4 :: ftext_gldas(nx,ny,3),fext_coarse(360,180)
      real*4 :: ftext_gldas_extend(nx,720,3),alpha

      character*80 title,title2,title3,title4,title5,
     &     name,nameout,nameout2,nameout3,nameout4,nameout5
      integer iu_SOIL,iuout,iuout2,iuout3,iuout4,iuout5,iuout6,i,k,
     &    iby2,jby2,iby3,jby3,iby4,jby4,icoarse,jcoarse
      character*1 :: kch,ich

      write(*,*) "ims,jms,nts",ims,jms,nts

      iu_SOIL=20
      name="S360X180_0098M.rep"
      
      open(iu_SOIL,FILE=name,FORM='unformatted', STATUS='old')
      
      read(iu_SOIL) dz,ftext,ftextk,sl

      close(iu_SOIL)
            
      iuout=iu_SOIL+1
      nameout="dz"
      open( iuout, FILE=nameout,
     &     FORM='unformatted', STATUS="UNKNOWN")
      

      do k=1,6
         write(kch,'(i1)') k
         title="thickness of layer " // kch
         write(*,*) title
        
         write(unit=iuout) title,dz(:,:,k)
      enddo

      close(iuout)
      
      iuout=iuout+1
      nameout="ftext"
      open( iuout, FILE=nameout,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do i=1,5
         do k=1,6
            write(kch,'(i1)') k
            write(ich,'(i1)') i
            title="texture " // ich // " layer " // kch
            write(*,*) title
            write(unit=iuout) title,ftext(:,:,i,k)
         enddo
      enddo
      close(iuout)

      iuout2=201
      nameout2="sumtext"
      open( iuout2, FILE=nameout2,
     &     FORM='unformatted', STATUS="UNKNOWN")
      
      do k=1,6
         sum(:,:,k)=0.
         write(kch,'(i1)') k
         do i=1,5
            sum(:,:,k)=sum(:,:,k)+ftext(:,:,i,k)
         enddo
         title2="sum textures for layer " // kch
         write(unit=iuout2) title2,sum(:,:,k)-1.0
      enddo
      close(iuout2)


      iuout3=202
      nameout3="sum-peat-bedrock"
      open( iuout3, FILE=nameout3,
     &     FORM='unformatted', STATUS="UNKNOWN")
      
      do k=1,6
         write(kch,'(i1)') k
         sum_rock_peat(:,:,k)= !ftext(:,:,4,k)+
     & ftext(:,:,5,k)
         title3="peat + bedrock for layer " // kch
         write(unit=iuout3) title3,sum_rock_peat(:,:,k)
      enddo
      close(iuout3)

      iuout=iuout+1
      nameout="ftextk"
      open( iuout, FILE=nameout,
     &     FORM='unformatted', STATUS="UNKNOWN")
      
      do i=1,5
         do k=1,6
            write(kch,'(i1)') k
            write(ich,'(i1)') i
            title="texturek " // ich // " layer " // kch
            write(*,*) title
            write(unit=iuout) title,ftextk(:,:,i,k)
         enddo
      enddo
         
      close(iuout)

      iuout=iuout+1
      nameout="sl"
      open( iuout, FILE=nameout,
     &     FORM='unformatted', STATUS="UNKNOWN")
      title="sl "
      write(*,*) title
      write(unit=iuout) title,sl(:,:)
      close(iuout)

      iuin=iuout+1
      open(unit=iuin,file='sand_tp1.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,ftext_fine(:,:,1,1)
      close(iuin)

      iuin=iuin+1
      open(unit=iuin,file='sand_sb1.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,ftext_fine(:,:,1,2)
      close(iuin)

      iuin=iuout+1
      open(unit=iuin,file='silt_tp1.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,ftext_fine(:,:,2,1)
      close(iuin)

      iuin=iuin+1
      open(unit=iuin,file='silt_sb1.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,ftext_fine(:,:,2,2)
      close(iuin)

      iuin=iuout+1
      open(unit=iuin,file='clay_tp1.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,ftext_fine(:,:,3,1)
      close(iuin)

      iuin=iuin+1
      open(unit=iuin,file='clay_sb1.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,ftext_fine(:,:,3,2)
      close(iuin)

      iuout4=203
      nameout4="sumtextfine"
      open( iuout4, FILE=nameout4,
     &     FORM='unformatted', STATUS="UNKNOWN")
            
      sumfine(:,:)=ftext_fine(:,:,1,1)
     &     + ftext_fine(:,:,2,1)
     &     + ftext_fine(:,:,3,1)

      title4="sum textures for top layer"
      write(unit=iuout4) title4,sumfine(:,:)
      close(iuout4)

c*     fix "no data" regions in Reynolds dataset
c*     depending on the land types in landclss.img
          
      iuin=iuin+1
      open(unit=iuin,file='landclss.giss',form='unformatted',
     &     status='old')
      read(unit=iuin) title,landclss
      close(iuin)

      iuin=iuin+1
      name="sandfao.1gd4r"
      open(iuin,file=name,form='unformatted',status='old',
     &    access='direct',recl=nx*ny*4)
      read(iuin,rec=1) ftext_gldas(:,:,1)
      close(iuin)

      iuin=iuin+1
      name="siltfao.1gd4r"
      open(iuin,file=name,form='unformatted',status='old', 
     &    access='direct',recl=nx*ny*4)
      read(iuin,rec=1) ftext_gldas(:,:,2)
      close(iuin)

      iuin=iuin+1
      name="clayfao.1gd4r"
      open(iuin,file=name,form='unformatted',status='old', 
     &    access='direct',recl=nx*ny*4)
      read(iuin,rec=1) ftext_gldas(:,:,3)
      close(iuin)

      ftext_gldas_extend=0.
      ftext_gldas_extend(:,121:720,:)=
     &     max(ftext_gldas(:,1:600,:),0.)

      do i=1,1440
         do j=1,720
            iby4=i/4
            jby4=j/4
            icoarse=min(iby4,1440)
            jcoarse=min(jby4,720)
            fext_coarse(icoarse,jcoarse)=ftext_gldas_extend(i,j,2)
         enddo
      enddo

      iuout=1000
      title="test"
      open(unit=iuout,file='test',form='unformatted',
     &     status='unknown')
      write(iuout) title,fext_coarse(:,:)
      close(iuout)

      do j=1,jmf
         do i=1,imf
               sand_top(i,j)=ftext_fine(i,j,1,1)
               silt_top(i,j)=ftext_fine(i,j,2,1)
               clay_top(i,j)=ftext_fine(i,j,3,1)

               sand_bot(i,j)=ftext_fine(i,j,1,2)
               silt_bot(i,j)=ftext_fine(i,j,2,2)
               clay_bot(i,j)=ftext_fine(i,j,3,2)

            if (landclss(i,j) .eq. 117.) then  !shifting dunes
             sand_top(i,j)=100.0
             silt_top(i,j)=0.0
             clay_top(i,j)=0.0

             sand_bot(i,j)=100.0
             silt_bot(i,j)=0.0
             clay_bot(i,j)=0.0
            endif
            if(landclss(i,j) .eq. 116.) then  !glaciers
             sand_top(i,j)=100./3.
             silt_top(i,j)=100./3.
             clay_top(i,j)=100./3.

             sand_bot(i,j)=100./3.
             silt_bot(i,j)=100./3.
             clay_bot(i,j)=100./3.
            endif
            if((landclss(i,j) .eq. 113) .or. (landclss(i,j) .eq. 114)
     &       .or. (landclss(i,j) .eq. 115)) then
             iby3=i/3
             jby3=j/3
             icoarse=min(iby3+1,1440)
             jcoarse=min(jby3+1,720)
             
             sand_top(i,j)=100.*ftext_gldas_extend(icoarse,jcoarse,1)
             silt_top(i,j)=100.*ftext_gldas_extend(icoarse,jcoarse,2)
             clay_top(i,j)=100.*ftext_gldas_extend(icoarse,jcoarse,3)

             sand_bot(i,j)=100.*ftext_gldas_extend(icoarse,jcoarse,1)
             silt_bot(i,j)=100.*ftext_gldas_extend(icoarse,jcoarse,2)
             clay_bot(i,j)=100.*ftext_gldas_extend(icoarse,jcoarse,3)
            endif
         enddo
      enddo

c*     transfer to 0.5 degree x 0.5 degree grid
      ratio = imf/720

      do i=1,720
         do j=1,360
            do k=1,ratio
               ires=(i-1)*ratio+k
               do m=1,ratio
                  jres=(j-1)*ratio+m
                  sand_top_coarse(i,j)=sand_top_coarse(i,j)
     &             +sand_top(ires,jres)
                  silt_top_coarse(i,j)=silt_top_coarse(i,j)
     &             +silt_top(ires,jres)
                  clay_top_coarse(i,j)=clay_top_coarse(i,j)
     &             +clay_top(ires,jres)

                  sand_bot_coarse(i,j)=sand_bot_coarse(i,j)
     &             +sand_bot(ires,jres)
                  silt_bot_coarse(i,j)=silt_bot_coarse(i,j)
     &             +silt_bot(ires,jres)
                  clay_bot_coarse(i,j)=clay_bot_coarse(i,j)
     &             +clay_bot(ires,jres)
               enddo
            enddo
           sand_top_coarse(i,j)=sand_top_coarse(i,j)/real(ratio*ratio)
           silt_top_coarse(i,j)=silt_top_coarse(i,j)/real(ratio*ratio)
           clay_top_coarse(i,j)=clay_top_coarse(i,j)/real(ratio*ratio)

           sand_bot_coarse(i,j)=sand_bot_coarse(i,j)/real(ratio*ratio)
           silt_bot_coarse(i,j)=silt_bot_coarse(i,j)/real(ratio*ratio)
           clay_bot_coarse(i,j)=clay_bot_coarse(i,j)/real(ratio*ratio)
         enddo
      enddo

c
c layer 1 = 100% top
c layer 2 = 100% top
c layer 3 = 10.17% top, 89.82% bottom
c layer 4 = 100% bottom
c layer 5 = 100% bottom
c layer 6 = 100% bottom 
c
      sand_lay(:,:,:)=0.
      silt_lay(:,:,:)=0.
      clay_lay(:,:,:)=0.
      peat_lay(:,:,:)=0.
      rock_lay(:,:,:)=0.

      sand_lay(:,:,1)= sand_top_coarse(:,:)
      sand_lay(:,:,2)= sand_top_coarse(:,:)
      sand_lay(:,:,3)= 0.8983*sand_top_coarse(:,:)
     &     + 0.1017*sand_bot_coarse(:,:)
      sand_lay(:,:,4)= sand_bot_coarse(:,:)
      sand_lay(:,:,5)= sand_bot_coarse(:,:)
      sand_lay(:,:,6)= sand_bot_coarse(:,:)

      silt_lay(:,:,1)= silt_top_coarse(:,:)
      silt_lay(:,:,2)= silt_top_coarse(:,:)
      silt_lay(:,:,3)= 0.8282*silt_top_coarse(:,:)
     &     + 0.1017*silt_bot_coarse(:,:)
      silt_lay(:,:,4)= silt_bot_coarse(:,:)
      silt_lay(:,:,5)= silt_bot_coarse(:,:)
      silt_lay(:,:,6)= silt_bot_coarse(:,:)

      clay_lay(:,:,1)= clay_top_coarse(:,:)
      clay_lay(:,:,2)= clay_top_coarse(:,:)
      clay_lay(:,:,3)= 0.8282*clay_top_coarse(:,:)
     &     + 0.1017*clay_bot_coarse(:,:)
      clay_lay(:,:,4)= clay_bot_coarse(:,:)
      clay_lay(:,:,5)= clay_bot_coarse(:,:)
      clay_lay(:,:,6)= clay_bot_coarse(:,:)

c
c if peat or bedrock exist in Cynthia's original file
c get their fractions and rescale the other fractions
c

      do j=1,360
         do i=1,720
            iby2=i/2
            jby2=j/2
            icoarse=min(iby2+1,360)
            jcoarse=min(jby2+1,180)
            do k=1,6
               if (sum_rock_peat(icoarse,jcoarse,k) .gt. 0.) then
                  peat_lay(i,j,k)=100.*ftext(icoarse,jcoarse,4,k)
                  rock_lay(i,j,k)=100.*ftext(icoarse,jcoarse,5,k)
                  if (sand_lay(i,j,k)+silt_lay(i,j,k)
     &                 +clay_lay(i,j,k) .gt. 0.) then
                     alpha=
     &               (100.-peat_lay(i,j,k)-rock_lay(i,j,k))/
     &               (sand_lay(i,j,k)+silt_lay(i,j,k)+clay_lay(i,j,k))
                     sand_lay(i,j,k)=alpha*sand_lay(i,j,k)
                     silt_lay(i,j,k)=alpha*silt_lay(i,j,k)
                     clay_lay(i,j,k)=alpha*clay_lay(i,j,k)
                  endif
               endif
               sand_lay(i,j,k)=sand_lay(i,j,k) / 100.
               silt_lay(i,j,k)=silt_lay(i,j,k) / 100.
               clay_lay(i,j,k)=clay_lay(i,j,k) / 100.
               peat_lay(i,j,k)=peat_lay(i,j,k) / 100.
               rock_lay(i,j,k)=rock_lay(i,j,k) / 100.
               
               dz_targ(i,j,k)=dz(icoarse,jcoarse,k)
            enddo               
            sl_targ(i,j)=sl(icoarse,jcoarse)
            
         enddo
      enddo

      iuout5=204
      nameout5="sumcoarse-fixed"
      open( iuout5, FILE=nameout5,
     &     FORM='unformatted', STATUS="UNKNOWN")
            
      sumfixed= sand_lay(:,:,1)+silt_lay(:,:,1)+clay_lay(:,:,1)
     &     +peat_lay(:,:,1)+rock_lay(:,:,1)

      title5="sum textures for top layer"
      write(unit=iuout5) title5,sumfixed
      close(iuout5)

      ftext_targ(:,:,1,:)=sand_lay(:,:,:)
      ftext_targ(:,:,2,:)=silt_lay(:,:,:)
      ftext_targ(:,:,3,:)=clay_lay(:,:,:)
      ftext_targ(:,:,4,:)=peat_lay(:,:,:)
      ftext_targ(:,:,5,:)=rock_lay(:,:,:)

      iuout6=205
      nameout5="SOIL720X360_Reynolds"
      open( iuout6, FILE=nameout5,
     &     FORM='unformatted', STATUS="UNKNOWN")
            
      write(unit=iuout6)  dz_targ,ftext_targ,ftext_targ,sl_targ 
      close(iuout6)


      end program convSOIL
