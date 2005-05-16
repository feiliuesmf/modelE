      subroutine inikpp
c
c --- hycom version 1.0
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kpp.h'
c
c -------------------------------------------------------------------
c --- initialize large, mc williams, doney kpp vertical mixing scheme
c -------------------------------------------------------------------
c
css   integer i,j
      integer    nzehat,nustar
      parameter (nzehat=890,nustar=192)
c
      real, dimension (0:nzehat+1,0:nustar+1) ::
     & wmt            ! momentum velocity scale table
     &,wst            ! scalar   velocity scale table
      common/kppltr/ wmt,wst
      save  /kppltr/
c
      real zehat,zeta,am,cm,c22,zetam,as,c33,zetas,usta,
     .     athird,ahalf,afourth
c
      data am,cm,c22,zetam/1.257,8.380,16.0,-0.2/
      data as,c33,zetas/-28.86,16.0,-1.0/
      data athird/0.3333333333333/,ahalf/0.5/,afourth/0.25/
c
c --- construct the velocity-scale lookup tables
c
      deltaz = (zmax-zmin)/(nzehat+1)
      deltau = (umax-umin)/(nustar+1)
c
      do 200 i=0,nzehat+1
      zehat=deltaz*float(i)+zmin
      do 200 j=0,nustar+1
      usta=deltau*float(j)+umin
      zeta=zehat/(usta**3+epsil)
      if (zehat.ge.0.) then
        wmt(i,j)=vonk*usta/(1.+c11*zeta)
        wst(i,j)=wmt(i,j)
      else
        if (zeta.gt.zetam) then
          wmt(i,j)=vonk*usta*(1.-c22*zeta)**afourth
        else
          wmt(i,j)=vonk*(am*usta**3-cm*zehat)**athird
        end if
        if (zeta.gt.zetas) then
          wst(i,j)=vonk*usta*(1.-c33*zeta)**ahalf
        else
          wst(i,j)=vonk*(as*usta**3-cs*zehat)**athird
        endif
      endif
 200  continue
c
c --- set constants
      vtc=cv*sqrt(.2/cs/epsilon)/vonk**2/ricr
      cg=cstar*vonk*(cs*vonk*epsilon)**athird
      dp0enh=2.0*tenm
c
c$OMP PARALLEL DO
      do 199 j=1,jj
      do 199 l=1,isp(j)
      do 199 i=ifp(j,l),ilp(j,l)
c --- map shallow depths to high jerlov numbers
      jerlov(i,j)=6-max(1,min(5,int(depths(i,j)/15.0)))
      jerlov(i,j)=max(jerlv0,jerlov(i,j))
 199  continue
c$OMP END PARALLEL DO
c
      return
      end
