      integer ia2o,ja2o,ka2o,io2a,jo2a,ko2a,katm,kocn
c --- accumulate fields in the coupler      
      real*8 ataux(iia,jja),atauy(iia,jja),aflxa2o(iia,jja)
     .      ,aemnp(iia,jja),aice(iia,jja),asalt(iia,jja)
     .      ,austar(iia,jja),aswflx(iia,jja),asst(iia,jja)
c
      character title*80
      integer afogcm,nsavea,nsaveo
      common /cpl/ataux,atauy,aflxa2o,aemnp,aice,asalt,austar
     .           ,aswflx,asst,afogcm,nsavea,nsaveo
