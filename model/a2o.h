c     integer iio,jjo,iia,jja,nwgta2o,nwgto2a
c     parameter (iio=181,jjo=180,iia=72,jja=46,nwgta2o=jjo
c    .                 ,nwgto2a=4*jja)
      integer nwgta2o,nwgto2a,nwgta2o2
      parameter (nwgta2o=18,nwgto2a=29,nwgta2o2=18)
c
      real*8 wlista2o(iio,jjo,nwgta2o),wtaua2o(iio,jjo,nwgta2o2)
     .    ,wlisto2a(iia,jja,nwgto2a),ofrac(iia,jja)
     .    ,cosa(iia,jja,2),sina(iia,jja,2),agcmgdsz(jja)
     .    ,ocellsz(iio,jjo),coso(iio,jjo),sino(iio,jjo)
      integer ilista2o(iio,jjo,nwgta2o),jlista2o(iio,jjo,nwgta2o)
     .                                 ,nlista2o(iio,jjo)
     .       ,itaua2o(iio,jjo,nwgta2o2), jtaua2o(iio,jjo,nwgta2o2)
     .                                 , ntaua2o(iio,jjo)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a(iia,jja,nwgto2a)
     .                                 ,nlisto2a(iia,jja)
      common /a2o/ilista2o,jlista2o,nlista2o,wlista2o
     .           ,itaua2o,jtaua2o,ntaua2o,wtaua2o
     .           ,coso,sino,ofrac
      common /o2a/ilisto2a,jlisto2a,nlisto2a,wlisto2a
     .           ,cosa,sina,ocellsz
