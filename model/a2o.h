      integer nwgta2o,nwgta2o2,nwgto2a
      parameter (nwgta2o=18,nwgta2o2=18,nwgto2a=39)
c
      real*8 wlista2o(iio,jjo,nwgta2o),wtaua2o(iio,jjo,nwgta2o2)
     .    ,wlisto2a(iia,jja,nwgto2a)
     .    ,wlisto2a_e(iia,jja,nwgto2a)
     .    ,wlisto2a_n(iia,jja,nwgto2a)
     .    ,ofrac(iia,jja),agcmgdsz(jja)
     .    ,ocellsz(iio,jjo),coso(iio,jjo),sino(iio,jjo)
      integer ilista2o(iio,jjo,nwgta2o),jlista2o(iio,jjo,nwgta2o)
     .                                 ,nlista2o(iio,jjo)
     .       ,itaua2o(iio,jjo,nwgta2o2), jtaua2o(iio,jjo,nwgta2o2)
     .                                 , ntaua2o(iio,jjo)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a(iia,jja,nwgto2a)
     .                                 ,nlisto2a(iia,jja)
     .     ,ilisto2a_e(iia,jja,nwgto2a),jlisto2a_e(iia,jja,nwgto2a)
     .                                 ,nlisto2a_e(iia,jja)
     .     ,ilisto2a_n(iia,jja,nwgto2a),jlisto2a_n(iia,jja,nwgto2a)
     .                                 ,nlisto2a_n(iia,jja)
      common /a2o/ilista2o,jlista2o,nlista2o,wlista2o
     .           ,itaua2o,jtaua2o,ntaua2o,wtaua2o
     .           ,coso,sino,ofrac
      common /o2a/ilisto2a,jlisto2a,nlisto2a,wlisto2a
     .     ,ilisto2a_e,jlisto2a_e,nlisto2a_e,wlisto2a_e
     .     ,ilisto2a_n,jlisto2a_n,nlisto2a_n,wlisto2a_n
     .     ,ocellsz
