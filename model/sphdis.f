      function sphdis(x1,y1,x2,y2)
c --- dist.(m) between 2 points on sphere, lat/lon (x1,y1) and lat/lon (x2,y2)
      implicit none
      real x1,y1,x2,y2,sphdis,ang,radius,radian
      data radius/6375.e3/,radian/57.2957795/
c
      ang=mod(y2-y1+540.,360.)-180.
      sphdis=radius*acos(min(1.,cosd(90.-x1)*cosd(90.-x2)
     .                         +sind(90.-x1)*sind(90.-x2)*cosd(ang)))
      if (sphdis.eq.0.) 
     .  sphdis=radius*sqrt((x2-x1)**2+(ang*cosd(.5*(x1+x2)))**2)/radian
cdiag if (sphdis.eq.0.) write (*,'(a,2f8.3,2x,2f8.3)')
cdiag.  'warning - zero distance between lat/lon points',x1,y1,x2,y2
      sphdis=max(sphdis,1.)
      return
      end
