c-----------------------------------------------------------------------------
      real c1,c2,c3,c4,c5,c6,c7,c8,c9
css  .  ,qttt,qtt,qt,qs,qst,qpt,qpst,qptt
css  .  ,pref,sclkap,reftem,refsal,refprs
     .  ,pref,sclkap
c
ccc      data c1,c2,c3,c4,c5,c6,c7/
c --- coefficients for sigma-0 (based on Brydon & Sun fit)
ccc     . -1.36471E-01, 4.68181E-02, 8.07004E-01,-7.45353E-03,-2.94418E-03,
ccc    .  3.43570E-05, 3.48658E-05/
ccc      data pref/0./
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
css     .  9.77093E+00,-2.26493E-02, 7.89879E-01,-6.43205E-03,-2.62983E-03,
css     .  2.75835E-05, 3.15235E-05/
c     data c1,c2,c3,c4,c5,c6,c7,c8,c9/  ! T:[-2,32],S:[16,38] c1-c9 fit for sig2
c    .  9.887030E+00,-1.571861E-02, 7.827997E-01,-6.575396E-03,
c    . -2.913760E-03, 2.974136E-05, 3.248280E-05, 1.062302E-04,
c    .  3.524143E-06/
      data c1,c2,c3,c4,c5,c6,c7,c8,c9/  ! T:[-2,30],S:[18,38] c1-c9 fit for sig2
     .  9.893322E+00,-1.500683E-02, 7.823698E-01,-6.630868E-03,
     . -2.936529E-03, 3.068372E-05, 3.343679E-05, 1.135806E-04,
     .  3.620535E-06/
      data pref/2.e7/
c
c --- sub-coefficients for locally referenced sigma
c --- a fit towards Jackett & McDougall (1995)
      real, parameter, dimension(7) ::
     &  alphap = (/ -0.1364705627213484   , 0.04681812123458564,
     &               0.80700383913187     ,-0.007453530323180844,
     &              -0.002944183249153631 , 0.00003435702568990446,
     &               0.0000348657661057688 /)
     & ,betap  = (/  0.05064226654169138  ,-0.0003571087848996894,
     &              -0.0000876148051892879, 5.252431910751829e-6,
     &               1.579762259448864e-6 ,-3.466867400295792e-8,
     &              -1.687643078774232e-8 /)
     & ,gammap = (/ -5.526396144304812e-6 , 4.885838128243163e-8,
     &               9.96026931578033e-9  ,-7.251389796582352e-10,
     &              -3.987360250058777e-11, 4.006307891935698e-12,
     &               8.26367520608008e-13 /)
c
c --- coefficients for kappa^(theta) (revised Sun et al. 1999):
c
      data sclkap/1.e-11/                                      !  SI units
c     data qttt,qtt,qt,qs,qst,qpt,qpst,qptt,reftem,refsal/
c --- reference: t= 1.0, s=34.0, p=0 bar,kap(4.5,34.5,1.e7)=  0.09265729
c    . -3.03869352E-05, 4.47509519E-03,-2.80147157E-01,-1.07490547E-01,
c    .  7.82551293E-04, 1.04449450E-09, 1.44433359E-11,-1.31383707E-11,
c    .   1.0, 34.0/
c --- reference: t= 0.0, s=34.5, p=0 pascal, kap_t(4.5,34.5,1.e7)=  0.11411243
c    . -3.03869354E-05, 4.56625601E-03,-2.88801209E-01,-1.08670290E-01,
c    .  7.90503772E-04, 1.07813750E-09, 1.41541548E-11,-1.31383708E-11,
c    .   0.0, 34.5/
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
c --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
      real, parameter, dimension(3) ::
     &  toff = (/  0.0,             3.0,            13.0 /)
     & ,soff = (/ 34.5,            35.0,            38.5 /)
     & ,qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /)
     & ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /)
     & ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /)
     & ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /)
     & ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /)
     & ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /)
     & ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /)
     & ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)
c
c> Revision history:
c>
c> Mar. 1999 - made thermobaric coefficients conform to Brydon & Sun
c> Nov. 2002 - updated thermobaric coefficients
c> May  2003 - further expansion of thermobaric options
c> July 2005 - added coeeficients for in-situ density (alphap,betap,gammap)
c-----------------------------------------------------------------------------
