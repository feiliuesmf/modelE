c-----------------------------------------------------------------------------
      real c1,c2,c3,c4,c5,c6,c7,c8,c9,qttt,qtt,qt,qs,qst,qpt,qpst,qptt
     .    ,pref,sclkap,qthref,reftem,refsal,refprs
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
c
      data pref/2000.e4/					!  SI units
c --- coefficients for sigma-4 (based on Brydon & Sun fit)
ccc     .  1.92362E+01,-8.82080E-02, 7.73552E-01,-5.46858E-03,-2.31866E-03,
ccc     .  2.11306E-05, 2.82474E-05/
ccc      data pref/4000.e4/					!  SI units
c
c --- coefficients for kappa^(theta):
c
      data sclkap/1.e-11/,qthref/1.e3/		      !  SI units
      data qttt,qtt,qt,qs,qst,qpt,qpst,qptt,reftem,refsal/
c --- reference: t= 0.0, s=34.0, p=0 bar
ccc     . -3.07410328E-05, 4.60243067E-03,-2.90445096E-01,-1.09102601E-01,
ccc     .  7.95549079E-04, 1.08918830E-09, 1.42421130E-11,-1.31927315E-11,
ccc     .   0.0, 34.0/
c --- reference: t= 0.0, s=35.0, p=0 bar
ccc     . -3.07410328E-05, 4.60243067E-03,-2.89649547E-01,-1.09102601E-01,
ccc     .  7.95549079E-04, 1.10343041E-09, 1.42421130E-11,-1.31927315E-11,
ccc     .   0.0, 35.0/
c --- reference: t= 3.0, s=35.0, p=0 bar
     . -3.07410328E-05, 4.32576137E-03,-2.62884494E-01,-1.05540980E-01,
     .  7.76025039E-04, 1.02498399E-09, 1.49520781E-11,-1.31927314E-11,
     .   3.0, 35.0/
c
c> Revision history:
c>
c> Mar. 1999 - made thermobaric coefficients conform to Brydon & Sun
c> Nov. 2002 - updated thermobaric coefficients
c-----------------------------------------------------------------------------
