C**** FFT72.S    Fast Fourier Transform for KM = 72     7/07/92
C****
      SUBROUTINE FFT (F,A,B)
C****
C**** FFT calculates a fast fourier transform of the input array F,
C**** producing the cosine and sine coefficients in the output arrays
C**** A and B.  F is dimensioned by 72 = KM; A and B are dimensioned
C**** 0:36.  Upon entering FFT, the total energy is:
C****   .5*sum(F(K)**2)
C**** with the sum being taken over all K from 1 to KM.  The
C**** Fourier coefficients are defined by:
C****   A(N)+i*B(N) = sum(F(K)*exp(-2*PI*i*N*K/KM))/KMH
C**** with the sum being taken over all K from 1 to KM.  KMH = KM
C**** when N = 0 or KM/2, and KMH = KM/2 otherwise.  In the
C**** program's notation, CPQN means:
C****   CPQN = sum(F(K)*cos(2*PI*N*K/KM))
C**** with the sum being taken over all K = Q mod(P).  SPQN has a
C**** similar definition, but cos is replaced by sin.  The notation
C**** A=10, B=11, etc. is used for P, Q and N.  The same total
C**** energy can be calculated from the spectral coefficients as:
C****   .5*sum(A(N)**2+B(N)**2)*KMH
C**** with the sum being taken over all wave numbers N from 0 to KM/2.
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KM=72, BYKM=1./72.D0, BYKMH=2./72.D0, BYKM2=1./144.D0,
     *  TWOPI=6.283185307179586477D0)
      REAL*8 F(KM),A(0:KM/2),B(0:KM/2),E(0:KM/2)
      REAL*8 SINN(17),SINNH(35),COSNH(35)
      LOGICAL*4 QENERG
      COMMON /FFTCOM/ SIN5,SIN10,SIN15,SIN20,SIN25,SIN30,SIN35,SIN40,
     *   SIN45,SIN50,SIN55,SIN60,SIN65,SIN70,SIN75,SIN80,SIN85,
     *   SINNH,COSNH, R2,R3,TSIN75,TSIN15
      EQUIVALENCE (SINN,SIN5)
      QENERG = .FALSE.
C**** Calculate expressions summed by increments of 24
   10 CO00 = F(72)+F(24)+F(48)
      CO10 = F( 1)+F(25)+F(49)
      CO20 = F( 2)+F(26)+F(50)
      CO30 = F( 3)+F(27)+F(51)
      CO40 = F( 4)+F(28)+F(52)
      CO50 = F( 5)+F(29)+F(53)
      CO60 = F( 6)+F(30)+F(54)
      CO70 = F( 7)+F(31)+F(55)
      CO80 = F( 8)+F(32)+F(56)
      CO90 = F( 9)+F(33)+F(57)
      COA0 = F(10)+F(34)+F(58)
      COB0 = F(11)+F(35)+F(59)
      COC0 = F(12)+F(36)+F(60)
      COD0 = F(13)+F(37)+F(61)
      COE0 = F(14)+F(38)+F(62)
      COF0 = F(15)+F(39)+F(63)
      COG0 = F(16)+F(40)+F(64)
      COH0 = F(17)+F(41)+F(65)
      COI0 = F(18)+F(42)+F(66)
      COJ0 = F(19)+F(43)+F(67)
      COK0 = F(20)+F(44)+F(68)
      COL0 = F(21)+F(45)+F(69)
      COM0 = F(22)+F(46)+F(70)
      CON0 = F(23)+F(47)+F(71)
      CO01 = F(72)-(F(24)+F(48))*SIN30
      CO11 = F(1)*SIN85-F(25)*SIN35-F(49)*SIN25
      CO21 = F(2)*SIN80-F(26)*SIN40-F(50)*SIN20
      CO31 = F(3)*SIN75-F(27)*SIN45-F(51)*SIN15
      CO41 = F(4)*SIN70-F(28)*SIN50-F(52)*SIN10
      CO51 = F(5)*SIN65-F(29)*SIN55-F(53)*SIN5
      CO61 = (F(6)-F(30))*SIN60
      CO71 = F(7)*SIN55-F(31)*SIN65+F(55)*SIN5
      CO81 = F(8)*SIN50-F(32)*SIN70+F(56)*SIN10
      CO91 = F(9)*SIN45-F(33)*SIN75+F(57)*SIN15
      COA1 = F(10)*SIN40-F(34)*SIN80+F(58)*SIN20
      COB1 = F(11)*SIN35-F(35)*SIN85+F(59)*SIN25
      COC1 = (F(12)+F(60))*SIN30-F(36)
      COD1 = F(13)*SIN25-F(37)*SIN85+F(61)*SIN35
      COE1 = F(14)*SIN20-F(38)*SIN80+F(62)*SIN40
      COF1 = F(15)*SIN15-F(39)*SIN75+F(63)*SIN45
      COG1 = F(16)*SIN10-F(40)*SIN70+F(64)*SIN50
      COH1 = F(17)*SIN5-F(41)*SIN65+F(65)*SIN55
      COI1 = (F(66)-F(42))*SIN60
      COJ1 = F(67)*SIN65-F(19)*SIN5-F(43)*SIN55
      COK1 = F(68)*SIN70-F(20)*SIN10-F(44)*SIN50
      COL1 = F(69)*SIN75-F(21)*SIN15-F(45)*SIN45
      COM1 = F(70)*SIN80-F(22)*SIN20-F(46)*SIN40
      CON1 = F(71)*SIN85-F(23)*SIN25-F(47)*SIN35
      SO01 = (F(24)-F(48))*SIN60
      SO11 = F(1)*SIN5+F(25)*SIN55-F(49)*SIN65
      SO21 = F(2)*SIN10+F(26)*SIN50-F(50)*SIN70
      SO31 = F(3)*SIN15+F(27)*SIN45-F(51)*SIN75
      SO41 = F(4)*SIN20+F(28)*SIN40-F(52)*SIN80
      SO51 = F(5)*SIN25+F(29)*SIN35-F(53)*SIN85
      SO61 = (F(6)+F(30))*SIN30-F(54)
      SO71 = F(7)*SIN35+F(31)*SIN25-F(55)*SIN85
      SO81 = F(8)*SIN40+F(32)*SIN20-F(56)*SIN80
      SO91 = F(9)*SIN45+F(33)*SIN15-F(57)*SIN75
      SOA1 = F(10)*SIN50+F(34)*SIN10-F(58)*SIN70
      SOB1 = F(11)*SIN55+F(35)*SIN5-F(59)*SIN65
      SOC1 = (F(12)-F(60))*SIN60
      SOD1 = F(13)*SIN65-F(37)*SIN5-F(61)*SIN55
      SOE1 = F(14)*SIN70-F(38)*SIN10-F(62)*SIN50
      SOF1 = F(15)*SIN75-F(39)*SIN15-F(63)*SIN45
      SOG1 = F(16)*SIN80-F(40)*SIN20-F(64)*SIN40
      SOH1 = F(17)*SIN85-F(41)*SIN25-F(65)*SIN35
      SOI1 = F(18)-(F(42)+F(66))*SIN30
      SOJ1 = F(19)*SIN85-F(43)*SIN35-F(67)*SIN25
      SOK1 = F(20)*SIN80-F(44)*SIN40-F(68)*SIN20
      SOL1 = F(21)*SIN75-F(45)*SIN45-F(69)*SIN15
      SOM1 = F(22)*SIN70-F(46)*SIN50-F(70)*SIN10
      SON1 = F(23)*SIN65-F(47)*SIN55-F(71)*SIN5
C**** Calculate expressions summed by increments of 8
      C800 = CO00+CO80+COG0
      C810 = CO10+CO90+COH0
      C820 = CO20+COA0+COI0
      C830 = CO30+COB0+COJ0
      C840 = CO40+COC0+COK0
      C850 = CO50+COD0+COL0
      C860 = CO60+COE0+COM0
      C870 = CO70+COF0+CON0
      C801 = CO01+CO81+COG1
      C811 = CO11+CO91+COH1
      C821 = CO21+COA1+COI1
      C831 = CO31+COB1+COJ1
      C841 = CO41+COC1+COK1
      C851 = CO51+COD1+COL1
      C861 = CO61+COE1+COM1
      C871 = CO71+COF1+CON1
      C802 = CO01-(CO81+COG1)*SIN30+(SO81-SOG1)*SIN60
      C812 = (CO11-SOH1)*SIN75+(SO91-CO91)*SIN45+(SO11-COH1)*SIN15
      C822 = (CO21-COA1)*SIN60+(SO21+SOA1)*SIN30-SOI1
      C832 = (CO31+SO31)*SIN45-(COB1+SOJ1)*SIN75+(COJ1+SOB1)*SIN15
      C842 = (CO41+COK1)*SIN30+(SO41-SOK1)*SIN60-COC1
      C852 = (CO51-SOD1)*SIN15+(COL1-SOL1)*SIN45-(COD1-SO51)*SIN75
      C862 = (COM1-COE1)*SIN60-(SOM1+SOE1)*SIN30+SO61
      C872 = (CON1+SO71)*SIN75-(SOF1+COF1)*SIN45-(SON1+CO71)*SIN15
      C803 = CO00-(CO80+COG0)*SIN30
      C813 = CO10*SIN75-CO90*SIN45-COH0*SIN15
      C823 = (CO20-COA0)*SIN60
      C833 = CO30*SIN45-COB0*SIN75+COJ0*SIN15
      C843 = (CO40+COK0)*SIN30-COC0
      C853 = CO50*SIN15-COD0*SIN75+COL0*SIN45
      C863 = (COM0-COE0)*SIN60
      C873 = CON0*SIN75-COF0*SIN45-CO70*SIN15
      C804 = CO01-(CO81+COG1)*SIN30-(SO81-SOG1)*SIN60
      C814 = (CO11+SOH1)*SIN75-(SO91+CO91)*SIN45-(SO11+COH1)*SIN15
      C824 = (CO21-COA1)*SIN60-(SO21+SOA1)*SIN30+SOI1
      C834 = (CO31-SO31)*SIN45-(COB1-SOJ1)*SIN75+(COJ1-SOB1)*SIN15
      C844 = (CO41+COK1)*SIN30-(SO41-SOK1)*SIN60-COC1
      C854 = (CO51+SOD1)*SIN15+(COL1+SOL1)*SIN45-(COD1+SO51)*SIN75
      C864 = (COM1-COE1)*SIN60+(SOM1+SOE1)*SIN30-SO61
      C874 = (CON1-SO71)*SIN75+(SOF1-COF1)*SIN45+(SON1-CO71)*SIN15
      S801 = SO01+SO81+SOG1
      S811 = SO11+SO91+SOH1
      S821 = SO21+SOA1+SOI1
      S831 = SO31+SOB1+SOJ1
      S841 = SO41+SOC1+SOK1
      S851 = SO51+SOD1+SOL1
      S861 = SO61+SOE1+SOM1
      S871 = SO71+SOF1+SON1
      S802 = (CO81-COG1)*SIN60+(SO81+SOG1)*SIN30-SO01
      S812 = (CO11+SOH1)*SIN15+(CO91+SO91)*SIN45-(COH1+SO11)*SIN75
      S822 = (COA1+CO21)*SIN30+(SOA1-SO21)*SIN60-COI1
      S832 = (CO31-SO31)*SIN45+(SOB1-COJ1)*SIN75+(COB1-SOJ1)*SIN15
      S842 = (CO41-COK1)*SIN60-(SO41+SOK1)*SIN30+SOC1
      S852 = (CO51+SOD1)*SIN75-(SOL1+COL1)*SIN45-(COD1+SO51)*SIN15
      S862 = CO61-(COE1+COM1)*SIN30+(SOE1-SOM1)*SIN60
      S872 = (CO71-SON1)*SIN75+(SOF1-COF1)*SIN45-(CON1-SO71)*SIN15
      S803 = (CO80-COG0)*SIN60
      S813 = CO10*SIN15+CO90*SIN45-COH0*SIN75
      S823 = (COA0+CO20)*SIN30-COI0
      S833 = CO30*SIN45-COJ0*SIN75+COB0*SIN15
      S843 = (CO40-COK0)*SIN60
      S853 = CO50*SIN75-COL0*SIN45-COD0*SIN15
      S863 = CO60-(COE0+COM0)*SIN30
      S873 = CO70*SIN75-COF0*SIN45-CON0*SIN15
      S804 = (CO81-COG1)*SIN60-(SO81+SOG1)*SIN30+SO01
      S814 = (CO11-SOH1)*SIN15+(CO91-SO91)*SIN45-(COH1-SO11)*SIN75
      S824 = (COA1+CO21)*SIN30-(SOA1-SO21)*SIN60-COI1
      S834 = (CO31+SO31)*SIN45-(SOB1+COJ1)*SIN75+(COB1+SOJ1)*SIN15
      S844 = (CO41-COK1)*SIN60+(SO41+SOK1)*SIN30-SOC1
      S854 = (CO51-SOD1)*SIN75+(SOL1-COL1)*SIN45-(COD1-SO51)*SIN15
      S864 = CO61-(COE1+COM1)*SIN30-(SOE1-SOM1)*SIN60
      S874 = (CO71+SON1)*SIN75-(SOF1+COF1)*SIN45-(CON1+SO71)*SIN15
C**** Calculate expressions summed by increments of 4
      C400 = C800+C840
      C401 = C801+C841
      C402 = C802+C842
      C403 = C803+C843
      C404 = C804+C844
      C405 = C804-C844
      C406 = C803-C843
      C407 = C802-C842
      C408 = C801-C841
      C409 = C800-C840
      C410 = C810+C850
      C411 = C811+C851
      C412 = C812+C852
      C413 = C813+C853
      C414 = C814+C854
      C415 = (C814-C854+S814-S854)*SIN45
      C416 = (C813-C853+S813-S853)*SIN45
      C417 = (C812-C852+S812-S852)*SIN45
      C418 = (C811-C851+S811-S851)*SIN45
      C419 = (C810-C850)*SIN45
      C420 = C820+C860
      C421 = C821+C861
      C422 = C822+C862
      C423 = C823+C863
      C424 = C824+C864
      C425 = S824-S864
      C426 = S823-S863
      C427 = S822-S862
      C428 = S821-S861
C     C429 = 0.
      C430 = C830+C870
      C431 = C831+C871
      C432 = C832+C872
      C433 = C833+C873
      C434 = C834+C874
      C435 = (C874-C834+S834-S874)*SIN45
      C436 = (C873-C833+S833-S873)*SIN45
      C437 = (C872-C832+S832-S872)*SIN45
      C438 = (C871-C831+S831-S871)*SIN45
      C439 = (C870-C830)*SIN45
      S401 = S841+S801
      S402 = S842+S802
      S403 = S843+S803
      S404 = S844+S804
      S405 = S844-S804
      S406 = S843-S803
      S407 = S842-S802
      S408 = S841-S801
C     S409 = 0.
      S411 = S851+S811
      S412 = S852+S812
      S413 = S853+S813
      S414 = S854+S814
      S415 = (S854-S814+C814-C854)*SIN45
      S416 = (S853-S813+C813-C853)*SIN45
      S417 = (S852-S812+C812-C852)*SIN45
      S418 = (S851-S811+C811-C851)*SIN45
      S419 = (C810-C850)*SIN45
      S421 = S821+S861
      S422 = S822+S862
      S423 = S823+S863
      S424 = S824+S864
      S425 = C824-C864
      S426 = C823-C863
      S427 = C822-C862
      S428 = C821-C861
      S429 = C820-C860
      S431 = S831+S871
      S432 = S832+S872
      S433 = S833+S873
      S434 = S834+S874
      S435 = (S834-S874+C834-C874)*SIN45
      S436 = (S833-S873+C833-C873)*SIN45
      S437 = (S832-S872+C832-C872)*SIN45
      S438 = (S831-S871+C831-C871)*SIN45
      S439 = (C830-C870)*SIN45
C**** Calculate expressions summed by increments of 2
      C200 = C400+C420
      C201 = C401+C421
      C202 = C402+C422
      C203 = C403+C423
      C204 = C404+C424
      C205 = C405+C425
      C206 = C406+C426
      C207 = C407+C427
      C208 = C408+C428
C     C209 = C409
      C20A = C408-C428
      C20B = C407-C427
      C20C = C406-C426
      C20D = C405-C425
      C20E = C404-C424
      C20F = C403-C423
      C20G = C402-C422
      C20H = C401-C421
C     C20I = C400-C420
      C210 = C410+C430
      C211 = C411+C431
      C212 = C412+C432
      C213 = C413+C433
      C214 = C414+C434
      C215 = C415+C435
      C216 = C416+C436
      C217 = C417+C437
      C218 = C418+C438
      C219 = S419-S439
      C21A = S418-S438
      C21B = S417-S437
      C21C = S416-S436
      C21D = S415-S435
      C21E = S414-S434
      C21F = S413-S433
      C21G = S412-S432
      C21H = S411-S431
C     C21I = 0.
C     S200 = 0.
      S201 = S421+S401
      S202 = S422+S402
      S203 = S423+S403
      S204 = S424+S404
      S205 = S425+S405
      S206 = S426+S406
      S207 = S427+S407
      S208 = S428+S408
C     S209 = S429
      S20A = S428-S408
      S20B = S427-S407
      S20C = S426-S406
      S20D = S425-S405
      S20E = S424-S404
      S20F = S423-S403
      S20G = S422-S402
      S20H = S421-S401
C     S20I = 0.
C     S210 = 0.
      S211 = S411+S431
      S212 = S412+S432
      S213 = S413+S433
      S214 = S414+S434
      S215 = S415+S435
      S216 = S416+S436
      S217 = S417+S437
      S218 = S418+S438
      S219 = C419-C439
      S21A = C418-C438
      S21B = C417-C437
      S21C = C416-C436
      S21D = C415-C435
      S21E = C414-C434
      S21F = C413-C433
      S21G = C412-C432
      S21H = C411-C431
C     S21I = C410-C430
      IF(QENERG)  GO TO 20
C**** Calculate final coefficients of fourier expansion
      A(0)  = (C200+C210)*BYKM
      A(1)  = (C201+C211)*BYKMH
      A(2)  = (C202+C212)*BYKMH
      A(3)  = (C203+C213)*BYKMH
      A(4)  = (C204+C214)*BYKMH
      A(5)  = (C205+C215)*BYKMH
      A(6)  = (C206+C216)*BYKMH
      A(7)  = (C207+C217)*BYKMH
      A(8)  = (C208+C218)*BYKMH
      A(9)  = (C409+C219)*BYKMH
      A(10) = (C20A+C21A)*BYKMH
      A(11) = (C20B+C21B)*BYKMH
      A(12) = (C20C+C21C)*BYKMH
      A(13) = (C20D+C21D)*BYKMH
      A(14) = (C20E+C21E)*BYKMH
      A(15) = (C20F+C21F)*BYKMH
      A(16) = (C20G+C21G)*BYKMH
      A(17) = (C20H+C21H)*BYKMH
      A(18) = (C400-C420)*BYKMH
      A(19) = (C20H-C21H)*BYKMH
      A(20) = (C20G-C21G)*BYKMH
      A(21) = (C20F-C21F)*BYKMH
      A(22) = (C20E-C21E)*BYKMH
      A(23) = (C20D-C21D)*BYKMH
      A(24) = (C20C-C21C)*BYKMH
      A(25) = (C20B-C21B)*BYKMH
      A(26) = (C20A-C21A)*BYKMH
      A(27) = (C409-C219)*BYKMH
      A(28) = (C208-C218)*BYKMH
      A(29) = (C207-C217)*BYKMH
      A(30) = (C206-C216)*BYKMH
      A(31) = (C205-C215)*BYKMH
      A(32) = (C204-C214)*BYKMH
      A(33) = (C203-C213)*BYKMH
      A(34) = (C202-C212)*BYKMH
      A(35) = (C201-C211)*BYKMH
      A(36) = (C200-C210)*BYKM
      B(0)  = 0.
      B(1)  = (S211+S201)*BYKMH
      B(2)  = (S212+S202)*BYKMH
      B(3)  = (S213+S203)*BYKMH
      B(4)  = (S214+S204)*BYKMH
      B(5)  = (S215+S205)*BYKMH
      B(6)  = (S216+S206)*BYKMH
      B(7)  = (S217+S207)*BYKMH
      B(8)  = (S218+S208)*BYKMH
      B(9)  = (S429+S219)*BYKMH
      B(10) = (S21A+S20A)*BYKMH
      B(11) = (S21B+S20B)*BYKMH
      B(12) = (S21C+S20C)*BYKMH
      B(13) = (S21D+S20D)*BYKMH
      B(14) = (S21E+S20E)*BYKMH
      B(15) = (S21F+S20F)*BYKMH
      B(16) = (S21G+S20G)*BYKMH
      B(17) = (S21H+S20H)*BYKMH
      B(18) = (C410-C430)*BYKMH
      B(19) = (S21H-S20H)*BYKMH
      B(20) = (S21G-S20G)*BYKMH
      B(21) = (S21F-S20F)*BYKMH
      B(22) = (S21E-S20E)*BYKMH
      B(23) = (S21D-S20D)*BYKMH
      B(24) = (S21C-S20C)*BYKMH
      B(25) = (S21B-S20B)*BYKMH
      B(26) = (S21A-S20A)*BYKMH
      B(27) = (S219-S429)*BYKMH
      B(28) = (S218-S208)*BYKMH
      B(29) = (S217-S207)*BYKMH
      B(30) = (S216-S206)*BYKMH
      B(31) = (S215-S205)*BYKMH
      B(32) = (S214-S204)*BYKMH
      B(33) = (S213-S203)*BYKMH
      B(34) = (S212-S202)*BYKMH
      B(35) = (S211-S201)*BYKMH
      B(36) = 0.
      RETURN
C****
C****
      ENTRY FFTI (A,B,F)
C****
C**** FFTI performs an inverse fast fourier transform
C**** Input: A,B = spectral coefficients
C**** Output:  F = grid point values
C****
C**** These have been divided by 18
      C200 = (A(0)+A(36))*2.
      C210 = (A(0)-A(36))*2.
      C201 = A( 1)+A(35)
      C211 = A( 1)-A(35)
      C202 = A( 2)+A(34)
      C212 = A( 2)-A(34)
      C203 = A( 3)+A(33)
      C213 = A( 3)-A(33)
      C204 = A( 4)+A(32)
      C214 = A( 4)-A(32)
      C205 = A( 5)+A(31)
      C215 = A( 5)-A(31)
      C206 = A( 6)+A(30)
      C216 = A( 6)-A(30)
      C207 = A( 7)+A(29)
      C217 = A( 7)-A(29)
      C208 = A( 8)+A(28)
      C218 = A( 8)-A(28)
      C209 = A( 9)+A(27)
      C219 = A( 9)-A(27)
      C20A = A(10)+A(26)
      C21A = A(10)-A(26)
      C20B = A(11)+A(25)
      C21B = A(11)-A(25)
      C20C = A(12)+A(24)
      C21C = A(12)-A(24)
      C20D = A(13)+A(23)
      C21D = A(13)-A(23)
      C20E = A(14)+A(22)
      C21E = A(14)-A(22)
      C20F = A(15)+A(21)
      C21F = A(15)-A(21)
      C20G = A(16)+A(20)
      C21G = A(16)-A(20)
      C20H = A(17)+A(19)
      C21H = A(17)-A(19)
      C20I = A(18)*2.
C     C21I = 0.
C     S200 = 0.
      S201 = B( 1)-B(35)
      S211 = B( 1)+B(35)
      S202 = B( 2)-B(34)
      S212 = B( 2)+B(34)
      S203 = B( 3)-B(33)
      S213 = B( 3)+B(33)
      S204 = B( 4)-B(32)
      S214 = B( 4)+B(32)
      S205 = B( 5)-B(31)
      S215 = B( 5)+B(31)
      S206 = B( 6)-B(30)
      S216 = B( 6)+B(30)
      S207 = B( 7)-B(29)
      S217 = B( 7)+B(29)
      S208 = B( 8)-B(28)
      S218 = B( 8)+B(28)
      S209 = B( 9)-B(27)
      S219 = B( 9)+B(27)
      S20A = B(10)-B(26)
      S21A = B(10)+B(26)
      S20B = B(11)-B(25)
      S21B = B(11)+B(25)
      S20C = B(12)-B(24)
      S21C = B(12)+B(24)
      S20D = B(13)-B(23)
      S21D = B(13)+B(23)
      S20E = B(14)-B(22)
      S21E = B(14)+B(22)
      S20F = B(15)-B(21)
      S21F = B(15)+B(21)
      S20G = B(16)-B(20)
      S21G = B(16)+B(20)
      S20H = B(17)-B(19)
      S21H = B(17)+B(19)
C     S20I = 0.
      S21I = B(18)*2.
C**** Now multiply by 2, so FACTOR = 2/18 = 9
      C400 = C200+C20I
      C420 = C200-C20I
      C401 = C201+C20H
      C421 = C201-C20H
      C402 = C202+C20G
      C422 = C202-C20G
      C403 = C203+C20F
      C423 = C203-C20F
      C404 = C204+C20E
      C424 = C204-C20E
      C405 = C205+C20D
      C425 = C205-C20D
      C406 = C206+C20C
      C426 = C206-C20C
      C407 = C207+C20B
      C427 = C207-C20B
      C408 = C208+C20A
      C428 = C208-C20A
      C409 = C209*2.
C     C429 = 0.
      C410 = C210+S21I
      C430 = C210-S21I
      C411 = C211+S21H
      C431 = C211-S21H
      C412 = C212+S21G
      C432 = C212-S21G
      C413 = C213+S21F
      C433 = C213-S21F
      C414 = C214+S21E
      C434 = C214-S21E
      C415 = C215+S21D
      C435 = C215-S21D
      C416 = C216+S21C
      C436 = C216-S21C
      C417 = C217+S21B
      C437 = C217-S21B
      C418 = C218+S21A
      C438 = C218-S21A
      S401 = S201-S20H
      S421 = S201+S20H
      S402 = S202-S20G
      S422 = S202+S20G
      S403 = S203-S20F
      S423 = S203+S20F
      S404 = S204-S20E
      S424 = S204+S20E
      S405 = S205-S20D
      S425 = S205+S20D
      S406 = S206-S20C
      S426 = S206+S20C
      S407 = S207-S20B
      S427 = S207+S20B
      S408 = S208-S20A
      S428 = S208+S20A
C     S409 = 0.
      S429 = S209*2.
      S411 = S211+C21H
      S431 = S211-C21H
      S412 = S212+C21G
      S432 = S212-C21G
      S413 = S213+C21F
      S433 = S213-C21F
      S414 = S214+C21E
      S434 = S214-C21E
      S415 = S215+C21D
      S435 = S215-C21D
      S416 = S216+C21C
      S436 = S216-C21C
      S417 = S217+C21B
      S437 = S217-C21B
      S418 = S218+C21A
      S438 = S218-C21A
      S419 = S219+C219
      S439 = S219-C219
C**** Multiply by 2 again, so FACTOR = 2*2/18 = 2/9
      C800 = C400+C409
      C840 = C400-C409
      C801 = C401+C408
      C841 = C401-C408
      C802 = C402+C407
      C842 = C402-C407
      C803 = C403+C406
      C843 = C403-C406
      C804 = C404+C405
      C844 = C404-C405
      C820 = C420+S429
      C860 = C420-S429
      C821 = C421+S428
      C861 = C421-S428
      C822 = C422+S427
      C862 = C422-S427
      C823 = C423+S426
      C863 = C423-S426
      C824 = C424+S425
      C864 = C424-S425
      S821 = S421+C428
      S861 = S421-C428
      S822 = S422+C427
      S862 = S422-C427
      S823 = S423+C426
      S863 = S423-C426
      S824 = S424+C425
      S864 = S424-C425
      S801 = S401-S408
      S841 = S401+S408
      S802 = S402-S407
      S842 = S402+S407
      S803 = S403-S406
      S843 = S403+S406
      S804 = S404-S405
      S844 = S404+S405
      C850 = C410-S419*R2
      C810 = C410+S419*R2
      S851 = S411+(S418-C418)*SIN45
      S811 = S411-(S418-C418)*SIN45
      C851 = C411-(S418+C418)*SIN45
      C811 = C411+(S418+C418)*SIN45
      S852 = S412+(S417-C417)*SIN45
      S812 = S412-(S417-C417)*SIN45
      C852 = C412-(S417+C417)*SIN45
      C812 = C412+(S417+C417)*SIN45
      S853 = S413+(S416-C416)*SIN45
      S813 = S413-(S416-C416)*SIN45
      C853 = C413-(S416+C416)*SIN45
      C813 = C413+(S416+C416)*SIN45
      S854 = S414+(S415-C415)*SIN45
      S814 = S414-(S415-C415)*SIN45
      C854 = C414-(S415+C415)*SIN45
      C814 = C414+(S415+C415)*SIN45
      C830 = C430+S439*R2
      C870 = C430-S439*R2
      S831 = S431+(S438+C438)*SIN45
      S871 = S431-(S438+C438)*SIN45
      C831 = C431+(S438-C438)*SIN45
      C871 = C431-(S438-C438)*SIN45
      S832 = S432+(S437+C437)*SIN45
      S872 = S432-(S437+C437)*SIN45
      C832 = C432+(S437-C437)*SIN45
      C872 = C432-(S437-C437)*SIN45
      S833 = S433+(S436+C436)*SIN45
      S873 = S433-(S436+C436)*SIN45
      C833 = C433+(S436-C436)*SIN45
      C873 = C433-(S436-C436)*SIN45
      S834 = S434+(S435+C435)*SIN45
      S874 = S434-(S435+C435)*SIN45
      C834 = C434+(S435-C435)*SIN45
      C874 = C434-(S435-C435)*SIN45
C**** Multiply by 3, so FACTOR = 3*2*2/18 = 3*2/9 = 2/3
      CO00 = C800+C803*2.
      COG0 = C800-C803-S803*R3
      CO80 = C800-C803+S803*R3
      COI0 = C820-S823*2.
      CO20 = C820+S823+C823*R3
      COA0 = C820+S823-C823*R3
      COC0 = C840-C843*2.
      CO40 = C840+C843+S843*R3
      COK0 = C840+C843-S843*R3
      CO60 = C860+S863*2.
      COE0 = C860-S863-C863*R3
      COM0 = C860-S863+C863*R3
      CO90 = C810-(C813-S813)*R2
      CO10 = C810+C813*TSIN75+S813*TSIN15
      COH0 = C810-C813*TSIN15-S813*TSIN75
      CO30 = C830+(C833+S833)*R2
      COB0 = C830-C833*TSIN75+S833*TSIN15
      COJ0 = C830+C833*TSIN15-S833*TSIN75
      COL0 = C850+(C853-S853)*R2
      COD0 = C850-C853*TSIN75-S853*TSIN15
      CO50 = C850+C853*TSIN15+S853*TSIN75
      COF0 = C870-(C873+S873)*R2
      CON0 = C870+C873*TSIN75-S873*TSIN15
      CO70 = C870-C873*TSIN15+S873*TSIN75
      CO01 =  C802+C804+C801
      SO01 = -S802+S804+S801
      CO81 =-(C802+C804)*SIN30+(S802+S804)*SIN60+C801
      COG1 =-(C802+C804)*SIN30-(S802+S804)*SIN60+C801
      SO81 = (C802-C804)*SIN60+(S802-S804)*SIN30+S801
      SOG1 =-(C802-C804)*SIN60+(S802-S804)*SIN30+S801
      COI1 = -S822-S824+C821
      SOI1 = -C822+C824+S821
      CO21 = (C822+C824)*SIN60+(S822+S824)*SIN30+C821
      COA1 =-(C822+C824)*SIN60+(S822+S824)*SIN30+C821
      SO21 = (C822-C824)*SIN30-(S822-S824)*SIN60+S821
      SOA1 = (C822-C824)*SIN30+(S822-S824)*SIN60+S821
      SOC1 =  S842-S844+S841
      COC1 = -C842-C844+C841
      SO41 = (C842-C844)*SIN60-(S842-S844)*SIN30+S841
      SOK1 =-(C842-C844)*SIN60-(S842-S844)*SIN30+S841
      CO41 = (C842+C844)*SIN30+(S842+S844)*SIN60+C841
      COK1 = (C842+C844)*SIN30-(S842+S844)*SIN60+C841
      CO61 =  S862+S864+C861
      SO61 =  C862-C864+S861
      COM1 = (C862+C864)*SIN60-(S862+S864)*SIN30+C861
      COE1 =-(C862+C864)*SIN60-(S862+S864)*SIN30+C861
      SOM1 =-(C862-C864)*SIN30-(S862-S864)*SIN60+S861
      SOE1 =-(C862-C864)*SIN30+(S862-S864)*SIN60+S861
      CO91 =(-C812-C814+S812+S814)*SIN45+C811
      SO91 = (C812-C814+S812-S814)*SIN45+S811
      CO11 = (C812+C814)*SIN75+(S812+S814)*SIN15+C811
      COH1 =-(C812+C814)*SIN15-(S812+S814)*SIN75+C811
      SO11 = (C812-C814)*SIN15-(S812-S814)*SIN75+S811
      SOH1 =-(C812-C814)*SIN75+(S812-S814)*SIN15+S811
      CO31 = (C832+C834+S832+S834)*SIN45+C831
      SO31 = (C832-C834-S832+S834)*SIN45+S831
      COB1 =-(C832+C834)*SIN75+(S832+S834)*SIN15+C831
      SOJ1 =-(C832-C834)*SIN75-(S832-S834)*SIN15+S831
      COJ1 = (C832+C834)*SIN15-(S832+S834)*SIN75+C831
      SOB1 = (C832-C834)*SIN15+(S832-S834)*SIN75+S831
      COL1 = (C852+C854-S852-S854)*SIN45+C851
      SOL1 =(-C852+C854-S852+S854)*SIN45+S851
      COD1 =-(C852+C854)*SIN75-(S852+S854)*SIN15+C851
      SOD1 =-(C852-C854)*SIN15+(S852-S854)*SIN75+S851
      CO51 = (C852+C854)*SIN15+(S852+S854)*SIN75+C851
      SO51 = (C852-C854)*SIN75-(S852-S854)*SIN15+S851
      COF1 =-(C872+C874+S872+S874)*SIN45+C871
      SOF1 =(-C872+C874+S872-S874)*SIN45+S871
      CON1 = (C872+C874)*SIN75-(S872+S874)*SIN15+C871
      CO71 =-(C872+C874)*SIN15+(S872+S874)*SIN75+C871
      SON1 =-(C872-C874)*SIN15-(S872-S874)*SIN75+S871
      SO71 = (C872-C874)*SIN75+(S872-S874)*SIN15+S871
C**** Multiply by FACTOR = 3/2
      F(27) =  CO30*SIN30-(CO31-SO31)*SIN45
      F(51) =  CO30*SIN30-CO31*SIN15-SO31*SIN75
      F(30) = (CO60+SO61)*SIN30-CO61*SIN60
      F(54) =  CO60*SIN30-SO61
      F(33) =  CO90*SIN30-CO91*SIN75+SO91*SIN15
      F(57) =  CO90*SIN30+CO91*SIN15-SO91*SIN75
      F(36) =  COC0*SIN30-COC1
      F(60) = (COC0+COC1)*SIN30-SOC1*SIN60
      F(39) =  COF0*SIN30-COF1*SIN75-SOF1*SIN15
      F(63) =  COF0*SIN30+(COF1-SOF1)*SIN45
      F(42) = (COI0-SOI1)*SIN30-COI1*SIN60
      F(66) = (COI0-SOI1)*SIN30+COI1*SIN60
      F(45) =  COL0*SIN30-(COL1+SOL1)*SIN45
      F(69) =  COL0*SIN30+COL1*SIN75-SOL1*SIN15
      F(48) = (CO00-CO01)*SIN30-SO01*SIN60
      F(72) =  CO00*SIN30+CO01
      F(25) = CO10*SIN30-CO11*SIN35+SO11*SIN55
      F(49) = CO10*SIN30-SO11*SIN65-CO11*SIN25
      F(26) = CO20*SIN30-CO21*SIN40+SO21*SIN50
      F(50) = CO20*SIN30-SO21*SIN70-CO21*SIN20
      F(28) = CO40*SIN30-CO41*SIN50+SO41*SIN40
      F(52) = CO40*SIN30-SO41*SIN80-CO41*SIN10
      F(29) = CO50*SIN30-CO51*SIN55+SO51*SIN35
      F(53) = CO50*SIN30-SO51*SIN85-CO51*SIN5
      F(31) = CO70*SIN30-CO71*SIN65+SO71*SIN25
      F(55) = CO70*SIN30-SO71*SIN85+CO71*SIN5
      F(32) = CO80*SIN30-CO81*SIN70+SO81*SIN20
      F(56) = CO80*SIN30-SO81*SIN80+CO81*SIN10
      F(34) = COA0*SIN30-COA1*SIN80+SOA1*SIN10
      F(58) = COA0*SIN30-SOA1*SIN70+COA1*SIN20
      F(35) = COB0*SIN30-COB1*SIN85+SOB1*SIN5
      F(59) = COB0*SIN30-SOB1*SIN65+COB1*SIN25
      F(37) = COD0*SIN30-COD1*SIN85-SOD1*SIN5
      F(61) = COD0*SIN30-SOD1*SIN55+COD1*SIN35
      F(38) = COE0*SIN30-COE1*SIN80-SOE1*SIN10
      F(62) = COE0*SIN30-SOE1*SIN50+COE1*SIN40
      F(40) = COG0*SIN30-COG1*SIN70-SOG1*SIN20
      F(64) = COG0*SIN30-SOG1*SIN40+COG1*SIN50
      F(41) = COH0*SIN30-COH1*SIN65-SOH1*SIN25
      F(65) = COH0*SIN30-SOH1*SIN35+COH1*SIN55
      F(43) = COJ0*SIN30-COJ1*SIN55-SOJ1*SIN35
      F(67) = COJ0*SIN30-SOJ1*SIN25+COJ1*SIN65
      F(44) = COK0*SIN30-COK1*SIN50-SOK1*SIN40
      F(68) = COK0*SIN30-SOK1*SIN20+COK1*SIN70
      F(46) = COM0*SIN30-COM1*SIN40-SOM1*SIN50
      F(70) = COM0*SIN30-SOM1*SIN10+COM1*SIN80
      F(47) = CON0*SIN30-CON1*SIN35-SON1*SIN55
      F(71) = CON0*SIN30-SON1*SIN5 +CON1*SIN85
      F( 1) = 1.5*CO10-F(25)-F(49)
      F( 2) = 1.5*CO20-F(26)-F(50)
      F( 3) = 1.5*CO30-F(27)-F(51)
      F( 4) = 1.5*CO40-F(28)-F(52)
      F( 5) = 1.5*CO50-F(29)-F(53)
      F( 6) = 1.5*CO60-F(30)-F(54)
      F( 7) = 1.5*CO70-F(31)-F(55)
      F( 8) = 1.5*CO80-F(32)-F(56)
      F( 9) = 1.5*CO90-F(33)-F(57)
      F(10) = 1.5*COA0-F(34)-F(58)
      F(11) = 1.5*COB0-F(35)-F(59)
      F(12) = 1.5*COC0-F(36)-F(60)
      F(13) = 1.5*COD0-F(37)-F(61)
      F(14) = 1.5*COE0-F(38)-F(62)
      F(15) = 1.5*COF0-F(39)-F(63)
      F(16) = 1.5*COG0-F(40)-F(64)
      F(17) = 1.5*COH0-F(41)-F(65)
      F(18) = 1.5*COI0-F(42)-F(66)
      F(19) = 1.5*COJ0-F(43)-F(67)
      F(20) = 1.5*COK0-F(44)-F(68)
      F(21) = 1.5*COL0-F(45)-F(69)
      F(22) = 1.5*COM0-F(46)-F(70)
      F(23) = 1.5*CON0-F(47)-F(71)
      F(24) = 1.5*CO00-F(48)-F(72)
      RETURN
C****
C****
C**** FFTE calculates the spectral energy E from the input grid
C**** point values F.
C****
      ENTRY FFTE (F,E)
      QENERG = .TRUE.
      GO TO 10
C****
   20 E(0)  =  (C200+C210)*(C200+C210)*BYKM2
      E(1)  = ((C201+C211)*(C201+C211)+(S211+S201)*(S211+S201))*BYKM
      E(2)  = ((C202+C212)*(C202+C212)+(S212+S202)*(S212+S202))*BYKM
      E(3)  = ((C203+C213)*(C203+C213)+(S213+S203)*(S213+S203))*BYKM
      E(4)  = ((C204+C214)*(C204+C214)+(S214+S204)*(S214+S204))*BYKM
      E(5)  = ((C205+C215)*(C205+C215)+(S215+S205)*(S215+S205))*BYKM
      E(6)  = ((C206+C216)*(C206+C216)+(S216+S206)*(S216+S206))*BYKM
      E(7)  = ((C207+C217)*(C207+C217)+(S217+S207)*(S217+S207))*BYKM
      E(8)  = ((C208+C218)*(C208+C218)+(S218+S208)*(S218+S208))*BYKM
      E(9)  = ((C409+C219)*(C409+C219)+(S429+S219)*(S429+S219))*BYKM
      E(10) = ((C20A+C21A)*(C20A+C21A)+(S21A+S20A)*(S21A+S20A))*BYKM
      E(11) = ((C20B+C21B)*(C20B+C21B)+(S21B+S20B)*(S21B+S20B))*BYKM
      E(12) = ((C20C+C21C)*(C20C+C21C)+(S21C+S20C)*(S21C+S20C))*BYKM
      E(13) = ((C20D+C21D)*(C20D+C21D)+(S21D+S20D)*(S21D+S20D))*BYKM
      E(14) = ((C20E+C21E)*(C20E+C21E)+(S21E+S20E)*(S21E+S20E))*BYKM
      E(15) = ((C20F+C21F)*(C20F+C21F)+(S21F+S20F)*(S21F+S20F))*BYKM
      E(16) = ((C20G+C21G)*(C20G+C21G)+(S21G+S20G)*(S21G+S20G))*BYKM
      E(17) = ((C20H+C21H)*(C20H+C21H)+(S21H+S20H)*(S21H+S20H))*BYKM
      E(18) = ((C400-C420)*(C400-C420)+(C410-C430)*(C410-C430))*BYKM
      E(19) = ((C20H-C21H)*(C20H-C21H)+(S21H-S20H)*(S21H-S20H))*BYKM
      E(20) = ((C20G-C21G)*(C20G-C21G)+(S21G-S20G)*(S21G-S20G))*BYKM
      E(21) = ((C20F-C21F)*(C20F-C21F)+(S21F-S20F)*(S21F-S20F))*BYKM
      E(22) = ((C20E-C21E)*(C20E-C21E)+(S21E-S20E)*(S21E-S20E))*BYKM
      E(23) = ((C20D-C21D)*(C20D-C21D)+(S21D-S20D)*(S21D-S20D))*BYKM
      E(24) = ((C20C-C21C)*(C20C-C21C)+(S21C-S20C)*(S21C-S20C))*BYKM
      E(25) = ((C20B-C21B)*(C20B-C21B)+(S21B-S20B)*(S21B-S20B))*BYKM
      E(26) = ((C20A-C21A)*(C20A-C21A)+(S21A-S20A)*(S21A-S20A))*BYKM
      E(27) = ((C409-C219)*(C409-C219)+(S219-S429)*(S219-S429))*BYKM
      E(28) = ((C208-C218)*(C208-C218)+(S218-S208)*(S218-S208))*BYKM
      E(29) = ((C207-C217)*(C207-C217)+(S217-S207)*(S217-S207))*BYKM
      E(30) = ((C206-C216)*(C206-C216)+(S216-S206)*(S216-S206))*BYKM
      E(31) = ((C205-C215)*(C205-C215)+(S215-S205)*(S215-S205))*BYKM
      E(32) = ((C204-C214)*(C204-C214)+(S214-S204)*(S214-S204))*BYKM
      E(33) = ((C203-C213)*(C203-C213)+(S213-S203)*(S213-S203))*BYKM
      E(34) = ((C202-C212)*(C202-C212)+(S212-S202)*(S212-S202))*BYKM
      E(35) = ((C201-C211)*(C201-C211)+(S211-S201)*(S211-S201))*BYKM
      E(36) =  (C200-C210)*(C200-C210)*BYKM2
C     QENERG =.FALSE.
      RETURN
C****
C****
      ENTRY FFT0 (IM)
C****
C**** Initialize sines and cosines used by FFT routines.
C****
      IF(IM.NE.72)  GO TO 130
      DO 110 N=1,17
  110 SINN(N) = SIN(TWOPI*N/72.D0)
      DO 120 N=1,35
      SINNH(N) = SIN(TWOPI*N/144.D0)
  120 COSNH(36-N) = SINNH(N)
      R2 = SQRT(2.D0)
      R3 = SQRT(3.D0)
      TSIN75 = 2.*SIN75
      TSIN15 = 2.*SIN15
      RETURN
  130 WRITE (6,*) ' This version of FFT is for 72.  IM =',IM
      STOP
      END

      SUBROUTINE FFT2 (F,A,B)
C****
C**** FFT2 calculates the Fourier coefficients for quantities
C**** defined on the secondary grid, i.e. the values in F(K) are
C**** located a half space to the right of those in FFT.
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(72),A(0:36),B(0:36)
      COMMON /FFTCOM/ SINN(17),SINNH(35),COSNH(35)
      CALL FFT (F,A,B)
      DO 10 N=1,35
      ANEW = A(N)*COSNH(N)-B(N)*SINNH(N)
      B(N) = A(N)*SINNH(N)+B(N)*COSNH(N)
   10 A(N) = ANEW
      B(36) = A(36)
      A(36) = 0.
      RETURN
      END
