C ----------------------------------------------------------------------
C                           run_maid_07.f
C                   5-fold diff. cross section
C    *****        Program Source: 2007  (18.05.07)         *****
C ----------------------------------------------------------------------
c

      program run_maid_07
      IMPLICIT REAL*8 (A-H,O-Z)
        EIMEV=4560
        EFMEV=4044.4376
        THE=7.66995
         TH=0
        PHI=0
      call DFTOT(EIMEV,EFMEV,THE,TH,PHI,S5FOLD)
      write(*,*) ' output'
      write(*,*) EIMEV,EFMEV,THE,TH,PHI,S5FOLD
      end
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      subroutine DFTOT(EIMEV,EFMEV,THE,TH,PHI,S5FOLD)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER VAR*11, UNIT*10, TEXT*64
      CHARACTER*10 MODEL(16)
      CHARACTER*20 SOLUTION
      COMPLEX*16 F1,F2,F3,F4,F5,F6,ACH,AIS,APN3,AMPL
      COMPLEX*16 A1,A2,A3,A4,A5,A6,H1,H2,H3,H4,H5,H6
      DIMENSION ARG(100),AMPL(100,7),ECM(100),ELAB(100)
      DIMENSION QPICM(100)
      DIMENSION ACH(10,8,4),AIS(10,8,3),APN3(10,8,3)
      COMMON/HELIC/ SIG1, SIG3, STL
      COMMON/QQMAX/ Q2MAX
      COMMON /SOLUTION/ XX(100),ISOL
c ***************************************************************+
      common/newres/ S31, P13, P31, D15, F35, F37
c *****************************************************************
      DATA MODEL/'Born ','Rho ','Omega ','P33(1232)', 'P11(1440)',
     * 'D13(1520) ','S11(1535) ','S31(1620)','S11(1650) ','F15(1680)'
     *,'D33(1700) ','P13(1720) ','P31(1910)','D15(1675) ','F35(1905)'
     *,'F37(1950) '/
c
c      OPEN(6,File='run_maid_07.out',form='formatted',status='unknown')
c
       Q2MAX=3.D0
	HQC=197.3285D0
        PI=3.1415926536D0
	AMP=938.2723D0/HQC
        AMN=939.5653D0/HQC
        AMPI0=134.9642D0/HQC
	AMPIP=139.5685D0/HQC
        LMAX=5
        IMULT=0
       Solution="maid07_final"
       ISOL=1
c       PRINT *, 'channel  1 - pi0 p; 2 - pi0 n; 3 - pi+ n; 4 - pi- p'
       ISO=1
C ******************************
        GO TO (31,32,33,34) ISO

  31     AMI=AMP
         AMF=AMP
         AM=AMP
         API=AMPI0
         GO TO 35
  32     AMI=AMN
         AMF=AMN
         AM=AMN
         API=AMPI0
         GO TO 35
  33     AMI=AMP
         AMF=AMN
         AM=(AMI+AMF)/2.
         API=AMPIP
         GO TO 35
  34     AMI=AMN
         AMF=AMP
         AM=(AMI+AMF)/2.
         API=AMPIP

  35    WTHR=AMF+API
        WTHR2=WTHR**2
c ***********************************************************************
       VEC=1
c ***********************************************************
       BORN=1
       RHO=1
       OMEGA=1
       P33=1
       P11=1
       D13=1
       S11F=1
       S31=1
       S11S=1
       F15=1
       D33=1
       P13=1
       P31=1
       D15=1
       F35=1
       F37=1
       XE=1
       XS=1
       XMIX=1
       X3P33=1
       X1P33=1
       XSP33=1
       X1S31=1
       XSS31=1
       X3D33=1
       X1D33=1
       XSD33=1 
       X1P31=1
       XSP31=1
       X3F35=1
       X1F35=1
       XSF35=1
       X3F37=1
       X1F37=1
       XSF37=1
c ***************************************************
       IF (ISO.EQ.1.OR.ISO.EQ.3) THEN
          X1P11p=1
          XSP11p=1
          X1S11p=1
          XSS11p=1
          X1S2p=1
          XSS2p=1 
          X3D13p=1
          X1D13p=1
          XSD13p=1
          X3F15p=1
          X1F15p=1
          XSF15p=1
          X3P13p=1
          X1P13p=1
          XSP13p=1
          X3D15p=1
          X1D15p=1
          XSD15p=1
      STD=BORN*RHO*OMEGA*P33*P11*D13*D33*S11F*S11S*F15*S31*
     &  P13*P31*D15*F35*F37*XE*XS*XMIX*
     &  X3P33*X1P33*XSP33*X1S31*XSS31*X3D33*X1D33*XSD33*
     &  X1P31*XSP31*X3F35*X1F35*XSF35*X3F37*X1F37*XSF37*
     &  X1P11p*XSP11p*X1S11p*XSS11p*X1S2p*XSS2p*
     &  X3D13p*X1D13p*XSD13p*X3F15p*X1F15p*XSF15p*
     &  X3P13p*X1P13p*XSP13p*X3D15p*X1D15p*XSD15p

       ELSE IF (ISO.EQ.2.OR.ISO.EQ.4) THEN
          X1P11n=1
          XSP11n=1
          X1S11n=1
          XSS11n=1
          X1S2n=1
          XSS2n=1
          X3D13n=1
          X1D13n=1
          XSD13n=1
          X3F15n=1
          X1F15n=1
          XSF15n=1
          X3P13n=1
          X1P13n=1
          XSP13n=1
          X3D15n=1
          X1D15n=1
          XSD15n=1
      STD=BORN*RHO*OMEGA*P33*P11*D13*D33*S11F*S11S*F15*S31*
     &  P13*P31*D15*F35*F37*XE*XS*XMIX*
     &  X3P33*X1P33*XSP33*X1S31*XSS31*X3D33*X1D33*XSD33*
     &  X1P31*XSP31*X3F35*X1F35*XSF35*X3F37*X1F37*XSF37*
     &  X1P11n*XSP11n*X1S11n*XSS11n*X1S2n*XSS2n*
     &  X3D13n*X1D13n*XSD13n*X3F15n*X1F15n*XSF15n*
     &  X3P13n*X1P13n*XSP13n*X3D15n*X1D15n*XSD15n

       END IF
c ****************************************************************
       CALL SOLUTIONS(SOLUTION,TEXT,XE,XS,XMIX,
     &  X3P33,X1P33,XSP33,X1S31,XSS31,X3D33,X1D33,XSD33,
     &  X1P31,XSP31,X3F35,X1F35,XSF35,X3F37,X1F37,XSF37,
     &  X1P11p,XSP11p,X1S11p,XSS11p,X1S2p,XSS2p,
     &  X3D13p,X1D13p,XSD13p,X3F15p,X1F15p,XSF15p,
     &  X3P13p,X1P13p,XSP13p,X3D15p,X1D15p,XSD15p,
     &  X1P11n,XSP11n,X1S11n,XSS11n,X1S2n,XSS2n,
     &  X3D13n,X1D13n,XSD13n,X3F15n,X1F15n,XSF15n,
     &  X3P13n,X1P13n,XSP13n,X3D15n,X1D15n,XSD15n)
c ***************************************************
      IUNI=IDINT(P33+P11+D13+D33+S11F+S11S+F15+S31+
     + P13+P31+D15+F35+F37  )
c ***************************************************
c*********************************************************************
	WMEVMAX=2000.0
c     	IF (ISOL .GT. 1) WMEVMAX=1790.0
	Q2GEVMAX=5.0



        HEL=0
        IVAR=2
        

c *********** Kinematics for the electrons *********
       IF (THE.LT.0.D0.OR.THE.GT.180.D0) GO TO 2000
       IF (EFMEV.GE.EIMEV) GO TO 2000
       EI=EIMEV/HQC
       EF=EFMEV/HQC
       EGL=EI-EF
       XE=DCOS(THE*PI/180.D0)
       TNE = DTAN(THE*PI/180./2.)
       QG2 = EI**2 + EF**2 -2.*EI*EF*XE
       Q2=2.*EI*EF*(1.-XE)
       Q2GEV=Q2*(HQC/1000.)**2
c *************************************************
       IF (Q2.LE.0.D0) GO TO 2000
       EPS=1./(1.+2.*QG2/Q2*TNE**2)
       WFM2 = 2.*AM*EGL + AM**2 - Q2
       IF (WFM2.LE.WTHR2) GO TO 2100
       WFM=DSQRT(WFM2)
       WMEV=WFM*HQC
       EGLEQ=(WFM**2-AM**2)/2./AM
       GMEV=(EF/EI)*(EGLEQ/Q2)/(2.*PI**2*137.*(1.-EPS))/HQC
C *******************************************
       IF (Q2GEV.GT.Q2GEVMAX.OR.WMEV.GT.WMEVMAX) GO TO 2200
      THPI=TH
      X=DCOS(THPI*PI/180.D0)
      PHPI=PHI
      CSF=DCOS(PHPI*PI/180.D0)
      SNF=DSIN(PHPI*PI/180.D0)
      CS2F=DCOS(2.*PHPI*PI/180.D0)
       CALL MAID(ISO,WFM,Q2,X,QPIAV,EGVCM,EGVLAB,
     & F1,F2,F3,F4,F5,F6,A1,A2,A3,A4,A5,A6,H1,H2,H3,H4,H5,H6,
     & IMULT,LMAX,ACH,AIS,APN3,
     & BORN,VEC,OMEGA,RHO,P33,P11,D13,S11F,S11S,F15,D33)
      EGEQ=(WFM**2-AMI**2)/2./WFM
      EGCM=(WFM**2-Q2-AMI**2)/2./WFM
      EGLAB=(WFM**2+Q2-AMI**2)/2./AMI
      EPI=(WFM**2+API**2-AMF**2)/2./WFM
      QPI=SQRT(EPI**2-API**2)

       CALL OBSERV(WFM,Q2,QPI,EGEQ,EGVCM,ST,SL,STL,STT,STLP,
     & H1,H2,H3,H4,H5,H6)
       SVIR=ST + EPS*SL + DSQRT(2.*EPS*(1.+EPS))*CSF*STL +
     + EPS*CS2F*STT + HEL*DSQRT(2.*EPS*(1.-EPS))*SNF*STLP

       S5FOLD=GMEV*SVIR

      GO TO 3000
c *************** Error messages **********************

2000  WRITE (*,2001)
2001  FORMAT(/,1X,
     * '******  W R O N G   K I N E M A T I C ! ******')
      GO TO 3000
2100  WRITE (*,2002)
2002  FORMAT(/,1X,
     * '******  B E L O W  T H R E S H O L D ! ******')
      GO TO 3000
2200  WRITE (*,2003) Q2GEVMAX,WTHR,WMEVMAX
2003  FORMAT(/,1X,'****** Q2 or W  out of limit (0-',F3.1,
     *         ' or ',F6.1,'-',F5.0,') ******')
      GO TO 3000
2300  WRITE (*,2004)
2004  FORMAT(/,1X,
     * '****** Electron helicity h is wrong ! ******')
3000  CONTINUE
c      close (6)
      return
      END


      SUBROUTINE OBSERV(WCM,Q2,QPI,EGEQ,EGCM,ST,SL,STL,STT,STLP,
     & H1,H2,H3,H4,H5,H6)
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON /F16H16/ F1,F2,F3,F4,F5,F6,H1,H2,H3,H4,H5,H6
      COMPLEX*16 F1,F2,F3,F4,F5,F6,H1,H2,H3,H4,H5,H6

C  *******  Unpolarized cross sections (in mcb/sr) *************

      FACT=10000.*QPI/EGEQ
      SQ2=DSQRT(2.D0)
c
      RT=DREAL(H1*DCONJG(H1)+H2*DCONJG(H2)+H3*DCONJG(H3)
     + +H4*DCONJG(H4))/2.
      ST=FACT*RT
c
      RL=DREAL(H5*DCONJG(H5)+H6*DCONJG(H6))
      SL=FACT*RL*Q2/(EGCM**2)
c
      RTL=DREAL((H1-H4)*DCONJG(H5)+(H2+H3)*DCONJG(H6))/SQ2
      STL=FACT*RTL*DSQRT(Q2)/EGCM
c
      RTLP=-DIMAG((H1-H4)*DCONJG(H5)+(H2+H3)*DCONJG(H6))/SQ2
      STLP=FACT*RTLP*DSQRT(Q2)/EGCM
c
      RTT=DREAL(H3*DCONJG(H2)-H4*DCONJG(H1))
      STT=FACT*RTT

      RETURN
      END

C!!!!!!!!!!!!! NEW with helicity amplitudes  !!!!!!!!!!!
C -------------------------- M A I D 2007  ----------------------
C          Pion Photo and Electroproduction on the Nucleon
C          with Born + Vector mesons: Omega, Rho
C          and Resonances: P33(M1+,E1+,L1+), P11(M1-,1-)
C                          D13, D33 (M2-,E2-), S11(E0+,L0+), F15(M3-,E3-)
C          D15 (E2+,M2+), P13(E1+,M1+), P31(M1-), F35(E3-,M3-), F37(E3+,M3+)
C    *****        Program Source: 2007  (18.05.07)         *****
C    *****              Version 01 (19.06.07)              *****
C ----------------------------------------------------------------------
      SUBROUTINE MAID(ISO,WCM,Q20,X,QPI,EGCM,EGLAB,
     & F1,F2,F3,F4,F5,F6,S1,S2,S3,S4,S5,S6,H1,H2,H3,H4,H5,H6,
     & IMULT,LMAX,ACH,AIS,APN3,
     & BORN,VEC,OMEGA,RHO,P33,P11,D13,S11F,S11S,F15,D33)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 F1,F2,F3,F4,F5,F6,MULTI,ACH,AIS,APN3,FPL,FMN,F0
      COMPLEX*16 H1,H2,H3,H4,H5,H6,A1,A2,A3,A4,A5,A6
	COMPLEX*16 S1,S2,S3,S4,S5,S6
      COMPLEX*16 F_back(6), MULT, AP0M(10,8,3), ATST(6,10,3)
      REAL*8 KAPPAP,KAPPAN,MPI,MPI0,MPIP,MP,MN,M,M1,M2,M1P,M1M
      INTEGER GAUGE, D13MODE
      DIMENSION AIS0(10,6,3)
      DIMENSION ACH(10,8,4),AIS(10,8,3),APN3(10,8,3)
      COMMON /MULTB/MULT(6,10,3), FPOLE(6,10,3), FBV(6,10,3)
      COMMON/COUPL/ IPVPS
      COMMON/GAUSS/EX(100),WX(100),NX
      COMMON/FORMF/ xlampi,xlamax,mode,gauge
      COMMON/bl/ PI,PIH,HQC,F2PI,E,KAPPAP,KAPPAN,CSUC,CT,COME,CRHO,
     & CBORN,CVEC,ALAMBDA,MPI,MPI0,MPIP,MP,MN,AVM,M,M1,M2,IPV,IFORM,
     & CP33,CP11,CD13,CD33,CS11F,CS11S,CF15,CS11T,CP11S
      common/newres/ S31, P13, P31, D15, F35, F37
      common/newpar/JRHOFORM,JOMEFORM,D13mode
      COMMON /E0CORR/ XE,XS,XMIX
      COMMON /SOLUTION/ XRES(100), ISOL
	COMMON /INVARI/ S,T,A1,A2,A3,A4,A5,A6
      DIMENSION OBS(10)
      DATA IFST/0/
      SAVE IFST
c *************  INPUT ***************
c  ISO :  1  for pi0 + proton   :  2  for pi0 + neutron
c      :  3  for pi+ + neutron  :  4  for pi- + proton
c  WCM :  Total pi N c.m. energy (in 1/fm) (WCM<1.7 GeV)
c  Q2  :  Square of virtual photon 4-mometum (in 1/fm/fm)
c      :  (0=< Q2 < 3 (GeV/c)^2)
c  X   : cos(Theta_pi) in the piN c.m. frame
c**************  OUTPUT  *************
c   F1, ....,F6  :  CGLN amplitudes (in fm). Convention as in
c                : ref. D. Drechsel and L. Tiator,
c                : J.Phys.G: Nucl. Phys. 18 (1992) 449.
c                : F5=F8(Berends)*w/q ,  F6=F7(Berends)*w/q
c    EGCM        : virtual photon energy (in 1/fm). It is useful
c                : to use this value for the Lorentz tranformations of
c                : of the longitudinal amplitudes F5, F6 in order
c                : to avoid spurious singularities
c ************   PARAMETERS ***************************************
c   CBORN,CVEC,COME,CRHO  :  1 or 0 :  with or without Born, Vect-mesons
c   ISO              :  1  for pi0 + proton   :  2  for pi0 + neutron
c                    :  3  for pi+ + neutron  :  4  for pi- + proton
c   MODE             :  0  use Sachs form factors in Dipole form (std)
c                    :  1  all Dirac form factors = 1
c                    :  3  all Dirac form factors in Dipole form (diff from 0)
c                    :  2  Fpi and Fax different F1 and F2 as in MODE=0
c  !!!!!!!!  Please, use option MODE=0  and Gauge=1 or 2  !!!!!!!!!!
c   GAUGE            :  1  Coulomb Gauge (only Jz is used)
c                    :  2  Lorentz Gauge
c                    :  3  Time (axial)  (only Rho is used)
c   IPVPS            :  0   PV coupling in the Born terms
c                    :  1   PS coupling in the Born terms
c                    :  2   mixed coupling (realistic)
c                    :  1   unitarized E0+ (realistic)
c CP33,....,CF15     :  1 or 0 with or without corresponding resonances
c *****************************************************************************
c   FTNORM           :  arbitrary renormalization, not used in MAID2000
C   IFORM=JFORM      :  0 or 1 without or with hadronic pion ff
C                    :  =0 in MAID2000
c*************************new March 2002***************************************
c   FRHOFORM         :  0 or 1 without or with hadronic rho ff
c   FOMEFORM         :  0 or 1 without or with hadronic omega ff
c   D13MODE          :  energy dependence for D13 gamma vertex
c                       1 or 2  for Maid2000 or Maid2002
c******************************************************************************

c     test with mode=2 (different form factors) (std:mode=0)
c    with MODE>0 also Lorentz gauge GAUGE=2 must be used!!!
c        MODE=0
c        GAUGE=2.
c  standard case mode=7
        MODE=7
        GAUGE=2.        
c        MODE=0
        IDISP=2

      IF (IDISP .EQ. 0) THEN
        IMAID=0
        IPOLE=1
      ELSE IF (IDISP .EQ. 1) THEN
        IMAID=1
        IPOLE=-1
      ELSE IF (IDISP .EQ. 2) THEN
        IMAID=1
        IPOLE=0
      ENDIF

c ****************************
      IF (IFST.EQ.1) GO TO 999
c	write(6,90) idisp,mode,gauge
c      IF (mode.gt.5) write(6,91) xlamax,xlampi
c  90  format('  Maid2007: disp =',I2,
c     &   ', e.m. form factors: mode =',i2,', gauge =',i2,/,1x,80(1H-))
c  91  format(10x,' Lambda_A =',F8.4,'  Lambda_pi =',F8.4)
c ***********  Interpolation of phases and inelasticities  *******************
      open (2, file = 'piN_2500.dat', status = 'old')
c      open (2, file = 'piN_2500fa02.dat', status = 'old')
      CALL SPLINES()
       IFST=1
c      CLOSE (2)
999   CONTINUE
C *********  Problem at Q2=0 ***********
       Q2=Q20
       IF (DABS(Q20).LE.1.E-6) Q2=1.E-6
c  *****  Q2 should not be much smaller than 1.E-6
C******************************************************************************

        IF (ISO .EQ. 0) RETURN
	  HQC=197.3285D0
	  MP=938.2723D0/HQC
        MN=939.5653D0/HQC
        AVM = (MP + MN) / 2.0
	  MPIP=139.5685D0/HQC
        MPI0=134.9642D0/HQC
	  F2PI=0.079D0
C         f2pi=0.079  yields  g2pi=14.28
c         PI=4.D0*ATAN(1.D0)
        PI=3.1415926536D0
        PIH=PI/2.D0
        E=DSQRT(4.D0*PI/137.D0)
	  KAPPAP= 1.7928
	  KAPPAN=-1.9130
c ****************************
        FTNORM=1.
        IPVPS=2
	IF (DABS(XMIX).LE.1.E-06) IPVPS=1
	IF (DABS(XMIX).GE.1000.) IPVPS=0
C ****************************
        COME=OMEGA
        CRHO=RHO
        CBORN=BORN
        CVEC=VEC
        CP33=P33
        CP11=P11
        P11S=0
c        P11S=1
	  CP11S=P11S
        CD13=D13
        CD33=D33
        CS11F=S11F
        CS11S=S11S
c	   S11T=S11S
c	   S11T=0
        S11T=1
	  CS11T=S11T
        CF15=F15
c std. MAID2007web is without 2nd P11 and with third S11: P11S=0, S11T=1
c	LMAX=3
        LMAX1=LMAX+1
c ************************
       IF (ISO.EQ.1) THEN
         M1=MP
         M2=MP
         MPI=MPI0
       ELSE IF (ISO.EQ.2) THEN
         M1=MN
         M2=MN
         MPI=MPI0
       ELSE IF (ISO.EQ.3) THEN
         M1=MP
         M2=MN
         MPI=MPIP
       ELSE IF (ISO.EQ.4) THEN
         M1=MN
         M2=MP
         MPI=MPIP
	   IF (WCM*HQC .LT. 1080) M1=MP
       ENDIF
C *********  N E W  20.09.98  **********
         M1=(M1+M2)/2.
         M2=M1
c  for consistency between multipoles and CGLN amplitudes (Dez 04/Jan 05)
           M1=MP
           M2=MP
           MPI=MPI0
C *********************************************************
      CALL NACHALO(WCM,Q2,M1,MPI,M2,CK0,CK,Q1,W1)
      EGCM=CK0
      QPI=Q1
      EGLAB=(WCM**2-M1**2+Q2)/2./M1
        IF (W1.EQ.0) THEN
         WTHR=M2+MPI
	 Q2GEV=Q2*(197.3285D0/1000.D0)**2
         WRITE(1,81) WCM*HQC,Q2GEV,WTHR*HQC
   81    FORMAT(' Kinematik: wcm,q2 ',F7.1,F7.4,' below threshold: ',
     &        ' wthr =',F7.2)
        ENDIF
C ************************
      F1=(0.D0,0.D0)
      F2=(0.D0,0.D0)
      F3=(0.D0,0.D0)
      F4=(0.D0,0.D0)
      F5=(0.D0,0.D0)
      F6=(0.D0,0.D0)
      H1=(0.D0,0.D0)
      H2=(0.D0,0.D0)
      H3=(0.D0,0.D0)
      H4=(0.D0,0.D0)
      H5=(0.D0,0.D0)
      H6=(0.D0,0.D0)

        DO 99  L1=1,10
        DO 99  K=1,8
        DO 98  I=1,3
        APN3(L1,K,I)=(0.D0,0.D0)
        AIS(L1,K,I)=(0.D0,0.D0)
 98     ACH(L1,K,I)=(0.D0,0.D0)
 99     ACH(L1,K,4)=(0.D0,0.D0)
      IF (W1.EQ.0) RETURN
C ****************************************************************
C   unitarization phases, CGLN amplitudes from background and res
c   additional terms for background unitarization
c ****************************************************************
      DO 101 K=1,6
 101  F_back(K)=(0.D0,0.D0)
      IF (CBORN.EQ.0 .AND. CRHO.EQ.0 .AND. COME.EQ.0) THEN
c ******************************  TEST !!!!!!!!!
	IBGR=1
*************************************************
	ELSE
	IBGR=1
	ENDIF
      IUNI=IDINT(P33+P11+D13+D33+S11F+S11S+F15+S31+P13+P31+
     + D15+F35+F37)
      IF (ISOL.EQ.1)
     & CALL MULT_TOT(WCM,Q2,M1,MPI,CK0,lmax1,AP0M,AIS,APN3,
     & ACH,AIS0,IUNI,IBGR)
      IF (IUNI.EQ.0) GO TO 113
c*****************************************************************
c       IF (ISOL.EQ.2)
c     & CALL INTAPRX(WCM,Q2,LMAX1,1,AP0M,AIS,APN3,ACH,AIS0,CK0)
c       IF (ISOL.EQ.3)
c     & CALL INTEXACT(WCM,Q2,LMAX1,1,AP0M,AIS,APN3,ACH,AIS0,CK0)
c ****************************************************************
      CALL NACHALO(WCM,Q2,M1,MPI,M2,CK0,CK,Q1,W1)
      CALL CGLN_DISP(ISO, WCM, Q2, X, m1, MPI, AIS, AIS0, F_back)
113   CALL CGLN(ISO,WCM,Q2,CK0,X,F1c,F2c,F3c,F4c,F5c,F6c)
      F1 = F1c + F_back(1)
      F2 = F2c + F_back(2)
      F3 = F3c + F_back(3)
      F4 = F4c + F_back(4)
      F5 = F5c + F_back(5)
      F6 = F6c + F_back(6)
c *****************  TEST LET **************
c       APN3(2,1,3) =  APN3(2,1,3) / CK
c       APN3(2,5,3) =  APN3(2,5,3) / CK     
c *******************************************
      DO 997 K=7,8
      DO 997 L1=1,LMAX1
      DO 998 I=1,3
      AIS(L1,K,I)=AIS(L1,K-2,I)*CK/CK0
      APN3(L1,K,I)=APN3(L1,K-2,I)*CK/CK0
998   ACH(L1,K,I)=ACH(L1,K-2,I)*CK/CK0
997   ACH(L1,K,4)=ACH(L1,K-2,4)*CK/CK0
C ******** extra line added for I=4 ***********
C ***********   Truncated Amplitudes  ***********

c      Ltr=0
c *******  amplitudes will only be changed if Lmax<6
c      Lt0=0
c      Call CGLN_TRUNC(ISO,Lt0,W,Q2,x,MN,MPI,AIS,F1,F2,F3,F4,F5,F6)

C ***********   HELICITY AMPLITUDES ***********
      XX=X
      YY=DSQRT(1-XX**2)
      X2Y=XX*XX-YY*YY
      SQ2=DSQRT(2.D0)
c *** old spin amplitudes of Knoechlein et al. used for total c.s. *****
c ******* recoil polarization along the incoming photon  *****
      S1=-YY*(F3+F4*XX)/SQ2
      S2=-(2.*F1-2.*F2*XX+F4*YY*YY)/SQ2
      S3=-YY**2*F4/SQ2
      S4=YY*(2.*F2+F3+F4*XX)/SQ2
      S5=F5+F6*XX
      S6=F6*YY
      XX2=DSQRT((1+XX)/2.D0)
      YY2=DSQRT((1-XX)/2.D0)
c ****** new helicity amplitudes (21.3.99)  ****************
c ****** recoil polarization along the outgoing pion  ******
      H1=-1/SQ2*YY*XX2*(F3+F4)
      H2=SQ2*XX2*(F2-F1+(F3-F4)*(1-XX)/2.D0)
      H3=1/SQ2*YY*YY2*(F3-F4)
      H4=SQ2*YY2*(F1+F2+(F3+F4)*(1+XX)/2.D0)
      H5= XX2*(F5+F6)
      H6=-YY2*(F5-F6)
c ***********************************************
c ***  invariant amplitudes (Pasquini notes, 2003)
      CALL InvFtoA(WCM,M1,M2,MPI,Q2,CK0,CK,QPI,X,S,T,F1,F2,F3,F4,F5,F6,  
     &             A1,A2,A3,A4,A5,A6)  
      RETURN
      END
c
c
      SUBROUTINE NACHALO(WCM,Q2,XM1,XMPI,XM2,CK0,CK,CKP,ECPI)
C******************************************************************************
C     pion cm momentum and energy
C******************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KAPPAP,KAPPAN,MPIP,MPI0,MPI,MP,MN,M,M1,M2
      COMMON/bl/ PI,PIH,HQC,F2PI,E,KAPPAP,KAPPAN,CSUC,CT,COME,CRHO,
     & CBORN,CVEC,ALAMBDA,MPI,MPI0,MPIP,MP,MN,M,AVM,M1,M2,IPV,IFORM,
     & CP33,CP11,CD13,CD33,CS11F,CS11S,CF15,CS11T,CP11S
        INTEGER IBORN,IVEC,IOME,IRHO,JFORM,IP33,IP11,ID13,IS11,IF15
        INTEGER GAUGE, IS11F, IS11S, ID33
      real*8 M1FM,M2FM,MPIFM
        COMMON /FORMF/ xlampi,xlamax,mode,gauge
        COMMON /LINK1/ M1FM,M2FM,MPIFM,IBORN,IVEC,IOME,IRHO,JFORM,
     &                IP33,IP11,ID13,ID33,IS11F,IS11S,IF15,IS11T,IP11S
C************************************
      IBORN=INT(CBORN+0.1)
      IVEC=INT(CVEC+0.1)
      IOME=INT(COME+0.1)
      IRHO=INT(CRHO+0.1)
      IP33=INT(CP33+0.1)
      IP11=INT(CP11+0.1)
      IP11S=INT(CP11S+0.1)
      ID13=INT(CD13+0.1)
      ID33=INT(CD33+0.1)
      IS11F=INT(CS11F+0.1)
      IS11S=INT(CS11S+0.1)
      IS11T=INT(CS11T+0.1)
      IF15=INT(CF15+0.1)
      JFORM = 0
c      JFORM=IFORM
      M1FM=M1
      M2FM=M2
      MPIFM=MPI
c *******************************

      S=WCM**2
      CK0=(S-Q2-XM1**2)/(2*WCM)
      CK=SQRT(CK0**2+Q2)
      IF (WCM.GT.XM2+XMPI) THEN
        ECPI=(S+XMPI*XMPI-XM2*XM2)/(2*WCM)
        CKP=SQRT(ECPI*ECPI-XMPI*XMPI)
      ELSE
        PRINT *,' below threshold '
        ECPI=0
        CKP=0
      ENDIF
      RETURN
      END

        SUBROUTINE CGLN(ISO,WCM,Q2FM,WGCM,XX,F1,F2,F3,F4,F5,F6)
C******************************************************************************
C   CGLN amplitudes for Born terms, Vector mesons
C******************************************************************************
      IMPLICIT NONE
      INTEGER ISO,IS11F,IS11S,IS11T,IS31,IP13,IP31,ID15,IF35,IF37,IP11S
      DOUBLE PRECISION WCM,Q2FM,XX,WGCM, F1,F2,F3,F4,F5,F6
      INTEGER IBORN,IVEC,IOME,IRHO,JFORM,IP33,IP11,ID13,ID33,IF15
      DOUBLE PRECISION M1FM,M2FM,MPIFM, MPIP
      COMMON /LINK1/ M1FM,M2FM,MPIFM,IBORN,IVEC,IOME,IRHO,JFORM,
     &                IP33,IP11,ID13,ID33,IS11F,IS11S,IF15,IS11T,IP11S
    	DOUBLE PRECISION OmegL, Q2, W, qcm, kcm, Q2G, wGacm
	COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
    	DOUBLE PRECISION mi, mf, mPi, FBV, CBV
	DIMENSION FBV(6,2), CBV(6)
	COMMON /Mass/ mi, mf, mPi
	INTEGER Reaction, k, l
        DOUBLE PRECISION F1b,F2b,F3b,F4b,F5b,F6b,F7b,F8b
        DOUBLE PRECISION F1v,F2v,F3v,F4v,F5v,F6v,kgcm,hqc,Pi
	DIMENSION F1b(2),F2b(2),F3b(2),F4b(2),F5b(2),
     * F6b(2),F7b(2),F8b(2)
	DIMENSION F1v(2),F2v(2),F3v(2),F4v(2),F5v(2),F6v(2)
	PARAMETER ( Pi = 3.1415926536D0, hqc = 197.3285)
c
	MPIP=139.5685D0/HQC
        Reaction=1+MOD(ISO+1,4)
        W=WCM*hqc
        Q2=Q2fm*hqc**2
c  same nucleon mass is necessary for gauge invariance
        mi=(m1fm+m2fm)/2.*hqc
        mf=mi
        mPi=mpifm*hqc
	OmegL = (W*W+Q2-mi*mi)/2./mi
	kcm = SQRT(((W*W+mPi*mPi-mf*mf)/(2*W))**2-mPi*mPi)
        kgcm = (W*W-mi*mi)/2./W
        wGacm=WGCM*hqc
        qcm=sqrt(wGacm**2+Q2)
c ********************************
	DO 99 k=1,6
	CBV(k)=0.
	DO 99 l=1,2
99	FBV(k,l)=0.

        DO 109 l=1,2
        F1b(l)=0
        F2b(l)=0
        F3b(l)=0
        F4b(l)=0
        F5b(l)=0
        F6b(l)=0
        F7b(l)=0
        F8b(l)=0
        F1v(l)=0
        F2v(l)=0
        F3v(l)=0
        F4v(l)=0
        F5v(l)=0
        F6v(l)=0
  109   CONTINUE
C ************************************************************************
C         standard Born terms (PS-PV) and vector meson contributions
C         without unitarization
C ************************************************************************
       IF(IBORN.EQ.1) CALL BORN2006(Reaction,xx,F1b,F2b,F3b,F4b,F5b,F6b)
       IF(IVEC.EQ.1) CALL VECTNEW(Reaction,xx,F1v,F2v,F3v,F4v,F5v,F6v)
C ************************************************************************

        DO 108 l=1,2
        FBV(1,l)=F1b(l)+F1v(l)
        FBV(2,l)=F2b(l)+F2v(l)
        FBV(3,l)=F3b(l)+F3v(l)
        FBV(4,l)=F4b(l)+F4v(l)
        FBV(5,l)=F5b(l)+F5v(l)
        FBV(6,l)=F6b(l)+F6v(l)
  108   CONTINUE
	DO 111 l = 1, 6
          IF (Reaction.EQ.1) THEN
	    CBV(l) = -FBV(l,2)/SQRT(3.D0)-FBV(l,1)*SQRT(2.D0/3.D0)
          ELSE IF (Reaction.EQ.2) THEN
	    CBV(l) = FBV(l,2)/SQRT(3.D0)-FBV(l,1)*SQRT(2.D0/3.D0)
          ELSE IF (Reaction.EQ.3) THEN
	    CBV(l) = FBV(l,2)*SQRT(2.D0/3.D0)-FBV(l,1)/SQRT(3.D0)
          ELSE IF (Reaction.EQ.4) THEN
	    CBV(l) = FBV(l,2)*SQRT(2.D0/3.D0)+FBV(l,1)/SQRT(3.D0)
	  END IF
 111	CONTINUE
        F1=CBV(1)/MPIP
        F2=CBV(2)/MPIP
        F3=CBV(3)/MPIP
        F4=CBV(4)/MPIP
        F5=CBV(5)/MPIP
        F6=CBV(6)/MPIP

        RETURN
    	END

 	SUBROUTINE BORN2006 (Reaction, CosTh, F1, F2, F3, F4, F5, F6)
c  *****************************************************************************        
c    ***    Born  contributions  pole terms + PV-PS mixing (June 2006)   ***
c               Ai(+) = Akp,   Ai(0) = Ak0,   Ai(-) = Akm
c    Result is given as   F1(1) = F1^p  and F1(2) = F1^Delta  for Reaction=1        
c                   and   F1(1) = F1^n  and F1(2) = F1^Delta  for Reaction=2
c    notation as in Hanstein, Diploma thesis 1993, page 13, eqs. 2.29-2.30
c    all energies are in MeV, Q2 in MeV^2 and CGLN amplitudes in 1/mPi  
c  *****************************************************************************        
      IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER Reaction
	REAL*8 mi,mf,mPi,mn,kcm,kqv
	REAL*8 Kopplung,mup,mun
	REAL*8 G(3,8),F1(2),F2(2),F3(2),F4(2),F5(2),F6(2)
	COMPLEX*16 A(6,2),F1c(2),F2c(2),F3c(2),F4c(2),F5c(2),F6c(2)
	COMMON /Mass/ mi, mf, mPi
	COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
      COMMON /COUPL/ IPVPS
      COMMON /E0CORR/ XA, XB, XMIX
	PARAMETER (mup = 2.7928, mun = -1.913, apvps0=450.)
	PARAMETER (Pi=3.1415926536D0,hqc=197.3285D0,
     * eLadung=.3028619D0,Kopplung=0.996365D0)
	mn = (mi+mf)/2.D0
      qcm=sqrt(wGacm**2+Q2)
	Ei = SQRT(mi*mi+qcm**2)
	wPicm = SQRT(mPi*mPi+kcm*kcm)
	kqv = kcm*qcm*CosTh
	u = mi*mi+mPi*mPi-2.*(Ei*wPicm+kqv)
	t = mPi*mPi-Q2-2.*(wGacm*wPicm-kqv)
      s = W**2
       amn=mn
	 ampi=mPi
       EEM=eLadung
	 GPI=(2*mn/mPi)*Kopplung * mPi
	 XKAPS=mup-1+mun
	 XKAPV=mup-1-mun
       Q2GEV=Q2*1.0D-6
	 CALL FORM (Q2GEV, D, F1S, F1V, F2S, F2V, FPI, FAX)
c ***************  test  for Dirac Currents only *********
c	F2S=0
c	F2V=0
c ***************  test  *********************************
  	a1s=EEM*GPI/2.D0*1./(s-amn**2)
  	a1u=EEM*GPI/2.D0*1./(u-amn**2)
          A10=(a1s+a1u)*F1S
          A1P=(a1s+a1u)*F1V
	    A1M=(a1s-a1u)*F1V
  	a2s=-EEM*GPI/(t-ampi**2)*1./(s-amn**2)
  	a2u=-EEM*GPI/(t-ampi**2)*1./(u-amn**2)
          A20=(a2s+a2u)*F1S
          A2P=(a2s+a2u)*F1V
          A2M=(a2s-a2u)*F1V
    	a3s=-EEM*GPI/(4*amn)*1./(s-amn**2)
    	a3u=-EEM*GPI/(4*amn)*1./(u-amn**2)
          A30=(a3s-a3u)*F2S
          A3P=(a3s-a3u)*F2V
          A3M=(a3s+a3u)*F2V
   	a4s=-EEM*GPI/(4*amn)*1./(s-amn**2)
   	a4u=-EEM*GPI/(4*amn)*1./(u-amn**2)
          A40=(a4s+a4u)*F2S
          A4P=(a4s+a4u)*F2V
          A4M=(a4s-a4u)*F2V
  	a5s=-EEM*GPI/2.D0/(t-ampi**2)*1./(s-amn**2)
  	a5u=-EEM*GPI/2.D0/(t-ampi**2)*1./(u-amn**2)
          A50=(a5s-a5u)*F1S
          A5P=(a5s-a5u)*F1V
          A5M=(a5s+a5u)*F1V + 2*EEM*GPI/(t-ampi**2)*(FPI-F1V)/Q2
          A60=0
          A6P=0
          A6M=0

c ******************  PV and MIXED couplings ************
      if (ipvps.ne.1) then
       if (ipvps.eq.0) then
	  fpvps=1.
	 else
	  apvps=apvps0*XMIX
        fpvps=apvps**2/(apvps**2+kcm**2)
       endif
        A10FFR=EEM*GPI/(4*amn**2)*F2S
        A1pFFR=EEM*GPI/(4*amn**2)*F2V
        A6mFFR=EEM*GPI/(2*amn)*(FAX-F1V)/Q2
	    A10=A10+fpvps*A10FFR 
	    A1p=A1p+fpvps*A1pFFR 
	    A6m=A6m+fpvps*A6mFFR 
      endif

C  the following part generates the old notation of Subroutine BORNEW
C  the units are in MeV^-n and all amplitudes are multiplied by mPi=MPI0
C  in CGLN sub. all amplitudes are divided by MPIP=139.5685 
      G(1,1) = A1p
	G(2,1) = A10
	G(3,1) = A1m
      G(1,2) = A2p
	G(2,2) = A20
	G(3,2) = A2m
      G(1,3) = A3p
	G(2,3) = A30
	G(3,3) = A3m
      G(1,4) = A4p
	G(2,4) = A40
	G(3,4) = A4m
      G(1,5) = A5p
	G(2,5) = A50
	G(3,5) = A5m
      G(1,6) = A6p
	G(2,6) = A60
	G(3,6) = A6m
	DO 200 i = 1, 6
	  IF (Reaction.EQ.1 .OR. Reaction.EQ.3) THEN
	    A(i,1) = -SQRT(1/3.D0)*(G(1,i)+2*G(3,i))-SQRT(3.D0)*G(2,i)
	    A(i,2) = SQRT(2/3.D0)*(G(1,i)-G(3,i))
	  ELSE IF (Reaction.EQ.2 .OR. Reaction.EQ.4) THEN
	    A(i,1) = SQRT(1/3.D0)*(G(1,i)+2*G(3,i))-SQRT(3.D0)*G(2,i)
	    A(i,2) = SQRT(2/3.D0)*(G(1,i)-G(3,i))
	  END IF
 200	CONTINUE
c
      Do k=1,2
         CALL InvAtoF(W,mi,mf,mPi,Q2,wGacm,qcm,kcm,CosTh,s,t,
     &        A(1,k),A(2,k),A(3,k),A(4,k),A(5,k),A(6,k),
     &        F1c(k),F2c(k),F3c(k),F4c(k),F5c(k),F6c(k) )
         F1(k)=DREAL(F1c(k))
         F2(k)=DREAL(F2c(k))
         F3(k)=DREAL(F3c(k))
         F4(k)=DREAL(F4c(k))
         F5(k)=DREAL(F5c(k))
         F6(k)=DREAL(F6c(k))
      End Do

      RETURN
      END
c
c

 	SUBROUTINE VECTNEW (Reaction, CosTh, F1, F2, F3, F4, F5, F6)
C******************************************************************************
C     vector mesons omega and rho  JFORM=0 : no hadronic ff
C******************************************************************************
	IMPLICIT NONE
	common/parvect/ xom1, xom2, xrho1, xrho2, xom, xrho
	DOUBLE PRECISION xom1, xom2, xrho1, xrho2, xom, xrho
      INTEGER IBORN,IVEC,IOME,IRHO,JFORM,IP33,IP11,ID13,IS11,IF15,IS11T
      DOUBLE PRECISION M1FM,M2FM,MPIFM
	INTEGER i, a, Reaction, MODE,GAUGE, IS11F, IS11S, ID33,IP11S
      COMMON /LINK1/ M1FM,M2FM,MPIFM,IBORN,IVEC,IOME,IRHO,JFORM,
     &              IP33,IP11,ID13,ID33,IS11F,IS11S,IF15,IS11T,IP11S
      DOUBLE PRECISION XLAMPI,XLAMAX, com,crho
      COMMON /FORMF/ xlampi,xlamax,mode,gauge
	DOUBLE PRECISION mi, mf, mPi, mn, corr
	COMMON /Mass/ mi, mf, mPi
    	DOUBLE PRECISION OmegL, Q2, W, qcm, kcm, kt, wGacm
	COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
    	DOUBLE PRECISION wPicm, CosTh, U, T, DU, DT,KQV, KQ
	DOUBLE PRECISION Ei, Ef, WPP, WPM, WMP, WMM, WM, WP
	DOUBLE PRECISION F(6,2),F1(2),F2(2),F3(2),F4(2)
      DOUBLE PRECISION F5(2),F6(2),G(3,6)
	DOUBLE PRECISION D,F1S,F1V,F2S,F2V,FPI,FAX,FPiNN,FOmNN,FRhoNN
	DOUBLE PRECISION Pi, eLadung, Kopplung, C1, C2, C3
	PARAMETER (Pi=3.1415926536D0,eLadung=.3028619D0,
     A     Kopplung=0.996365D0)
	DOUBLE PRECISION TMP1, TMP2, TMP3, VEK1, VEK2, Q2GEV
       DOUBLE PRECISION ggpo, ggpr, gom1,gom2, grho1, grho2
	DOUBLE PRECISION mOm,mRho,DTOm,DTRho,LamPi,LamOm, LamRho
	PARAMETER(mOm=783.0D0,mRho=770.D0,LamPi=0.6D3,
     * LamOm=1.2D3,LamRho=1.5D3)
	INTEGER JRHOFORM,JOMEFORM,D13mode
	common/newpar/JRHOFORM,JOMEFORM,D13mode

	mn = (mi+mf)/2.D0
c ****************   NEW **************
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	wPicm = SQRT(mPi*mPi+kcm*kcm)
	Ei = SQRT(mi*mi+qcm**2)
	Ef = SQRT(mf*mf+kcm**2)
	KQV = kcm*qcm*CosTh
	KQ = 2.*(wGacm*wPicm-KQV)
	kt = kcm*kcm+qcm*qcm-2.*KQV
	U = mi*mi+mPi*mPi-2.*(Ei*wPicm+KQV)
	DU = 1./(U-mf*mf)
	T = mPi*mPi-Q2-2.*(wGacm*wPicm-KQV)
	DTOm = 1./(T-mOm*mOm)
	DTRho = 1./(T-mRho*mRho)

c      VZ geaendert am 19.5.95
c ****************************************
      ggpo = -IVEC*IOME*0.314*eLadung
      FOmNN = LamOm*LamOm*xom**2/(LamOm*LamOm*xom**2+kt)
      IF (JOMEFORM .EQ. 0) FOmNN=1
      gom1 = 21.*FOmNN*xom1
      gom2= -12.*FOmNN*xom2
c **************   rho ******************
	ggpr = -IVEC*IRHO*0.103*eLadung
      grho1 = SQRT(4.*Pi*0.84)
      FRhoNN =LamRho*LamRho*xrho**2/(LamRho*LamRho*xrho**2+kt)
      IF (JRHOFORM .EQ. 0) FRhoNN=1
      grho1 =  2.*FRhoNN*xrho1 
      grho2 = 13.*FRhoNN*xrho2
c *************************************
      WP = W+mn
	WM = W-mn
	WPP = SQRT((Ef+mf)*(Ei+mi))/(2.*W)
	WPM = SQRT((Ef+mf)*(Ei-mi))/(2.*W)
	WMP = SQRT((Ef-mf)*(Ei+mi))/(2.*W)
	WMM = SQRT((Ef-mf)*(Ei-mi))/(2.*W)

      Q2GEV=Q2/1.0D6
	CALL FORM (Q2GEV, D, F1S, F1V, F2S, F2V, FPI, FAX)
c	CALL FORM (D, F1S, F1V, F2S, F2V, FPI, FAX)

	C2 = ggpo/4./Pi
	C3 = ggpr/4./Pi

	corr=1.0

       VEK1 = WM*WM-KQ/2.
	VEK2 = VEK1+  corr*W/mn*((Ef-mf)*(Ei-mi)-KQV)
        G(1,1) = C2*WPP*D*DTOm*(gom1*VEK1+gom2*VEK2)
        G(2,1) = C3*WPP*D*DTRho*(grho1*VEK1+grho2*VEK2)
        G(3,1) = 0

	VEK1 = WP*WP-KQ/2.
	VEK2 = VEK1+  corr*W/mn*(KQV-(Ef+mf)*(Ei+mi))
       G(1,2) = C2*WMM*D*DTOm*(gom1*VEK1+gom2*VEK2)
	 G(2,2) = C3*WMM*D*DTRho*(grho1*VEK1+grho2*VEK2)
       G(3,2) = 0

	VEK1 = -(Ef+mn)*WP
	VEK2 = VEK1+  corr*W/mn*(Ef+mf)*(Ei+mi)
       G(1,3) = C2*WMM*D*DTOm*(gom1*VEK1+gom2*VEK2)
	 G(2,3) = C3*WMM*D*DTRho*(grho1*VEK1+grho2*VEK2)
       G(3,3) = 0

	VEK1 = -WM*(Ef-mn)
	VEK2 = VEK1-  corr*W/mn*(Ef-mf)*(Ei-mi)
       G(1,4) = C2*WPP*D*DTOm*(gom1*VEK1+gom2*VEK2)
	 G(2,4) = C3*WPP*D*DTRho*(grho1*VEK1+grho2*VEK2)
       G(3,4) = 0

	VEK1 = wGacm*WM-KQ/2.-WP*KQV/(Ei+mn)
	VEK2 = VEK1
       G(1,5) = C2*WPP*D*DTOm*(gom1*VEK1+gom2*VEK2)
	 G(2,5) = C3*WPP*D*DTRho*(grho1*VEK1+grho2*VEK2)
       G(3,5) = 0

	VEK1 = wGacm*WP-KQ/2.-WM*KQV/(Ei-mn)
	VEK2 = VEK1
       G(1,6) = C2*WMM*D*DTOm*(gom1*VEK1+gom2*VEK2)
	 G(2,6) = C3*WMM*D*DTRho*(grho1*VEK1+grho2*VEK2)
       G(3,6) = 0

	DO 200 i = 1, 6

	  IF (Reaction.EQ.1 .OR. Reaction.EQ.3) THEN
	    F(i,1) = -SQRT(1./3.)*(G(1,i)+2.*G(3,i))-SQRT(3.)*G(2,i)
	    F(i,2) = SQRT(2./3.)*(G(1,i)-G(3,i))
	  ELSE IF (Reaction.EQ.2 .OR. Reaction.EQ.4) THEN
	    F(i,1) = SQRT(1./3.)*(G(1,i)+2.*G(3,i))-SQRT(3.)*G(2,i)
	    F(i,2) = SQRT(2./3.)*(G(1,i)-G(3,i))
	  END IF

 200	CONTINUE

	DO a = 1, 2
        F1(a) = F(1,a)
	  F2(a) = F(2,a)
	  F3(a) = F(3,a)
	  F4(a) = F(4,a)
        F5(a) = F(5,a)
	  F6(a) = F(6,a)
	END DO
        return
        END

	SUBROUTINE FORM (Q2GEV, F, F1S, F1V, F2S, F2V, FPI, FAX)
C******************************************************************************
C     electromagnetic form factors
C     all masses, Q2 and Lambda parameters in GeV
C     cut-off Lamdas are now explicitly given here
C******************************************************************************
	IMPLICIT NONE
      INTEGER MODE,GAUGE
      DOUBLE PRECISION XLAMPI,XLAMAX, FPI,FAX, XFPI,XFAX
	DOUBLE PRECISION Q2GEV,F,mg2,tau,GEP,GEN,GMN,GMP,F1P,F1N,F2P,F2N
	DOUBLE PRECISION F1S, F1V, F2S, F2V
	DOUBLE PRECISION mup, mun, mn
c	DOUBLE PRECISION mi, mf, mPi, mn
c	COMMON /Mass/ mi, mf, mPi
c	DOUBLE PRECISION OmegL, Q2, W, wGacm, kcm
c	COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
      COMMON /FORMF/ xlampi,xlamax,mode,gauge
c       COMMON /EPLCORR/ XFPI,XFAX
	PARAMETER (mn = 0.9382723, mg2 = .71, mup = 2.793, mun = -1.913)
c      mn proton mass
c      mn = (mi+mf)/2.
	tau = Q2GEV/(4.*mn*mn)
	F = 1./(1.+Q2GEV/mg2)**2
	GEP = F
	GMN = mun*F
	GEN = (-tau)/(1.+4.*tau)*GMN
	GMP = mup*F
	F1P = (GEP+tau*GMP)/(1.+tau)
	F1N = (GEN+tau*GMN)/(1.+tau)
	F2P = (GMP-GEP)/((mup-1.)*(1.+tau))
	F2N = (GMN-GEN)/(mun*(1.+tau))
	F1S = F1P+F1N
	F1V = F1P-F1N
	F2S = (mup-1.)*F2P+mun*F2N
	F2V = (mup-1.)*F2P-mun*F2N
      FPI = F1V
      FAX = F1V
c   MODE  :  0  Sachs form factors in Dipole form and Fpi=Fax=F1V  (std)
c         :  1  all Dirac form factors and Fpi, Fax = 1
c         :  2  as MODE=0, but Fpi and Fax different from Dipole
c         :  3  all Dirac form factors and Fpi, Fax in Dipole form (diff from 0)
c         :  4  parametrization by J.J. Kelly, PRC70, 068202 (2004) and Fpi, Fax diff. from Dipole
c         :  5  Kelly plus new Fpi, GA
c         :  6  Kelly plus ad hoc Fpi, GA
c         :  7  Kelly plus Fpi, GA fitted (Maid2007std)
c         :  9  Test  Fpi=Q2GeV all others = 0

        IF (MODE.EQ.0) RETURN
        IF (MODE.EQ.1) GOTO 100
        IF (MODE.EQ.2) GOTO 200
        IF (MODE.EQ.3) GOTO 300
        IF (MODE.EQ.4) GOTO 400
        IF (MODE.EQ.5) GOTO 500
        IF (MODE.EQ.6) GOTO 600
        IF (MODE.EQ.7) GOTO 700
        IF (MODE.EQ.9) GOTO 900
        F=0
  100   F=1
  300   F1S=F
        F1V=F
        F2S=(mup-1+mun)*F
        F2V=(mup-1-mun)*F
        FPI = F
        FAX = F
        RETURN
  200   XLAMPI=0.643
        XLAMAX=1.010
        GOTO 899
  400   XLAMPI=0.680
        XLAMAX=1.030        
        GOTO 899
  500   XLAMPI=0.730
        XLAMAX=1.026
C   new (July 2006), pion ff Amendolia, axial ff (Liesenfeld, Sirca)
        GOTO 899
  600   CONTINUE
c   xlampi and xlamax are defined in the calling program
        GOTO 899
c   xlampi and xlamax are fitted and modified by solution file (Maid2007)
  700   CONTINUE
c        XLAMAX=1.030*FAX         
c        XLAMPI=0.680*XFPI 
  899   CONTINUE
        FPI = 1.0/(1.0+Q2GEV/XLAMPI**2)
        FAX = 1.0/(1.0+Q2GEV/XLAMAX**2)**2
	  GEP =     (1-0.24*tau)/(1+10.98*tau+12.82*tau**2+21.97*tau**3)
	  GMP = mup*(1+0.12*tau)/(1+10.97*tau+18.86*tau**2+ 6.55*tau**3)
	  GMN = mun*(1+2.33*tau)/(1+14.72*tau+24.20*tau**2+84.10*tau**3)
	  GEN = (1.70*tau)/(1.+3.30*tau)*F
	  F1P = (GEP+tau*GMP)/(1.+tau)
	  F1N = (GEN+tau*GMN)/(1.+tau)
	  F2P = (GMP-GEP)/((mup-1.)*(1.+tau))
	  F2N = (GMN-GEN)/(mun*(1.+tau))
	  F1S = F1P+F1N
	  F1V = F1P-F1N
	  F2S = (mup-1.)*F2P+mun*F2N
	  F2V = (mup-1.)*F2P-mun*F2N
        RETURN
  900 FPI = Q2GEV
      FAX = 0
      F=0
      F1S=0
	F1V=0
	F2S=0
	F2V=0
        RETURN
	END



	SUBROUTINE FORM_old (Q2GEV, F, F1S, F1V, F2S, F2V, FPI, FAX)
C******************************************************************************
C     electromagnetic form factors
C     all masses, Q2 and Lambda parameters in GeV
C     cut-off Lamdas are now explicitly given here
C******************************************************************************
	IMPLICIT NONE
      INTEGER MODE,GAUGE
      DOUBLE PRECISION XLAMPI,XLAMAX, FPI,FAX
	DOUBLE PRECISION Q2GEV,F,mg2,tau,GEP,GEN,GMN,GMP,F1P,F1N,F2P,F2N
	DOUBLE PRECISION F1S, F1V, F2S, F2V
	DOUBLE PRECISION mup, mun, mn
c	DOUBLE PRECISION mi, mf, mPi, mn
c	COMMON /Mass/ mi, mf, mPi
c	DOUBLE PRECISION OmegL, Q2, W, wGacm, kcm
c	COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
      COMMON /FORMF/ xlampi,xlamax,mode,gauge
	PARAMETER (mn = 0.9382723, mg2 = .71, mup = 2.793, mun = -1.913)
c      mn proton mass
c      mn = (mi+mf)/2.
	tau = Q2GEV/(4.*mn*mn)
	F = 1./(1.+Q2GEV/mg2)**2
	GEP = F
	GMN = mun*F
	GEN = (-tau)/(1.+4.*tau)*GMN
	GMP = mup*F
	F1P = (GEP+tau*GMP)/(1.+tau)
	F1N = (GEN+tau*GMN)/(1.+tau)
	F2P = (GMP-GEP)/((mup-1.)*(1.+tau))
	F2N = (GMN-GEN)/(mun*(1.+tau))
	F1S = F1P+F1N
	F1V = F1P-F1N
	F2S = (mup-1.)*F2P+mun*F2N
	F2V = (mup-1.)*F2P-mun*F2N
      FPI = F1V
      FAX = F1V
c   MODE  :  0  use Sachs form factors in Dipole form (std)
c         :  1  all Dirac form factors = 1
c         :  3  all Dirac form factors in Dipole form (diff from 0)
c         :  4  parametrization by J.J. Kelly, PRC70, 068202 (2004)
c         :  2  Fpi and Fax different F1 and F2 as in MODE=0
c         :  9  Test  Fpi=Q2GeV all others = 0

        IF (MODE.EQ.0) RETURN
        IF (MODE.EQ.1) GOTO 100
        IF (MODE.EQ.2) GOTO 200
        IF (MODE.EQ.3) GOTO 300
        IF (MODE.EQ.4) GOTO 400
        IF (MODE.EQ.9) GOTO 900
        F=0
  100   F=1
  300   F1S=F
        F1V=F
        F2S=(mup-1+mun)*F
        F2V=(mup-1-mun)*F
        FPI = F
        FAX = F
        RETURN
  200   XLAMPI=0.643
        XLAMAX=1.010
        FPI = 1.0/(1.0+Q2GEV/XLAMPI**2)
        FAX = 1.0/(1.0+Q2GEV/XLAMAX**2)**2
        RETURN
  400   XLAMPI=0.680
        XLAMAX=1.030
        FPI = 1.0/(1.0+Q2GEV/XLAMPI**2)
        FAX = 1.0/(1.0+Q2GEV/XLAMAX**2)**2
	  GEP =     (1-0.24*tau)/(1+10.98*tau+12.82*tau**2+21.97*tau**3)
	  GMP = mup*(1+0.12*tau)/(1+10.97*tau+18.86*tau**2+ 6.55*tau**3)
	  GMN = mun*(1+2.33*tau)/(1+14.72*tau+24.20*tau**2+84.10*tau**3)
	  GEN = (1.70*tau)/(1.+3.30*tau)*F
	  F1P = (GEP+tau*GMP)/(1.+tau)
	  F1N = (GEN+tau*GMN)/(1.+tau)
	  F2P = (GMP-GEP)/((mup-1.)*(1.+tau))
	  F2N = (GMN-GEN)/(mun*(1.+tau))
	  F1S = F1P+F1N
	  F1V = F1P-F1N
	  F2S = (mup-1.)*F2P+mun*F2N
	  F2V = (mup-1.)*F2P-mun*F2N
        RETURN
  900 FPI = Q2GEV
      FAX = 0
      F=0
      F1S=0
	F1V=0
	F2S=0
	F2V=0
        RETURN
	END



       SUBROUTINE P33NEW(ReM1P,ImM1P)
C******************************************************************************
       IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fke,fpi,DE,DM,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION W0,Ga0,X,fkm,kglb,FIM,FIE,FIL
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM1P(3),ImM1P(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM1(3),ImM1(3),WX
       DOUBLE PRECISION FMQ,FEQ,FLQ,Fq,GMED,GNR,WRMED,XCM
       COMMON /MEDIUM/ GMED,WRMED,XCM
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37
       DOUBLE PRECISION RM,RE,RL, DS,A1,A3,S1
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
       COMMON /PARP33/ RM,RE,RL, DM, DE, DS
       INTEGER i
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1232.,Ga0=130.)
c ***************10.06.02 *****************************
c	X=570.
	X=540.
c *********************************************
	mn = (mi+mf)/2.
        mPi=139.5685
	  kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        qcm0 = SQRT(((W0*W0+mi*mi)/(2*W0))**2-mi*mi)
        kgcm = (W*W-mi*mi)/2./W
        kglb=(W*W-mi*mi)/2./mi
        qcm=sqrt(wGacm**2+Q2)
	Q2G = Q2/1.D6
C **********************************************
	fkm = (kgcm/qcm0)**2*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
	fke = (qcm0/kgcm)*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
	fpi =(kcm/kcm0)**3*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))
c *************************************************
         CALL HP33(Q2G,qcm,kgcm,DE,DM,DS,A1,A3,S1,0)
c *************************************************

	GaGa = kgcm*qcm/Pi*mn/W0/4./1.D9
        GaPi=Ga0*fpi*W0/W
	GaT =GaPi

 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)

	  ReM1(1) = Fak*(W0*W0-W*W)*DM*mPi*fkm 
	  ImM1(1) = Fak*W0*GaT*DM*mPi*fkm 
	  ReM1(2) = Fak*(W0*W0-W*W)*DE*mPi*fke 
	  ImM1(2) = FaK*W0*GaT*DE*mPi*fke 
	  ReM1(3) = wGacm/qcm*Fak*(W0*W0-W*W)*DS*mPi*fke
	  ImM1(3) = wGacm/qcm*FaK*W0*GaT*DS*mPi*fke
C******************************************************************************
C                (1), (2), (3)  for  M1+, E1+, L1+
c******************************************************************************
c resonance phase is the same for all 3 multipoles
c it must be calculated for an amplitude .ne. 0.0
c******************************************************************************
          if (ReM1(1).ne.0.0) then
           del_R = datan2(ImM1(1),ReM1(1))
          elseif (ReM1(2).ne.0.0) then
           del_R = datan2(ImM1(2),ReM1(2))
          elseif (ReM1(3).ne.0.0) then
           del_R = datan2(ImM1(3),ReM1(3))
          else
	    del_R = pi/2.
          endif
          IF (del_R.LT.0.0) del_R = del_R + pi
          phi_c = phi_p33 - del_R
          Phi(1) = phi_c
          Phi(2) = phi_c
          Phi(3) = phi_c
c	 write (6,145) w,phi_p33*180./pi,del_R*180./pi,phi_c*180./pi
c145      format(1x,f8.2,2x,2x,3f12.6)
c*************************************************************
          DO 111 I=1,3
       ReM1P(i) = ReM1(i)*DCOS(Phi(i))-ImM1(i)*DSIN(Phi(i))
       ImM1P(i) = ReM1(i)*DSIN(Phi(i))+ImM1(i)*DCOS(Phi(i))
 111      CONTINUE
       RETURN
	END


	SUBROUTINE P11NEW(ISO,ReM1M,ImM1M)
C******************************************************************************
C     multipoles  M1- and L1-  for P11(1440)  with fit parameters RM and RL
C******************************************************************************
	IMPLICIT NONE
       DOUBLE PRECISION fk,fpi,DE,DM,DS,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION W0, Ga0, X, S1, FQP, FQN, FQ, WX
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM1M(3),ImM1M(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM1(3),ImM1(3)
       DOUBLE PRECISION alp,q2EVF0,q2EVF,gam0,gam,EXP0
	DOUBLE PRECISION k2pi0,k2pi,w2thr,finel, PhiR
       INTEGER i, ISO
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37, Phi_R, CS0, A1
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71, alp=370.)
       PARAMETER(W0=1440.,Ga0=350., X=470.)

	mn = (mi+mf)/2.
c **************************************
        MPI=139.5685
        qcm=sqrt(wGacm**2+Q2)
c *******************************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
C ******************************************************
	fk = (kgcm/qcm0)**0*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**3*(kcm0*kcm0+X*X)/(kcm*kcm+X*X)
	finel=(k2pi/k2pi0)**6*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**3
c ******************************************************
	GaGa = kgcm*qcm/2./Pi*mn/W0/1.D9
        GaPi=0.70*Ga0*fpi*W0/W
	GaT =GaPi+0.30*Ga0*finel
        Fak=W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c ******************************************************
       	 Q2G = Q2/1.D6
         CALL HP11(ISO,PhiR,Q2G,DM,DS,A1,S1,0)
c ******************************************************
	  IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DM=DM*fk
        DS=DS*fk

	  ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DM=DM*fk *(kgcm/qcm0)**(-1) 
        DS=DS*fk
	  END IF

c*************************************************************
	  ReM1M(1) = 0
	  ImM1M(1) = 0
	  ReM1(2) = Fak* (W0*W0-W*W)*DM*mPi
	  ImM1(2) = Fak* W0*GaT*DM*mPi
	  ReM1(3) = Fak* (W0*W0-W*W)*DS*mPi*wGacm/qcm
	  ImM1(3) = Fak* W0*GaT*mPi*DS*wGacm/qcm

c******************************************************************************
c resonance phase is the same for all 2 multipoles
c it must be calculated for an amplitude .ne. 0.0
c******************************************************************************
          if (ReM1(2).ne.0.0) then
           del_R = datan2(ImM1(2),ReM1(2))
          elseif (ReM1(3).ne.0.0) then
           del_R = datan2(ImM1(3),ReM1(3))
          else
	    del_R = pi/2.
          endif
          IF (del_R.LT.0.0) del_R = del_R + pi
          phi_c = phi_p11 - del_R
          Phi(1) = phi_c
          Phi(2) = phi_c
          Phi(3) = phi_c
c	 write (6,145) w,iso,phi_p11*180./pi,del_R*180./pi,phi_c*180./pi
c145      format(1x,f8.2,2xI3,2x,3f12.6)
c*************************************************************
          ReM1M(1)= 0.
          ImM1M(1)= 0.
       ReM1M(2)=ReM1(2)*DCOS(Phi(2))-ImM1(2)*DSIN(Phi(2))
       ImM1M(2)=ReM1(2)*DSIN(Phi(2))+ImM1(2)*DCOS(Phi(2))
       ReM1M(3)=ReM1(3)*DCOS(Phi(3))-ImM1(3)*DSIN(Phi(3))
       ImM1M(3)=ReM1(3)*DSIN(Phi(3))+ImM1(3)*DCOS(Phi(3))
       return
       END

	SUBROUTINE P11sec(ISO,ReM1M,ImM1M)
C******************************************************************************
C     multipoles  M1- and L1-  for P11(1440)  with fit parameters RM and RL
C******************************************************************************
	IMPLICIT NONE
       DOUBLE PRECISION fk,fpi,DE,DM,DS,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION W0, Ga0, X, S1, FQP, FQN, FQ, WX
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM1M(3),ImM1M(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM1(3),ImM1(3)
       DOUBLE PRECISION alp,q2EVF0,q2EVF,gam0,gam,EXP0
	DOUBLE PRECISION k2pi0,k2pi,w2thr,finel, PhiR
       INTEGER i, ISO
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37, Phi_R, CS0, A1
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71, alp=370.)
       PARAMETER(W0=1700.,Ga0=30., X=500.)

	mn = (mi+mf)/2.
c **************************************
        MPI=139.5685
        qcm=sqrt(wGacm**2+Q2)
c *******************************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
C ******************************************************
	fk = (kgcm/qcm0)**0*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**3*(kcm0*kcm0+X*X)/(kcm*kcm+X*X)
	finel=(k2pi/k2pi0)**6*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**3
c ******************************************************
	GaGa = kgcm*qcm/2./Pi*mn/W0/1.D9
        GaPi=0.1*Ga0*fpi*W0/W
	GaT =GaPi+0.9*Ga0*finel
        Fak=W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c *****************************************************
       	 Q2G = Q2/1.D6
         CALL HP11sec(ISO,PhiR,Q2G,DM,DS,A1,S1,0)
c ******************************************************
	  IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DM=DM*fk
        DS=DS*fk

	  ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DM=DM*fk *(kgcm/qcm0)**(-1) 
        DS=DS*fk
	  END IF

c*************************************************************
	  ReM1M(1) = 0
	  ImM1M(1) = 0
	  ReM1(2) = Fak* (W0*W0-W*W)*DM*mPi
	  ImM1(2) = Fak* W0*GaT*DM*mPi
	  ReM1(3) = Fak* (W0*W0-W*W)*DS*mPi*wGacm/qcm
	  ImM1(3) = Fak* W0*GaT*mPi*DS*wGacm/qcm

c******************************************************************************
c resonance phase is the same for all 2 multipoles
c it must be calculated for an amplitude .ne. 0.0
c******************************************************************************
          if (ReM1(2).ne.0.0) then
           del_R = datan2(ImM1(2),ReM1(2))
          elseif (ReM1(3).ne.0.0) then
           del_R = datan2(ImM1(3),ReM1(3))
          else
	    del_R = pi/2.
          endif
          IF (del_R.LT.0.0) del_R = del_R + pi
          phi_c = phi_p11 - del_R
          Phi(1) = phi_c
          Phi(2) = phi_c
          Phi(3) = phi_c
c	 write (6,145) w,iso,phi_c*180./pi
c145      format(1x,f8.2,2xI3,2x,e10.4)
c*************************************************************
          ReM1M(1)= 0.
          ImM1M(1)= 0.
       ReM1M(2)=ReM1(2)*DCOS(Phi(2))-ImM1(2)*DSIN(Phi(2))
       ImM1M(2)=ReM1(2)*DSIN(Phi(2))+ImM1(2)*DCOS(Phi(2))
       ReM1M(3)=ReM1(3)*DCOS(Phi(3))-ImM1(3)*DSIN(Phi(3))
       ImM1M(3)=ReM1(3)*DSIN(Phi(3))+ImM1(3)*DCOS(Phi(3))
       return
       END
	

	SUBROUTINE D13NEW(ISO,ReM2M,ImM2M)
C*****************************************************************
C      multipoles M2- , E2- and L2- for D13(1520)
C*****************************************************************
	IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT,Ga0,X,W0
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM2M(3),ImM2M(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM2(3),ImM2(3)
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM, A1,A3,S1
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel,DS,PhiR
       INTEGER i, ISO
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
c       PARAMETER(W0=1520.,Ga0=130., X=500.)
       PARAMETER(W0=1530.,Ga0=130., X=500.)

	mn = (mi+mf)/2.
c **************************************
        MPI=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
c *********************************************************
	fk = (kgcm/qcm0)**2*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**5*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**2
	finel=(k2pi/k2pi0)**8*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**4
        GaPi=0.60*Ga0*fpi*W0/W
	GaT =GaPi+0.40*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/4./1.D9
 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c *********************************************************
	Q2G = Q2/1.D6
	CALL HD13(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR
c ******************************************************
	  IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DE=DE*fk*(kgcm/qcm0)
        DM=DM*fk*(kgcm/qcm0)**2.2
        DS=DS*fk*(kgcm/qcm0)

	  ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE=DE*fk 
        DM=DM*fk*(kgcm/qcm0)**5  
        DS=DS*fk

	  END IF
c	 write (6,145) w,iso,PhiR*180./pi
c145      format(1x,f8.2,2xI3,2x,3f12.6)
c **************************************************

        ReM2(1) =Fak* (W0*W0-W*W)*DE*mPi
	ImM2(1) =Fak* W0*GaT*DE*mPi
	ReM2(2) =Fak*(W0*W0-W*W)*DM*mPi
	ImM2(2) =Fak* W0*GaT*DM*mPi
	ReM2(3) =wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
	ImM2(3) =wGacm/qcm*Fak* W0*GaT*DS*mPi
c *****************************************************
          DO 111 i=1,3
       ReM2M(i) = ReM2(i)*DCOS(Phi(i))-ImM2(i)*DSIN(Phi(i))
       ImM2M(i) = ReM2(i)*DSIN(Phi(i))+ImM2(i)*DCOS(Phi(i))
 111      CONTINUE

       RETURN
	END

C
       SUBROUTINE D33NEW(ReD33,ImD33)
C*************************************************
C      multipoles M2- , E2- and L2- for D33(1700)
C*************************************************
       IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION W0,Ga0,X,k2pi0,k2pi,w2thr,finel, PhiR
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReD33(3),ImD33(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReD3(3),ImD3(3)
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,DS,A1,A3,S1
       INTEGER i
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1740.,Ga0=450., X=700.)

	mn = (mi+mf)/2.
        mPi=139.5685
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        qcm0 = SQRT(((W0*W0+mi*mi)/(2*W0))**2-mi*mi)
        kgcm = (W*W-mi*mi)/2./W
        qcm=sqrt(wGacm**2+Q2)
	Q2G = Q2/1.D6

	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
C **********************************************
	fk = (kgcm/qcm0)**2*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
	fpi = (kcm/kcm0)**5*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**2
	finel=(k2pi/k2pi0)**8*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**4
        GaPi=0.15*Ga0*fpi*W0/W
	GaT =GaPi+0.85*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/4./1.D9
 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c **************************************************
	Q2G = Q2/1.D6
	CALL HD33(Q2G,PhiR,DE,DM,DS,A1,A3,S1,0)
	Phi(1) = PhiR
	Phi(2) = PhiR
	Phi(3) = PhiR
c *********************************************************
        DE=DE*fk*(kgcm/qcm0)**2
        DM=DM*fk*(kgcm/qcm0)**3.6
        DS=DS*fk*(kgcm/qcm0)**2
c *************************************************

	  ReD3(1) = Fak* (W0*W0-W*W)*DE*mPi
	  ImD3(1) = Fak* W0*GaT*DE*mPi
	  ReD3(2) = Fak* (W0*W0-W*W)*DM*mPi
	  ImD3(2) = Fak* W0*GaT*DM*mPi
          ReD3(3)= wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
          ImD3(3)= wGacm/qcm*Fak* W0*GaT*DS*mPi

          DO 111 I=1,3
       ReD33(i) = ReD3(i)*DCOS(Phi(i))-ImD3(i)*DSIN(Phi(i))
       ImD33(i) = ReD3(i)*DSIN(Phi(i))+ImD3(i)*DCOS(Phi(i))
  111     CONTINUE

       RETURN
	END
c
       SUBROUTINE S11fst(ISO,ReM0P,ImM0P)
C********************************************
C      multipoles E0+ and L0+  for S11(1535)
C********************************************
       IMPLICIT NONE
c ********************************************************************
       DOUBLE PRECISION fk,fpi,DE,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION ALPHE, ALPHM, W0, Ga0, X
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM0P(3),ImM0P(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,Fak,kx,ReM0(3),ImM0(3)
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel, PhiR
       DOUBLE PRECISION qj, qj0, GaEta, mEta, Dis, FE
       DOUBLE PRECISION TE,TM,SE,SM,DS,A1,S1
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
       INTEGER i, ISO
       COMMON /Mass/ mi, mf, mPION
       COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
c       PARAMETER(W0=1520.,Ga0=100., X=500.,mEta=549.)
       PARAMETER(W0=1535.,Ga0=100., X=500.,mEta=549.)       
c
        mn = (mi+mf)/2.
        mPi=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
        kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        kgcm = (W*W-mi*mi)/2./W
        qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
        k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *  k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
        finel=(k2pi/k2pi0)**4*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**2
c ********************************************************
        fk = (kgcm/qcm0)**2*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)

        GaPi = 0.40*Ga0*kcm/kcm0*W0/W

        Dis = ((W*W-mEta*mEta+mf*mf)/2./W)**2-mf*mf
        IF (Dis.LE.0.) THEN
          qj = 0.
        ELSE IF (Dis.GT.0.) THEN
          qj = SQRT(Dis)
        END IF
        qj0 = SQRT(((W0*W0-mEta*mEta+mf*mf)/2./W0)**2-mf*mf)
        GaEta = 0.50*Ga0*qj/qj0*W0/W
        GaT = GaPi+GaEta+0.10*Ga0*finel

        GaGa = kgcm*qcm/Pi*mn/W0/2.D9
        Fak=W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     *  *SQRT(GaGa*GaPi/qcm/kcm)
c *************************************************************
	Q2G = Q2/1.D6
        CALL HS11f(ISO,PhiR,Q2G,DE,DS,A1,S1,0)
c **********************************************
         FE=1.
        IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DE=DE*fk
	DS=DS*fk
c		IF (W.LT.1480.) FE=dexp(25.*(W-1480.)/1480.)
        ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE=DE*fk 
	DS=DS*fk
	END IF
c *************************************************************
          ReM0(1) = Fak*(W0*W0-W*W)*DE*mPi*FE
          ImM0(1) = Fak*W0*GaT*DE*mPi*FE
          ReM0(3) = (wGacm/qcm)*Fak*(W0*W0-W*W)*DS*mPi
          ImM0(3) = (wGacm/qcm)*Fak*W0*GaT*DS*mPi
c *************************************************************
          del_R=pi/2.
          IF (DABS(ReM0(1)).GT.1.E-6)
     &     del_R = datan2(ImM0(1),ReM0(1))
          IF (del_R.LT.0.0) del_R = del_R + pi
          Phi(1) = phi_s11 - del_R
          phi_c = phi_s11 - del_R
         IF (W.GT.1300.) Phi(1)=0.143
c	   Phi(1)=0
c	 write (6,145) w,iso,phi_s11*180./pi,del_R*180./pi,phi_c*180./pi,
c     &               Phi(1)*180./pi
c145      format(1x,f8.2,2xI3,2x,4f12.6)
c*************************************************************
          ReM0P(1)=DCOS(Phi(1))*ReM0(1)-DSIN(Phi(1))*ImM0(1)
          ImM0P(1)=DSIN(Phi(1))*ReM0(1)+DCOS(Phi(1))*ImM0(1)
          ReM0P(2) = 0.
          ImM0P(2) = 0.
          ReM0P(3)=DCOS(Phi(1))*ReM0(3)-DSIN(Phi(1))*ImM0(3)
          ImM0P(3)=DSIN(Phi(1))*ReM0(3)+DCOS(Phi(1))*ImM0(3)

        RETURN
        END
c
       SUBROUTINE S11sec(ISO,ReM0P,ImM0P)
C********************************************
C      multipoles E0+ and L0+ for S11(1650)
C********************************************
       IMPLICIT NONE
c *****************************************************
       DOUBLE PRECISION fk,fpi,DE,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION ALPHE, ALPHM, W0, Ga0, X
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM0P(3),ImM0P(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM0(3),ImM0(3)
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel
       DOUBLE PRECISION TE,TM,SE,SM,DS,A1,S1
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37, PhiR
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
       INTEGER i, ISO
       COMMON /Mass/ mi, mf, mPION
       COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
c       PARAMETER(W0=1670.,Ga0=110., X=500.)
       PARAMETER(W0=1690.,Ga0=100., X=500.)
c
        mn = (mi+mf)/2.
        mPi=139.5685
        qcm=sqrt(wGacm**2+Q2)
        kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
        k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *  k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
        kgcm = (W*W-mi*mi)/2./W
c ************************************************************
        fk = (kgcm/qcm0)**4*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
        finel=(k2pi/k2pi0)**4*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**2
c
        GaPi = 0.85*Ga0*kcm/kcm0*W0/W
        GaT =GaPi+0.15*Ga0*finel
        GaGa = kgcm*qcm/Pi*mn/W0/2.D9
        Fak=W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     *  *SQRT(GaGa*GaPi/qcm/kcm)
c *************************************************************
	Q2G = Q2/1.D6
	CALL HS11s(ISO,PhiR,Q2G,DE,DS,A1,S1,0)
c ****************************************************
          IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DE=DE*fk
        DS=DS*fk

          ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE=DE*fk*(kgcm/qcm0)**(0)
	DS=DS*fk
          END IF

c
c *************************************************************
          ReM0(1) = Fak*(W0*W0-W*W)*DE*mPi
          ImM0(1) = Fak*W0*GaT*DE*mPi
          ReM0(3) = (wGacm/qcm)*Fak*(W0*W0-W*W)*DS*mPi
          ImM0(3) = (wGacm/qcm)*Fak*W0*GaT*DS*mPi

c*************************************************************
          del_R=pi/2.
          IF (DABS(ReM0(1)).GT.1.E-6)
     &     del_R = datan2(ImM0(1),ReM0(1))
          IF (del_R.LT.0.0) del_R = del_R + pi
          Phi(1) = phi_s11 - del_R
          phi_c = phi_s11 - del_R
         IF (W.GT.1300.) Phi(1)=0.1215
c	   Phi(1)=0
c	 write (6,145) w,iso,phi_s11*180./pi,del_R*180./pi,phi_c*180./pi,
c     &               Phi(1)*180./pi
c145      format(1x,f8.2,2xI3,2x,4f12.6)
c*************************************************************
          ReM0P(1)=DCOS(Phi(1))*ReM0(1)-DSIN(Phi(1))*ImM0(1)
          ImM0P(1)=DSIN(Phi(1))*ReM0(1)+DCOS(Phi(1))*ImM0(1)
          ReM0P(2) = 0.
          ImM0P(2) = 0.
          ReM0P(3)=DCOS(Phi(1))*ReM0(3)-DSIN(Phi(1))*ImM0(3)
          ImM0P(3)=DSIN(Phi(1))*ReM0(3)+DCOS(Phi(1))*ImM0(3)
c ***************************************************************
        RETURN
        END
c
       SUBROUTINE S11trd(ISO,ReM0P,ImM0P)
C********************************************
C      multipoles E0+ and L0+ for S11(1650)
C********************************************
       IMPLICIT NONE
c *****************************************************
       DOUBLE PRECISION fk,fpi,DE,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION ALPHE, ALPHM, W0, Ga0, X
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM0P(3),ImM0P(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM0(3),ImM0(3)
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel
       DOUBLE PRECISION TE,TM,SE,SM,DS,A1,S1
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37, PhiR
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
       INTEGER i, ISO
       COMMON /Mass/ mi, mf, mPION
       COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1950.,Ga0=200., X=500.)
c
        mn = (mi+mf)/2.
        mPi=139.5685
        qcm=sqrt(wGacm**2+Q2)
        kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
        k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *  k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
        kgcm = (W*W-mi*mi)/2./W
c ************************************************************
        fk = (kgcm/qcm0)**5*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
        finel=(k2pi/k2pi0)**4*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**2
c
        GaPi = 0.4*Ga0*kcm/kcm0*W0/W
        GaT =GaPi+0.60*Ga0*finel
        GaGa = kgcm*qcm/Pi*mn/W0/2.D9
        Fak=W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     *  *SQRT(GaGa*GaPi/qcm/kcm)
c *************************************************************
C ****************************************************
	Q2G = Q2/1.D6
	CALL HS11trd(ISO,PhiR,Q2G,DE,DS,A1,S1,0)
c ****************************************************
          IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DE=DE*fk
        DS=DS*fk

          ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE=DE*fk
	DS=DS*fk
          END IF

c
c *************************************************************
          ReM0(1) = Fak*(W0*W0-W*W)*DE*mPi
          ImM0(1) = Fak*W0*GaT*DE*mPi
          ReM0(3) = (wGacm/qcm)*Fak*(W0*W0-W*W)*DS*mPi
          ImM0(3) = (wGacm/qcm)*Fak*W0*GaT*DS*mPi

c*************************************************************
          del_R=pi/2.
          IF (DABS(ReM0(1)).GT.1.E-6)
     &     del_R = datan2(ImM0(1),ReM0(1))
          IF (del_R.LT.0.0) del_R = del_R + pi
          Phi(1) = phi_s11 - del_R
          IF (W.GT.1300.) Phi(1)=0.1215
c	 write (6,145) iso,phi(1),Del_R
c145      format(1x,I3,2x,2(e10.4,2x))
c*************************************************************
          ReM0P(1)=DCOS(Phi(1))*ReM0(1)-DSIN(Phi(1))*ImM0(1)
          ImM0P(1)=DSIN(Phi(1))*ReM0(1)+DCOS(Phi(1))*ImM0(1)
          ReM0P(2) = 0.
          ImM0P(2) = 0.
          ReM0P(3)=DCOS(Phi(1))*ReM0(3)-DSIN(Phi(1))*ImM0(3)
          ImM0P(3)=DSIN(Phi(1))*ReM0(3)+DCOS(Phi(1))*ImM0(3)
c ***************************************************************
        RETURN
        END

C
       SUBROUTINE S31NEW(ReS31,ImS31)
       IMPLICIT NONE
c ****************************************************************
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT, XE
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION W0,Ga0,X,k2pi0,k2pi,w2thr,finel,DS, PhiR
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReS31(3),ImS31(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReS3(3),ImS3(3)
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,TE0,TM0,SE0,SM0
       INTEGER i
       DOUBLE PRECISION del_R, phi_c, phi_s11,phi_s31,phi_p11,phi_p31
       DOUBLE PRECISION phi_p13, phi_p33, phi_d13, phi_d33
       DOUBLE PRECISION phi_d15, phi_d35, phi_f15, phi_f35
       DOUBLE PRECISION phi_f17, phi_f37, Phi_R, CS0,A1,S1
       COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1620.,Ga0=150., X=470.)

	mn = (mi+mf)/2.
        mPi=139.5685
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        qcm0 = SQRT(((W0*W0+mi*mi)/(2*W0))**2-mi*mi)
        kgcm = (W*W-mi*mi)/2./W
        qcm=sqrt(wGacm**2+Q2)
	Q2G = Q2/1.D6

	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
C **********************************************
	fk = (kgcm/qcm0)**5*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
	fpi = kcm/kcm0
	finel=(k2pi/k2pi0)**4*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**2
        GaPi=0.25*Ga0*fpi*W0/W
	GaT =GaPi+0.75*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/2./1.D9
 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c **************************************************
	Q2G = Q2/1.D6
	CALL HS31(Q2G,PhiR,DE,DS,A1,S1,0)
c **************************************************
        DE =DE*fk
        DS =DS*fk
c *************************************************

	  ReS3(1) = Fak* (W0*W0-W*W)*DE*mPi
	  ImS3(1) = Fak* W0*GaT*DE*mPi
          ReS3(2)=0.D0
          ImS3(2)=0.D0
          ReS3(3)=(wGacm/qcm)*Fak*(W0*W0-W*W)*DS*mPi
          ImS3(3)=(wGacm/qcm)*Fak*W0*GaT*DS*mPi

c*************************************************************
          del_R = atan(ImS3(1)/ReS3(1))
          IF (del_R.LT.0.0) del_R = del_R + pi
          phi_c = phi_S31 - del_R
          Phi(1) = phi_c
c	   Phi(1)=0
c	 write (6,145) w,phi_s31*180./pi,del_R*180./pi,phi_c*180./pi,
c     &               Phi(1)*180./pi
c145      format(1x,f8.2,2x,2x,4f12.6)
c*************************************************************
          ReS31(2)=0.D0
          ImS31(2)=0.D0
       ReS31(1)=ReS3(1)*DCOS(Phi(1))-ImS3(1)*DSIN(Phi(1))
       ImS31(1)=ReS3(1)*DSIN(Phi(1))+ImS3(1)*DCOS(Phi(1))
       ReS31(3)=ReS3(3)*DCOS(Phi(1))-ImS3(3)*DSIN(Phi(1))
       ImS31(3)=ReS3(3)*DSIN(Phi(1))+ImS3(3)*DCOS(Phi(1))
c *****************************************************

       RETURN
	END

	SUBROUTINE F15NEW(ISO,ReM3M,ImM3M)
C**************************************************
C      multipoles E3- , M3-  and L3- for F15(1680)
C**************************************************
	IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT,W0, Ga0, X
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM3M(3),ImM3M(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM3(3),ImM3(3)
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel,PhiR
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,DS,A1,A3,S1
       INTEGER i, ISO
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1680.,Ga0=135., X=500.)

	mn = (mi+mf)/2.
c **************************************
        mPi=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
c ***************************************************************
	fk = (kgcm/qcm0)**3*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**7*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**3
	finel=(k2pi/k2pi0)**10*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**5
        GaPi=0.70*Ga0*fpi*W0/W
	GaT =GaPi+0.30*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/6./1.D9
	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c ***************************************************************
	Q2G = Q2/1.D6
	CALL HF15(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR
c ***************************************************************
	IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN

        DE=DE*fk
        DM=DM*fk*(kgcm/qcm0)**0       
        DS=DS*fk

	ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE=DE*fk*(kgcm/qcm0)**(-1) 
        DM=DM*fk*(kgcm/qcm0)**(-1) 
	DS=DS*fk

	END IF
c ******************************************************

	ReM3(1)=Fak* (W0*W0-W*W)*DE*mPi
	ImM3(1)=Fak* W0*GaT*DE*mPi
	ReM3(2)=Fak* (W0*W0-W*W)*DM*mPi
	ImM3(2)=Fak* W0*GaT*DM*mPi
	ReM3(3) =wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
	ImM3(3) =wGacm/qcm*Fak* W0*GaT*DS*mPi
c *****************************************************

          DO 111 i=1,3
       ReM3M(i) = ReM3(i)*DCOS(Phi(i))-ImM3(i)*DSIN(Phi(i))
       ImM3M(i) = ReM3(i)*DSIN(Phi(i))+ImM3(i)*DCOS(Phi(i))
 111      CONTINUE

        RETURN
	END


       SUBROUTINE P31NEW(ReP31,ImP31)
C******************************************************************************
C      multipoles M1- and L1- for P31(1910)
C******************************************************************************
       IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION W0,Ga0,X,k2pi0,k2pi,w2thr,finel, PhiR
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReP31(3),ImP31(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReD3(3),ImD3(3)
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,DS,A1,S1
       INTEGER i
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1910.,Ga0=250., X=500.)

	mn = (mi+mf)/2.
        mPi=139.5685
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
        qcm0 = SQRT(((W0*W0+mi*mi)/(2*W0))**2-mi*mi)
        kgcm = (W*W-mi*mi)/2./W
        qcm=sqrt(wGacm**2+Q2)
	Q2G = Q2/1.D6

	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
C **********************************************
	fk = (kgcm/qcm0)**1*(qcm0*qcm0+X*X)/(kgcm*kgcm+X*X)
	fpi = (kcm/kcm0)**3*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**2
	finel=(k2pi/k2pi0)**6*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**4
        GaPi=0.25*Ga0*fpi*W0/W
	GaT =GaPi+0.75*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/2./1.D9
 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c **************************************************
	Q2G = Q2/1.D6
        CALL HP31(Q2G,PhiR,DM,DS,A1,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR
c *********************************************************
        DM=DM*fk
        DS=DS*fk
c *************************************************
	  ReD3(1) = Fak* (W0*W0-W*W)*DM*mPi
	  ImD3(1) = Fak* W0*GaT*DM*mPi
	  ReD3(2) = 0.
	  ImD3(2) = 0.
          ReD3(3) = wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
          ImD3(3) = wGacm/qcm*Fak* W0*GaT*DS*mPi
          DO 111 I=1,3
       ReP31(i) = ReD3(i)*DCOS(Phi(i))-ImD3(i)*DSIN(Phi(i))
       ImP31(i) = ReD3(i)*DSIN(Phi(i))+ImD3(i)*DCOS(Phi(i))
  111     CONTINUE

       RETURN
	END
c

	SUBROUTINE D15NEW(ISO,ReM2P,ImM2P)
C*****************************************************
C      multipoles M2+, E2+ and L2+ for D15(1675)
C*****************************************************
	IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT,Ga0,X,W0
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM2P(3),ImM2P(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM2(3),ImM2(3)
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,DS,A1,A3,S1
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel, PhiR
       INTEGER i, ISO
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1675.,Ga0=150., X=500.)


	mn = (mi+mf)/2.
c **************************************
        MPI=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
c *********************************************************
	fk = (kgcm/qcm0)**3*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**5*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**2
	finel=(k2pi/k2pi0)**8*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**4
        GaPi=0.45*Ga0*fpi*W0/W
	GaT =GaPi+0.55*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/6./1.D9
 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c *********************************************************
	Q2G = Q2/1.D6
	CALL HD15(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR
c ******************************************************
	  IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DE=DE*fk
        DM=DM*(kgcm/qcm0)**2.5*fk
        DS=DS*fk

	  ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE= DE*fk
        DM= DM*(kgcm/qcm0)*fk
  	DS= DS*fk
	  END IF
c **************************************************

        ReM2(1) =Fak* (W0*W0-W*W)*DE*mPi
	ImM2(1) =Fak* W0*GaT*DE*mPi
	ReM2(2) =Fak*(W0*W0-W*W)*DM*mPi
	ImM2(2) =Fak* W0*GaT*DM*mPi
	ReM2(3) =wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
	ImM2(3) =wGacm/qcm*Fak* W0*GaT*DS*mPi

          DO 111 i=1,3
       ReM2P(i) = ReM2(i)*DCOS(Phi(i))-ImM2(i)*DSIN(Phi(i))
       ImM2P(i) = ReM2(i)*DSIN(Phi(i))+ImM2(i)*DCOS(Phi(i))
 111      CONTINUE
       RETURN
	END

	SUBROUTINE F35NEW(ReM35,ImM35)
C*************************************************************
C      multipoles M3-, E3- and L3-  for F35(1905)
C*************************************************************
	IMPLICIT NONE
c ***************************************************************+
c       COMMON/parF3/ XEF35,XMF35,XSF35,XEF37,XMF37,XSF37
c       DOUBLE PRECISION XEF35,XMF35,XSF35,XEF37,XMF37,XSF37
c *****************************************************************
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT,W0, Ga0, X
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM35(3),ImM35(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM3(3),ImM3(3)
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel,DS,PhiR
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,A1,A3,S1
       INTEGER i
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1905.,Ga0=350., X=500.)

	mn = (mi+mf)/2.
c **************************************
        mPi=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
c ***************************************************************
	fk = (kgcm/qcm0)**5*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**7*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**3
	finel=(k2pi/k2pi0)**10*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**5
        GaPi=0.10*Ga0*fpi*W0/W
	GaT =GaPi+0.90*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/6./1.D9
	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c ***************************************************************
	Q2G = Q2/1.D6
        CALL HF35(Q2G,PhiR,DE,DM,DS,A1,A3,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR

c ***************************************************************

        DM=DM*fk
	DE=DE*fk*(qcm0/kgcm)
	DS=DS*fk*(qcm0/kgcm)
c ******************************************************

	ReM3(1)=Fak* (W0*W0-W*W)*DM*mPi
	ImM3(1)=Fak* W0*GaT*DM*mPi
	ReM3(2)=Fak* (W0*W0-W*W)*DE*mPi
	ImM3(2)=Fak* W0*GaT*DE*mPi
	ReM3(3) =wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
	ImM3(3) =wGacm/qcm*Fak* W0*GaT*DS*mPi

          DO 111 i=1,3
       ReM35(i) = ReM3(i)*DCOS(Phi(i))-ImM3(i)*DSIN(Phi(i))
       ImM35(i) = ReM3(i)*DSIN(Phi(i))+ImM3(i)*DCOS(Phi(i))
 111      CONTINUE

        RETURN
	END

	SUBROUTINE F37NEW(ReM37,ImM37)
C*****************************************************
C      multipoles M3+, E3+ and L3+  for F37(1950)
C*****************************************************
	IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT,W0, Ga0, X
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReM37(3),ImM37(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM3(3),ImM3(3)
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel,DS,PhiR
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,A1,A3,S1
       INTEGER i
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1945.,Ga0=280., X=500.)

	mn = (mi+mf)/2.
c **************************************
        mPi=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
c ***************************************************************
	fk = (kgcm/qcm0)**6*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**7*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))**3
	finel=(k2pi/k2pi0)**10*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**5
        GaPi=0.40*Ga0*fpi*W0/W
	GaT =GaPi+0.60*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/8./1.D9
	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c ***************************************************************
	Q2G = Q2/1.D6
        CALL HF37(Q2G,PhiR,DE,DM,DS,A1,A3,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR

c ***************************************************************

        DM=DM*fk
	DE=DE*fk
	DS=DS*fk
c ******************************************************

	ReM3(1)=Fak* (W0*W0-W*W)*DM*mPi
	ImM3(1)=Fak* W0*GaT*DM*mPi
	ReM3(2)=Fak* (W0*W0-W*W)*DE*mPi
	ImM3(2)=Fak* W0*GaT*DE*mPi
	ReM3(3) =wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
	ImM3(3) =wGacm/qcm*Fak* W0*GaT*DS*mPi

          DO 111 i=1,3
       ReM37(i) = ReM3(i)*DCOS(Phi(i))-ImM3(i)*DSIN(Phi(i))
       ImM37(i) = ReM3(i)*DSIN(Phi(i))+ImM3(i)*DCOS(Phi(i))
 111      CONTINUE

        RETURN
	END

	SUBROUTINE P13NEW(ISO,ReE1P,ImE1P)
C*******************************************
C      multipoles E1+ and L1+ for P13(1720)
C*******************************************
	IMPLICIT NONE
c ***************************************************************+
       DOUBLE PRECISION fk,fpi,DE,DM,Fak,GaGa,GaPi,GaT,Ga0,X,W0
       DOUBLE PRECISION mi,mf,mPION,OmegL,Q2,Q2G,W,qcm,kcm,qcm0,kgcm
       DOUBLE PRECISION Phi(3),kcm0,Ga,Gpi,ReE1P(3),ImE1P(3),mPi
       DOUBLE PRECISION Pi,mn,mdip,wGacm,kx,ReM2(3),ImM2(3)
       DOUBLE PRECISION TE,TM,SE,SM,FQE,FQM,A1,A3,S1
       DOUBLE PRECISION k2pi0,k2pi,w2thr,finel,DS,PhiR
       INTEGER i, ISO
          COMMON /Mass/ mi, mf, mPION
          COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
       PARAMETER (Pi=3.14159265358979D0, mdip=0.71)
       PARAMETER(W0=1740.,Ga0=250., X=500.)


	mn = (mi+mf)/2.
c **************************************
        MPI=139.5685
        qcm=sqrt(wGacm**2+Q2)
c **************************************
	kcm0 = SQRT(((W0*W0+mPi*mPi-mf*mf)/(2*W0))**2-mPi*mPi)
	kgcm = (W*W-mi*mi)/2./W
	qcm0 = SQRT(((W0*W0+mi*mi)/2./W0)**2-mi*mi)
	k2pi0=SQRT(((W0*W0+4.*mPi*mPi-mf*mf)/(2*W0))**2-4.*mPi*mPi)
        w2thr=mf+2.*mPi
        k2pi=0.
        if (W.gt.w2thr)
     *	k2pi=SQRT(((W*W+4.*mPi*mPi-mf*mf)/(2*W))**2-4.*mPi*mPi)
c *********************************************************
	fk = (kgcm/qcm0)**3*((qcm0*qcm0+X*X)/(kgcm*kgcm+X*X))
	fpi = (kcm/kcm0)**3*((kcm0*kcm0+X*X)/(kcm*kcm+X*X))
	finel=(k2pi/k2pi0)**6*((k2pi0*k2pi0+X*X)/
     /  (k2pi*k2pi+X*X))**3
        GaPi=0.20*Ga0*fpi*W0/W
	GaT =GaPi+0.80*Ga0*finel
	GaGa = kgcm*qcm/Pi*mn/W0/4./1.D9
 	Fak = W0/((W0*W0-W*W)**2+W0*W0*GaT*GaT)
     * *SQRT(GaGa*GaPi/qcm/kcm)
c *********************************************************
	Q2G = Q2/1.D6
        CALL HP13(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,0)
        Phi(1) = PhiR
        Phi(2) = PhiR
        Phi(3) = PhiR
c ******************************************************
	  IF (ISO.EQ.1 .OR. ISO.EQ.3) THEN
        DE=DE*fk
        DM=DM*fk
	DS=DS*fk

	  ELSE IF (ISO.EQ.2 .OR. ISO.EQ.4) THEN
        DE=DE*fk 
        DM=DM*fk *(kgcm/qcm0)**0
	DS=DS*fk
	  END IF
c **************************************************

        ReM2(1) =Fak* (W0*W0-W*W)*DE*mPi
	ImM2(1) =Fak* W0*GaT*DE*mPi
	ReM2(2) =Fak*(W0*W0-W*W)*DM*mPi
	ImM2(2) =Fak* W0*GaT*DM*mPi
	ReM2(3) =wGacm/qcm*Fak* (W0*W0-W*W)*DS*mPi
	ImM2(3) =wGacm/qcm*Fak* W0*GaT*DS*mPi

          DO 111 i=1,3
       ReE1P(i) = ReM2(i)*DCOS(Phi(i))-ImM2(i)*DSIN(Phi(i))
       ImE1P(i) = ReM2(i)*DSIN(Phi(i))+ImM2(i)*DCOS(Phi(i))
 111      CONTINUE
       RETURN
	END


      SUBROUTINE GAUSSP(N,E,W)
      IMPLICIT REAL*8 (A-H,O-Z)
C     GAUSSP BERECHNET ABSZISSEN UND GEWICHTE FUER GAUSSINTEGRATION
C     IM INTERVALL "-1,1@ IN AUFSTEIGENDER REIHENFOLGE
C     REFERENZ:   ABRAMOWITZ AND STEGUN
C                 HANDBOOK OF MATHEMATICAL FUNCTIONS
C                 IM FOLGENDEN ABGEKUERZT A+S
C     DEFINITION VON ABSZISSEN UND GEWICHTEN: A+S(25.4.29)
C     METHODE:
C             BESTIMMUNG DER NULLSTELLEN VON PN DURCH NEWTONSCHE
C             ITERATION
C     GENAUIGKEIT:     1.E-15
C     N:     ANZAHL DER PUNKTE
C     E      RESULTAT:   ABSZISSEN
C     W      RESULTAT:   GEWICHTE
      DIMENSION E(1),W(1)
      DATA PI/3.141592653589793238462643D0/,EPS/1.D-16/
      M=(N+1)/2
      DN=N
      DO 1000 I=1,M
      DI=I
C     STARTWERTE A+S(22.16.6)
      X=PI*(4.D0*(DN-DI)+3.D0)/(4.D0*DN+2.D0)
      XN=(1.D0-(DN-1.D0)/(8.D0*DN*DN*DN))*DCOS(X)
      IF(I.GT.N/2) XN=0
      DO 100  ITER=1,10
      X=XN
C     BERECHNUNG VON PN(X) DURCH ITERATION A+S(8.5.3)
      Y1=1.D0
      Y=X
      IF(N.LT.2) GOTO 250
      DO 200 J=2,N
      DJ=J
      Y2=Y1
      Y1=Y
200   Y=((2.D0*DJ-1.D0)*X*Y1-(DJ-1.D0)*Y2)/DJ
250   CONTINUE
C     BERECHNUNG DER ABLEITUNG A+S(8.5.4)
      YS=DN*(X*Y-Y1)/(X*X-1.D0)
C     NEWTONITERATION
      H=-Y/YS
      XN=X+H
      IF (DABS(H) .LT. EPS) GOTO 110
  100 CONTINUE
  110 E(I)=X
      E(N-I+1)=-X
      GEW=2.D0/((1.D0-X*X)*YS*YS)
      W(I)=GEW
      W(N-I+1)=GEW
1000  CONTINUE
      RETURN
      END

      SUBROUTINE LET_CORR(W,Q2,MM,MMPI,Epl,Spl,E0_corr,S0_corr,xpi)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MM, MMPI, kcm
      COMMON /KinVar/ OmegL, Q2MEV, WWW, wGacm, kcm
      COMMON /E0CORR/ XE, XS, XMIX
c ***********************************
      HQC = 197.3285
      ampi=139.4999/hqc
      amn=939.5/hqc
c here we use the same mass for all channels
      apiN=0.124/ampi
      qpi_fm = 0.d0
      W_MEV=W*HQC
      W0=1090./hqc
      qpi0 = SQRT(((W0*W0+ampi**2-amn**2)/(2*W0))**2-ampi**2)
      IF (W_MeV.GE.1079.)
     & qpi_fm = SQRT(((W*W+ampi**2-amn**2)/(2*W))**2-ampi**2)
      xpi=qpi_fm/qpi0
c ******  from fit of Fucs data  ********
c      Ae=0.135E-4
c      Be=0.68
c      Ce=-0.14668
c *******  from fit of new Schmidt data *****
      Ae=0.14415E-4
      Be=0.70868
      Ce=2.5E-5
c *******************************************
      Wres=1232./HQC
      qres = sqrt(((Wres**2+MMPI**2-MM**2)/(2*Wres))**2-MMPI**2)
      qpi = SQRT(((W*W+MMPI**2-MM**2)/(2*W))**2-MMPI**2)
        	egcm=(W**2-MM**2-Q2)/2./W
                qcm=sqrt(egcm**2+Q2)
	        Q2GEV=Q2*(197.3285/1000.)**2
	        FDIP = 1./(1.+Q2GEV/0.71)**2
            CALL FORM (Q2GEV,D,F1S,F1V,F2S,F2V,FPI,FAX)
    
      E0_corr = Ae/((1.+Be**2*qpi**2)**2)*F2V/3.7058 * XE
      Q2pt = -(W-MM)**2*HQC**2*1.D-6
      S0_corr = E0_corr * egcm/(W-MM)*XS      
     $ * dexp(-10.*Q2GEV)/dexp(-10.*Q2pt) 
c      ccc = dexp(-36.*Q2GEV)/dexp(-36.*Q2pt)*XS        
c      write (6,167) W, ccc
c 167  format(1x,2(e10.4,2x))     
c **********  cusp effects (same for proton and neutron) ***********
       qpi_eff=0.
      IF (W_MeV.LT.1079.)
     & qpi_eff = SQRT(-((W*W+ampi**2-amn**2)/(2*W))**2+ampi**2)
      IF (W_MeV.LT.1079.) E0_corr = E0_corr - Epl*apiN*qpi_eff
      IF (W_MeV.LT.1079.) S0_corr = S0_corr - Spl*apiN*qpi_eff
c *******************************************************************

      RETURN
      END

      SUBROUTINE SPLINES()
      IMPLICIT NONE
      INTEGER i
      REAL*8 a
      REAL*8 W_sp(143)
      REAL*8 del_s11(143),del_s11_ss(143),del_s31(143),del_s31_ss(143)
      REAL*8 del_p11(143),del_p11_ss(143),del_p31(143),del_p31_ss(143)
      REAL*8 del_p13(143),del_p13_ss(143),del_p33(143),del_p33_ss(143)
      REAL*8 del_d13(143),del_d13_ss(143),del_d33(143),del_d33_ss(143)
      REAL*8 del_d15(143),del_d15_ss(143),del_d35(143),del_d35_ss(143)
      REAL*8 del_f15(143),del_f15_ss(143),del_f35(143),del_f35_ss(143)
      REAL*8 del_f17(143),del_f17_ss(143),del_f37(143),del_f37_ss(143)
      REAL*8 eta_s11(143),eta_s11_ss(143),eta_s31(143),eta_s31_ss(143)
      REAL*8 eta_p11(143),eta_p11_ss(143),eta_p31(143),eta_p31_ss(143)
      REAL*8 eta_p13(143),eta_p13_ss(143),eta_p33(143),eta_p33_ss(143)
      REAL*8 eta_d13(143),eta_d13_ss(143),eta_d33(143),eta_d33_ss(143)
      REAL*8 eta_d15(143),eta_d15_ss(143),eta_d35(143),eta_d35_ss(143)
      REAL*8 eta_f15(143),eta_f15_ss(143),eta_f35(143),eta_f35_ss(143)
      REAL*8 eta_f17(143),eta_f17_ss(143),eta_f37(143),eta_f37_ss(143)
      COMMON /SPL/ W_sp,del_s11,del_s11_ss,del_s31,del_s31_ss,
     *      del_p11,del_p11_ss,del_p31,del_p31_ss,
     *      del_p13,del_p13_ss,del_p33,del_p33_ss,
     *      del_d13,del_d13_ss,del_d33,del_d33_ss,
     *      del_d15,del_d15_ss,del_d35,del_d35_ss,
     *      del_f15,del_f15_ss,del_f35,del_f35_ss,
     *      del_f17,del_f17_ss,del_f37,del_f37_ss,
     *      eta_s11,eta_s11_ss,eta_s31,eta_s31_ss,
     *      eta_p11,eta_p11_ss,eta_p31,eta_p31_ss,
     *      eta_p13,eta_p13_ss,eta_p33,eta_p33_ss,
     *      eta_d13,eta_d13_ss,eta_d33,eta_d33_ss,
     *      eta_d15,eta_d15_ss,eta_d35,eta_d35_ss,
     *      eta_f15,eta_f15_ss,eta_f35,eta_f35_ss,
     *      eta_f17,eta_f17_ss,eta_f37,eta_f37_ss
      REAL*8 yp1, ypn,pi
      PARAMETER(yp1=1.0E30, ypn = 1.0E30,pi = 3.1415926536D0)

            READ(2, 2111)
            READ(2, 2111)
 2111        FORMAT(20X)
      DO 2112 i = 1, 143
         READ(2, *) W_sp(i), del_s11(i),a,eta_s11(i)
         del_s11(i) = del_s11(i) * pi / 180.0
         eta_s11(i) = SQRT(1.0 - eta_s11(i))
 2112  CONTINUE
      CALL SPLINE(W_sp, del_s11, 143, yp1, ypn, del_s11_ss)
      CALL SPLINE(W_sp, eta_s11, 143, yp1, ypn, eta_s11_ss)

            READ(2, 2211)
            READ(2, 2211)
            READ(2, 2211)
 2211        FORMAT(20X)
      DO 2212 i = 1, 143
         READ(2, *) W_sp(i), del_s31(i),a,eta_s31(i)
         del_s31(i) = del_s31(i) * pi / 180.0
         eta_s31(i) = SQRT(1.0 - eta_s31(i))
 2212  CONTINUE
      CALL SPLINE(W_sp, del_s31, 143, yp1, ypn, del_s31_ss)
      CALL SPLINE(W_sp, eta_s31, 143, yp1, ypn, eta_s31_ss)

            READ(2, 2311)
            READ(2, 2311)
            READ(2, 2311)
 2311        FORMAT(20X)
      DO 2312 i = 1, 143
         READ(2, *) W_sp(i), del_p11(i),a,eta_p11(i)
         del_p11(i) = del_p11(i) * pi / 180.0
         eta_p11(i) = SQRT(1.0 - eta_p11(i))
 2312  CONTINUE
      CALL SPLINE(W_sp, del_p11, 143, yp1, ypn, del_p11_ss)
      CALL SPLINE(W_sp, eta_p11, 143, yp1, ypn, eta_p11_ss)

            READ(2, 2411)
            READ(2, 2411)
            READ(2, 2411)
 2411        FORMAT(20X)
      DO 2412 i = 1, 143
         READ(2, *) W_sp(i), del_p31(i),a,eta_p31(i)
         del_p31(i) = del_p31(i) * pi / 180.0
         eta_p31(i) = SQRT(1.0 - eta_p31(i))
 2412  CONTINUE
      CALL SPLINE(W_sp, del_p31, 143, yp1, ypn, del_p31_ss)
      CALL SPLINE(W_sp, eta_p31, 143, yp1, ypn, eta_p31_ss)

            READ(2, 2511)
            READ(2, 2511)
            READ(2, 2511)
 2511        FORMAT(20X)
      DO 2512 i = 1, 143
         READ(2, *) W_sp(i), del_p13(i),a,eta_p13(i)
         del_p13(i) = del_p13(i) * pi / 180.0
         eta_p13(i) = SQRT(1.0 - eta_p13(i))
 2512  CONTINUE
      CALL SPLINE(W_sp, del_p13, 143, yp1, ypn, del_p13_ss)
      CALL SPLINE(W_sp, eta_p13, 143, yp1, ypn, eta_p13_ss)

            READ(2, 2611)
            READ(2, 2611)
            READ(2, 2611)
 2611        FORMAT(20X)
      DO 2612 i = 1, 143
         READ(2, *) W_sp(i), del_p33(i),a,eta_p33(i)
         del_p33(i) = del_p33(i) * pi / 180.0
         IF (del_p33(i).LT.0.0.OR.W_sp(i).GT.1800.0)
     *       del_p33(i) = del_p33(i) + pi
         eta_p33(i) = SQRT(1.0 - eta_p33(i))
 2612  CONTINUE
      CALL SPLINE(W_sp, del_p33, 143, yp1, ypn, del_p33_ss)
      CALL SPLINE(W_sp, eta_p33, 143, yp1, ypn, eta_p33_ss)

            READ(2, 2811)
            READ(2, 2811)
            READ(2, 2811)
 2811        FORMAT(20X)
      DO 2812 i = 1, 143
         READ(2, *) W_sp(i), del_d13(i),a,eta_d13(i)
         del_d13(i) = del_d13(i) * pi / 180.0
         eta_d13(i) = SQRT(1.0 - eta_d13(i))
 2812  CONTINUE
      CALL SPLINE(W_sp, del_d13, 143, yp1, ypn, del_d13_ss)
      CALL SPLINE(W_sp, eta_d13, 143, yp1, ypn, eta_d13_ss)

            READ(2, 2911)
            READ(2, 2911)
            READ(2, 2911)
 2911        FORMAT(20X)
      DO 2912 i = 1, 143
         READ(2, *) W_sp(i), del_d33(i),a,eta_d33(i)
         del_d33(i) = del_d33(i) * pi / 180.0
         eta_d33(i) = SQRT(1.0 - eta_d33(i))
 2912  CONTINUE
      CALL SPLINE(W_sp, del_d33, 143, yp1, ypn, del_d33_ss)
      CALL SPLINE(W_sp, eta_d33, 143, yp1, ypn, eta_d33_ss)

            READ(2, 3011)
            READ(2, 3011)
            READ(2, 3011)
 3011        FORMAT(20X)
      DO 3012 i = 1, 143
         READ(2, *) W_sp(i), del_d15(i),a,eta_d15(i)
         del_d15(i) = del_d15(i) * pi / 180.0
         eta_d15(i) = SQRT(1.0 - eta_d15(i))
 3012  CONTINUE
      CALL SPLINE(W_sp, del_d15, 143, yp1, ypn, del_d15_ss)
      CALL SPLINE(W_sp, eta_d15, 143, yp1, ypn, eta_d15_ss)

            READ(2, 3111)
            READ(2, 3111)
            READ(2, 3111)
 3111        FORMAT(20X)
      DO 3112 i = 1, 143
         READ(2, *) W_sp(i), del_d35(i),a,eta_d35(i)
         del_d35(i) = del_d35(i) * pi / 180.0
         eta_d35(i) = SQRT(1.0 - eta_d35(i))
 3112  CONTINUE
      CALL SPLINE(W_sp, del_d35, 143, yp1, ypn, del_d35_ss)
      CALL SPLINE(W_sp, eta_d35, 143, yp1, ypn, eta_d35_ss)

            READ(2, 3211)
            READ(2, 3211)
            READ(2, 3211)
 3211        FORMAT(20X)
      DO 3212 i = 1, 143
         READ(2, *) W_sp(i), del_f15(i),a,eta_f15(i)
         del_f15(i) = del_f15(i) * pi / 180.0
         eta_f15(i) = SQRT(1.0 - eta_f15(i))
 3212  CONTINUE
      CALL SPLINE(W_sp, del_f15, 143, yp1, ypn, del_f15_ss)
      CALL SPLINE(W_sp, eta_f15, 143, yp1, ypn, eta_f15_ss)

            READ(2, 3311)
            READ(2, 3311)
            READ(2, 3311)
 3311        FORMAT(20X)
      DO 3312 i = 1, 143
         READ(2, *) W_sp(i), del_f35(i),a,eta_f35(i)
         del_f35(i) = del_f35(i) * pi / 180.0
         eta_f35(i) = SQRT(1.0 - eta_f35(i))
 3312  CONTINUE
      CALL SPLINE(W_sp, del_f35, 143, yp1, ypn, del_f35_ss)
      CALL SPLINE(W_sp, eta_f35, 143, yp1, ypn, eta_f35_ss)

            READ(2, 3411)
            READ(2, 3411)
            READ(2, 3411)
 3411        FORMAT(20X)
      DO 3412 i = 1, 143
         READ(2, *) W_sp(i), del_f17(i),a,eta_f17(i)
         del_f17(i) = del_f17(i) * pi / 180.0
         eta_f17(i) = SQRT(1.0 - eta_f17(i))
 3412  CONTINUE
      CALL SPLINE(W_sp, del_f17, 143, yp1, ypn, del_f17_ss)
      CALL SPLINE(W_sp, eta_f17, 143, yp1, ypn, eta_f17_ss)

            READ(2, 3511)
            READ(2, 3511)
            READ(2, 3511)
 3511        FORMAT(20X)
      DO 3512 i = 1, 143
         READ(2, *) W_sp(i), del_f37(i),a,eta_f37(i)
         del_f37(i) = del_f37(i) * pi / 180.0
         eta_f37(i) = SQRT(1.0 - eta_f37(i))
 3512  CONTINUE
      CALL SPLINE(W_sp, del_f37, 143, yp1, ypn, del_f37_ss)
      CALL SPLINE(W_sp, eta_f37, 143, yp1, ypn, eta_f37_ss)

      END

      SUBROUTINE BORN_MUL(IPVPS,W,Q2,M,MPI,LMAX,FB)
c *****************************************************************************
c     Born multipoles l=0,..,lmax
c *****************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 M,MPI,MPIP,K0
      DIMENSION FPOL(6,10,3),AEP(3),AEM(3),AMP(3),AMM(3)
      DIMENSION FB(6,10,3),ALP(3),ALM(3)
c   IPVPS            :  0   PV coupling in the Born terms
c                    :  1   PS coupling in the Born terms
c                    :  2   mixed coupling (realistic)
      LMAX1=LMAX+1
      HQC=197.3285D0
      MPIP=139.5685D0/HQC
      K0=(W**2-M**2-Q2)/2.D0/W
      
c      Q2GEV=Q2*(197.337/1000.)**2
      Q2GEV=Q2*(HQC/1000.D0)**2
      CALL FORM (Q2GEV,D,F1S,F1V,F2S,F2V,FPI,FAX)
c  e.m. form factors are also calculated with FORM inside FPOLE 
c   FPOLE calculates multipoles from the pole (PS) terms
      FFA=(FAX-F1V)/Q2
      CALL FPOLE(W,Q2,M,MPI,LMAX1,FPOL)

      DO 1001 L1=1,LMAX1
      L=L1-1

      DO 99 ISO=1,3
      AEP(ISO)=0.
      AEM(ISO)=0.
      AMP(ISO)=0.
      AMM(ISO)=0.
      ALP(ISO)=0.
 99   ALM(ISO)=0.


      CL1=2*L1/MPIP/1000.
      CL=2*L/MPIP/1000.

      DO 1 ISO=1,3
      AEP(ISO)=FPOL(1,1,ISO)/CL1
      ALP(ISO)=K0*FPOL(5,L1,ISO)*MPIP*1000./L1
      IF (L1.EQ.1) GO TO 1
      AEP(ISO)=(FPOL(1,L1,ISO)-L*FPOL(3,L1,ISO))/CL1
      AMP(ISO)=(FPOL(1,L1,ISO)+(L+2)*FPOL(3,L1,ISO))/CL1
      AMM(ISO)=-FPOL(2,1,ISO)/CL1
      ALM(ISO)=K0*FPOL(6,L,ISO)*MPIP*1000./L
      IF (L.LT.1) GO TO 1
      AMM(ISO)=(-FPOL(2,L,ISO)+(L-1)*FPOL(4,L,ISO))/CL
       IF (L.LT.2) GO TO 1
      AEM(ISO)=(FPOL(2,L,ISO)+(L+1)*FPOL(4,L,ISO))/CL
 1    CONTINUE

c *********  PV-PS couplings *************
      DF1=0.
      DF2=0.
      DF5=0.
      DF6=0.
       DF1m=0.D0
       DF2m=0.D0
       DF5m=0.D0
       DF6m=0.D0

      IF (L.GT.1) GO TO 111
c  changes on 16. Nov 2006 (new FFR terms)
      CALL DFPVPS(IPVPS,W,Q2,M,MPI,DF1,DF2,DF5,DF6,DF1m,DF2m,DF5m,DF6m)   
      IF (L.EQ.1) GO TO 222
      AEP(1)=AEP(1)+F2V*DF1*MPIP*1000.
      AEP(2)=AEP(2)+F2S*DF1*MPIP*1000.
      AEP(3)=AEP(3)+FFA*DF1m*MPIP*1000.

      ALP(1)=ALP(1)+F2V*DF5*MPIP*1000.
      ALP(2)=ALP(2)+F2S*DF5*MPIP*1000.
      ALP(3)=ALP(3)+FFA*DF5m*MPIP*1000.
      GO TO 111
 222  AMM(1)=AMM(1)+F2V*DF2*MPIP*1000.
      AMM(2)=AMM(2)+F2S*DF2*MPIP*1000.
      AMM(3)=AMM(3)+FFA*DF2m*MPIP*1000.

      ALM(1)=ALM(1)+F2V*DF6*MPIP*1000.
      ALM(2)=ALM(2)+F2S*DF6*MPIP*1000.
      ALM(3)=ALM(3)+FFA*DF6m*MPIP*1000.
 111  CONTINUE
C **************************************
C  AEP(1),AEP(2),AEP(3)  = E0+(+,0,-)
C  AEM=EL-, AMP=ML+, ALP=LL+, etc
c    (1/2) multipoles:
        FB(1,L1,1)=AEP(1)+2.*AEP(3)
        FB(2,L1,1)=AEM(1)+2.*AEM(3)
        FB(3,L1,1)=AMP(1)+2.*AMP(3)
        FB(4,L1,1)=AMM(1)+2.*AMM(3)
        FB(5,L1,1)=ALP(1)+2.*ALP(3)
        FB(6,L1,1)=ALM(1)+2.*ALM(3)
c    (3/2) multipoles:
        FB(1,L1,3)=AEP(1)-AEP(3)
        FB(2,L1,3)=AEM(1)-AEM(3)
        FB(3,L1,3)=AMP(1)-AMP(3)
        FB(4,L1,3)=AMM(1)-AMM(3)
        FB(5,L1,3)=ALP(1)-ALP(3)
        FB(6,L1,3)=ALM(1)-ALM(3)
c    (0) multipoles:
        FB(1,L1,2)=AEP(2)
        FB(2,L1,2)=AEM(2)
        FB(3,L1,2)=AMP(2)
        FB(4,L1,2)=AMM(2)
        FB(5,L1,2)=ALP(2)
        FB(6,L1,2)=ALM(2)

 1001  CONTINUE
       RETURN
       END


      SUBROUTINE DFPVPS(IPVPS,WFM,LAM2,M,MPI,DF1,DF2,DF5,DF6,
     &   DF1m,DF2m,DF5m,DF6m)
c *******************  PS-PV mixing  ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /E0CORR/ XA, XB, XMIX
      REAL*8 M,MPI,LAM2,K0,K2,KF,NU,Kopplung,MPIP
      PARAMETER (Pi=3.1415926536D0,eLadung=.3028619D0,
     * Kopplung=0.996365D0, apvps0=450.,HQC=197.3285 )
C  here the FFR multipoles are calculated for E0+, L0+, M1-, L1-
C  Kopplung = f = g*MPI/(2*M)
C  see electroproduction paper by B.P.,D.D.,L.T. EPJA 2007

      MPIP=139.5685D0/HQC
       DF1=0.D0
       DF2=0.D0
       DF5=0.D0
       DF6=0.D0
       DF1m=0.D0
       DF2m=0.D0
       DF5m=0.D0
       DF6m=0.D0
       
      if (ipvps.eq.1) return

      apvps=apvps0*XMIX
      W=WFM
      SIG=W**2-M**2
      BET=SIG+LAM2/2.
      K0=(SIG-LAM2)/2./W
      K2=K0**2+LAM2
      KF=DSQRT(K2)
      Q0=(SIG+MPI**2)/2./W
      Q2=Q0**2-MPI**2
      Q=DSQRT(Q2)
C  Q is pion and KF is photon momentum, LAM2=QQ**2
      E1=DSQRT(M**2+K2)
      E2=DSQRT(M**2+Q2)

        fpvps=apvps**2/(apvps**2+Q2*hqc**2)
        IF(IPVPS.EQ.0) FPVPS=1.

        EFIP=DSQRT((E1+M)*(E2+M))
        EFIM=DSQRT((E1-M)/(E2+M))
c	C0 =fpvps*eLadung*Kopplung/(8.*Pi*W*M)
	C1 =fpvps*eLadung*Kopplung/(8.*Pi*W*MPIP)
	C0 =C1/(2*M)
        DF1= EFIP*C0*(W-M)
        DF2=-EFIM*C0*(W+M)*Q
        DF5= EFIP*C0*K0
        DF6=-EFIM*C0*Q*K0
        DF1m= EFIP*C1*LAM2
        DF2m= EFIM*C1*Q*LAM2
        DF5m=-EFIP*C1*(W-M)*K0
        DF6m=-EFIM*C1*(W+M)*Q*K0

	RETURN
	END


      SUBROUTINE FPOLE(WFM,LAM2,M,MPI,LMAX1,FPOL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 M,MPI,LAM2
      DIMENSION HH1(12,3),HH3(12,3),HH5(12,3),FPOL(6,10,3)
      DIMENSION HH1M(12,3),HH3M(12,3),HH5M(12,3)

      CALL ODDH(WFM,LAM2,M,MPI,LMAX1+1,HH1,HH3,HH5)
      CALL ODDH(-WFM,LAM2,M,MPI,LMAX1+1,HH1M,HH3M,HH5M)

      DO 1 ISO=1,3
      DO 1 L1=1,LMAX1

      FPOL(1,L1,ISO) = HH1(L1,ISO) + HH1M(L1+1,ISO)
      FPOL(2,L1,ISO) = HH1M(L1,ISO) + HH1(L1+1,ISO)
      FPOL(3,L1,ISO) = HH3(L1,ISO) - HH3M(L1+1,ISO)
      FPOL(4,L1,ISO) =- HH3M(L1,ISO) + HH3(L1+1,ISO)
      FPOL(5,L1,ISO) = HH5(L1,ISO) + HH5M(L1+1,ISO)
      FPOL(6,L1,ISO) = HH5M(L1,ISO) + HH5(L1+1,ISO)
 1    CONTINUE

      RETURN
      END

      SUBROUTINE ODDH(WFM,LAM2,M,MPI,LMAX1,HH1,HH3,HH5)
c ***************  helicity multipoles  ****************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 M,MPI,LAM2,K0,K2,KF,NU,MPIP
      DIMENSION HM(6),FIJ(6,12),NEPS(3)
      DIMENSION QLZ(12),QLY(12),HH1(12,3),HH3(12,3),HH5(12,3)
      DATA NEPS/1,1,-1/
      PARAMETER (Pi=3.1415926536D0,er=.3028619D0,fr=0.996365D0)

      W=WFM
      COF=8.*PI*W
      SIG=W**2-M**2
      BET=SIG+LAM2/2.
      K0=(SIG-LAM2)/2./W
      K2=K0**2+LAM2
      KF=DSQRT(K2)
      Q0=(SIG+MPI**2)/2./W
      Q2=Q0**2-MPI**2
      Q=DSQRT(Q2)
      NU=-(W+M)
      OM=W-M
      E1=W+M-K0
      E2=W+M-Q0
      G1=-W+M+K0
      G2=-W+M+Q0

      MPIP=139.5685/197.3285
      GPIN=2.*M/MPIP*fr
      QKF2=2*Q*KF
      Z=(2.*Q0*K0+LAM2)/QKF2
      Y=(2.*Q0*K0-SIG)/QKF2

      Q2GEV=LAM2*(197.3285/1000.)**2
      CALL FORM (Q2GEV,D,F1S,F1V,F2S,F2V,FPI,FAX)
      CALL HFMATR(W,LAM2,M,MPI,LMAX1+1,HM,FIJ)
      CALL QLEG(Z,LMAX1+1,QLZ)
      CALL QLEG(Y,LMAX1+1,QLY)

      DO 1 ISO=1,3

      F1I=F1V*er
      F2I=F2V*er/2./M
      FPII=0.5*er*(1.-NEPS(ISO))*FPI
      IF (ISO.NE.2) GO TO 111
      F1I=F1S*er
      F2I=F2S*er/2./M
 111  CONTINUE

      A1 = GPIN/K2 * FPII * (LAM2+2.*NU*(Q0-K0)) / QKF2
      B1 = GPIN/2/K2 * NEPS(ISO) * ( (SIG+2*NU*Q0)*F1I +
     + ((NU+2*K0)*SIG + 2*Q0*LAM2) * F2I ) / QKF2
      C1 = GPIN/2/K2* ((F1I-OM*F2I)*2*E1/NU + 2*FPII +
     + NEPS(ISO)*(F1I + (NU+2*K0)*F2I))

      A3 = -2*GPIN * FPII/QKF2
      B3 = -GPIN * NEPS(ISO)*(F1I-NU*F2I)/QKF2

      A5 = GPIN/K2 *FPII*(2*Q0-K0)/QKF2
      B5 = - GPIN/K2/2. * NEPS(ISO) * ((OM-2*Q0)*F1I +
     + ((2*Q0-K0)*NU+SIG+LAM2)*F2I)/QKF2
      C5 =  GPIN/K2/2.*((2*FPII-(1-NEPS(ISO))*F1I)*K0/LAM2 +
     + (F1I+(E1-NEPS(ISO)*NU)*F2I)/NU)

      DL=1.
      RLZ=0.
      RLY=0.

      DO 2 L1=1,LMAX1
      IF (L1.EQ.1) GO TO 222
      DL=0.
      RLZ=(QLZ(L1-1)-QLZ(L1+1))/(2*L1-1)
      RLY=(QLY(L1-1)-QLY(L1+1))/(2*L1-1)
 222  HH1(L1,ISO)=HM(1)*(A1*QLZ(L1)+B1*QLY(L1)+C1*DL)
      HH3(L1,ISO)=HM(3)*(A3*RLZ + B3*RLY)
      HH5(L1,ISO)=HM(5)*(A5*QLZ(L1)+B5*QLY(L1)+C5*DL)

  2   CONTINUE
  1   CONTINUE

      RETURN
      END


      SUBROUTINE QLEG(Z,LMAX1,QL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION QL(1000)


      DO 1 L1=1,LMAX1
      L=L1-1
1     QL(L1)=QLEGY(Z,L)

      RETURN
      END

      FUNCTION QLEGY(ZZ,L)
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DURALT(11),DALT(11),DNEU(11),C(20)
C     WRITE(6,2)L,Z
    2 FORMAT(' INPUT OF QLEG L,Z=',I2,D20.8)
      Z=ZZ
      SIGN=Z/DABS(Z)
      IF(Z .LT. 0.0) Z=-Z
      IF((Z-1.)*(L+1).LT.2.) GO TO 100
      MMAX=20
      EPS=2.D-7
      LG=L+1
      X=1./Z
      XX=X*X
      SUM=0
      TERM=1.0
      M=0
      C(1)=1.0
      IF(L.EQ.0) GO TO 5
      DO 4 N=2,LG
      XN=N-2
4     C(N)=(XN+1.)/(2.*XN+3.)*C(N-1)
5     SUM=SUM+C(L+1)*TERM
      NALT=1
      DALT(1)=SUM
      GO TO 14
10    NURALT=NALT
      DO 11 N=1,NURALT
11    DURALT(N)=DALT(N)
      NALT=NNEU
      DO 12 N=1,NALT
12    DALT(N)=DNEU(N)
14    M=M+1
      XM=M-1
      TERM=XX*TERM
      C(1)=1./(2.*XM+3.)
      IF(L.EQ.0) GO TO 16
      DO 15 N=2,LG
      XN=N-1
15    C(N)=((2.*XM+XN+1.)*C(N)+XN*C(N-1))/(2.*(XN+XM)+3.)
16    SUM=SUM+C(L+1)*TERM
      NNEU=M/2+1
      DNEU(1)=SUM
      IF(M.LE.1) GO TO 10
      Y1=DNEU(1)-DALT(1)
      Y2=DURALT(1)-DALT(1)
      Z1=1.D99
      IF(Y1.NE.0) Z1=1./Y1
      Z2=1.D99
      IF(Y2.NE.0) Z2=1./Y2
      DNEU(2)=DALT(1)+1./(Z1+Z2)
      IF(NNEU.LE.2) GO TO 20
      DO 18 N=3,NNEU
      Y1=DNEU(N-1)-DALT(N-1)
      Y2=DURALT(N-1)-DALT(N-1)
      Y3=DURALT(N-2)-DALT(N-1)
      Z1=1.D99
      IF(Y1.NE.0) Z1=1./Y1
      Z2=1.D99
      IF(Y2.NE.0) Z2=1./Y2
      Z3=1.D99
      IF(Y3.NE.0) Z3=1./Y3
      DNEU(N)=DALT(N-1)+1./(Z1+Z2-Z3)
18    CONTINUE
20    CONTINUE
      IF(M.GE.MMAX) GO TO 30
      IF(ABS(DNEU(NNEU)-DNEU(NNEU-1)).LT.EPS*ABS(DNEU(NNEU))) GO TO 40
      GO TO 10
30    WRITE(6,31) MMAX,EPS,L,Z
31    FORMAT(////10X,'WARNUNG :'/10X,'IN QLEG KONNTE NACH',I3/10X,
     *	     'SCHRITTEN DIE GENAUIGKEIT VON ',E8.2,/10X,'FUER Q',I2,
     *	     1H(,E12.6,') NICHT ERREICHT WERDEN'/)
40    QLEGY=X**(L+1)*DNEU(NNEU)
      QLEGY=SIGN**(L+1)*QLEGY
C     PRINT 2,L,Z,QLEG
      RETURN
100   Y=0.5*DLOG((Z+1.D0)/(Z-1.D0))
      IF(L.EQ.0) GO TO 120
      YALT=Y
      Y=Z*Y-1.0
      IF(L.EQ.1) GO TO 120
      DO 110 N=2,L
      XN=N-1
      YURALT=YALT
      YALT=Y
110   Y=((2.*XN+1.)*Z*YALT-XN*YURALT)/(XN+1.)
120   QLEGY=Y
      QLEGY=SIGN**(L+1)*QLEGY

      RETURN
      END


      SUBROUTINE FNFPI(lam2,MFM,F1S,F1V,F2S,F2V,FPI)
c  ****************** e.m. form  factors *********************+*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 mn,lam2,mg2,mup,mun,MFM
      PARAMETER (mg2 =.71D6,mup=2.7928,mun= -1.913,HQC=197.3285)

        Q2=lam2*HQC**2
        mn=MFM*HQC

	tau = Q2/(4.*mn*mn)
	F = 1./(1.+Q2/mg2)**2
	GEP = F
	GMN = mun*F
	GEN = (-tau)/(1.+4.*tau)*GMN
	GMP = mup*F
	F1P = (GEP+tau*GMP)/(1.+tau)
	F1N = (GEN+tau*GMN)/(1.+tau)
	F2P = (GMP-GEP)/((mup-1.)*(1.+tau))
	F2N = (GMN-GEN)/(mun*(1.+tau))
	F1S =F1P+F1N
	F1V =F1P-F1N
	F2S = (mup-1.)*F2P+mun*F2N
	F2V = (mup-1.)*F2P-mun*F2N
        FPI = F1V

      RETURN
      END


      SUBROUTINE HFMATR(WFM,LAM2,M,MPI,LMAX1,HM,FIJ)
c  ***************  helicity multipoles  *****************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 M,MPI,LAM2,K0,K2,KF,NU
      DIMENSION HM(6),FIJ(6,12)
      DATA PI/3.1415926536/

      W=WFM
      COF=8.*PI*W
      SIG=W**2-M**2
      BET=SIG+LAM2/2.
      K0=(SIG-LAM2)/2./W
      K2=K0**2+LAM2
      KF=DSQRT(K2)
      Q0=(SIG+MPI**2)/2./W
      Q2=Q0**2-MPI**2
      Q=DSQRT(Q2)
      NU=-(W+M)
      OM=W-M
      E1=W+M-K0
      E2=W+M-Q0
      G1=-W+M+K0
      G2=-W+M+Q0

      HM(1)=G1*DSQRT(E1*E2)/COF
      HM(2)=-E1*DSQRT(G1*G2)/COF
      HM(3)=E2*DSQRT(G1*G2)/COF
      HM(4)=G2*DSQRT(E1*E2)/COF
      HM(5)=-HM(1)
      HM(6)=-HM(2)

      DO 1 L1=1,LMAX1
      L=L1-1
      COF1=(Q*KF)**L*W
      COF2=(Q*KF)**(L-1)
      DO 1 I=1,6
      C=COF1
      IF (I.EQ.3) C=COF2
      IF (I.EQ.4) C=COF2
      FIJ(I,L1)=C*HM(I)
 1    CONTINUE

      RETURN
      END

      SUBROUTINE EL(W, Q2_MeV, m_N, m_pi, l, spin, iso, Mul)
      IMPLICIT REAL*8(A-H,O-Z)
	common/parvect/ xom1, xom2, xrho1, xrho2, xom, xrho
      INTEGER spin,D13mode
	common/newpar/JRHOFORM,JOMEFORM,D13mode
      REAL*8 Mul
      DIMENSION Q_l(1000)
      REAL*8 m_N, m_pi, m_V
      REAL*8 l_V, l_om, l_rho
      REAL*8 k, l_1, l_2, l_3
      PARAMETER (pi = 3.1415926536D0, mg2 = 0.71D6)
      PARAMETER (m_om = 783.0, m_rho = 770.0, l_om = 0.314,
     *           gV_om = 21.0, gT_om = -0.5714, cut_om = 1.2E3,
     *           l_rho = 0.103, gV_rho = 2.0, gT_rho = 6.5,
     *           cut_rho = 1.5E3, el_ch = 0.3028619)

      dip = mg2

      W_m = W - m_N
      W_p = W + m_N
      sigma = W_m * W_p
      om_gam = (sigma - Q2_MeV) / 2.0 / W
      q = SQRT(om_gam**2 + Q2_MeV)
      om_pi = (sigma + m_pi**2) / 2.0 / W
      k = SQRT(om_pi**2 - m_pi**2)
      E_i = W - om_gam
      D_i = SQRT(E_i + m_N)
      E_f = W - om_pi
      D_f = SQRT(E_f + m_N)

      IF (iso.EQ.1) THEN
      m_V = m_om
      l_V = el_ch * l_om / m_pi
      gV_V0 = gV_om
      gV_V = gV_om*xom1
      gT_V = gT_om * gV_V0 / 2.0 / m_N*xom2
      cut_V = cut_om*xom
	F_VNN=cut_V**2 / (cut_V**2 + k * k + q * q)
        IF (JOMEFORM .EQ. 0) F_VNN=1
      ELSE IF (iso.EQ.2) THEN
      m_V = m_rho
      l_V = el_ch * l_rho / m_pi
      gV_V0 = gV_rho
      gV_V = gV_rho*xrho1
      gT_V = gT_rho * gV_V0 / 2.0 / m_N*xrho2
      cut_V = cut_rho*xrho
	F_VNN=cut_V**2 / (cut_V**2 + k * k + q * q)
        IF (JRHOFORM .EQ. 0) F_VNN=1
      ENDIF
      arg_Q = (m_V**2 - m_pi**2 + 2.0 * om_gam * om_pi
     *      + Q2_MeV) / 2.0 / k / q
      CALL QLEG(arg_Q, l + 5, Q_l)

      IF (l.LT.2.AND.spin.EQ.2) GO TO 333

      IF (spin.EQ.1) THEN
         l_1 = 1.0 / (real(l) + 1)
         l_2 = real(l)
         l_3 = (real(l) + 1.0) / (2.0 * real(l) + 3.0)
         Q_2 = Q_l(l + 2)
         Q_4 = Q_l(l + 1) - Q_l(l + 3)

      ELSE IF (spin.EQ.2) THEN
         l_1 = 1.0 / real(l)
         l_2 = - (real(l) + 1.0)
         l_3 = real(l) / (2.0 * real(l) - 1.0)
         Q_2 = Q_l(l)
         Q_4 = Q_l(l + 1) - Q_l(l - 1)
      ENDIF
      delta = 0.0
      IF (l.EQ.0) delta = 1.0

      F_V = F_VNN / (1.0 + Q2_MeV / mg2)**2

c     *    * cut_V**2 / (cut_V**2 + k * k + q * q)
C ********* above the hadronic form factor appears  ******************
      Q_1 = Q_l(l + 1)
      IF (l.GT.0) THEN
         Q_3 = (Q_l(l) - Q_l(l + 2)) / (2.0 * real(l) + 1.0)
      ELSE IF (l.EQ.0) THEN
         Q_3 = 0
      END IF

      Mul = l_1 * l_V * F_V / 16.0 / pi / W * D_f / D_i
     *     * (D_i**2 / k / q * Q_1 * ((W_m**2
     *     + (m_V**2 - m_pi**2  + Q2_MeV) / 2.0) * gV_V
     *     - W_m * m_V**2 * gT_V)
     *     - Q_2 / D_f**2 * ((W_p**2
     *     + (m_V**2 - m_pi**2 + Q2_MeV) / 2.0) * gV_V
     *     + W_p * m_V**2 * gT_V)
     *     - l_2 * Q_3 * (W_p * gV_V - (W_p * W_m + Q2_MeV) * gT_V)
     *     - l_3 * Q_4 * k / q * D_i**2 / D_f**2
     *     * (W_m * gV_V + (W_p * W_m + Q2_MeV) * gT_V)
     *     + delta * (-gV_V + 2.0 * W_m * gT_V) * D_i**2)
c '-' in last line before gV_V needs further confirmation

      GO TO 334

 333  Mul = 0.0
 334  CONTINUE

      END

      SUBROUTINE ML(W, Q2_MeV, m_N, m_pi, l, spin, iso, Mul)
      IMPLICIT REAL*8(A-H,O-Z)
	common/parvect/ xom1, xom2, xrho1, xrho2, xom, xrho
      INTEGER spin,D13mode
	common/newpar/JRHOFORM,JOMEFORM,D13mode
      REAL*8 Mul
      DIMENSION Q_l(20)
      REAL*8 m_N, m_pi, m_V
      REAL*8 l_V, l_om, l_rho
      REAL*8 k, l_1
      PARAMETER (pi = 3.1415926536D0, mg2 = 0.71D6)
      PARAMETER (m_om = 783.0, m_rho = 770.0, l_om = 0.314,
     *           gV_om = 21.0, gT_om = -0.5714, cut_om = 1.2E3,
     *           l_rho = 0.103, gV_rho = 2.0, gT_rho = 6.5,
     *           cut_rho = 1.5E3, el_ch = 0.3028619)

      dip = mg2

      W_m = W - m_N
      W_p = W + m_N
      sigma = W_m * W_p
      om_gam = (sigma - Q2_MeV) / 2.0 / W
      q = SQRT(om_gam**2 + Q2_MeV)
      om_pi = (sigma + m_pi**2) / 2.0 / W
      k = SQRT(om_pi**2 - m_pi**2)
      E_i = W - om_gam
      D_i = SQRT(E_i + m_N)
      E_f = W - om_pi
      D_f = SQRT(E_f + m_N)

      IF (iso.EQ.1) THEN
      m_V = m_om
      l_V = el_ch * l_om / m_pi
      gV_V0 = gV_om
      gV_V = gV_om*xom1
      gT_V = gT_om * gV_V0 / 2.0 / m_N*xom2
      cut_V = cut_om*xom
	F_VNN=cut_V**2 / (cut_V**2 + k * k + q * q)
        IF (JOMEFORM .EQ. 0) F_VNN=1
      ELSE IF (iso.EQ.2) THEN
      m_V = m_rho
      l_V = el_ch * l_rho / m_pi
      gV_V0 = gV_rho
      gV_V = gV_rho*xrho1
      gT_V = gT_rho * gV_V0 / 2.0 / m_N*xrho2
      cut_V = cut_rho*xrho
	F_VNN=cut_V**2 / (cut_V**2 + k * k + q * q)
        IF (JRHOFORM .EQ. 0) F_VNN=1
      ENDIF

      arg_Q = (m_V**2 - m_pi**2 + 2.0 * om_gam * om_pi
     *      + Q2_MeV) / 2.0 / k / q
      CALL QLEG(arg_Q, l + 5, Q_l)

      IF (l.LT.1) GO TO 335

      delta = 0.0

      IF (spin.EQ.1) THEN
         l_1 = 1.0 / (real(l) + 1)
         Q_2 = Q_l(l + 2)

      ELSE IF (spin.EQ.2) THEN
         l_1 = - 1.0 / real(l)
         Q_2 = Q_l(l)
         IF (l.EQ.1) delta = 1.0
      ENDIF

      F_V = F_VNN / (1.0 + Q2_MeV / mg2)**2

c     *    * cut_V**2 / (cut_V**2 + k * k + q * q)
C ********* above the hadronic form factor appears  ******************
      Q_1 = Q_l(l + 1)
      Q_3 = (Q_l(l) - Q_l(l+2)) / (2.0 * real(l) + 1.0)

      Mul = l_1 * l_V * F_V / 16.0 / pi / W * D_f / D_i *
     *    (Q_1 * D_i**2 / k / q * ((W_m**2
     *  + (m_V**2 - m_pi**2 + Q2_MeV) / 2.0) * gV_V
     *  - W_m * m_V**2 * gT_V)
     *  - Q_2 / D_f**2 * ((W_p**2
     *  + (m_V**2 - m_pi**2 + Q2_MeV) / 2.0) * gV_V
     *  + W_p * m_V**2 * gT_V)
     *  + Q_3 * (W_p * gV_V - (W_m * W_p + Q2_MeV) * gT_V)
     *  + delta * k * q / D_f**2 * (gV_V + 2.0 * W_p * gT_V))

      GO TO 336

 335  Mul = 0.0

 336  CONTINUE

      END

      SUBROUTINE LL(W, Q2_MeV, m_N, m_pi, l, spin, iso, Mul)
      IMPLICIT REAL*8(A-H,O-Z)
	common/parvect/ xom1, xom2, xrho1, xrho2, xom, xrho
      INTEGER spin,D13mode
	common/newpar/JRHOFORM,JOMEFORM,D13mode
      REAL*8 Mul
      DIMENSION Q_l(20)
      REAL*8 m_N, m_pi, m_V
      REAL*8 l_V, l_om, l_rho
      REAL*8 k, l_1
      PARAMETER (pi = 3.1415926536D0, mg2 = 0.71D6)
      PARAMETER (m_om = 783.0, m_rho = 770.0, l_om = 0.314,
     *           gV_om = 21.0, gT_om = -0.5714, cut_om = 1.2E3,
     *           l_rho = 0.103, gV_rho = 2.0, gT_rho = 6.5,
     *           cut_rho = 1.5E3, el_ch = 0.3028619)


      dip = mg2

      W_m = W - m_N
      W_p = W + m_N
      sigma = W_m * W_p
      om_gam = (sigma - Q2_MeV) / 2.0 / W
      q = SQRT(om_gam**2 + Q2_MeV)
      om_pi = (sigma + m_pi**2) / 2.0 / W
      k = SQRT(om_pi**2 - m_pi**2)
      E_i = W - om_gam
      D_i = SQRT(E_i + m_N)
      E_f = W - om_pi
      D_f = SQRT(E_f + m_N)

      IF (iso.EQ.1) THEN
      m_V = m_om
      l_V = el_ch * l_om / m_pi
      gV_V0 = gV_om
      gV_V = gV_om*xom1
      gT_V = gT_om * gV_V0 / 2.0 / m_N*xom2
      cut_V = cut_om*xom
	F_VNN=cut_V**2 / (cut_V**2 + k * k + q * q)
        IF (JOMEFORM .EQ. 0) F_VNN=1
      ELSE IF (iso.EQ.2) THEN
      m_V = m_rho
      l_V = el_ch * l_rho / m_pi
      gV_V0 = gV_rho
      gV_V = gV_rho*xrho1
      gT_V = gT_rho * gV_V0 / 2.0 / m_N*xrho2
      cut_V = cut_rho*xrho
	F_VNN=cut_V**2 / (cut_V**2 + k * k + q * q)
        IF (JRHOFORM .EQ. 0) F_VNN=1
      ENDIF

      arg_Q = (m_V**2 - m_pi**2 + 2.0 * om_gam * om_pi
     *      + Q2_MeV) / 2.0 / k / q
      CALL QLEG(arg_Q, l + 5, Q_l)

      delta_1 = 0.0
      delta_2 = 0.0
      IF (l.LT.1.AND.spin.EQ.2) GO TO 337

      IF (spin.EQ.1) THEN
         l_1 = 1.0 / (real(l) + 1)
         Q_2 = Q_l(l + 2)
         Q_4 = Q_l(l + 2)
         IF (l.EQ.0) delta_1 = 1.0

      ELSE IF (spin.EQ.2) THEN
         l_1 = 1.0 / real(l)
         Q_2 = Q_l(l)
         Q_4 = Q_l(l)
         IF (l.EQ.1) delta_2 = 1.0
      ENDIF

      F_V = F_VNN / (1.0 + Q2_MeV / mg2)**2
c      F_V = 1.0 / (1.0 + Q2_MeV / dip)**2
c     *    * cut_V**2 / (cut_V**2 + k * K + q * q)
C ********* above the hadronic form factor appears  ******************
      Q_1 = Q_l(l + 1)
      Q_3 = Q_l(l + 1)

      Mul = l_1 * l_V * F_V / 16.0 / pi / W * D_f / D_i
     *    * (((-D_i**2 * (m_pi**2 - Q2_MeV - 2.0 * om_gam * om_pi)
     *    - 2.0 * E_f * (om_gam**2 + Q2_MeV))* gT_V
     *    + (E_f - m_N) * (E_i + m_N) * gV_V) * Q_1 / k
     *    + D_i**2 / D_f**2 * (((E_i - m_N)
     *    * (m_pi**2 - Q2_MeV - 2.0 * om_gam * om_pi)
     *    + 2.0 * E_f * (om_gam**2 + Q2_MeV)) * gT_V
     *    + (E_i - m_N) * (E_f + m_N) * gV_V) * Q_2 / q
     *    + (2.0 * m_N * gT_V + gV_V) * (q * delta_1
     *    + k * delta_2 * D_i**2 / D_f **2
     *    - arg_Q * (q * Q_3 + k * D_i**2 / D_f **2 * Q_4)))
     *    * om_gam / q

      GO TO 338

 337  Mul = 0.0

 338  CONTINUE

      END


      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x, y, y2
      dimension x(200),y(200),y2(200)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u
      dimension u(500)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
C        PRINT *, k, y2(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software VsXz&2+L.0(9p+.

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
       write(6,'(''bad xa input in splint'')')
       stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software VsXz&2+L.0(9p+.



      SUBROUTINE ODDSR(WFM,WFMP,LAM2,M,MPI,SIK,RIK)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 M,MPI,LAM2,K0,K2,NU,K0P,K2P,NUP
      DIMENSION SIK(3,3,3),RIK(3,3,4)

      W=WFM
      SIG=W**2-M**2
      BET=SIG+LAM2/2.
      K0=(SIG-LAM2)/2/W
      K2=K0**2+LAM2
      Q0=(SIG+MPI**2)/2./W
      Q2=Q0**2-MPI**2
      NU=-(W+M)
      OM=W-M
      E1=W+M-K0
      E2=W+M-Q0
      G1=-W+M+K0
      G2=-W+M+Q0

      WP=WFMP
      SIGP=WP**2-M**2
      BETP=SIGP+LAM2/2.
      K0P=(SIGP-LAM2)/2/WP
      K2P=K0P**2+LAM2
      Q0P=(SIGP+MPI**2)/2./WP
      Q2P=Q0P**2-MPI**2
      NUP=-(WP+M)
      OMP=WP-M
      E1P=WP+M-K0P
      E2P=WP+M-Q0P
      G1P=-WP+M+K0P
      G2P=-WP+M+Q0P

      GAM1=-LAM2*(BET+BETP)+2.*E1*LAM2*(2.*OM+NUP)
      GAM2=NU*LAM2*(OMP-2.*Q0P)+(BET+BETP)*(LAM2+2.*OMP*K0P)-
     - 4.*WP*K2P*OMP
      GAM3=2.*W*K2-(BET+BETP)*K0+E1*(OM*NUP-2.*LAM2)
      GAM4=2.*(Q0P-K0P)*(G1P*NUP+LAM2)+K0P*LAM2+2.*OM*K2P
      GAM5=2.*G1P*(Q0P-K0P)
      D1=Q0-K0
      D2=2.*Q0-K0

      SIK(1,1,1)=GAM1+4*WP*K0P*E1*OM+OMP*LAM2*(2*D1-NU)
      SIK(1,2,1)=(-4*WP*E1*(OM*GAM4-LAM2*GAM5)+GAM1*(LAM2+2*OMP*Q0P)-
     - GAM2*(LAM2+2*NU*D1)+4*E1*LAM2*W*(OMP-2.*Q0P)*D1)/2/K2P
      SIK(1,3,1)=2*OMP*GAM1-8*WP*LAM2*E1*OM-2*LAM2**2*(2*D1-NU)
      SIK(2,1,1)=2*K2*NU*OMP
      SIK(2,2,1)=K2*GAM2/K2P
      SIK(2,3,1)=-4*K2*NU*LAM2
      SIK(3,1,1)=GAM3-2*WP*K0P*E1-OMP*NU*D2
      SIK(3,2,1)=(2*WP*E1*(GAM4+OM*GAM5)+GAM3*(LAM2+2*OMP*Q0P)-
     - GAM2*D2)/2/K2P
      SIK(3,3,1)=2*OMP*GAM3+4*WP*LAM2*E1+2*NU*LAM2*D2

      SIK(1,1,2)=OMP*(E1-K0)
      SIK(1,2,2)=(GAM2-LAM2*(OMP-2*K0P)*(2*D1-NU)+4*E1*(OM*(M*K0P-
     - G1P*NUP)-LAM2*(WP+Q0P-K0P)))/2/K2P
      SIK(1,3,2)=-2*LAM2*(E1-K0)
      SIK(2,1,2)=0.
      SIK(2,2,2)=K2*NU*(2*K0P-OMP)/K2P
      SIK(2,3,2)=0.
      SIK(3,1,2)=OMP
      SIK(3,2,2)=(LAM2*(OMP-2*Q0P)-GAM3-NU*D2*(G1P+K0P)+2*E1*(K0P*(OM+
     + OMP-M)-2*WP*OM)-4*OMP*(W*K0P+WP*K0-1.5*K0*K0P))/2/K2P
      SIK(3,3,2)=-2*LAM2

      SIK(1,1,3)=0.
      SIK(1,2,3)=(G1P+K0P)*(E1-K0)/2/K2P
      SIK(1,3,3)=0.
      SIK(2,1,3)=0.
      SIK(2,2,3)=0.
      SIK(2,3,3)=0.
      SIK(3,1,3)=0.
      SIK(3,2,3)=(G1P+K0P)/2/K2P
      SIK(3,3,3)=0.

      RIK(1,1,1)=-LAM2-2*NU*D1
      RIK(1,2,1)=-(LAM2+2*OMP*Q0P)*(LAM2+2*NU*D1)/2/K2P
      RIK(1,3,1)=-2*OMP*(LAM2+2*NU*D1)
      RIK(2,1,1)=2*K2
      RIK(2,2,1)=K2*(LAM2+2*OMP*Q0P)/K2P
      RIK(2,3,1)=4*OMP*K2
      RIK(3,1,1)=-D2
      RIK(3,2,1)=-(LAM2+2*OMP*Q0P)*D2/2/K2P
      RIK(3,3,1)=-2*OMP*D2

      RIK(1,1,4)=LAM2
      RIK(1,2,4)=(LAM2+2*OMP*Q0P)*LAM2/2/K2P
      RIK(1,3,4)=2*OMP*LAM2
      RIK(2,1,4)=0.
      RIK(2,2,4)=0.
      RIK(2,3,4)=0.
      RIK(3,1,4)=K0
      RIK(3,2,4)=(LAM2+2*OMP*Q0P)*K0/2/K2P
      RIK(3,3,4)=2*K0*OMP

      RIK(1,1,2)=LAM2+2*OM*E1
      RIK(1,2,2)=LAM2*(G1P*(NU+2*Q0)-E1*(OMP-2*Q0P)-
     - NU*D1+OMP*Q0P)/K2P
      RIK(1,3,2)=2*LAM2*(OMP+2*E1)
      RIK(2,1,2)=0.
      RIK(2,2,2)=-K2*(LAM2+2*OMP*K0P)/K2P
      RIK(2,3,2)=0.
      RIK(3,1,2)=NU+2*K0
      RIK(3,2,2)=((NU+2*K0)*(LAM2+2*OMP*Q0P)+2*G1P*(WP+W)*
     * (NU+2*Q0)-(LAM2+2*OMP*G1P)*D2)/2./K2P
      RIK(3,3,2)=2*OMP*K0+2*OM*E1

      RIK(1,1,3)=0.
      RIK(1,2,3)=(LAM2+2*G1P*NU-2*E1*OMP)/2./K2P
      RIK(1,3,3)=0.
      RIK(2,1,3)=0.
      RIK(2,2,3)=0.
      RIK(2,3,3)=0.
      RIK(3,1,3)=0.
      RIK(3,2,3)=(NU+2*K0+2*G1P)/2./K2P
      RIK(3,3,3)=0.

      DO 1 I=1,3
      DO 1 K=1,3
      RIK(I,K,1)=RIK(I,K,1)*LAM2/8/WP/WP/K2
      RIK(I,K,4)=RIK(I,K,4)/8/WP/WP/K2
      SIK(I,K,1)=SIK(I,K,1)/8/WP/WP/K2
      DO 1 J=2,3
      RIK(I,K,J)=RIK(I,K,J)/8/WP/WP/K2
 1    SIK(I,K,J)=SIK(I,K,J)/8/WP/WP/K2

      RETURN
      END

      SUBROUTINE VIJL(LMAX1,X,VIJ)
c    V_{ij} matrix defined by Eq. (20)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VIJ(6,6,10),PLD(10),PLDD(10)

      DO 1 L1=1,10
      DO 1 I=1,6
      DO 1 J=1,6
  1   VIJ(I,J,L1)=0.D0

      CALL PLEGD(X,LMAX1+1,PLD)
      CALL PLEGDD(X,LMAX1+1,PLDD)

      DO 2 L1=1,LMAX1
      VIJ(1,1,L1)=PLD(L1+1)
      VIJ(1,2,L1)=-PLD(L1)
      VIJ(2,2,L1)=PLD(L1+1)
      VIJ(2,1,L1)=-PLD(L1)
      VIJ(5,5,L1)=PLD(L1+1)
      VIJ(5,6,L1)=-PLD(L1)
      VIJ(6,6,L1)=PLD(L1+1)
      VIJ(6,5,L1)=-PLD(L1)
      VIJ(3,3,L1)=PLDD(L1+1)
      VIJ(3,4,L1)=-PLDD(L1)
      VIJ(4,4,L1)=PLDD(L1+1)
      VIJ(4,3,L1)=-PLDD(L1)
  2   CONTINUE

      RETURN
      END

      SUBROUTINE WIJL(LMAX1,X,WIJ)
c  W_{ij) matrix defined by Eq. (28)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WIJ(6,6,10),PL(10),RL(10)

      DO 1 L1=1,10
      DO 1 I=1,6
      DO 1 J=1,6
  1   WIJ(I,J,L1)=0.D0

      LMAX12=LMAX1+2
      CALL PLEG(X,LMAX12,PL)

      RL(1)=0.
      LMAX11=LMAX1+1
      DO 22 L1=2,LMAX11
 22   RL(L1)=(PL(L1-1)-PL(L1+1))/(2*L1-1)

      DL=0.
      DO 2 L1=1,LMAX1
      IF (L1.GT.1) DL=1.D0
      WIJ(1,1,L1)=PL(L1)
      WIJ(1,2,L1)=PL(L1+1)
      WIJ(2,2,L1)=PL(L1)
      WIJ(2,1,L1)=PL(L1+1)
      WIJ(5,5,L1)=PL(L1)
      WIJ(5,6,L1)=PL(L1+1)
      WIJ(6,6,L1)=PL(L1)
      WIJ(6,5,L1)=PL(L1+1)
      WIJ(3,3,L1)=DL*RL(L1)
      WIJ(3,4,L1)=DL*RL(L1+1)
      WIJ(4,4,L1)=DL*RL(L1)
      WIJ(4,3,L1)=DL*RL(L1+1)
  2   CONTINUE

      RETURN
      END

      SUBROUTINE PLEG(X,LMAX1,PL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PL(10)

      X2=X*X
      PL(1)=1.D0
      PL(2)=X
      PL(3)=(3*X2-1.)/2.
      PL(4)=X*(5*X2-3)/2.
      PL(5)=(35*X2*X2-30*X2+3)/8.

      IF (LMAX1.LE.5) RETURN

      DO 1 L1=5,LMAX1
      L=L1-1
 1    PL(L1+1)=((2*L+1)*X*PL(L1)-L*PL(L1-1))/(L+1)

      RETURN
      END

      SUBROUTINE PLEGD(X,LMAX1,PLD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PLD(10),PL(10)

      X2=X*X
      PLD(1)=0.D0
      PLD(2)=1.D0
      PLD(3)=3*X
      PLD(4)=(15*X2-3.)/2.
      PLD(5)=X*(140*X2-60.)/8.

      IF (LMAX1.LE.5) RETURN

      CALL PLEG(X,LMAX1,PL)

      DO 1 L1=5,LMAX1
      L=L1-1
 1    PLD(L1+1)=((2*L+1)*(PL(L1)+X*PLD(L1))-L*PLD(L1-1))/(L+1)
      RETURN
      END


      SUBROUTINE PLEGDD(X,LMAX1,PLDD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PLDD(10),PLD(10)

      X2=X*X
      PLDD(1)=0.D0
      PLDD(2)=0.D0
      PLDD(3)=3
      PLDD(4)=15*X
      PLDD(5)=(105*X2-15.)/2.

      IF (LMAX1.LE.5) RETURN

      CALL PLEGD(X,LMAX1,PLD)

      DO 1 L1=5,LMAX1
      L=L1-1
 1    PLDD(L1+1)=((2*L+1)*(2*PLD(L1)+X*PLDD(L1))-L*PLDD(L1-1))/(L+1)

      RETURN
      END

      SUBROUTINE GSET(AX,BX,NX,Z,W)
      IMPLICIT REAL*8(A-H,O-Z)
C     N-POINT GAUSS ZEROS AND WEIGHTS FOR THE INTERVAL (AX,BX) ARE
C           STORED IN  ARRAYS Z AND W RESPECTIVELY.
C
      COMMON /GQCOM/A(273),X(273),KTAB(96)
      DIMENSION Z(200),W(200)
C
      DATA IBD/0/
C     IF(IBD.EQ.1)CALL D106BD
C     THIS IS A CALL TO A SUBROUTINE WHICH ONLY GIVES VALUES TO THE
C     ARRAYS IN THE COMMON BLOCK GQCOM BY DATA STATEMENTS
C     THE SUBROUTINE D106BD ONLY HAS TO BE LOADED, NOT EXECUTED
C     IF SUBROUTINE D106BD IS CHANGED TO BLOCK DATA THIS CALL SHOULD
C     BE OMITTED
C
C-----TEST N
      N=NX
      ALPHA=0.5*(BX+AX)
      BETA=0.5*(BX-AX)
      IF(N.LT.1) GO TO 100
      IF(N.NE.1) GO TO 1
      Z(1)=ALPHA
      W(1)=BX-AX
      RETURN
C
    1 IF(N.LE.16) GO TO 2
      IF(N.EQ.20) GO TO 2
      IF(N.EQ.24) GO TO 2
      IF(N.EQ.32) GO TO 2
      IF(N.EQ.40) GO TO 2

      GO TO 100
C
C----- SET K EQUAL TO INITIAL SUBSCRIPT AND STORE RESULTS
    2 K=KTAB(N)
      M=N/2
C
      DO 3 J=1,M
      JTAB=K-1+J
      WTEMP=BETA*A(JTAB)
      DELTA=BETA*X(JTAB)
      Z(J)=ALPHA-DELTA
      W(J)=WTEMP
      JP=N+1-J
      Z(JP)=ALPHA+DELTA
      W(JP)=WTEMP
    3 CONTINUE
C
      IF((N-M-M).EQ.0) RETURN
      Z(M+1)=ALPHA
      JMID=K+M
      W(M+1)=BETA*A(JMID)
      RETURN
C
  100 ZN=N
      PRINT 200,ZN
      RETURN
C
  200 FORMAT(  41H GSET ... N HAS THE NON-PERMISSIBLE VALUE E11.3 )
      END

      BLOCK DATA
      IMPLICIT REAL*8(A-H,O-Z)
C     SUBROUTINE D106BD
C     THIS SUBROUTINE GIVES VALUES TO THE ARRAYS IN THE COMMON BLOCK
C     GQCOM AND COULD BE CHANGED TO A BLOCK DATA
C     THE CALLS IN GQUAD AND GSET SHOULD THEN BE REMOVED
C     (DATA BLOCK FOR GQUAD AND GSET)
      COMMON /GQCOM/A(273),X(273),KTAB(96)
C
C-----TABLE OF INITIAL SUBSCRIPTS FOR N=2(1)16(4)96
      DATA KTAB(2)/1/
      DATA KTAB(3)/2/
      DATA KTAB(4)/4/
      DATA KTAB(5)/6/
      DATA KTAB(6)/9/
      DATA KTAB(7)/12/
      DATA KTAB(8)/16/
      DATA KTAB(9)/20/
      DATA KTAB(10)/25/
      DATA KTAB(11)/30/
      DATA KTAB(12)/36/
      DATA KTAB(13)/42/
      DATA KTAB(14)/49/
      DATA KTAB(15)/56/
      DATA KTAB(16)/64/
      DATA KTAB(20)/72/
      DATA KTAB(24)/82/
      DATA KTAB(28)/82/
      DATA KTAB(32)/94/
      DATA KTAB(36)/94/
      DATA KTAB(40)/110/
C
C-----TABLE OF ABSCISSAE (X) AND WEIGHTS (A) FOR INTERVAL (-1,+1).
C
C-----N=2
      DATA X(1)/0.577350269189626  /, A(1)/1.000000000000000  /
C-----N=3
      DATA X(2)/0.774596669241483  /, A(2)/0.555555555555556  /
      DATA X(3)/0.000000000000000  /, A(3)/0.888888888888889  /
C-----N=4
      DATA X(4)/0.861136311594053  /, A(4)/0.347854845137454  /
      DATA X(5)/0.339981043584856  /, A(5)/0.652145154862546  /
C-----N=5
      DATA X(6)/0.906179845938664  /, A(6)/0.236926885056189  /
      DATA X(7)/0.538469310105683  /, A(7)/0.478628670499366  /
      DATA X(8)/0.000000000000000  /, A(8)/0.568888888888889  /
C-----N=6
      DATA X(9)/0.932469514203152  /, A(9)/0.171324492379170  /
      DATA X(10)/0.661209386466265 /, A(10)/0.360761573048139 /
      DATA X(11)/0.238619186083197 /, A(11)/0.467913934572691 /
C-----N=7
      DATA X(12)/0.949107912342759 /, A(12)/0.129484966168870 /
      DATA X(13)/0.741531185599394 /, A(13)/0.279705391489277 /
      DATA X(14)/0.405845151377397 /, A(14)/0.381830050505119 /
      DATA X(15)/0.000000000000000 /, A(15)/0.417959183673469 /
C-----N=8
      DATA X(16)/0.960289856497536 /, A(16)/0.101228536290376 /
      DATA X(17)/0.796666477413627 /, A(17)/0.222381034453374 /
      DATA X(18)/0.525532409916329 /, A(18)/0.313706645877887 /
      DATA X(19)/0.183434642495650 /, A(19)/0.362683783378362 /
C-----N=9
      DATA X(20)/0.968160239507626 /, A(20)/0.081274388361574 /
      DATA X(21)/0.836031107326636 /, A(21)/0.180648160694857 /
      DATA X(22)/0.613371432700590 /, A(22)/0.260610696402935 /
      DATA X(23)/0.324253423403809 /, A(23)/0.312347077040003 /
      DATA X(24)/0.000000000000000 /, A(24)/0.330239355001260 /
C-----N=10
      DATA X(25)/0.973906528517172 /, A(25)/0.066671344308688 /
      DATA X(26)/0.865063366688985 /, A(26)/0.149451349150581 /
      DATA X(27)/0.679409568299024 /, A(27)/0.219086362515982 /
      DATA X(28)/0.433395394129247 /, A(28)/0.269266719309996 /
      DATA X(29)/0.148874338981631 /, A(29)/0.295524224714753 /
C-----N=11
      DATA X(30)/0.978228658146057 /, A(30)/0.055668567116174 /
      DATA X(31)/0.887062599768095 /, A(31)/0.125580369464905 /
      DATA X(32)/0.730152005574049 /, A(32)/0.186290210927734 /
      DATA X(33)/0.519096129206812 /, A(33)/0.233193764591990 /
      DATA X(34)/0.269543155952345 /, A(34)/0.262804544510247 /
      DATA X(35)/0.000000000000000 /, A(35)/0.272925086777901 /
C-----N=12
      DATA X(36)/0.981560634246719 /, A(36)/0.047175336386512 /
      DATA X(37)/0.904117256370475 /, A(37)/0.106939325995318 /
      DATA X(38)/0.769902674194305 /, A(38)/0.160078328543346 /
      DATA X(39)/0.587317954286617 /, A(39)/0.203167426723066 /
      DATA X(40)/0.367831498998180 /, A(40)/0.233492536538355 /
      DATA X(41)/0.125233408511469 /, A(41)/0.249147045813403 /
C-----N=13
      DATA X(42)/0.984183054718588 /, A(42)/0.040484004765316 /
      DATA X(43)/0.917598399222978 /, A(43)/0.092121499837728 /
      DATA X(44)/0.801578090733310 /, A(44)/0.138873510219787 /
      DATA X(45)/0.642349339440340 /, A(45)/0.178145980761946 /
      DATA X(46)/0.448492751036447 /, A(46)/0.207816047536889 /
      DATA X(47)/0.230458315955135 /, A(47)/0.226283180262897 /
      DATA X(48)/0.000000000000000 /, A(48)/0.232551553230874 /
C-----N=14
      DATA X(49)/0.986283808696812 /, A(49)/0.035119460331752 /
      DATA X(50)/0.928434883663574 /, A(50)/0.080158087159760 /
      DATA X(51)/0.827201315069765 /, A(51)/0.121518570687903 /
      DATA X(52)/0.687292904811685 /, A(52)/0.157203167158194 /
      DATA X(53)/0.515248636358154 /, A(53)/0.185538397477938 /
      DATA X(54)/0.319112368927890 /, A(54)/0.205198463721296 /
      DATA X(55)/0.108054948707344 /, A(55)/0.215263853463158 /
C-----N=15
      DATA X(56)/0.987992518020485 /, A(56)/0.030753241996117 /
      DATA X(57)/0.937273392400706 /, A(57)/0.070366047488108 /
      DATA X(58)/0.848206583410427 /, A(58)/0.107159220467172 /
      DATA X(59)/0.724417731360170 /, A(59)/0.139570677926154 /
      DATA X(60)/0.570972172608539 /, A(60)/0.166269205816994 /
      DATA X(61)/0.394151347077563 /, A(61)/0.186161000015562 /
      DATA X(62)/0.201194093997435 /, A(62)/0.198431485327111 /
      DATA X(63)/0.000000000000000 /, A(63)/0.202578241925561 /
C-----N=16
      DATA X(64)/0.989400934991650 /, A(64)/0.027152459411754 /
      DATA X(65)/0.944575023073233 /, A(65)/0.062253523938648 /
      DATA X(66)/0.865631202387832 /, A(66)/0.095158511682493 /
      DATA X(67)/0.755404408355003 /, A(67)/0.124628971255534 /
      DATA X(68)/0.617876244402644 /, A(68)/0.149595988816577 /
      DATA X(69)/0.458016777657227 /, A(69)/0.169156519395003 /
      DATA X(70)/0.281603550779259 /, A(70)/0.182603415044924 /
      DATA X(71)/0.095012509837637 /, A(71)/0.189450610455069 /
C-----N=20
      DATA X(72)/0.993128599185094 /, A(72)/0.017614007139152 /
      DATA X(73)/0.963971927277913 /, A(73)/0.040601429800386 /
      DATA X(74)/0.912234428251325 /, A(74)/0.062672048334109 /
      DATA X(75)/0.839116971822218 /, A(75)/0.083276741576704 /
      DATA X(76)/0.746331906460150 /, A(76)/0.101930119817240 /
      DATA X(77)/0.636053680726515 /, A(77)/0.118194531961518 /
      DATA X(78)/0.510867001950827 /, A(78)/0.131688638449176 /
      DATA X(79)/0.373706088715419 /, A(79)/0.142096109318382 /
      DATA X(80)/0.227785851141645 /, A(80)/0.149172986472603 /
      DATA X(81)/0.076526521133497 /, A(81)/0.152753387130725 /
C-----N=24
      DATA X(82)/0.995187219997021 /, A(82)/0.012341229799987 /
      DATA X(83)/0.974728555971309 /, A(83)/0.028531388628933 /
      DATA X(84)/0.938274552002732 /, A(84)/0.044277438817419 /
      DATA X(85)/0.886415527004401 /, A(85)/0.059298584915436 /
      DATA X(86)/0.820001985973902 /, A(86)/0.073346481411080 /
      DATA X(87)/0.740124191578554 /, A(87)/0.086190161531953 /
      DATA X(88)/0.648093651936975 /, A(88)/0.097618652104113 /
      DATA X(89)/0.545421471388839 /, A(89)/0.107444270115965 /
      DATA X(90)/0.433793507626045 /, A(90)/0.115505668053725 /
      DATA X(91)/0.315042679696163 /, A(91)/0.121670472927803 /
      DATA X(92)/0.191118867473616 /, A(92)/0.125837456346828 /
      DATA X(93)/0.064056892862605 /, A(93)/0.127938195346752 /
C-----N=32
      DATA X(94)/0.997263861849481 /, A(94)/0.007018610009470 /
      DATA X(95)/0.985611511545268 /, A(95)/0.016274394730905 /
      DATA X(96)/0.964762255587506 /, A(96)/0.025392065309262 /
      DATA X(97)/0.934906075937739 /, A(97)/0.034273862913021 /
      DATA X(98)/0.896321155766052 /, A(98)/0.042835898022226 /
      DATA X(99)/0.849367613732569 /, A(99)/0.050998059262376 /
      DATA X(100)/0.794483795967942/, A(100)/0.058684093478535/
      DATA X(101)/0.732182118740289/, A(101)/0.065822222776361/
      DATA X(102)/0.663044266930215/, A(102)/0.072345794108848/
      DATA X(103)/0.587715757240762/, A(103)/0.078193895787070/
      DATA X(104)/0.506899908932229/, A(104)/0.083311924226946/
      DATA X(105)/0.421351276130635/, A(105)/0.087652093004403/
      DATA X(106)/0.331868602282127/, A(106)/0.091173878695763/
      DATA X(107)/0.239287362252137/, A(107)/0.093844399080804/
      DATA X(108)/0.144471961582796/, A(108)/0.095638720079274/
      DATA X(109)/0.048307665687738/, A(109)/0.096540088514727/
C-----N=40
      DATA X(110)/0.998237709710559/, A(110)/0.004521277098533/
      DATA X(111)/0.990726238699457/, A(111)/0.010498284531152/
      DATA X(112)/0.977259949983770/, A(112)/0.016421058381907/
      DATA X(113)/0.957916819213791/, A(113)/0.022245849194166/
      DATA X(114)/0.932812808278676/, A(114)/0.027937006980023/
      DATA X(115)/0.902098806968874/, A(115)/0.033460195282547/
      DATA X(116)/0.865959503212259/, A(116)/0.038782167974472/
      DATA X(117)/0.824612230833311/, A(117)/0.043870908185673/
      DATA X(118)/0.778305651426519/, A(118)/0.048695807635072/
      DATA X(119)/0.727318255189927/, A(119)/0.053227846983936/
      DATA X(120)/0.671956684614179/, A(120)/0.057439769099391/
      DATA X(121)/0.612553889667980/, A(121)/0.061306242492928/
      DATA X(122)/0.549467125095128/, A(122)/0.064804013456601/
      DATA X(123)/0.483075801686178/, A(123)/0.067912045815233/
      DATA X(124)/0.413779204371605/, A(124)/0.070611647391286/
      DATA X(125)/0.341994090825758/, A(125)/0.072886582395804/
      DATA X(126)/0.268152185007253/, A(126)/0.074723169057968/
      DATA X(127)/0.192697580701371/, A(127)/0.076110361900626/
      DATA X(128)/0.116084070675255/, A(128)/0.077039818164247/
      DATA X(129)/0.038772417506050/, A(129)/0.077505947978424/

C     RETURN
      END

      SUBROUTINE BACKGR(W,Q2,MM,MMPI,l_max,IUNI)
c *****************************************************************************
c     calculates unitarization phases and unitarized background multipoles
c     l=0, 1, 2, 3
c     1+iT = (PiN + 1)/2   or  (piN-1)/2 when non-unitarized bg is subtracted
c     S11, P11 prepared for alternative unitarization (commented)
c     no background unitarization for D13, D33 and F15
c *****************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MULTBR/MULT(6,10,3), FBV(6,10,3)
      REAL*8 KAPPAP,KAPPAN,MPI,MPI0,MPIP,MP,MN,M,M1,M2,M1P,M1M
      REAL*8 m_N, m_pi, MM, MMPI, m_pip
      REAL*8 MV_p, MV_m, LV_p, LV_m
      COMPLEX*16 piN, tmp, phase
      COMPLEX*16 MULT, tmp11, tmp31
      COMMON/COUPL/ IPVPS
      DIMENSION FB(6,10,3), FV(6, 10, 2)
      DIMENSION W_sp(143)
      DIMENSION del_s11(143),del_s11_ss(143),del_s31(143),
     & del_s31_ss(143)
      DIMENSION del_p11(143),del_p11_ss(143),del_p31(143),
     & del_p31_ss(143)
      DIMENSION del_p13(143),del_p13_ss(143),del_p33(143),
     & del_p33_ss(143)
      DIMENSION del_d13(143),del_d13_ss(143),del_d33(143),
     & del_d33_ss(143)
      DIMENSION del_d15(143),del_d15_ss(143),del_d35(143),
     & del_d35_ss(143)
      DIMENSION del_f15(143),del_f15_ss(143),del_f35(143),
     & del_f35_ss(143)
      DIMENSION del_f17(143),del_f17_ss(143),del_f37(143),
     & del_f37_ss(143)
      DIMENSION eta_s11(143),eta_s11_ss(143),eta_s31(143),
     & eta_s31_ss(143)
      DIMENSION eta_p11(143),eta_p11_ss(143),eta_p31(143),
     & eta_p31_ss(143)
      DIMENSION eta_p13(143),eta_p13_ss(143),eta_p33(143),
     & eta_p33_ss(143)
      DIMENSION eta_d13(143),eta_d13_ss(143),eta_d33(143),
     & eta_d33_ss(143)
      DIMENSION eta_d15(143),eta_d15_ss(143),eta_d35(143),
     & eta_d35_ss(143)
      DIMENSION eta_f15(143),eta_f15_ss(143),eta_f35(143),
     & eta_f35_ss(143)
      DIMENSION eta_f17(143),eta_f17_ss(143),eta_f37(143),
     & eta_f37_ss(143)
      COMMON /SPL/ W_sp,del_s11,del_s11_ss,del_s31,del_s31_ss,
     *      del_p11,del_p11_ss,del_p31,del_p31_ss,
     *      del_p13,del_p13_ss,del_p33,del_p33_ss,
     *      del_d13,del_d13_ss,del_d33,del_d33_ss,
     *      del_d15,del_d15_ss,del_d35,del_d35_ss,
     *      del_f15,del_f15_ss,del_f35,del_f35_ss,
     *      del_f17,del_f17_ss,del_f37,del_f37_ss,
     *      eta_s11,eta_s11_ss,eta_s31,eta_s31_ss,
     *      eta_p11,eta_p11_ss,eta_p31,eta_p31_ss,
     *      eta_p13,eta_p13_ss,eta_p33,eta_p33_ss,
     *      eta_d13,eta_d13_ss,eta_d33,eta_d33_ss,
     *      eta_d15,eta_d15_ss,eta_d35,eta_d35_ss,
     *      eta_f15,eta_f15_ss,eta_f35,eta_f35_ss,
     *      eta_f17,eta_f17_ss,eta_f37,eta_f37_ss
      COMMON /PHASES/ phi_s11, phi_s31, phi_p11, phi_p31,
     *                phi_p13, phi_p33, phi_d13, phi_d33,
     *                phi_d15, phi_d35, phi_f15, phi_f35,
     *                phi_f17, phi_f37
      COMMON/bl/ PI,PIH,HQC,F2PI,E,KAPPAP,KAPPAN,CSUC,CT,COME,CRHO,
     & CBORN,CVEC,ALAMBDA,MPI,MPI0,MPIP,MP,MN,AVM,M,M1,M2,IPV,IFORM,
     & CP33,CP11,CD13,CD33,CS11F,CS11S,CF15,CS11T,CP11S

      HQC = 197.3285
      W_MeV = HQC * W
      Q2_MeV = HQC * HQC * Q2
      m_N = HQC * MM
      m_pi = HQC * MMPI
      m_pip=139.5685
      CNORM=m_pi/m_pip

      DO 99 I=1,6
      DO 99 J=1,l_max
      DO 99 K=1,3
      FB(I,J,K)=(0.D0,0.D0)
      FBV(I,J,K)=0.D0
99    MULT(I,J,K)=(0.D0,0.D0)


c ******************************************************

      IF (CBORN.EQ.1) THEN
c ****************************  Born multipoles  ******************************
         CALL BORN_MUL(IPVPS,W,Q2,MM,MMPI,l_max,FB)

         DO 1234 i = 1, 6
            DO 1235 j = 1, l_max
               DO 1236 k = 1, 3
                  FB(i, j, k) = FB(i, j, k) / 1000.0 /m_pip
 1236          CONTINUE
 1235       CONTINUE
 1234    CONTINUE
      ELSE
         DO 1237 i = 1, 6
            DO 1238 j = 1, l_max
               DO 1239 k = 1, 3
                  FB(i, j, k) = 0.0
 1239          CONTINUE
 1238       CONTINUE
 1237    CONTINUE
      END IF
c ************************  vector mesons: omega and rho  *********************
      IF (COME.EQ.1) THEN
         DO 2100 l = 1, l_max + 1
            CALL EL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 1, 1, EV_p)
            CALL EL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 2, 1, EV_m)
            CALL ML(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 1, 1, MV_p)
            CALL ML(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 2, 1, MV_m)
            CALL LL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 1, 1, LV_p)
            CALL LL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 2, 1, LV_m)
            FV(1, l, 1) = EV_p*CNORM
            FV(2, l, 1) = EV_m*CNORM
            FV(3, l, 1) = MV_p*CNORM
            FV(4, l, 1) = MV_m*CNORM
            FV(5, l, 1) = LV_p*CNORM
            FV(6, l, 1) = LV_m*CNORM
 2100    CONTINUE
      ELSE
         DO 2300 l = 1, l_max + 1
            FV(1, l, 1) = 0.0
            FV(2, l, 1) = 0.0
            FV(3, l, 1) = 0.0
            FV(4, l, 1) = 0.0
            FV(5, l, 1) = 0.0
            FV(6, l, 1) = 0.0
 2300    CONTINUE
      END IF

      IF (CRHO.EQ.1) THEN
         DO 2200 l = 1, l_max + 1
            CALL EL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 1, 2, EV_p)
            CALL EL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 2, 2, EV_m)
            CALL ML(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 1, 2, MV_p)
            CALL ML(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 2, 2, MV_m)
            CALL LL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 1, 2, LV_p)
            CALL LL(W_MeV, Q2_MeV, m_N, m_pi, l - 1, 2, 2, LV_m)
            FV(1, l, 2) = EV_p*CNORM
            FV(2, l, 2) = EV_m*CNORM
            FV(3, l, 2) = MV_p*CNORM
            FV(4, l, 2) = MV_m*CNORM
            FV(5, l, 2) = LV_p*CNORM
            FV(6, l, 2) = LV_m*CNORM
 2200    CONTINUE
      ELSE
         DO 2400 l = 1, l_max + 1
            FV(1, l, 2) = 0.0
            FV(2, l, 2) = 0.0
            FV(3, l, 2) = 0.0
            FV(4, l, 2) = 0.0
            FV(5, l, 2) = 0.0
            FV(6, l, 2) = 0.0
 2400    CONTINUE
      END IF

C -----------------------------------  S11  -----------------------------------
c ********  LET correction, 22.10.2000 ***********
      EplB = sqrt(2.)*(FB(1,1,2)+(FB(1,1,1)-FB(1,1,3))/3)
      EplV = sqrt(2.)*FV(1,1,2)
      Epl = EplB + EplV
      SplB = sqrt(2.)*(FB(5,1,2)+(FB(5,1,1)-FB(5,1,3))/3)
      SplV = sqrt(2.)*FV(5,1,2)
      Spl = SplB + SplV

      CALL LET_CORR(W,Q2,MM,MMPI,Epl,Spl,E0_corr,S0_corr,xpi)
C ************************************************
       FCORR=0.
       Q2GEV=Q2_MeV/1.D6
       Fdip = 1./(1.+Q2GEV/0.71)**2       
       if (W_MEV.GT.1690.) 
     * FCORR=0.0006*(W_MEV-1690.)/W_MEV*dexp(-0.8*(W_MEV/1690.-1.))*Fdip
c *****************************************************************     
      CALL SPLINT(W_sp, del_s11, del_s11_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_s11, eta_s11_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1090.) del=3.83*PI*xpi/180.
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0
      tmp11=tmp
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
      phi_s11 = DATAN((1.0 - eta * COS(2.0 * del))
     *        / eta / SIN(2.0 * del))
c      phi_s11=phi_tmp
      IF (phi_s11.LT.0.0) phi_s11 = phi_s11 + pi
      del_phi = phi_s11 - phi_tmp
      phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)

      FBV(1, 1, 1) = (FB(1, 1, 1) + FV(1, 1, 1))
      MULT(1, 1, 1) = (FB(1, 1, 1) + FV(1, 1, 1)) *tmp + E0_corr*tmp
     * + FCORR * tmp

      FBV(1, 1, 2) = (FB(1, 1, 2) + FV(1, 1, 2))
      MULT(1, 1, 2) = (FB(1, 1, 2) + FV(1, 1, 2)) * tmp
     *  - FCORR/3. *tmp
     
      FBV(5, 1, 1) = (FB(5, 1, 1) + FV(5, 1, 1))
      MULT(5, 1, 1) = (FB(5, 1, 1) + FV(5, 1, 1))*tmp  + S0_corr*tmp

      FBV(5, 1, 2) = (FB(5, 1, 2) + FV(5, 1, 2))
      MULT(5, 1, 2) = (FB(5, 1, 2) + FV(5, 1, 2)) * tmp
C
C -----------------------------------  S31  -----------------------------------
      CALL SPLINT(W_sp, del_s31, del_s31_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_s31, eta_s31_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1090.) del=-2.08*PI*xpi/180.
      IF (W_MeV.LT.1079.) del=-0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
      phi_s31 = DATAN((1.0 - eta * COS(2.0 * del))
     *        / eta / SIN(2.0 * del))
c      phi_s31=phi_tmp
      IF (phi_s31.LT.0.0) phi_s31 = phi_s31 + pi
      del_phi = phi_s31 - phi_tmp
      phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)

      FBV(1, 1, 3) = (FB(1, 1, 3) + FV(1, 1, 1))
      MULT(1, 1, 3) = (FB(1, 1, 3) + FV(1, 1, 1))*tmp + E0_corr*tmp

      FBV(5, 1, 3) = (FB(5, 1, 3) + FV(5, 1, 1))
      MULT(5, 1, 3) = (FB(5, 1, 3) + FV(5, 1, 1))*tmp + S0_corr*tmp
C
C -----------------------------------  P11  -----------------------------------
      CALL SPLINT(W_sp, del_p11, del_p11_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_p11, eta_p11_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
      phi_p11 = DATAN((1.0 - eta * COS(2.0 * del))
     *        / eta / SIN(2.0 * del))
      IF (W_MeV.GT.1400.AND.phi_p11.LT.0.0) phi_p11 = phi_p11 + pi
c      IF (phi_p11.LT.0.0) phi_p11 = phi_p11 + pi
      del_phi = phi_p11 - phi_tmp
      phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)

      FBV(4, 2, 1) = (FB(4, 2, 1) + FV(4, 2, 1))
      MULT(4, 2, 1) = (FB(4, 2, 1) + FV(4, 2, 1)) * tmp

      FBV(4, 2, 2) = (FB(4, 2, 2) + FV(4, 2, 2))
      MULT(4, 2, 2) = (FB(4, 2, 2) + FV(4, 2, 2)) * tmp

      FBV(6, 2, 1) = (FB(6, 2, 1) + FV(6, 2, 1))
      MULT(6, 2, 1) = (FB(6, 2, 1) + FV(6, 2, 1)) * tmp

      FBV(6, 2, 2) = (FB(6, 2, 2) + FV(6, 2, 2))
      MULT(6, 2, 2) = (FB(6, 2, 2) + FV(6, 2, 2)) * tmp
C
C -----------------------------------  P31  -----------------------------------
      CALL SPLINT(W_sp, del_p31, del_p31_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_p31, eta_p31_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0
c ******** correction 30.07.2005 ****************
      FCORR=1.
      IF (W_MEV. GT.1600.) 
     * FCORR = dexp(-4.*(W_MEV/1600.-1.))
c ***************************************************     
      FBV(4, 2, 3) = (FB(4, 2, 3) + FV(4, 2, 1))
      MULT(4, 2, 3) = (FB(4, 2, 3) + FV(4, 2, 1)) * tmp * FCORR

      FBV(6, 2, 3) = (FB(6, 2, 3) + FV(6, 2, 1))
      MULT(6, 2, 3) = (FB(6, 2, 3) + FV(6, 2, 1)) * tmp
C
C -----------------------------------  P13  -----------------------------------
C ************************************************
       FCORR=0.
       Q2GEV=Q2_MeV/1.D6
       Fdip = 1./(1.+Q2GEV/0.71)**2       
       if (W_MEV.GT.1690.) 
     * FCORR=-0.00007*(W_MEV-1690.)/W_MEV*dexp(-1*(W_MEV/1690.-1.))
     * *Fdip
c *****************************************************************     
      CALL SPLINT(W_sp, del_p13, del_p13_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_p13, eta_p13_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0

      FBV(1, 2, 1) = (FB(1, 2, 1) + FV(1, 2, 1))
      MULT(1, 2, 1) = (FB(1, 2, 1) + FV(1, 2, 1)) * tmp

      FBV(1, 2, 2) = (FB(1, 2, 2) + FV(1, 2, 2))
      MULT(1, 2, 2) = (FB(1, 2, 2) + FV(1, 2, 2)) * tmp

      FBV(3, 2, 1) = (FB(3, 2, 1) + FV(3, 2, 1))
      MULT(3, 2, 1) = (FB(3, 2, 1) + FV(3, 2, 1)) * tmp + FCORR*tmp

      FBV(3, 2, 2) = (FB(3, 2, 2) + FV(3, 2, 2))
      MULT(3, 2, 2) = (FB(3, 2, 2) + FV(3, 2, 2)) * tmp - FCORR/3*tmp

      FBV(5, 2, 1) = (FB(5, 2, 1) + FV(5, 2, 1))
      MULT(5, 2, 1) = (FB(5, 2, 1) + FV(5, 2, 1)) * tmp

      FBV(5, 2, 2) = (FB(5, 2, 2) + FV(5, 2, 2))
      MULT(5, 2, 2) = (FB(5, 2, 2) + FV(5, 2, 2)) * tmp

C
C -----------------------------------  P33  -----------------------------------
      CALL SPLINT(W_sp, del_p33, del_p33_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_p33, eta_p33_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_p33 = DATAN(ARG)
      IF (phi_p33.LT.0.0) phi_p33 = phi_p33 + pi
      IF (W_MeV.GT.1400.AND.phi_p33.LT.PI/2.0) phi_p33 = phi_p33 + pi
      IF (W_MeV.LT.1300.) phi_p33=del
c ************  correction 30.07.2005 ***************
      FCORR = 1.
      if (W_MEV.GT.1500.) 
     * FCORR=(1+3.*(W_MEV-1500.)/W_MEV)*dexp(-4.3*(W_MEV/1500.-1.))
c *************************************************
      FBV(1, 2, 3) = (FB(1, 2, 3) + FV(1, 2, 1))
      MULT(1, 2, 3) = (FB(1, 2, 3) + FV(1, 2, 1)) * tmp * FCORR
c ************  correction 30.07.2005 ***************
      FCORR = 1.
      if (W_MEV.GT.1650.)
     * FCORR=(1+1.2*(W_MEV-1650.)/W_MEV)*dexp(-1.85*(W_MEV/1650.-1.))      
c     * FCORR=(1+1.5*(W_MEV-1650.)/W_MEV)*dexp(-1.85*(W_MEV/1650.-1.))
c *************************************************
      FBV(3, 2, 3) = (FB(3, 2, 3) + FV(3, 2, 1))
      MULT(3, 2, 3) = (FB(3, 2, 3) + FV(3, 2, 1)) * tmp * FCORR

      FBV(5, 2, 3) = (FB(5, 2, 3) + FV(5, 2, 1))
      MULT(5, 2, 3) = (FB(5, 2, 3) + FV(5, 2, 1)) * tmp

C
C -----------------------------------  D13  -----------------------------------
      CALL SPLINT(W_sp, del_d13, del_d13_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_d13, eta_d13_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      IF (del.EQ.0.0) THEN
         phi_d13 = 0.0
         phase = (1.0, 0.0)
         tmp=(1.0D0,0.0D0)
      ELSE
         tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
         phi_d13 =DATAN((1.0 - eta * COS(2.0 * del))
     *           / eta / SIN(2.0 * del))
c         IF (W_MeV.GT.1400.AND.phi_d13.LT.0.0) phi_d13 = phi_d13 + pi
         IF (phi_d13.LT.0.0) phi_d13 = phi_d13 + pi
         del_phi = phi_d13 - phi_tmp
         phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)
      END IF
      FBV(2, 3, 1) = (FB(2, 3, 1) + FV(2, 3, 1))
      MULT(2, 3, 1) = (FB(2, 3, 1) + FV(2, 3, 1))*tmp

      FBV(2, 3, 2) = (FB(2, 3, 2) + FV(2, 3, 2))
      MULT(2, 3, 2) = (FB(2, 3, 2) + FV(2, 3, 2))*tmp

      FBV(4, 3, 1) = (FB(4, 3, 1) + FV(4, 3, 1))
      MULT(4, 3, 1) = (FB(4, 3, 1) + FV(4, 3, 1))*tmp

      FBV(4, 3, 2) = (FB(4, 3, 2) + FV(4, 3, 2))
      MULT(4, 3, 2) = (FB(4, 3, 2) + FV(4, 3, 2))*tmp

      FBV(6, 3, 1) = (FB(6, 3, 1) + FV(6, 3, 1))
      MULT(6, 3, 1) = (FB(6, 3, 1) + FV(6, 3, 1))*tmp

      FBV(6, 3, 2) = (FB(6, 3, 2) + FV(6, 3, 2))
      MULT(6, 3, 2) = (FB(6, 3, 2) + FV(6, 3, 2))*tmp
C
C -----------------------------------  D33  -----------------------------------
      CALL SPLINT(W_sp, del_d33, del_d33_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_d33, eta_d33_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      IF (del.EQ.0.0) THEN
         phi_d33 = 0.0
         phase = (1.0, 0.0)
         tmp=(1.0D0,0.0D0)
      ELSE
         tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
         phi_d33 = DATAN((1.0 - eta * COS(2.0 * del))
     *           / eta / SIN(2.0 * del))
c         IF (W_MeV.GT.1400.AND.phi_d33.LT.0.0) phi_d33 = phi_d33 + pi
         IF (phi_d33.LT.0.0) phi_d33 = phi_d33 + pi
         del_phi = phi_d33 - phi_tmp
         phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)
      END IF

      FBV(2, 3, 3) = (FB(2, 3, 3) + FV(2, 3, 1))
      MULT(2, 3, 3) = (FB(2, 3, 3) + FV(2, 3, 1))*tmp
c ************  correction 01.08.2005 ***************
      FCORR = 1.
      if (W_MEV.GT.1690.) 
     * FCORR=(1+8*(W_MEV-1690.)/W_MEV)*dexp(-4.*(W_MEV/1690.-1.))
c *************************************************
      FBV(4, 3, 3) = (FB(4, 3, 3) + FV(4, 3, 1))
      MULT(4, 3, 3) = (FB(4, 3, 3) + FV(4, 3, 1))*tmp * FCORR

      FBV(6, 3, 3) = (FB(6, 3, 3) + FV(6, 3, 1))
      MULT(6, 3, 3) = (FB(6, 3, 3) + FV(6, 3, 1))*tmp
C
C -----------------------------------  D15  -----------------------------------
      CALL SPLINT(W_sp, del_d15, del_d15_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_d15, eta_d15_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
         tmp = (piN + 1.0) / 2.0
      IF (del.EQ.0.0) THEN
         phi_d15 = 0.0
         phase = (1.0, 0.0)
      ELSE
         tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
         phi_d15 = DATAN((1.0 - eta * COS(2.0 * del))
     *           / eta / SIN(2.0 * del))
c         IF (W_MeV.GT.1400.AND.phi_d15.LT.0.0) phi_d15 = phi_d15 + pi
         IF (phi_d15.LT.0.0) phi_d15 = phi_d15 + pi
         del_phi = phi_d15 - phi_tmp
         phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)
      END IF

      FBV(1, 3, 1) = (FB(1, 3, 1) + FV(1, 3, 1))
      MULT(1, 3, 1) = (FB(1, 3, 1) + FV(1, 3, 1)) * tmp

      FBV(1, 3, 2) = (FB(1, 3, 2) + FV(1, 3, 2))
      MULT(1, 3, 2) = (FB(1, 3, 2) + FV(1, 3, 2)) * tmp

      FBV(3, 3, 1) = (FB(3, 3, 1) + FV(3, 3, 1))
      MULT(3, 3, 1) = (FB(3, 3, 1) + FV(3, 3, 1)) * tmp

      FBV(3, 3, 2) = (FB(3, 3, 2) + FV(3, 3, 2))
      MULT(3, 3, 2) = (FB(3, 3, 2) + FV(3, 3, 2)) * tmp

      FBV(5, 3, 1) = (FB(5, 3, 1) + FV(5, 3, 1))
      MULT(5, 3, 1) = (FB(5, 3, 1) + FV(5, 3, 1)) * tmp

      FBV(5, 3, 2) = (FB(5, 3, 2) + FV(5, 3, 2))
      MULT(5, 3, 2) = (FB(5, 3, 2) + FV(5, 3, 2)) * tmp
C
C -----------------------------------  D35  -----------------------------------
      CALL SPLINT(W_sp, del_d35, del_d35_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_d35, eta_d35_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
         tmp = (piN + 1.0) / 2.0

      FBV(1, 3, 3) = (FB(1, 3, 3) + FV(1, 3, 1))
      MULT(1, 3, 3) = (FB(1, 3, 3) + FV(1, 3, 1)) * tmp

      FBV(3, 3, 3) = (FB(3, 3, 3) + FV(3, 3, 1))
      MULT(3, 3, 3) = (FB(3, 3, 3) + FV(3, 3, 1)) * tmp

      FBV(5, 3, 3) = (FB(5, 3, 3) + FV(5, 3, 1))
      MULT(5, 3, 3) = (FB(5, 3, 3) + FV(5, 3, 1)) * tmp
C
C -----------------------------------  F15  -----------------------------------
      CALL SPLINT(W_sp, del_f15, del_f15_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_f15, eta_f15_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
      tmp = (piN + 1.0) / 2.0
      IF (del.EQ.0.0) THEN
         phi_f15 = 0.0
         phase = (1.0, 0.0)
      ElSE
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
         phi_f15 = DATAN((1.0 - eta * COS(2.0 * del))
     *           / eta / SIN(2.0 * del))
c         IF (W_MeV.GT.1400.AND.phi_15.LT.0.0) phi_f15 = phi_f15 + pi
         IF (phi_f15.LT.0.0) phi_f15 = phi_f15 + pi
         del_phi = phi_f15 - phi_tmp
         phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)
      END IF
      FBV(2, 4, 1) = (FB(2, 4, 1) + FV(2, 4, 1))
      MULT(2, 4, 1) = (FB(2, 4, 1) + FV(2, 4, 1))*tmp

      FBV(2, 4, 2) = (FB(2, 4, 2) + FV(2, 4, 2))
      MULT(2, 4, 2) = (FB(2, 4, 2) + FV(2, 4, 2))*tmp

      FBV(4, 4, 1) = (FB(4, 4, 1) + FV(4, 4, 1))
      MULT(4, 4, 1) = (FB(4, 4, 1) + FV(4, 4, 1))*tmp

      FBV(4, 4, 2) = (FB(4, 4, 2) + FV(4, 4, 2))
      MULT(4, 4, 2) = (FB(4, 4, 2) + FV(4, 4, 2))*tmp

      FBV(6, 4, 1) = (FB(6, 4, 1) + FV(6, 4, 1))
      MULT(6, 4, 1) = (FB(6, 4, 1) + FV(6, 4, 1))*tmp

      FBV(6, 4, 2) = (FB(6, 4, 2) + FV(6, 4, 2))
      MULT(6, 4, 2) = (FB(6, 4, 2) + FV(6, 4, 2))*tmp
C
C -----------------------------------  F35  -----------------------------------
      CALL SPLINT(W_sp, del_f35, del_f35_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_f35, eta_f35_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
         tmp = (piN + 1.0) / 2.0
      IF (del.EQ.0.0) THEN
         phi_f35 = 0.0
         phase = (1.0, 0.0)
      ELSE
         tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
         phi_f35 = DATAN((1.0 - eta * COS(2.0 * del))
     *           / eta / SIN(2.0 * del))
c         IF (W_MeV.GT.1400.AND.phi_f35.LT.0.0) phi_f35 = phi_f35 + pi
         IF (phi_f35.LT.0.0) phi_f35 = phi_f35 + pi
         del_phi = phi_f35 - phi_tmp
         phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)
       END IF
      FBV(2, 4, 3) = (FB(2, 4, 3) + FV(2, 4, 1))
      MULT(2, 4, 3) = (FB(2, 4, 3) + FV(2, 4, 1)) * tmp

      FBV(4, 4, 3) = (FB(4, 4, 3) + FV(4, 4, 1))
      MULT(4, 4, 3) = (FB(4, 4, 3) + FV(4, 4, 1)) * tmp

      FBV(6, 4, 3) = (FB(6, 4, 3) + FV(6, 4, 1))
      MULT(6, 4, 3) = (FB(6, 4, 3) + FV(6, 4, 1)) * tmp
C
C -----------------------------------  F17  -----------------------------------
      CALL SPLINT(W_sp, del_f17, del_f17_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_f17, eta_f17_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
         tmp = (piN + 1.0) / 2.0

      FBV(1, 4, 1) = (FB(1, 4, 1) + FV(1, 4, 1))
      MULT(1, 4, 1) = (FB(1, 4, 1) + FV(1, 4, 1)) * tmp

      FBV(1, 4, 2) = (FB(1, 4, 2) + FV(1, 4, 2))
      MULT(1, 4, 2) = (FB(1, 4, 2) + FV(1, 4, 2)) * tmp

      FBV(3, 4, 1) = (FB(3, 4, 1) + FV(3, 4, 1))
      MULT(3, 4, 1) = (FB(3, 4, 1) + FV(3, 4, 1)) * tmp

      FBV(3, 4, 2) = (FB(3, 4, 2) + FV(3, 4, 2))
      MULT(3, 4, 2) = (FB(3, 4, 2) + FV(3, 4, 2)) * tmp

      FBV(5, 4, 1) = (FB(5, 4, 1) + FV(5, 4, 1))
      MULT(5, 4, 1) = (FB(5, 4, 1) + FV(5, 4, 1)) * tmp

      FBV(5, 4, 2) = (FB(5, 4, 2) + FV(5, 4, 2))
      MULT(5, 4, 2) = (FB(5, 4, 2) + FV(5, 4, 2)) * tmp
C
C -----------------------------------  F37  -----------------------------------
      CALL SPLINT(W_sp, del_f37, del_f37_ss, 143, W_MeV, del)
      CALL SPLINT(W_sp, eta_f37, eta_f37_ss, 143, W_MeV, eta)
      IF (W_MeV.LT.1079.) del=0.00001
      IF (W_MeV.LT.1079.) eta=1.
      piN = eta * (COS(2.0 * del) + (0.0, 1.0) * SIN(2.0 * del))
         tmp = (piN + 1.0) / 2.0
      IF (del.EQ.0.0) THEN
         phi_f37 = 0.0
         phase = (1.0, 0.0)
      ELSE
         tmp = (piN + 1.0) / 2.0
      ARG= DIMAG(tmp) /DREAL(tmp)
      phi_tmp = DATAN(ARG)
         phi_f37 = DATAN((1.0 - eta * COS(2.0 * del))
     *           / eta / SIN(2.0 * del))
c         IF (W_MeV.GT.1400.AND.phi_f37.LT.0.0) phi_f37 = phi_f37 + pi
         IF (phi_f37.LT.0.0) phi_f37 = phi_f37 + pi
         del_phi = phi_f37 - phi_tmp
         phase = COS(del_phi) + (0.0, 1.0) * SIN(del_phi)
       END IF

      FBV(1, 4, 3) = (FB(1, 4, 3) + FV(1, 4, 1))
      MULT(1, 4, 3) = (FB(1, 4, 3) + FV(1, 4, 1)) * tmp

      FBV(3, 4, 3) = (FB(3, 4, 3) + FV(3, 4, 1))
      MULT(3, 4, 3) = (FB(3, 4, 3) + FV(3, 4, 1)) * tmp

      FBV(5, 4, 3) = (FB(5, 4, 3) + FV(5, 4, 1))
      MULT(5, 4, 3) = (FB(5, 4, 3) + FV(5, 4, 1)) * tmp
c ********************************************************
       cof=1000.*mpip*HQC

       DO 2222 l = 1, l_max + 1
       DO 2222 k = 1, 6
       DO 2222 i = 1, 3
       iv=i
       if (i.eq.3) iv=1
       MULT(k,l,i) = MULT(k,l,i) * cof
       FBV(k,l,i) = FBV(k,l,i) * cof
       IF (IUNI.EQ.0) MULT(k,l,i)=FBV(k,l,i)
       IF (l.GT.4) FBV(k,l,i)=(FB(k,l,i)+FV(k,l,iv))*cof
       IF (l.GT.4) MULT(k,l,i)=FBV(k,l,i)
2222   CONTINUE

c       write (6,156) w*HQC, (MULT(1,2,k), k=1,3)
c       write (6,156) w*HQC, (MULT(2,2,k), k=1,3)
c       write (6,156) w*HQC, (MULT(3,2,k), k=1,3)
c       write (6,156) w*HQC, (MULT(4,2,k), k=1,3)
c156    format(1x,f10.2,2x, 6(E10.4,2x))

c ********************************************************

      RETURN
      END

      SUBROUTINE CGLN_DISP(ISO, W, Q2, x, MN, MPI, AIS, AIS0, F_back)
c *****************************************************************************
c     CGLN amplitudes F_background(i) from multipoles
c     F_unitarized = F_Born+Vector  +  F_background
c *****************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AIS0(10,6,3)
      REAL*8 MN, MPI , MPIP
      DIMENSION P_S(10), P_SS(10)
      COMPLEX*16 AIS(10,8,3), F_iso(6,3), F_back(6), MNEW(6)
      PARAMETER (pi = 3.1415926536D0)

	MPIP=139.5685D0/197.3285
      l_max = 4
c      l_max=2

      HQC = 197.3285


      P_S(1) = 0.0
      P_S(2) = 1.0
      P_S(3) = 3.0 * x
      P_S(4) = (15.0 * x * x - 3.0) / 2.0
      P_S(5) = (35.0 * x * x * x - 15.0 * x) / 2.0
      P_SS(1) = 0.0
      P_SS(2) = 0.0
      P_SS(3) = 3.0
      P_SS(4) = 15.0 * x
      P_SS(5) = 105.0 / 2.0 * x * x - 7.5

      DO 6000 k = 1, 6
      f_back(k)=(0.d0,0.d0)
         DO 7000 j = 1, 3
            F_iso(k, j) = (0.0, 0.0)
 7000    CONTINUE
 6000 CONTINUE

      DO 4000 l = 1, l_max
          DO 5000 j = 1, 3
      DO 5001 k = 1, 6
5001  MNEW(k) = AIS(l,k,j) - AIS0(l,k,j)
             F_iso(1, j) = F_iso(1, j)
     *                + ((REAL(l) - 1.0) * MNEW(3)
     *                     + MNEW(1)) * P_S(l + 1)
             F_iso(5, j) = F_iso(5, j)
     *                + REAL(l) * MNEW(5) * P_S(l + 1)
             IF (l.GT.1) THEN
                F_iso(2, j) = F_iso(2, j)
     *                   + (REAL(l) * MNEW(3)
     *                   + (REAL(l) - 1.0) * MNEW(4)) * P_S(l)
                F_iso(3, j) = F_iso(3, j)
     *                   + (MNEW(1)
     *                       - MNEW(3)) * P_SS(l + 1)
     *                   + (MNEW(2)
     *                       + MNEW(4)) * P_SS(l - 1)
                F_iso(6, j) = F_iso(6, j)
     *                   + ((REAL(l) - 1.0) * MNEW(6)
     *                       - REAL(l) * MNEW(5)) * P_S(l)
             END IF
             IF (l.GT.2) THEN
                F_iso(1, j) = F_iso(1, j)
     *                   + (REAL(l) * MNEW(4)
     *                   + MNEW(2)) * P_S(l - 1)
                F_iso(4, j) = F_iso(4, j)
     *                   + (MNEW(3) - MNEW(1)
     *                   - MNEW(4) - MNEW(2)) * P_SS(l)
                F_iso(5, j) = F_iso(5, j)
     *                 - (REAL(l) - 1.0) * MNEW(6) * P_S(l - 1)
             END IF
 5000     CONTINUE
 4000 CONTINUE

         DO 4500 k = 1, 6
      IF (ISO.EQ.1) THEN
            F_back(k) = F_iso(k, 2) + F_iso(k, 1) / 3.0
     *                 + 2.0 * F_iso(k, 3) / 3.0
      ELSE IF (ISO.EQ.2) THEN
            F_back(k) = - F_iso(k, 2) + F_iso(k, 1) / 3.0
     *                 + 2.0 * F_iso(k, 3) / 3.0
      ELSE IF (ISO.EQ.3) THEN
            F_back(k) = SQRT(2.0) * (F_iso(k, 2) + F_iso(k, 1) / 3.0
     *                 - F_iso(k, 3) / 3.0)
      ELSE IF (ISO.EQ.4) THEN
            F_back(k) = SQRT(2.0) * (F_iso(k, 2) - F_iso(k, 1) / 3.0
     *                 + F_iso(k, 3) / 3.0)
      END IF
 4500       CONTINUE

      DO 4900 k = 1, 6
         F_back(k) = F_back(k)/1000./MPIP
 4900 CONTINUE

      RETURN
      END

        SUBROUTINE MRES(ISO,WCM,Q2FM,WGCM,AR1,AR3)
C******************************************************************************
        IMPLICIT NONE
        INTEGER ISO,MODE,GAUGE, I33, IS11F, IS11S, IS31, IS11T,IP11S
        INTEGER IP13, IP31, ID15, IF35, IF37
	COMMON/newres/ S31, P13, P31, D15, F35, F37
	DOUBLE PRECISION  P13, P31, D15, F35, F37, S31
        DOUBLE PRECISION WCM,Q2FM,XX,WGCM
        COMPLEX*16 AR1(10,8), AR3(10,8)
        INTEGER IBORN,IVEC,IOME,IRHO,JFORM,IP33,IP11,ID13,ID33,IS11,IF15
        DOUBLE PRECISION M1FM,M2FM,MPIFM
        COMMON /LINK1/ M1FM,M2FM,MPIFM,IBORN,IVEC,IOME,IRHO,JFORM,
     &                IP33,IP11,ID13,ID33,IS11F,IS11S,IF15,IS11T,IP11S
    	DOUBLE PRECISION OmegL,Q2,W,qcm,kcm,Q2G,wGacm,kgcm
	COMMON /KinVar/ OmegL, Q2, W, wGacm, kcm
    	DOUBLE PRECISION mi, mf, mPi , Pi, hqc
	COMMON /Mass/ mi, mf, mPi
	INTEGER  k, l
	DIMENSION ReM3MN(3), ReM2MN(3),ImM2MN(3),ImM3MN(3)
       DOUBLE PRECISION ReM1P, ImM1P, ReM1M, ImM1M, ReM2M, ImM2M,
     * ReM1MN,ImM1MN,ReM0PL,ImM0PL,ReS31(3),ImS31(3),ReM2MN,
     * ReF37(3),ImF37(3), ReE1P(3), ImE1P(3), P4S, P4SS,ImM2MN,
     * ReM2P(3),ImM2P(3), ReP31(3),ImP31(3), ReM35(3), ImM35(3),
     * ReM1Ms, ImM1Ms
       DIMENSION ReM1P(3),ImM1P(3),ReM1M(3),ImM1M(3),
     * ReM1Ms(3),ImM1Ms(3),
     * ReM2M(3),ImM2M(3),ReM1MN(3), ImM1MN(3), ReM0PL(3),ImM0pl(3)
	DOUBLE PRECISION ReM0P,ImM0P,ReM3M,ImM3M,ReM3MN,ImM3MN
	DIMENSION ReM0P(3),ImM0P(3),ReM3M(3),ImM3M(3)
	DIMENSION ReM0PL2(3),ImM0PL2(3),ReM1PL(3),ImM1PL(3)
	DIMENSION ReM0PL3(3),ImM0PL3(3)	
	DOUBLE PRECISION ReM0PL2,ImM0PL2,ReM1PL,ImM1PL,ReM0PL3,ImM0PL3	
        DOUBLE PRECISION ReD33M(3), ImD33M(3),ReD33(3), ImD33(3)
	PARAMETER ( Pi = 3.1415926536D0, hqc = 197.3285)
c
        W=WCM*hqc
        Q2=Q2fm*hqc**2
c  same nucleon mass is necessary for gauge invariance
        mi=(m1fm+m2fm)/2.*hqc
        mf=mi
        mPi=mpifm*hqc
	OmegL = (W*W+Q2-mi*mi)/2./mi
	kcm = SQRT(((W*W+mPi*mPi-mf*mf)/(2*W))**2-mPi*mPi)
        kgcm = (W*W-mi*mi)/2./W
        wGacm=WGCM*hqc
        qcm=sqrt(wGacm**2+Q2)
c ********************************
	DO 99 L=1,10
	DO 99 K=1,8
	AR1(L,K)=(0.D0,0.D0)
	AR3(L,K)=(0.D0,0.D0)
 99     CONTINUE
c *****************************
        DO 110 l=1,3
        ReM0P(l)=0
        ImM0P(l)=0
        ReM1P(l)=0
        ImM1P(l)=0
        ReM1M(l)=0
        ImM1M(l)=0
        
        ReM1Ms(l)=0
        ImM1Ms(l)=0
        
        ReM2M(l)=0
        ImM2M(l)=0
        ReM3M(l)=0
        ImM3M(l)=0
        ReD33M(l)=0
        ImD33M(l)=0
        ReM0PL(l)=0
        ImM0PL(l)=0
        ReM0PL2(l)=0
        ImM0PL2(l)=0
        ReM0PL3(l)=0
        ImM0PL3(l)=0
        ReS31(l)=0.
	ImS31(l)=0.
	ReF37(l)=0.
	ImF37(l)=0.
	ReE1P(l)=0.
	ImE1P(l)=0.
	ReM2P(l)=0.
	ImM2P(l)=0.
	ReP31(l)=0.
	ImP31(l)=0.
	ReM35(l)=0.
	ImM35(l)=0.
  110   CONTINUE
c *****************************************************************
c *******     resonance contributions from 14 resonances     *******
c *****************************************************************
        IF (IP33.EQ.1) CALL P33NEW(ReM1P,ImM1P)
        IF (IP11.EQ.1) CALL P11NEW(ISO,ReM1M,ImM1M)
        IF (IP11S.EQ.1) CALL P11sec(ISO,ReM1Ms,ImM1Ms)
        
        IF (ID13.EQ.1) CALL D13NEW(ISO,ReM2M,ImM2M)
        IF (ID33.EQ.1) CALL D33NEW(ReD33M,ImD33M)
        IF (IS11F.EQ.1) CALL S11FST(ISO,ReM0PL,ImM0PL)
        IF (IS11S.EQ.1) CALL S11SEC(ISO,ReM0PL2,ImM0PL2)
        IF (IS11T.EQ.1) CALL S11trd(ISO,ReM0PL3,ImM0PL3)        
        IF (IF15.EQ.1) CALL F15NEW(ISO,ReM3M,ImM3M)
	IS31=INT(S31+0.1)
        IF (IS31.EQ.1) CALL S31NEW(ReS31,ImS31)
c *************************************************
	IP31=INT(P31+0.1)
        IF (IP31.EQ.1) CALL P31NEW(ReP31,ImP31)
	IF35=INT(F35+0.1)
        IF (IF35.EQ.1) CALL F35NEW(ReM35,ImM35)
	IF37=INT(F37+0.1)
        IF (IF37.EQ.1) CALL F37NEW(ReF37,ImF37)
	IP13=INT(P13+0.1)
        IF (IP13.EQ.1) CALL P13NEW(ISO,ReE1P,ImE1P)
	ID15=INT(D15+0.1)
        IF (ID15.EQ.1) CALL D15NEW(ISO,ReM2P,ImM2P)
c *************************************************
        DO 156 L=1,3
        ReM1M(L)=ReM1M(L)+ReM1Ms(L)
        ImM1M(L)=ImM1M(L)+ImM1Ms(L)        
        REM0P(L)=REM0PL(L)+REM0PL2(L)+REM0PL3(L)
156     IMM0P(L)=IMM0PL(L)+IMM0PL2(L)+IMM0PL3(L)
c ******************  E0+  and L0+ *****************
	  AR1(1,1) = CMPLX(ReM0P(1),ImM0P(1))
	  AR3(1,1) = CMPLX(ReS31(1),ImS31(1))
	  AR1(1,5) = CMPLX(ReM0P(3),ImM0P(3))
	  AR3(1,5) = CMPLX(ReS31(3),ImS31(3))
c ******************  E1+  and L1+ *****************
	  AR1(2,1) = CMPLX(ReE1P(1),ImE1P(1))
	  AR3(2,1) = CMPLX(ReM1P(2),ImM1P(2))
	  AR1(2,5) = CMPLX(ReE1P(3),ImE1P(3))
	  AR3(2,5) = CMPLX(ReM1P(3),ImM1P(3))
c ******************  E2+  and L2+ *****************
	  AR1(3,1) = CMPLX(ReM2P(1),ImM2P(1))
	  AR1(3,5) = CMPLX(ReM2P(3),ImM2P(3))
c ******************  E3+  and L3+ *****************
	  AR3(4,1) = CMPLX(ReF37(2),ImF37(2))
	  AR3(4,5) = CMPLX(ReF37(3),ImF37(3))
c ******************  E2-  and L2- *****************
	  AR1(3,2) = CMPLX(ReM2M(1),ImM2M(1))
	  AR1(3,6) = CMPLX(ReM2M(3),ImM2M(3))
	  AR3(3,2) = CMPLX(ReD33M(1),ImD33M(1))
	  AR3(3,6) = CMPLX(ReD33M(3),ImD33M(3))
c ******************  E3-  and L3- *****************
	  AR1(4,2) = CMPLX(ReM3M(1),ImM3M(1))
	  AR1(4,6) = CMPLX(ReM3M(3),ImM3M(3))
	  AR3(4,2) = CMPLX(ReM35(2),ImM35(2))
	  AR3(4,6) = CMPLX(ReM35(3),ImM35(3))
c ******************  M1+ ,  M1-  and L1-*****************
	  AR1(2,3) = CMPLX(ReE1P(2),ImE1P(2))
	  AR1(2,4) = CMPLX(ReM1M(2),ImM1M(2))
	  AR1(2,6) = CMPLX(ReM1M(3),ImM1M(3))
	  AR3(2,3) = CMPLX(ReM1P(1),ImM1P(1))
	  AR3(2,4) = CMPLX(ReP31(1),ImP31(1))
	  AR3(2,6) = CMPLX(ReP31(3),ImP31(3))
c ******************  M2+  and  M3+ *****************
	  AR1(3,3) = CMPLX(ReM2P(2),ImM2P(2))
	  AR1(3,5) = CMPLX(ReM2P(3),ImM2P(3))
	  AR3(4,3) = CMPLX(ReF37(1),ImF37(1))
	  AR3(4,5) = CMPLX(ReF37(3),ImF37(3))
c ******************  M2- and M3- *****************
	  AR1(3,4) = CMPLX(ReM2M(2),ImM2M(2))
	  AR3(3,4) = CMPLX(ReD33M(2),ImD33M(2))
	  AR1(4,4) = CMPLX(ReM3M(2),ImM3M(2))
	  AR3(4,4) = CMPLX(ReM35(1),ImM35(1))
c *************************************************

        RETURN
    	END


      SUBROUTINE MULT_TOT(W, Q2, MN, MPI, CK0, lmax1,
     & AP0M, AIS, APN3, ACH, AIS0,IUNI,IBGR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AIS0(10,6,3)
      COMMON /MULTBR/MULT(6,10,3), FBV(6,10,3)
      REAL*8 MN, MPI
      COMPLEX*16 MULT, ARP1(10,8), ARN1(10,8), AR3(10,8)
      COMPLEX*16 AP0M(10,8,3),AIS(10,8,3),APN3(10,8,3)
      COMPLEX*16 ARP, ARM, AR0, ABP, ABM, AB0, ACH(10,8,4)

        DO 99  L1=1,10
        DO 99  K=1,6
        DO 98  I=1,3
	AIS0(L1,K,I)=0.D0
        MULT(K,L1,I)=(0.D0,0.D0)
        APN3(L1,K,I)=(0.D0,0.D0)
        AIS(L1,K,I)=(0.D0,0.D0)
        AP0M(L1,K,I)=(0.D0,0.D0)
 98     ACH(L1,K,I)=(0.D0,0.D0)
 99     ACH(L1,K,4)=(0.D0,0.D0)
c *********************************************************
      IF (IBGR .EQ. 1) CALL BACKGR(W, Q2, MN, MPI, lmax1,IUNI)
	IF (IUNI.EQ.0) GO TO 100
        CALL MRES(2,W,Q2,CK0,ARN1,AR3)
        CALL MRES(1,W,Q2,CK0,ARP1,AR3)
100     CONTINUE
c *********************************************************


      DO 1 L1=1,lmax1
      DO 1 K=1,6

      ABP = ( MULT(K,L1,1) + 2.*MULT(K,L1,3) )/3.
      ABM = ( MULT(K,L1,1) - MULT(K,L1,3) )/3.
      AB0 = MULT(K,L1,2)

      ARP =( SQRT(2.)*AR3(L1,K) + (ARN1(L1,K)-ARP1(L1,K))/2.)/SQRT(3.)
      ARM =(-AR3(L1,K)/SQRT(2.) + (ARN1(L1,K)-ARP1(L1,K))/2.)/SQRT(3.)
      AR0 = - (ARN1(L1,K) + ARP1(L1,K) )/2./SQRT(3.)

      AP0M(L1,K,1) = ABP + ARP*1000.
      AP0M(L1,K,2) = AB0 + AR0*1000.
      AP0M(L1,K,3) = ABM + ARM*1000.

      AIS(L1,K,1) = AP0M(L1,K,1) + 2.*AP0M(L1,K,3)
      AIS(L1,K,2) = AP0M(L1,K,2)
      AIS(L1,K,3) = AP0M(L1,K,1) - AP0M(L1,K,3)

      APN3(L1,K,1)= AIS(L1,K,2) + AIS(L1,K,1)/3.
      APN3(L1,K,2)= AIS(L1,K,2) - AIS(L1,K,1)/3.
      APN3(L1,K,3)= AIS(L1,K,3)

      ACH(L1,K,1) = AP0M(L1,K,1) + AP0M(L1,K,2)
      ACH(L1,K,2) = AP0M(L1,K,1) - AP0M(L1,K,2)
      ACH(L1,K,3) = SQRT(2.)*(AP0M(L1,K,2) + AP0M(L1,K,3))
      ACH(L1,K,4) = SQRT(2.)*(AP0M(L1,K,2) - AP0M(L1,K,3))

c ***************  only Born + omega + rho  *************
      AIS0(L1,K,1)= FBV(K,L1,1)
      AIS0(L1,K,2)= FBV(K,L1,2)
      AIS0(L1,K,3)= FBV(K,L1,3)


1     CONTINUE


      RETURN
      END

      SUBROUTINE SOLUTIONS(SOLUTION,TEXT,XE0,XS0,XMIX0,
     &  X3P330,X1P330,XSP330,X1S310,XSS310,X3D330,X1D330,XSD330,
     &  X1P310,XSP310,X3F350,X1F350,XSF350,X3F370,X1F370,XSF370,
     &  X1P11p0,XSP11p0,X1S11p0,XSS11p0,X1S2p0,XSS2p0,
     &  X3D13p0,X1D13p0,XSD13p0,X3F15p0,X1F15p0,XSF15p0,
     &  X3P13p0,X1P13p0,XSP13p0,X3D15p0,X1D15p0,XSD15p0,
     &  X1P11n0,XSP11n0,X1S11n0,XSS11n0,X1S2n0,XSS2n0,
     &  X3D13n0,X1D13n0,XSD13n0,X3F15n0,X1F15n0,XSF15n0,
     &  X3P13n0,X1P13n0,XSP13n0,X3D15n0,X1D15n0,XSD15n0)

      IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION X(100)
       COMMON /SOLUTION/ XX(100),ISOL
       COMMON /E0CORR/ XE,XS,XMIX
       COMMON /EPLCORR/ XEPL,XSPL
       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
       COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
       COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
       COMMON/secS11/ X1S2p,XSS2p,X1S2n,XSS2n
       COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
       COMMON/parPPn/ X1P13n,X3P13n,XSP13n
       COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n
       COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37
       common/parvect/ xom1, xom2, xrho1, xrho2, xom, xrho
	common/newpar/JRHOFORM,JOMEFORM,D13mode
      COMMON/FORMF/ xlampi,xlamax,mode,gauge
c ***************************************************************

	INTEGER D13MODE,Dispersion,GAUGE
	CHARACTER*20 SOLUTION, filename
      CHARACTER*64 text, text1, text2, tnew(3)
c ************************************************************
      DATA tnew/' ******* standard MAID calculation  *********',
     & ' ******** DR MAID calculation (approximate) **********',
     & '     ****** DR MAID calculation (exact) *******'/
c ************************************************************
        NPAR=65
c ************************************************************
       I = 1
       DO WHILE (Solution(I:I) .NE. ' ')
       I = I + 1
       END DO
	 filename=Solution(1:I-1)//'.dat'
       OPEN(21,file=filename,STATUS='OLD')
c ************************************************************
      read (21,111) text
 111  format(A64)
      read (21,*) JRHOFORM,JOMEFORM,D13mode,Dispersion
 112  format(3I3)
      DO 1 I=1,NPAR
      read (21,*) text1, text2, X(I)
1     continue
      CLOSE (21)
c ************************************************************
c      write (6,113) text,JRHOFORM,JOMEFORM,D13mode,Dispersion
c 113  format(/,8x,A64,/,55x,4I2)
      If (Dispersion .eq. 0 .and. isol .gt. 1) then
	 isol=1
	 write (6,114)
 114	 format(5x,'The current solution does not allow the use of',
     &         ' dispersion relations.')
	endif
c      write (6,115) tnew(ISOL)
c 115  format(15x,A64,/)
c	 ISOL1=ISOL
	 xom10=1
	 xom20=1
	 xom0 =1
c	 if (JOMEFORM .eq. 0) xom0=10000
	 xrho10=1
	 xrho20=1
	 xrho0=1
c	 if (JRHOFORM .eq. 0) xrho0=10000
	 X3P33    =  Param(X3P330 ,X(1))
	 X1P33    =  Param(X1P330 ,X(2))
	 X1P11p   =  Param(X1P11p0,X(3))
	 X1S11p   =  Param(X1S11p0,X(4))
	 X3D13p   =  Param(X3D13p0,X(5))
	 X1D13p   =  Param(X1D13p0,X(6))
	 X1S2p    =  Param(X1S2p0 ,X(7))
	 X3F15p   =  Param(X3F15p0,X(8))
	 X1F15p   =  Param(X1F15p0,X(9))
	 X1S31    =  Param(X1S310 ,X(10))
	 X3D33    =  Param(X3D330 ,X(11))
	 X1D33    =  Param(X1D330 ,X(12))
	 X1P13p   =  Param(X1P13p0,X(13))
	 X3P13p   =  Param(X3P13p0,X(14))
	 X1D15p   =  Param(X1D15p0,X(15))
	 X3D15p   =  Param(X3D15p0,X(16))
	 X1P31    =  Param(X1P310 ,X(17))
	 X1F35    =  Param(X1F350 ,X(18))
	 X3F35    =  Param(X3F350 ,X(19))
	 X1F37    =  Param(X1F370 ,X(20))
	 X3F37    =  Param(X3F370 ,X(21))
	 X1P11n   =  Param(X1P11n0,X(22))
	 X1S11n   =  Param(X1S11n0,X(23))
	 X3D13n   =  Param(X3D13n0,X(24))
	 X1D13n   =  Param(X1D13n0,X(25))
	 X1S2n    =  Param(X1S2n0 ,X(26))
	 X3F15n   =  Param(X3F15n0,X(27))
	 X1F15n   =  Param(X1F15n0,X(28))
	 X1P13n   =  Param(X1P13n0,X(29))
	 X3P13n   =  Param(X3P13n0,X(30))
	 X1D15n   =  Param(X1D15n0,X(31))
	 X3D15n   =  Param(X3D15n0,X(32))
       XSP33    =  Param(XSP330 ,X(33))
       XSP11p   =  Param(XSP11p0,X(34))
	 XSS11p   =  Param(XSS11p0,X(35))
	 XSD13p   =  Param(XSD13p0,X(36))
	 XSS2p    =  Param(XSS2p0 ,X(37))
	 XSF15p   =  Param(XSF15p0,X(38))
	 XSS31    =  Param(XSS310 ,X(39))
	 XSD33    =  Param(XSD330 ,X(40))
	 XSP13p   =  Param(XSP13p0,X(41))
	 XSD15p   =  Param(XSD15p0,X(42))
	 XSP31    =  Param(XSP310 ,X(43))
	 XSF35    =  Param(XSF350 ,X(44))
	 XSF37    =  Param(XSF370 ,X(45))
	 XSP11n   =  Param(XSP11n0,X(46))
	 XSS11n   =  Param(XSS11n0,X(47))
	 XSD13n   =  Param(XSD13n0,X(48))
	 XSS2n    =  Param(XSS2n0 ,X(49))
	 XSF15n   =  Param(XSF15n0,X(50))
	 XSP13n   =  Param(XSP13n0,X(51))
	 XSD15n   =  Param(XSD15n0,X(52))
	 XE       =  Param(XE0    ,X(53))
       XS       =  Param(XS0    ,X(54))
       XMIX     =  Param(XMIX0  ,X(55))
c  the following parameters are presently not changed by the input
	 xom1     =   xom10  *  X(56)
	 xom2     =   xom20  *  X(57)
	 xom      =    xom0  *  X(58)
	 xrho1    =  xrho10  *  X(59)
	 xrho2    =  xrho20  *  X(60)
	 xrho     =   xrho0  *  X(61)
	 XEPL     =             X(62)
	 XSPL     =             X(63)
       XLAMAX   =             X(64)
       XLAMPI   =             X(65)
c  XEPL, XSPL are no longer used in Maid2007
      RETURN
      END
c
	Double Precision Function Param(Xinp,Xfit)
	Real*8 Xinp,Xfit,eps
	eps=2.E-4
	IF (Xfit .EQ. 0) THEN
	 IF (ABS(Xinp-1) .LT. eps) THEN
	  Param=0
	 ELSE
	  Param=Xinp
	 ENDIF
	ELSE
	 Param=Xinp*Xfit
	ENDIF
	RETURN
	END

      subroutine HP33(Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
       COMMON/param3/ XMP33,XEP33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
       COMMON /Mass/ ami, amf, amPION
       COMMON /KinVar/ OmegL, Q2, W, wGacm, akcm
c ***************************************************************
        Fq=dexp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)
c *****************************************************************
      Phi_R=0.
c      X1P33=1
c      X3P33=1
C ******  fit parameters changed to XE and XM, 19 April 2005
C ***     changes onlyn inside of HP33
C ***     outside XMP33 is still calles X3P33 and XEP33 -> X1P33
c      A10=-140.385 *X1P33
c      A30=-265.220 *X3P33      
c      AE0 = -6.36992  * XEP33
c      AM0 = 299.880   * XMP33

      AE0 = -6.36992
      AM0 = 299.880
 
      AE1= AE0*(1.-0.0205657*Q2G)*Fq*dexp(0.05*Q2G)
c      AM= AM0*(1.-0.012*Q2G)*Fq
      AM1= AM0*(1.+0.0095*Q2G)*Fq*DEXP(-0.02*Q2G) 
      A11=-(1./2.)*(3*AE1+AM1)
      A31= (1./2.)*Sqrt(3.)*(AE1-AM1)
      
      AE=AE1 * XEP33
      AM=AM1 * XMP33
      A1=-(1./2.)*(3*AE+AM)
      A3= (1./2.)*Sqrt(3.)*(AE-AM)

      X1P33=A1/A11
      X3P33=A3/A31
c ***********  parametrization based on Siegert theorem *******************
       CNORM=1.
      IF (IPRN.NE.1) GO TO 111
       ami=938.2723
       W = 1232.
       CNORM=1000.
111    Q2pt= -(W-ami)**2*1.D-6
c
c      S11=AE*(1.+ 15.168012*Q2G**2)/(1.+15.168012*Q2pt**2)
c     & *dexp(-4.444704*(Q2G-Q2pt))*CNORM
c
c       ami_GEV=ami*1.D-3     
c      S12 = 108.302*(Q2pt-Q2G)/(1+4.740117*Q2G)**2
c     & *DEXP(-Q2G*0.375903)/ami_GEV**2  
c   
c      S1 = -sqrt(2.)*(S11 + S12*XSP33)*qcm/(W-ami) 
c
c ---------  Buchman parametrization -----------
c       Ab=1.2 
c       Db=4.9
c       tau=Q2/(4.*ami*ami)
c       CMR=Ab/(1.+Db*tau)
c       AMb= AM0*(1.+0.12*Q2G)*Fq*DEXP(-0.02*Q2G)       
c       S1 = sqrt(2.)*qcm*CMR*AMb/(8*ami) * XSP33
c --------  our new parametrization -------------
       Db=4.9
       Q2MEV=Q2G*1.E+6
       qcm_D=(1232**2-ami**2)/2./1232. 
       tau=Q2MEV/(4.*ami*ami)
       AS0=12.403
       S1=sqrt(2.)*AS0*(1.+0.12*Q2G)/(1+Db*tau)*Fq*DEXP(-0.02*Q2G)
     &   *qcm/qcm_D * CNORM* XSP33
         
c       Rasym=-(AS0*0.12*2*ami**2)/(AM0*0.01*1232*qcm_D*Db)
c       write (6,156) Q2G, qcm, S1, S1tst, AS0, Rasym
c156    format(1x,7(E12.6,2x))   
c ***************************************************        
c      A1=-(1./2.)*(3*AE+AM)
c      A3= (1./2.)*Sqrt(3.)*(AE-AM)

       AS=-(1./2.)*Sqrt(2.0)*S1
      
      IF (IPRN.EQ.0) return
      write (6,123) Phi_R,  X1P33, A1,X3P33, A3, XSP33, S1
123   format( 5x,'P33(1232):',1x,F6.2,3(1x,2(F8.3,1x)))
      return
	end

      subroutine HP11(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
       COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
        PI=3.1415926536D0
	  Phi_R = -15.48
        PhiR=Phi_R*Pi/180.
	  COSR=DCOS(PhiR)
        IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
c ***************************************************************
      A10 = -59.134/COSR  *X1P11p
      S10 = 4.02/COSR     *XSP11p 
      
      A1=A10*(1.-1.221691*Q2G-0.55*Q2G**4)*exp(-1.51*Q2G)
      S1= S10*(1.+ 41.001*Q2G+1.5*Q2G**4)*exp(-1.75*Q2G)            

      GO TO 20
c ***************************************
10    A10 = 52.137/COSR *X1P11n
      S10 = -40./COSR   *XSP11n
      A1= A10*(1.+0.9450*Q2G)*exp(-1.77*Q2G)
      S1=S10*(1.+ 2.97843*Q2G)*exp(-1.55*Q2G)
c **************************************************
20    AM=A1
      AS=-sqrt(2.)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123)  Phi_R, X1P11p, A1, XSP11p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123   format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

        return
	end

      subroutine HD13(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
c ***************************************************************
        PI=3.1415926536D0
      IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
	Phi_R = 32.
        PhiR=Phi_R*Pi/180.
	COSR=DCOS(PhiR)

        A10=-23.2016/COSR *X1D13p
        A30=136.2258/COSR *X3D13p
        S10=-53.9019/COSR *XSD13p  

c      A1= A10*(1.+8.067698*Q2G)*exp(-1.08576*Q2G)  
      A1= A10*(1.+7.7698*Q2G)*exp(-1.08576*Q2G)          

      A3= A30*(1.+0.69263*Q2G)*exp(-2.104*Q2G)      
      S1= S10*(1.+4.19237*Q2G)*exp(-3.4*Q2G)
 

      GO TO 20
c ***************************************************
 10   CONTINUE
	Phi_R = 19.
        PhiR=Phi_R*Pi/180.
	COSR=DCOS(PhiR)

        A10=-72.362/COSR   *X1D13n
        A30=-145.620/COSR  *X3D13n
        S10= 12.85/COSR    *XSD13n

      A1= A10*(1.-0.533924*Q2G)*exp(-1.55*Q2G)
      A3= A30*(1.+0.578587*Q2G)*exp(-1.75*Q2G)
      S1= S10*(1.+15.74199*Q2G)*exp(-1.5738*Q2G)
c ***************************************************
20    AE=-(1./2.)*(sqrt(3.)*A3+A1)
      AM=-(1./2.)*(A3/Sqrt(3.0)-A1)
      AS=-(1./2.)*Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123)  Phi_R, X1D13p, A1, X3D13p, A3, XSD13p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123)  Phi_R, X1D13n, A1, X3D13n, A3, XSD13n, S1
123   format( 5x,'D13(1520):',1x,F6.2,3(1x,2(F8.3,1x)))
        return
	end

      subroutine HD33(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
c*****************************************************************
        PI=3.1415926536D0
	  Phi_R = 61.
        PhiR = Phi_R*Pi/180.
	  COSR=DCOS(PhiR)

        A10=109.631/COSR *X1D33        
        A30=101.914/COSR *X3D33
	  S10=1./COSR *XSD33
c *************************************************

      A1= A10*(1.+1.906628*Q2G)*exp(-1.77207*Q2G)
      A3= A30*(1.+1.9722*Q2G)*exp(-2.2*Q2G)
      S1= S10*exp(-2.0*Q2G)
      
      AE=-(1./2.)*(A3*sqrt(3.)+A1)
      AM=-(1./2.)*(A3/sqrt(3.)-A1)
      AS=-(1./2.)*Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      write (6,123) Phi_R, X1D33, A1, X3D33, A3, XSD33, S1
123   format( 5x,'D33(1700):',1x,F6.2,3(1x,2(F8.3,1x)))

      Return
	End

      subroutine HS11f(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
       COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
c *****************************************************************
      PI=3.1415926536D0
	Phi_R=8.193
	PhiR=Phi_R*PI/180.
	COSR=DCOS(PhiR)

      IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

      A10=65.751/COSR   *X1S11p
      S10=-2.0/COSR     *XSS11p  

 
      A1=A10*(1.+1.6083226*Q2G)*exp(-0.70*Q2G)
      S1=S10*(1.+23.90148*Q2G)*exp(-0.81*Q2G)       
      
      GO TO 20
c **************************************************

10      A10=-50.148/COSR *X1S11n
	S10=28.18/COSR   *XSS11n

      A1=A10*(1.+4.746117*Q2G)*exp(-1.68723*Q2G)
      S1=S10*(1.+0.35874*Q2G)*exp(-1.55*Q2G)

c ****************************************************
20    AE=-A1
      AS=-sqrt(2.)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123)  Phi_R, X1S11p, A1, XSS11p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123   format( 5x,'S11(1535):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

        return
	end

      subroutine HS11s(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n
c *****************************************************************
      PI=3.1415926536D0
	Phi_R=6.961
	PhiR=Phi_R*PI/180.
	COSR=DCOS(PhiR)

      IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

      A10=33.0210/COSR *X1S11p
      S10=-3.489/COSR  *XSS11p

      A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
      S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

      GO TO 20
c **************************************************

10    A10=9.186/COSR *X1S11n
      S10=10./COSR   *XSS11n

      A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
      S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

c ********************************************
20    AE=-A1
      AS=-Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123) Phi_R, X1S11p, A1, XSS11p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123   format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

        return
	end
c

      subroutine HS31(Q2G,PhiR,AE,AS,A1,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
c *****************************************************************
        PI=3.1415926536D0
	Phi_R=22.54
	PhiR=Phi_R*PI/180.
	COSR=DCOS(PhiR)

        A10 = 60.6258/COSR *X1S31
        S10 =15./COSR      *XSS31

      A1=A10*(1.+1.858125*Q2G)*exp(-2.5*Q2G)
      S1=S10*(1.+2.82996*Q2G)*exp(-2.*Q2G)      

      AE=-A1
      AS=-Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      write (6,123) Phi_R,  X1S31, A1, XSS31, S1
123   format( 5x,'S31(1620):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))
        return
	end
c
c
      subroutine HF15(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
c ***************************************************************
      PI=3.1415926536D0
      IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
	Phi_R = 10.
      PhiR = Phi_R*Pi/180.
	COSR=DCOS(PhiR)

        A10=-24.7443/COSR  *X1F15p
        A30=132.2624/COSR  *X3F15p 
        S10=-43.2904/COSR  *XSF15p   


      A1= A10*(1.+3.978924*Q2G )*exp(-1.2*Q2G)
      A3= A30*(1.+0.996276*Q2G )*exp(-2.22357*Q2G)      
      S1= S10*(1.+3.138554*Q2G )*exp(-1.58*Q2G)
   

      GO TO 20
c ***************************************************************
10    CONTINUE
	Phi_R = 15.
      PhiR = Phi_R*Pi/180.
	COSR=DCOS(PhiR)

        A10=26.94/COSR   *X1F15n
        A30=-37.07/COSR  *X3F15n
        S10=1./COSR*XSF15n

      A1= A10*(1.+0.001*Q2G)*exp(-1.2*Q2G)
      A3= A30*(1.+4.09308*Q2G)*exp(-1.75*Q2G)
      S1= S10 *exp(-1.55*Q2G)    
c *******************************************************
20    AE=-(1./3.)*( A3*sqrt(2.)+A1)
      AM=-(1./3.)*( A3/sqrt(2.)-A1)
      AS=-(1./3.)*Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123) Phi_R,  X1F15p, A1, X3F15p, A3, XSF15p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123) Phi_R,  X1F15n, A1, X3F15n, A3, XSF15n, S1
123   format( 5x,'F15(1680):',1x,F6.2,3(1x,2(F8.3,1x)))

      return
	end

      subroutine HP31(Q2G,PhiR,AM,AS,A1,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31

      PI=3.1415926536D0
	Phi_R=35.
      PhiR = Phi_R*PI/180.
	COSR=DCOS(PhiR)
        A10=14.786/COSR*X1P31
        S10=1./COSR*XSP31

      A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
      S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)   

      AM= A1
      AS=-Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      write (6,123) Phi_R, X1P31, A1, XSP31, S1
123   format( 5x,'P31(1910):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

      Return
	End

      subroutine HD15(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n

      PI=3.1415926536D0
      IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
      Phi_R = 20.
      PhiR = Phi_R*Pi/180.
	COSR=DCOS(PhiR)
      A10=14.356/COSR  *X1D15p
      A30=20.322/COSR  *X3D15p
      S10=1./COSR*XSD15p

      
      A1= A10*(1.+0.1*Q2G)*exp(-2.*Q2G)
      A3= A30*(1.+0.1*Q2G)*exp(-2.*Q2G)
      S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)   
      GO TO 20
c ********************************************
10    CONTINUE
      Phi_R = 0.
      PhiR = Phi_R*Pi/180.
	COSR=DCOS(PhiR)
      A10=-61.738/COSR *X1D15n
      A30=-83.868/COSR *X3D15n
	S10=-1./COSR   *XSD15n

      A1= A10*(1.+0.01*Q2G)*exp(-2.*Q2G)
      A3= A30*(1.+0.01*Q2G)*exp(-2.*Q2G)
      S1= S10*(1.+0.01*Q2G)*exp(-2.*Q2G)


c *************************************************
20    AE= (1./3.)*(A3/sqrt(2.)-A1)
      AM=-(1./3.)*(A3*sqrt(2.)+A1)
      AS=-(1./3.)*Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123) Phi_R,  X1D15p, A1, X3D15p, A3, XSD15p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123) Phi_R,  X1D15n, A1, X3D15n, A3, XSD15n, S1
123   format( 5x,'D15(1675):',1x,F6.2,3(1x,2(F8.3,1x)))

      return
	end


      subroutine HF35(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

      PI=3.1415926536D0
      Phi_R =40.
      PhiR =Phi_R*Pi/180.
	COSR=DCOS(PhiR)
	A10=13.934/COSR*X1F35
      A30=-21.427/COSR*X3F35
	S10=1./COSR*XSF35

      A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
      A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
      S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)   

      AE=-(1./3.)*(A3*sqrt(2.)+A1)
      AM=-(1./3.)*(A3/sqrt(2.)-A1)
      AS=-(1./3.)*S1*sqrt(2.)

      IF (IPRN.EQ.0) return
      write (6,123) Phi_R, X1F35, A1, X3F35, A3, XSF35, S1
123   format( 5x,'F35(1905):',1x,F6.2,3(1x,2(F8.3,1x)))

      return
	end

      subroutine HF37(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

        PI=3.1415926536D0
        Phi_R = 30.
        PhiR = Phi_R*Pi/180.
	  COSR=DCOS(PhiR)

        A10=-81.06/COSR*X1F37
        A30=-104.65/COSR*X3F37
	S10=1./COSR*XSF37
      
      A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
      A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
      S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)   

      AE= (1./4.)*(A3*3./sqrt(15.)-A1)
      AM=-(1./4.)*(A3*5./sqrt(15.)+A1) 
      AS=-(1./4.)*Sqrt(2.)*S1

      IF (IPRN.EQ.0) return
      write (6,123) Phi_R, X1F37, A1, X3F37, A3, XSF37, S1
123   format( 5x,'F37(1950):',1x,F6.2,3(1x,2(F8.3,1x)))

      return
	end

      subroutine HP13(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
      COMMON/parPPn/ X1P13n,X3P13n,XSP13n

        PI=3.1415926536D0
        IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
        Phi_R = 0.
        PhiR = Phi_R*PI/180.
	  COSR = DCOS(PhiR)

        A10=73.002/COSR   *X1P13p
        A30=-11.465/COSR  *X3P13p
	S10=-53.03/COSR   *XSP13p

      A1= A10*(1.+1.891651*Q2G)*DEXP(-1.55*Q2G)
      A3= A30*(1.+15.9896*Q2G)*DEXP(-1.55*Q2G)
      S1= S10*(1.+2.45819*Q2G)*DEXP(-1.55*Q2G)
      GO TO 20
c ***************************************************
10    CONTINUE
        Phi_R = 0.
        PhiR = Phi_R*PI/180.
	  COSR = DCOS(PhiR)

        A10=-2.904/COSR  *X1P13n
        A30=-30.972/COSR *X3P13n
	  S10=-1./COSR   *XSP13n

      A1= A10*(1.+12.72411*Q2G)*DEXP(-1.55*Q2G)
      A3= A30*(1.+4.987*Q2G)*DEXP(-1.55*Q2G)
      S1= S10*DEXP(-1.55*Q2G)
c ****************************************************

20    AE= (1./2.)*(A3/sqrt(3.)-A1)
      AM=-(1./2.)*(A3*sqrt(3.)+A1)
      AS=-(1./2.)*Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123) Phi_R, X1P13p, A1, X3P13p, A3, XSP13p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123) Phi_R, X1P13n, A1, X3P13n, A3, XSP13n, S1
123   format( 5x,'P13(1720):',1x,F6.2,3(1x,2(F8.3,1x)))

        return
	end

      subroutine HS11trd(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
      IMPLICIT REAL*8 (A-H,O-Z)
c ***************************************************************
       COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n
c *****************************************************************
      PI=3.1415926536D0
	Phi_R=0.
	PhiR=Phi_R*PI/180.
	COSR=DCOS(PhiR)

      IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

      A10=28.0/COSR !  *X1S11p
      S10=-3.489/COSR  !*XSS11p

      A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
      S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

      GO TO 20
c **************************************************

10    A10=6.32/COSR *0 ! *X1S11n
      S10=10./COSR    !*XSS11n

      A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
      S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

c ********************************************
20    AE=-A1
      AS=-Sqrt(2.0)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123) Phi_R, X1S11p, A1, XSS11p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123   format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

        return
	end
c	

      subroutine HP11sec(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
       COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
        PI=3.1415926536D0
	  Phi_R = 0
        PhiR=Phi_R*Pi/180.
	  COSR=DCOS(PhiR)
        IF (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
c ***************************************************************
      A10 = -30.134/COSR  ! *X1P11p
      S10 =  0.1/COSR     ! *XSP11p 
      
      A1=A10*(1.+ 0.001*Q2G)*exp(-1.51*Q2G)
      S1= S10*(1.+ 0.001*Q2G)*exp(-1.75*Q2G)            

      GO TO 20
c ***************************************
10    A10 = 0.1/COSR ! *X1P11n
      S10 = 0.1/COSR   ! *XSP11n
      A1= A10*(1.+0.001*Q2G)*exp(-1.77*Q2G)
      S1=S10*(1.+ 0.001*Q2G)*exp(-1.55*Q2G)
c **************************************************
20    AM=A1
      AS=-sqrt(2.)*S1

      IF (IPRN.EQ.0) return
      IF (ISO.EQ.1.OR.ISO.EQ.3)
     & write (6,123)  Phi_R, X1P11p, A1, XSP11p, S1
      IF (ISO.EQ.2.OR.ISO.EQ.4)
     & write (6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123   format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

        return
	end

      SUBROUTINE HEL_OUT(ISO,Q2G)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION mi, kgcm0
       IF (ISO.EQ.1.OR.ISO.EQ.3) write (6,2) Q2G
 2    FORMAT(5x,'proton e.m. helicity amplitudes at Q^2=',F7.3,2x,
     &'in units 10^-3/sqrt(GeV)')
       IF (ISO.EQ.2.OR.ISO.EQ.4) write (6,3) Q2G
 3    FORMAT(5x,'neutron e.m. helicity amplitudes at Q^2=',F7.3,2x,
     &'in units 10^-3/sqrt(GeV)')
      write (6,4)
 4    format(17X,'Phi_R',6x,'X1',5x,'A1/2',8x,'X3',5x,'A3/2',8x,'XS',
     & 5x,'S1/2')

	W0=1.232
	mi=0.9382723
        kgcm0 = (W0*W0-mi*mi)/2./W0
        egcm = (W0*W0-Q2G-mi*mi)/2./W0
	qcm=dsqrt(egcm**2+Q2G)

      CALL HP33(Q2G,qcm,kgcm0,DE,DM,DS,A1,A3,S1,1)
      CALL HP11(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
      CALL HD13(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
      CALL HS11f(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
      CALL HS31(Q2G,PhiR,DM,DS,A1,S1,1)
      CALL HS11s(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
      CALL HD15(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
      CALL HF15(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
      CALL HD33(Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
      CALL HP13(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
      CALL HF35(Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
      CALL HP31(Q2G,PhiR,DM,DS,A1,S1,1)
      CALL HF37(Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)

      RETURN
      END
C
      SUBROUTINE InvFtoA(WCM,XM1,XM2,XMPI,Q2,CK0,CK,QPI,X,
     &  S,T,F1,F2,F3,F4,F5,F6,A1,A2,A3,A4,A5,A6)
c     &  F1,F2,F3,F4,F5,F6)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	REAL*8 k,k0,M,mpi
	COMPLEX*16 F1,F2,F3,F4,F5,F6,A1,A2,A3,A4,A5,A6
	COMPLEX*16 A10,A20,A30,A40,A50,A60
	COMPLEX*16 F1t,F2t,F3t,F4t,F5t,F6t,F3test,F4test
c	COMMON /INVARI/ S,T,A1,A2,A3,A4,A5,A6
C ***************************************************************
C Definition of the invariant amplitudes is the same as in
c      CGLN                  Phys. Rev.  106, 1345 (1957)
c      DENNERY P.            Phys. Rev. 124, 2000 (1961)
c      BERENDS F.A.  at al   Nucl. Phys. B4, 1 (1967)
c      GEHLEN G.v            Nucl. Phys. B9, 17 (1969)
c  Expressions used below are from Barbara Pasquini, 16 Dec 2003
c ****************************************************************
      PI=3.141592654D0
      XMN=(XM1+XM2)/2.
      ECPI=SQRT(QPI**2+XMPI**2)
	s=WCM**2
	t=2*CK*QPI*X-2*CK0*ECPI+XMPI**2-Q2
	q=QPI
	q0=ECPI
	k=CK
	k0=CK0
	Ei=Wcm-k0
	Ef=Wcm-q0
	Wp=wcm+xmn
	Wm=wcm-xmn
	Eip=Ei+xmn
	Eim=Ei-xmn
	Efp=Ef+xmn
	Efm=Ef-xmn
	W=wcm
	M=xmn
	mpi=xmpi

      F1t=8.*PI*W/(Wm*sqrt(Eip*Efp))*F1
      F2t=8.*PI*W*sqrt(Eip*Efp)/(q*k*Wp)*F2
      F3t=8.*PI*W*sqrt(Eip/Efp)/(q*k*Wp)*F3
      F4t=8.*PI*W*sqrt(Efp/Eip)/(q**2*Wm)*F4
      F5t=8.*PI*W*sqrt(Eip/Efp)/(k**2)*F5
      F6t=8.*PI*W*sqrt(Efp/Eip)/(q*k)*F6
  
      alpha=k0*(t-mpi**2+Q2)-2*q0*Q2
      beta=s-M**2 + Q2/2.0
  
      A1=( 2*Wm*(W+M*Wm/Eim)*F1t - Wp/Eip*(s-M**2+Q2)*F2t
     & + M/k**2*alpha*(Wp*F3t + Wm*F4t) 
     & + 2*M/k0*Q2*(F5t + F6t) )/(4*W**2)

      A2=( -Wm/Eim*F1t + Wp/Eip*F2t
     & + (s-M**2)/(2*k**2*Q2)*alpha*(F3t - F4t)
     & + Wm/k0*F5t - Wp/k0*F6t )*Q2/(2*W**2*(t-mpi**2))

      A3=( Wm**2/Eim*F1t + Wp**2/Eip*F2t 
     & + (alpha/k**2+4*W)/2*(Wp*F3t + Wm*F4t) 
     & + Q2/k0*(F5t + F6t) )/(4*W**2)

      A4=( Wm**2/Eim*F1t + Wp**2/Eip*F2t 
     & + alpha/(2*k**2)*(Wp*F3t + Wm*F4t) 
     & + Q2/k0*(F5t + F6t) )/(4*W**2)

      A5=( -Wm/Eim*beta*F1t + Wp/Eip*beta*F2t
     & + (s-M**2)/(2*k**2)*((t-mpi**2+Q2)*(3*k0/2.D0-2*W)-2*q0*beta
     &  +2*k**2*W)*(F3t - F4t) 
     & + beta/k0*(Wm*F5t - Wp*F6t) )/(2*W**2*(t-mpi**2))

      A6=( F1t/Eim + F2t/Eip  + ((t-mpi**2+Q2)*(2*W-k0)
     &  +(W**2-M**2)*2*q0)/(2*k**2)*(F3t/Wm + F4t/Wp)
     & - 1/k0*(F5t + F6t) )*(s-M**2)/(4*W**2)

      return
	end
C
      SUBROUTINE InvAtoF(WCM,XM1,XM2,XMPI,Q2,CK0,CK,QPI,X,
     &  S,T,A1,A2,A3,A4,A5,A6,F1,F2,F3,F4,F5,F6)
c     &  F1,F2,F3,F4,F5,F6)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	REAL*8 k,k0,M,mpi
	COMPLEX*16 F1,F2,F3,F4,F5,F6,A1,A2,A3,A4,A5,A6
	COMPLEX*16 F1t,F2t,F3t,F4t,F5t,F6t
	COMPLEX*16 F1test,F2test,F3test,F4test,F5test,F6test
c	COMMON /INVARI/ S,T,A1,A2,A3,A4,A5,A6
C ***************************************************************
C Definition of the invariant amplitudes is the same as in
c      CGLN                  Phys. Rev.  106, 1345 (1957)
c      DENNERY P.            Phys. Rev. 124, 2000 (1961)
c      BERENDS F.A.  at al   Nucl. Phys. B4, 1 (1967)
c      GEHLEN G.v            Nucl. Phys. B9, 17 (1969)
c  Expressions used below are from D. Drechsel, 12 February 2005
c ****************************************************************
      PI=3.141592654D0
      XMN=(XM1+XM2)/2.
      ECPI=SQRT(QPI**2+XMPI**2)
	s=WCM**2
	t=2*CK*QPI*X-2*CK0*ECPI+XMPI**2-Q2
	q=QPI
	q0=ECPI
	k=CK
	k0=CK0
	Ei=Wcm-k0
	Ef=Wcm-q0
	Wp=wcm+xmn
	Wm=wcm-xmn
	Eip=Ei+xmn
	Eim=Ei-xmn
	Efp=Ef+xmn
	Efm=Ef-xmn
	W=wcm
	M=xmn
	mpi=xmpi
  
      alpha=t-mpi**2+Q2
      beta=s-M**2+Q2/2.D0
  
      F1t= A1 + Wm*A4 - alpha/2.D0/Wm*(A3-A4) + Q2/Wm*A6 
      F2t=-A1 + Wp*A4 - alpha/2.D0/Wp*(A3-A4) + Q2/Wp*A6 
      F3t= beta/Wp*A2 + A3 - A4 - Q2/Wp*A5
      F4t=-beta/Wm*A2 + A3 - A4 + Q2/Wm*A5

	F5t= Eip*A1 + (alpha*(W-3.D0/4.D0*k0)-k**2*W+q0*beta)*A2
     &    + (q0*Wp+alpha/2)*A3 + (Eip*Wm-q0*Wp-alpha/2)*A4
     &    + (k0/2*alpha-q0*Q2)*A5 - Eip*Wm*A6

	F6t=-Eim*A1 - (alpha*(W-3.D0/4.D0*k0)-k**2*W+q0*beta)*A2
     &    + (q0*Wm+alpha/2)*A3 + (Eim*Wp-q0*Wm-alpha/2)*A4
     &    - (k0/2*alpha-q0*Q2)*A5 - Eim*Wp*A6

      F1test=Wm/(8*PI*W)*sqrt(Eip*Efp)*F1t
      F2test=Wp/(8*PI*W)*sqrt(Eim/Efp)*F2t*q
      F3test=Wp/(8*PI*W)*sqrt(Eim*Efp)*F3t*q
      F4test=Wm/(8*PI*W)*sqrt(Eip/Efp)*F4t*q**2
      F5test=k0/(8*PI*W)*sqrt(Efp/Eip)*F5t
      F6test=k0/(8*PI*W)/sqrt(Efp*Eim)*F6t*q

      F1=F1test
      F2=F2test
      F3=F3test
      F4=F4test
      F5=F5test
      F6=F6test

      return
	end

c
      SUBROUTINE DispPole(WCM,XM1,XM2,XMPI,Q2,CK0,CK,QPI,X,
     &  F1,F2,F3,F4,F5,F6,A5xp)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	REAL*8 k,k0,M,mpi
	COMPLEX*16 F1,F2,F3,F4,F5,F6,A5xp
	COMPLEX*16 F1t,F2t,F3t,F4t,F5t,F6t
      PI=3.141592654D0
      XMN=(XM1+XM2)/2.
      ECPI=SQRT(QPI**2+XMPI**2)
	s=WCM**2
	t=2*CK*QPI*X-2*CK0*ECPI+XMPI**2-Q2
	q=QPI
	q0=ECPI
	k=CK
	k0=CK0
	Ei=Wcm-k0
	Ef=Wcm-q0
	Wp=wcm+xmn
	Wm=wcm-xmn
	Eip=Ei+xmn
	Eim=Ei-xmn
	Efp=Ef+xmn
	Efm=Ef-xmn
	W=wcm
	M=xmn
	mpi=xmpi

      F1t=8.*PI*W/(Wm*sqrt(Eip*Efp))*F1
      F2t=8.*PI*W*sqrt(Eip*Efp)/(q*k*Wp)*F2
      F3t=8.*PI*W*sqrt(Eip/Efp)/(q*k*Wp)*F3
      F4t=8.*PI*W*sqrt(Efp/Eip)/(q**2*Wm)*F4
      F5t=8.*PI*W*sqrt(Eip/Efp)/(k**2)*F5
      F6t=8.*PI*W*sqrt(Efp/Eip)/(q*k)*F6
  
      beta=s-M**2 + Q2/2.0
  
      A5xp=( -Wm/Eim*beta*F1t + Wp/Eip*beta*F2t
     & + (s-M**2)/(2*k**2)*(Q2*(3*k0/2.D0-2*W)-2*q0*beta
     &  +2*k**2*W)*(F3t - F4t) 
     & + beta/k0*(Wm*F5t - Wp*F6t) )/(2*W**2)

      return
	end
C	

